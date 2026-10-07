# ----------------------------------------------------------------------------
# R/models.R: the five model families of v1.0 with hand-rolled inner tuning
#
# Purpose: fit the random forest (ranger), ridge, lasso and elastic net
#   (glmnet) and the linear SVM (e1071) for regression and classification, with
#   the declared grids and the inner resampling of sections 4.2, 6 and 7 of
#   docs/v1_spec.md, and predict on a held-out matrix.
# Inputs: a standardised training matrix (rows: cell lines, columns: the 200
#   genes selected on the training fold), the training response (numeric AUC
#   or a factor with levels resistant and sensitive), the inner fold vector
#   foldid shared by every family, a grid from model_grid() and the split seed.
# Outputs: fit_model() returns the fitted object, the chosen hyperparameters
#   (best), the search table (tuning, one inner score per candidate) and the
#   elapsed seconds; predict_model() returns a numeric vector (regression) or
#   list(prob, class, score) (classification).
# Transformations: the SVM regression centres and scales the training
#   response with the training mean and standard deviation and transforms the
#   predictions back (section 5, point 4). Nothing else is rescaled here: the
#   matrices arrive standardised from preprocess_fold(), so glmnet is called
#   with standardize = FALSE and svm with scale = FALSE.
# Seeds (section 4.3): seed is the split seed S + 1000 * t + j. ranger uses
#   seed + 1; set.seed(seed + 5) precedes the first glmnet and the first svm
#   fit of each family; the random forest draw of model_grid() uses the seed
#   passed to it (seed + 100 in the runner).
# Errata fixed here: 2 (the forest is fitted with the tuned values), 3 (alpha
#   tuned on a grid, lambda.min everywhere), 4 (declared grids, same inner
#   folds for every family), 7 (every family receives the same matrix).
# Dependencies: ranger, glmnet, e1071 and the functions rmse() and roc_auc()
#   of R/metrics.R, which must be sourced before this file.
# ----------------------------------------------------------------------------

MODEL_NAMES <- c("rf", "ridge", "lasso", "enet", "svm_linear")
TASKS <- c("regression", "classification")

model_names <- function() MODEL_NAMES

check_task <- function(task) match.arg(task, TASKS)

check_model <- function(model) match.arg(model, MODEL_NAMES)

# Declared grids of section 7. The random forest grid (288 configurations)
# is sampled without replacement with rf_budget draws after set.seed(seed);
# the other families return their full grid.
model_grid <- function(model, task, rf_budget = 12, seed = NULL, num_trees = 500) {
  model <- check_model(model)
  task <- check_task(task)
  if (model == "rf") {
    full <- expand.grid(
      mtry = c(5, 10, 14, 20, 35, 50),
      min.node.size = c(1, 3, 5, 10),
      sample.fraction = c(0.55, 0.632, 0.75),
      max.depth = c(0, 20, 30, 40),
      KEEP.OUT.ATTRS = FALSE)
    full$num.trees <- num_trees
    if (is.null(rf_budget) || rf_budget >= nrow(full)) {
      grid <- full
    } else {
      if (rf_budget < 1) stop("rf_budget must be at least 1")
      if (!is.null(seed)) set.seed(seed)
      grid <- full[sort(sample.int(nrow(full), rf_budget, replace = FALSE)), , drop = FALSE]
    }
    rownames(grid) <- NULL
    return(grid)
  }
  if (model == "ridge") return(data.frame(alpha = 0))
  if (model == "lasso") return(data.frame(alpha = 1))
  if (model == "enet") return(data.frame(alpha = seq(0.1, 0.9, by = 0.1)))
  # svm_linear
  if (task == "regression") {
    grid <- expand.grid(cost = c(0.01, 0.1, 1, 10), epsilon = c(0.05, 0.1, 0.2),
                        KEEP.OUT.ATTRS = FALSE)
    rownames(grid) <- NULL
    return(grid)
  }
  data.frame(cost = c(0.01, 0.1, 1, 10, 100))
}

# Inner selection metric (section 7): RMSE for regression (lower is better)
# and ROC AUC for classification (higher is better), where predicted is the
# score oriented towards sensitive.
inner_score <- function(task, observed, predicted) {
  task <- check_task(task)
  if (task == "regression") return(rmse(observed, predicted))
  roc_auc(observed, predicted)
}

score_is_lower_better <- function(task) check_task(task) == "regression"

best_candidate <- function(scores, task) {
  if (all(is.na(scores))) stop("every candidate has a missing inner score")
  if (score_is_lower_better(task)) which.min(scores) else which.max(scores)
}

check_training_inputs <- function(x_train, y_train, task, foldid, grid) {
  if (!is.matrix(x_train) || !is.numeric(x_train)) {
    stop("x_train must be a numeric matrix")
  }
  if (is.null(colnames(x_train))) stop("x_train must have column names (genes)")
  if (anyNA(x_train)) stop("x_train must not contain missing values")
  n <- nrow(x_train)
  if (length(y_train) != n) stop("y_train must have one value per row of x_train")
  if (task == "regression") {
    if (!is.numeric(y_train) || anyNA(y_train)) {
      stop("y_train must be a numeric vector without missing values for regression")
    }
  } else {
    if (!is.factor(y_train) || !identical(levels(y_train), CLASS_LEVELS) || anyNA(y_train)) {
      stop("y_train must be a factor with levels resistant, sensitive for classification")
    }
    if (length(unique(y_train)) < 2) stop("y_train has a single class")
  }
  if (length(foldid) != n) stop("foldid must have one value per row of x_train")
  foldid <- as.integer(foldid)
  if (anyNA(foldid) || any(foldid < 1)) stop("foldid must be positive integers")
  if (length(unique(foldid)) < 2) stop("foldid must define at least two inner folds")
  if (!is.data.frame(grid) || nrow(grid) == 0) stop("grid must be a non-empty data frame")
  foldid
}

# Mean of the inner validation scores over the inner folds for one function
# that fits on rows and predicts a score on other rows.
inner_cv_score <- function(fit_predict, x_train, y_train, task, foldid) {
  folds <- sort(unique(foldid))
  scores <- vapply(folds, function(v) {
    tr <- which(foldid != v)
    va <- which(foldid == v)
    pred <- fit_predict(tr, va)
    inner_score(task, y_train[va], pred)
  }, numeric(1))
  if (all(is.na(scores))) return(NA_real_)
  mean(scores, na.rm = TRUE)
}

# --- random forest ---------------------------------------------------------

ranger_args <- function(candidate, task, seed, threads, p) {
  mtry <- as.integer(candidate$mtry)
  if (mtry > p) {
    warning("mtry = ", mtry, " exceeds the number of genes (", p,
            "); using mtry = ", p)
    mtry <- p
  }
  args <- list(
    num.trees = as.integer(candidate$num.trees),
    mtry = mtry,
    min.node.size = as.integer(candidate$min.node.size),
    sample.fraction = as.numeric(candidate$sample.fraction),
    max.depth = as.integer(candidate$max.depth),
    replace = FALSE,
    importance = "none",
    num.threads = as.integer(threads),
    seed = as.integer(seed))
  if (task == "regression") {
    args$splitrule <- "variance"
  } else {
    args$probability <- TRUE
    args$splitrule <- "gini"
  }
  args
}

fit_ranger <- function(x, y, args) {
  do.call(ranger::ranger, c(list(x = x, y = y), args))
}

predict_ranger_score <- function(fit, x, task, threads = 1) {
  pred <- stats::predict(fit, data = x, num.threads = as.integer(threads))$predictions
  if (task == "regression") return(as.numeric(pred))
  as.numeric(pred[, POSITIVE_CLASS])
}

fit_rf <- function(x_train, y_train, task, foldid, grid, seed, threads) {
  needed <- c("mtry", "min.node.size", "sample.fraction", "max.depth", "num.trees")
  missing_cols <- setdiff(needed, names(grid))
  if (length(missing_cols) > 0) {
    stop("rf grid lacks column(s): ", paste(missing_cols, collapse = ", "))
  }
  p <- ncol(x_train)
  tuning <- grid[needed]
  tuning$inner_score <- NA_real_
  for (i in seq_len(nrow(grid))) {
    args <- ranger_args(grid[i, ], task, seed + 1, threads, p)
    tuning$inner_score[i] <- inner_cv_score(function(tr, va) {
      fit <- fit_ranger(x_train[tr, , drop = FALSE], y_train[tr], args)
      predict_ranger_score(fit, x_train[va, , drop = FALSE], task, threads)
    }, x_train, y_train, task, foldid)
  }
  i_best <- best_candidate(tuning$inner_score, task)
  args <- ranger_args(grid[i_best, ], task, seed + 1, threads, p)
  fit <- fit_ranger(x_train, y_train, args)
  best <- list(mtry = args$mtry, min.node.size = args$min.node.size,
               sample.fraction = args$sample.fraction, max.depth = args$max.depth,
               num.trees = args$num.trees)
  list(fit = fit, best = best, tuning = tuning)
}

# --- glmnet: ridge, lasso, elastic net -------------------------------------

glmnet_measure <- function(task) if (task == "regression") "mse" else "auc"

glmnet_family <- function(task) if (task == "regression") "gaussian" else "binomial"

# cv.glmnet with the shared foldid; the inner score is the cross-validated
# measure at lambda.min, reported as RMSE for regression and AUC for
# classification. If glmnet replaces the AUC measure (too few rows per fold),
# the fallback measure is reported and named in the tuning table.
fit_cv_glmnet <- function(x, y, alpha, foldid, task) {
  cv <- glmnet::cv.glmnet(x = x, y = y, alpha = alpha, foldid = foldid,
                          nlambda = 100, standardize = FALSE,
                          family = glmnet_family(task),
                          type.measure = glmnet_measure(task))
  idx <- which(cv$lambda == cv$lambda.min)[1]
  measure <- names(cv$name)
  score <- cv$cvm[idx]
  if (task == "regression" && identical(measure, "mse")) score <- sqrt(score)
  list(cv = cv, score = score, measure = measure)
}

fit_glmnet_family <- function(x_train, y_train, task, foldid, grid, seed) {
  if (!"alpha" %in% names(grid)) stop("glmnet grid lacks the column alpha")
  if (any(grid$alpha < 0 | grid$alpha > 1)) stop("alpha must lie in [0, 1]")
  set.seed(seed + 5)
  tuning <- data.frame(alpha = as.numeric(grid$alpha), lambda = NA_real_,
                       measure = NA_character_, inner_score = NA_real_,
                       stringsAsFactors = FALSE)
  fits <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    res <- fit_cv_glmnet(x_train, y_train, grid$alpha[i], foldid, task)
    fits[[i]] <- res$cv
    tuning$lambda[i] <- res$cv$lambda.min
    tuning$measure[i] <- res$measure
    tuning$inner_score[i] <- res$score
  }
  lower_better <- if (task == "regression") TRUE else !all(tuning$measure == "auc")
  i_best <- if (lower_better) which.min(tuning$inner_score) else which.max(tuning$inner_score)
  fit <- fits[[i_best]]
  best <- list(alpha = tuning$alpha[i_best], lambda = fit$lambda.min,
               lambda_rule = "lambda.min")
  list(fit = fit, best = best, tuning = tuning)
}

predict_glmnet_score <- function(fit, x, task) {
  if (task == "regression") {
    return(as.numeric(stats::predict(fit, newx = x, s = "lambda.min")))
  }
  as.numeric(stats::predict(fit, newx = x, s = "lambda.min", type = "response"))
}

# --- linear SVM ------------------------------------------------------------

# Score oriented towards sensitive (section 7.4): e1071 returns decision
# values positive for the class named first in the column name.
svm_oriented_score <- function(pred) {
  dv <- attr(pred, "decision.values")
  if (is.null(dv)) stop("the SVM prediction carries no decision values")
  score <- as.numeric(dv[, 1])
  first <- colnames(dv)[1]
  if (!is.null(first) && startsWith(first, CLASS_LEVELS[1])) score <- -score
  score
}

fit_svm_regression <- function(x, y, cost, epsilon) {
  centre <- mean(y)
  scale <- stats::sd(y)
  if (is.na(scale) || scale == 0) scale <- 1
  fit <- e1071::svm(x = x, y = (y - centre) / scale, kernel = "linear",
                    scale = FALSE, type = "eps-regression",
                    cost = cost, epsilon = epsilon)
  list(fit = fit, y_center = centre, y_scale = scale)
}

predict_svm_regression <- function(obj, x) {
  as.numeric(stats::predict(obj$fit, newdata = x)) * obj$y_scale + obj$y_center
}

fit_svm_classification <- function(x, y, cost) {
  e1071::svm(x = x, y = y, kernel = "linear", scale = FALSE,
             type = "C-classification", cost = cost)
}

predict_svm_classification <- function(fit, x) {
  pred <- stats::predict(fit, newdata = x, decision.values = TRUE)
  list(class = factor(as.character(pred), levels = CLASS_LEVELS),
       score = svm_oriented_score(pred))
}

fit_svm_linear <- function(x_train, y_train, task, foldid, grid, seed) {
  set.seed(seed + 5)
  if (task == "regression") {
    needed <- c("cost", "epsilon")
  } else {
    needed <- "cost"
  }
  missing_cols <- setdiff(needed, names(grid))
  if (length(missing_cols) > 0) {
    stop("svm grid lacks column(s): ", paste(missing_cols, collapse = ", "))
  }
  tuning <- grid[needed]
  tuning$inner_score <- NA_real_
  for (i in seq_len(nrow(grid))) {
    tuning$inner_score[i] <- inner_cv_score(function(tr, va) {
      if (task == "regression") {
        obj <- fit_svm_regression(x_train[tr, , drop = FALSE], y_train[tr],
                                  grid$cost[i], grid$epsilon[i])
        predict_svm_regression(obj, x_train[va, , drop = FALSE])
      } else {
        fit <- fit_svm_classification(x_train[tr, , drop = FALSE], y_train[tr],
                                      grid$cost[i])
        predict_svm_classification(fit, x_train[va, , drop = FALSE])$score
      }
    }, x_train, y_train, task, foldid)
  }
  i_best <- best_candidate(tuning$inner_score, task)
  if (task == "regression") {
    obj <- fit_svm_regression(x_train, y_train, grid$cost[i_best], grid$epsilon[i_best])
    best <- list(cost = grid$cost[i_best], epsilon = grid$epsilon[i_best])
    return(list(fit = obj$fit, best = best, tuning = tuning,
                y_center = obj$y_center, y_scale = obj$y_scale))
  }
  fit <- fit_svm_classification(x_train, y_train, grid$cost[i_best])
  list(fit = fit, best = list(cost = grid$cost[i_best]), tuning = tuning)
}

# --- public interface ------------------------------------------------------

# Run the inner loop over the grid, refit the best candidate on all training
# rows and return the fit, the chosen hyperparameters and the search table.
fit_model <- function(model, x_train, y_train, task, foldid, grid, seed, threads = 1) {
  model <- check_model(model)
  task <- check_task(task)
  foldid <- check_training_inputs(x_train, y_train, task, foldid, grid)
  if (!is.numeric(seed) || length(seed) != 1 || is.na(seed)) {
    stop("seed must be a single number")
  }
  seed <- as.integer(seed)
  t0 <- proc.time()[["elapsed"]]
  res <- switch(model,
    rf = fit_rf(x_train, y_train, task, foldid, grid, seed, threads),
    ridge = fit_glmnet_family(x_train, y_train, task, foldid, grid, seed),
    lasso = fit_glmnet_family(x_train, y_train, task, foldid, grid, seed),
    enet = fit_glmnet_family(x_train, y_train, task, foldid, grid, seed),
    svm_linear = fit_svm_linear(x_train, y_train, task, foldid, grid, seed))
  elapsed <- proc.time()[["elapsed"]] - t0
  out <- list(model = model, task = task, fit = res$fit, best = res$best,
              tuning = res$tuning, elapsed = elapsed,
              genes = colnames(x_train), n_train = nrow(x_train),
              foldid = foldid, seed = seed, threads = as.integer(threads))
  if (!is.null(res$y_center)) {
    out$y_center <- res$y_center
    out$y_scale <- res$y_scale
  }
  out
}

# Align the held-out matrix with the training genes.
align_columns <- function(fitted, x_test) {
  if (!is.matrix(x_test) || !is.numeric(x_test)) stop("x_test must be a numeric matrix")
  if (is.null(colnames(x_test))) stop("x_test must have column names (genes)")
  missing_genes <- setdiff(fitted$genes, colnames(x_test))
  if (length(missing_genes) > 0) {
    stop("x_test lacks ", length(missing_genes), " training gene(s), for example ",
         missing_genes[1])
  }
  x_test[, fitted$genes, drop = FALSE]
}

# Regression: numeric vector. Classification: list(prob, class, score) with
# score oriented towards sensitive; prob is NA for the SVM, which has no
# probability model (section 7.4 uses its decision value).
predict_model <- function(fitted, x_test) {
  x_test <- align_columns(fitted, x_test)
  model <- fitted$model
  task <- fitted$task
  if (task == "regression") {
    pred <- switch(model,
      rf = predict_ranger_score(fitted$fit, x_test, task, fitted$threads),
      svm_linear = predict_svm_regression(fitted, x_test),
      predict_glmnet_score(fitted$fit, x_test, task))
    return(as.numeric(pred))
  }
  if (model == "svm_linear") {
    res <- predict_svm_classification(fitted$fit, x_test)
    return(list(prob = rep(NA_real_, nrow(x_test)), class = res$class,
                score = res$score))
  }
  prob <- if (model == "rf") {
    predict_ranger_score(fitted$fit, x_test, task, fitted$threads)
  } else {
    predict_glmnet_score(fitted$fit, x_test, task)
  }
  cls <- factor(ifelse(prob >= 0.5, CLASS_LEVELS[2], CLASS_LEVELS[1]),
                levels = CLASS_LEVELS)
  list(prob = prob, class = cls, score = prob)
}
