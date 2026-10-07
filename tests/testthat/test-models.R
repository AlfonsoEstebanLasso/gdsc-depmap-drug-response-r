# The five model families on both tasks with tiny grids: tuning tables,
# chosen hyperparameters, predictions and the shared inner folds.

m <- make_synthetic_matrix()
train <- 1:45
test <- 46:60

prepared <- list(
  regression = preprocess_fold(m$x, m$auc, train, test, "regression", n_genes = 200),
  classification = preprocess_fold(m$x, m$auc, train, test, "classification", n_genes = 200)
)
foldids <- list(
  regression = make_inner_folds(prepared$regression$y_train, 3, "regression", seed = 3),
  classification = make_inner_folds(prepared$classification$y_train, 3, "classification", seed = 3)
)

fit_one <- function(model, task, foldid = foldids[[task]], seed = 11) {
  pp <- prepared[[task]]
  grid <- tiny_grids(task)[[model]]
  list(grid = grid,
       fitted = quiet_small_classes(
         fit_model(model, pp$x_train, pp$y_train, task, foldid, grid, seed = seed, threads = 1)))
}

# A value of the ranger call, read from the object when stored there and
# from the recorded call otherwise.
ranger_setting <- function(fit, name) {
  if (!is.null(fit[[name]])) return(fit[[name]])
  value <- as.list(fit$call)[[name]]
  if (is.null(value)) return(NULL)
  if (is.numeric(value)) return(value)
  tryCatch(eval(value), error = function(e) NULL)
}

test_that("model_names and model_grid declare the five families and their grids", {
  expect_equal(model_names(), c("rf", "ridge", "lasso", "enet", "svm_linear"))
  for (task in c("regression", "classification")) {
    rf <- model_grid("rf", task, rf_budget = 4, seed = 1, num_trees = 50)
    expect_s3_class(rf, "data.frame")
    expect_equal(nrow(rf), 4)
    expect_true(all(c("mtry", "min.node.size", "sample.fraction", "max.depth", "num.trees") %in% names(rf)))
    expect_true(all(rf$num.trees == 50))
    expect_true(all(rf$mtry %in% c(5, 10, 14, 20, 35, 50)))
    expect_true(all(rf$min.node.size %in% c(1, 3, 5, 10)))
    expect_true(all(rf$sample.fraction %in% c(0.55, 0.632, 0.75)))
    expect_true(all(rf$max.depth %in% c(0, 20, 30, 40)))
    expect_identical(rf, model_grid("rf", task, rf_budget = 4, seed = 1, num_trees = 50))
    expect_equal(model_grid("ridge", task)$alpha, 0)
    expect_equal(model_grid("lasso", task)$alpha, 1)
    expect_equal(model_grid("enet", task)$alpha, seq(0.1, 0.9, by = 0.1))
    svm <- model_grid("svm_linear", task)
    if (task == "regression") {
      expect_equal(nrow(svm), 12)
      expect_setequal(unique(svm$cost), c(0.01, 0.1, 1, 10))
      expect_setequal(unique(svm$epsilon), c(0.05, 0.1, 0.2))
    } else {
      expect_equal(svm$cost, c(0.01, 0.1, 1, 10, 100))
    }
  }
})

test_that("every family fits, tunes within its grid and predicts the right type", {
  for (task in c("regression", "classification")) {
    pp <- prepared[[task]]
    for (model in model_names()) {
      res <- fit_one(model, task)
      fitted <- res$fitted
      grid <- res$grid
      expect_false(is.null(fitted$fit), info = paste(model, task))
      expect_equal(fitted$model, model)
      expect_equal(fitted$task, task)
      expect_type(fitted$best, "list")
      expect_s3_class(fitted$tuning, "data.frame")
      expect_true("inner_score" %in% names(fitted$tuning), info = paste(model, task))
      expect_equal(nrow(fitted$tuning), nrow(grid), info = paste(model, task))
      expect_true(all(is.finite(fitted$tuning$inner_score)), info = paste(model, task))
      for (p in intersect(names(fitted$best), names(grid))) {
        expect_true(fitted$best[[p]] %in% grid[[p]], info = paste(model, task, p))
      }
      expect_true(is.numeric(fitted$elapsed))

      pred <- predict_model(fitted, pp$x_test)
      if (task == "regression") {
        expect_true(is.numeric(pred), info = paste(model, task))
        expect_length(pred, length(test))
        expect_true(all(is.finite(pred)))
      } else {
        expect_type(pred, "list")
        expect_true(all(c("prob", "class", "score") %in% names(pred)))
        expect_length(pred$prob, length(test))
        expect_length(pred$class, length(test))
        expect_length(pred$score, length(test))
        if (model == "svm_linear") {
          # The linear SVM has no probability model: its ROC score is the
          # oriented decision value (section 7.4) and prob may be NA.
          expect_true(all(is.na(pred$prob) | (pred$prob >= 0 & pred$prob <= 1)),
                      info = paste(model, task))
        } else {
          expect_true(all(pred$prob >= 0 & pred$prob <= 1), info = paste(model, task))
        }
        expect_true(is.factor(pred$class))
        expect_equal(levels(pred$class), c("resistant", "sensitive"))
        expect_true(is.numeric(pred$score))
      }
    }
  }
})

test_that("the random forest refit carries exactly the tuned hyperparameters", {
  for (task in c("regression", "classification")) {
    fitted <- fit_one("rf", task)$fitted
    fit <- fitted$fit
    expect_s3_class(fit, "ranger")
    best <- fitted$best
    expect_equal(fit$mtry, best$mtry)
    expect_equal(fit$min.node.size, best$min.node.size)
    expect_equal(fit$num.trees, 50)
    expect_equal(fit$num.trees, best$num.trees)
    expect_equal(ranger_setting(fit, "max.depth"), best$max.depth)
    sample_fraction <- ranger_setting(fit, "sample.fraction")
    if (is.null(sample_fraction)) {
      skip("sample.fraction is not recoverable from the ranger object or its call")
    }
    expect_equal(sample_fraction, best$sample.fraction)
    expect_false(fit$replace)
    if (task == "classification") expect_equal(fit$treetype, "Probability estimation")
    if (task == "regression") expect_equal(fit$treetype, "Regression")
  }
})

test_that("the glmnet models use lambda.min and the elastic net alpha comes from the grid", {
  for (task in c("regression", "classification")) {
    for (model in c("ridge", "lasso", "enet")) {
      fitted <- fit_one(model, task)$fitted
      expect_s3_class(fitted$fit, "cv.glmnet")
      expect_true("alpha" %in% names(fitted$best))
      expect_true(fitted$best$alpha %in% tiny_grids(task)[[model]]$alpha)
      if (!is.null(fitted$best$lambda)) {
        expect_equal(fitted$best$lambda, fitted$fit$lambda.min)
      }
      if (model == "enet") {
        expect_equal(nrow(fitted$tuning), 2)
      }
    }
  }
})

test_that("the SVM uses the 200-gene matrix", {
  for (task in c("regression", "classification")) {
    fitted <- fit_one("svm_linear", task)$fitted
    expect_s3_class(fitted$fit, "svm")
    expect_equal(ncol(fitted$fit$SV), 200)
    expect_true(fitted$best$cost %in% tiny_grids(task)$svm_linear$cost)
    if (task == "regression") {
      expect_true(fitted$best$epsilon %in% tiny_grids(task)$svm_linear$epsilon)
    }
  }
})

test_that("the inner folds drive the tuning of every family", {
  for (task in c("regression", "classification")) {
    foldid <- foldids[[task]]
    other <- rev(foldid)
    for (model in model_names()) {
      same <- fit_one(model, task, foldid = foldid)$fitted$tuning
      again <- fit_one(model, task, foldid = foldid)$fitted$tuning
      different <- fit_one(model, task, foldid = other)$fitted$tuning
      expect_equal(same$inner_score, again$inner_score, info = paste(model, task))
      expect_false(isTRUE(all.equal(same$inner_score, different$inner_score)),
                   info = paste(model, task))
    }
  }
})
