# ----------------------------------------------------------------------------
# R/metrics.R: evaluation metrics and fold summaries of v1.0
#
# Purpose: compute the held-out metrics of section 9 of docs/v1_spec.md and
#   aggregate them over the outer folds with t-based confidence intervals.
# Inputs: observed and predicted values of one held-out fold (numeric for
#   regression, factors with levels resistant and sensitive for
#   classification) and the metrics_by_fold table written by run_drug().
# Outputs: named numeric vectors (one fold) and a summary data frame.
# Transformations: none of the inputs is modified; the ROC AUC uses pROC with
#   fixed levels and direction (section 7.4 of the specification), so the score
#   must increase with the probability of the sensitive class.
# Erratum 1 of docs/errata.md is fixed here: R squared is 1 minus SSres over
#   SStot on the held-out fold, never a squared correlation.
# Dependencies: pROC. No other file of R/ is required.
# ----------------------------------------------------------------------------

CLASS_LEVELS <- c("resistant", "sensitive")
POSITIVE_CLASS <- "sensitive"

regression_metric_names <- function() c("r2", "rmse", "mae")

classification_metric_names <- function() {
  c("accuracy", "balanced_accuracy", "auc", "sensitivity", "specificity")
}

metric_names <- function(task) {
  task <- match.arg(task, c("regression", "classification"))
  if (task == "regression") regression_metric_names() else classification_metric_names()
}

# Check that two numeric vectors are usable together.
check_numeric_pair <- function(observed, predicted) {
  if (!is.numeric(observed) || !is.numeric(predicted)) {
    stop("observed and predicted must be numeric")
  }
  if (length(observed) != length(predicted)) {
    stop("observed and predicted must have the same length")
  }
  if (length(observed) == 0) stop("observed and predicted are empty")
  if (anyNA(observed) || anyNA(predicted)) {
    stop("observed and predicted must not contain missing values")
  }
  invisible(TRUE)
}

# Coerce a class vector to the factor with the fixed levels.
as_class_factor <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  bad <- setdiff(unique(x[!is.na(x)]), CLASS_LEVELS)
  if (length(bad) > 0) {
    stop("unknown class label(s): ", paste(bad, collapse = ", "),
         "; expected ", paste(CLASS_LEVELS, collapse = " or "))
  }
  factor(x, levels = CLASS_LEVELS)
}

# Coefficient of determination: 1 minus SSres over SStot, SStot around the
# mean of the observed (held-out) values. Can be negative; NA when the
# observed values are constant.
r_squared <- function(observed, predicted) {
  check_numeric_pair(observed, predicted)
  ss_res <- sum((observed - predicted)^2)
  ss_tot <- sum((observed - mean(observed))^2)
  if (ss_tot == 0) return(NA_real_)
  1 - ss_res / ss_tot
}

rmse <- function(observed, predicted) {
  check_numeric_pair(observed, predicted)
  sqrt(mean((observed - predicted)^2))
}

mae <- function(observed, predicted) {
  check_numeric_pair(observed, predicted)
  mean(abs(observed - predicted))
}

# ROC AUC with the sensitive class as positive and the score increasing with
# the probability of sensitive (section 7.4). Returns NA when one class is
# absent from observed or when the score contains missing values.
roc_auc <- function(observed, score) {
  observed <- as_class_factor(observed)
  if (length(observed) != length(score)) {
    stop("observed and score must have the same length")
  }
  if (!is.numeric(score)) stop("score must be numeric")
  if (anyNA(observed) || anyNA(score)) return(NA_real_)
  if (length(unique(observed)) < 2) return(NA_real_)
  r <- pROC::roc(response = observed, predictor = score,
                 levels = CLASS_LEVELS, direction = "<", quiet = TRUE)
  as.numeric(pROC::auc(r))
}

# Accuracy, balanced accuracy, ROC AUC, sensitivity (recall of sensitive) and
# specificity (recall of resistant) of one held-out fold. When a class is
# absent from observed, its recall, the balanced accuracy and the AUC are NA.
classification_metrics <- function(observed, predicted_class, score = NULL) {
  observed <- as_class_factor(observed)
  predicted_class <- as_class_factor(predicted_class)
  if (length(observed) != length(predicted_class)) {
    stop("observed and predicted_class must have the same length")
  }
  if (length(observed) == 0) stop("observed is empty")
  if (anyNA(observed) || anyNA(predicted_class)) {
    stop("observed and predicted_class must not contain missing values")
  }
  is_pos <- observed == POSITIVE_CLASS
  n_pos <- sum(is_pos)
  n_neg <- sum(!is_pos)
  accuracy <- mean(observed == predicted_class)
  sensitivity <- if (n_pos > 0) mean(predicted_class[is_pos] == POSITIVE_CLASS) else NA_real_
  specificity <- if (n_neg > 0) mean(predicted_class[!is_pos] != POSITIVE_CLASS) else NA_real_
  balanced <- (sensitivity + specificity) / 2
  auc <- if (is.null(score)) NA_real_ else roc_auc(observed, score)
  c(accuracy = accuracy, balanced_accuracy = balanced, auc = auc,
    sensitivity = sensitivity, specificity = specificity)
}

# Named numeric vector with the metrics of the task (section 9). For
# classification, predicted is the predicted class and score the ROC score.
evaluate <- function(task, observed, predicted, score = NULL) {
  task <- match.arg(task, c("regression", "classification"))
  if (task == "regression") {
    out <- c(r2 = r_squared(observed, predicted),
             rmse = rmse(observed, predicted),
             mae = mae(observed, predicted))
  } else {
    out <- classification_metrics(observed, predicted, score)
  }
  out[metric_names(task)]
}

# Two-sided t-based interval of a vector of fold values; NA values are
# dropped and counted out of n_folds.
t_interval <- function(values, conf_level = 0.95) {
  values <- values[!is.na(values)]
  n <- length(values)
  if (n == 0) {
    return(c(n_folds = 0, mean = NA_real_, sd = NA_real_,
             ci_low = NA_real_, ci_high = NA_real_))
  }
  m <- mean(values)
  if (n < 2) {
    return(c(n_folds = n, mean = m, sd = NA_real_,
             ci_low = NA_real_, ci_high = NA_real_))
  }
  s <- stats::sd(values)
  half <- stats::qt(1 - (1 - conf_level) / 2, df = n - 1) * s / sqrt(n)
  c(n_folds = n, mean = m, sd = s, ci_low = m - half, ci_high = m + half)
}

# Aggregate the fold metrics per drug, task, model and metric. The input is
# the metrics_by_fold table (columns drug, task, model, repeat, fold, n_test
# and one column per metric); any numeric column that is not an identifier is
# treated as a metric. Folds with NA are dropped from that metric and
# n_folds reports the folds used.
summarise_metrics <- function(metrics_by_fold, conf_level = 0.95) {
  required <- c("drug", "task", "model")
  missing_cols <- setdiff(required, names(metrics_by_fold))
  if (length(missing_cols) > 0) {
    stop("metrics_by_fold lacks column(s): ", paste(missing_cols, collapse = ", "))
  }
  id_cols <- c("drug", "task", "model", "repeat", "fold", "n_test", "n_train", "seconds")
  metric_cols <- setdiff(names(metrics_by_fold), id_cols)
  metric_cols <- metric_cols[vapply(metrics_by_fold[metric_cols], is.numeric, logical(1))]
  if (length(metric_cols) == 0) stop("metrics_by_fold has no metric column")
  known <- c(regression_metric_names(), classification_metric_names())
  metric_cols <- c(intersect(known, metric_cols), setdiff(metric_cols, known))

  groups <- unique(metrics_by_fold[required])
  rows <- vector("list", nrow(groups) * length(metric_cols))
  i <- 0
  for (g in seq_len(nrow(groups))) {
    sel <- metrics_by_fold$drug == groups$drug[g] &
      metrics_by_fold$task == groups$task[g] &
      metrics_by_fold$model == groups$model[g]
    for (metric in metric_cols) {
      stats <- t_interval(metrics_by_fold[[metric]][sel], conf_level)
      i <- i + 1
      rows[[i]] <- data.frame(
        drug = groups$drug[g], task = groups$task[g], model = groups$model[g],
        metric = metric, n_folds = as.integer(stats[["n_folds"]]),
        mean = stats[["mean"]], sd = stats[["sd"]],
        ci_low = stats[["ci_low"]], ci_high = stats[["ci_high"]],
        ci_method = "t", stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
