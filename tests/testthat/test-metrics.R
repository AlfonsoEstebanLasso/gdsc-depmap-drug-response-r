# Metric formulas, classification metrics on a hand-built confusion matrix,
# ROC orientation and the fold summary with t intervals.

test_that("r_squared is 1 minus SSres over SStot on held-out data", {
  observed <- c(0.2, 0.4, 0.6, 0.8, 1.0)
  expect_equal(r_squared(observed, observed), 1)
  expect_equal(r_squared(observed, rep(mean(observed), 5)), 0)
  expect_lt(r_squared(observed, rev(observed)), 0)
  predicted <- c(0.3, 0.3, 0.7, 0.7, 0.9)
  ss_res <- sum((observed - predicted)^2)
  ss_tot <- sum((observed - mean(observed))^2)
  expect_equal(r_squared(observed, predicted), 1 - ss_res / ss_tot)
  # The 2019 formula, cor(y, y - prediction)^2, differs from the coefficient
  # of determination for a non-trivial prediction.
  expect_false(isTRUE(all.equal(r_squared(observed, predicted),
                                stats::cor(observed, observed - predicted)^2)))
})

test_that("rmse and mae on known values", {
  observed <- c(1, 2, 3, 4)
  predicted <- c(1, 2, 5, 0)
  expect_equal(rmse(observed, predicted), sqrt((0 + 0 + 4 + 16) / 4))
  expect_equal(mae(observed, predicted), (0 + 0 + 2 + 4) / 4)
  expect_equal(rmse(observed, observed), 0)
  expect_equal(mae(observed, observed), 0)
})

test_that("classification metrics on a hand-built confusion matrix", {
  observed <- factor(c("resistant", "resistant", "resistant", "resistant", "sensitive", "sensitive"),
                     levels = c("resistant", "sensitive"))
  predicted <- factor(c("resistant", "resistant", "resistant", "sensitive", "sensitive", "resistant"),
                      levels = c("resistant", "sensitive"))
  score <- c(0.1, 0.2, 0.3, 0.6, 0.9, 0.4)
  cm <- classification_metrics(observed, predicted, score)
  expect_equal(unname(cm[["accuracy"]]), 4 / 6)
  expect_equal(unname(cm[["sensitivity"]]), 1 / 2)
  expect_equal(unname(cm[["specificity"]]), 3 / 4)
  expect_equal(unname(cm[["balanced_accuracy"]]), (1 / 2 + 3 / 4) / 2)
  expect_equal(unname(cm[["auc"]]), 7 / 8)
  expect_equal(unname(roc_auc(observed, score)), 7 / 8)
  ev <- evaluate("classification", observed, predicted, score)
  expect_setequal(names(ev), c("accuracy", "balanced_accuracy", "auc", "sensitivity", "specificity"))
  expect_equal(unname(ev[["accuracy"]]), 4 / 6)
  reg <- evaluate("regression", c(1, 2, 3), c(1, 2, 4))
  expect_setequal(names(reg), c("r2", "rmse", "mae"))
  expect_equal(unname(reg[["mae"]]), 1 / 3)
})

test_that("metrics are NA when a class is absent from the held-out fold", {
  observed <- factor(rep("resistant", 4), levels = c("resistant", "sensitive"))
  predicted <- factor(c("resistant", "resistant", "sensitive", "resistant"),
                      levels = c("resistant", "sensitive"))
  cm <- classification_metrics(observed, predicted, c(0.1, 0.2, 0.8, 0.3))
  expect_true(is.na(cm[["sensitivity"]]))
  expect_true(is.na(cm[["auc"]]))
  expect_equal(unname(cm[["specificity"]]), 3 / 4)
  expect_equal(unname(cm[["accuracy"]]), 3 / 4)
})

test_that("roc_auc is above 0.5 for a score increasing with sensitive and the SVM score is oriented", {
  set.seed(21)
  observed <- factor(sample(c("resistant", "sensitive"), 100, replace = TRUE, prob = c(0.75, 0.25)),
                     levels = c("resistant", "sensitive"))
  score <- ifelse(observed == "sensitive", 1, 0) + stats::rnorm(100, sd = 0.5)
  expect_gt(roc_auc(observed, score), 0.5)
  expect_lt(roc_auc(observed, -score), 0.5)

  m <- make_synthetic_matrix(n = 120, seed = 4)
  train <- 1:90
  test <- 91:120
  pp <- preprocess_fold(m$x, m$auc, train, test, "classification", n_genes = 200)
  foldid <- make_inner_folds(pp$y_train, 3, "classification", seed = 3)
  grid <- model_grid("svm_linear", "classification")
  fitted <- fit_model("svm_linear", pp$x_train, pp$y_train, "classification", foldid, grid,
                      seed = 11, threads = 1)
  pred <- predict_model(fitted, pp$x_test)
  expect_gt(roc_auc(pp$y_test, pred$score), 0.6)
  expect_gt(roc_auc(pp$y_train, predict_model(fitted, pp$x_train)$score), 0.8)
})

test_that("summarise_metrics reproduces a hand-computed t interval and handles NA", {
  metrics_by_fold <- data.frame(
    drug = "synthetic", task = "regression", model = "ridge",
    `repeat` = c(1, 1, 1), fold = 1:3, n_test = 20,
    r2 = c(0.1, 0.2, 0.3), rmse = c(0.5, NA, 0.7), mae = c(0.4, 0.4, 0.4),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  s <- summarise_metrics(metrics_by_fold)
  expect_equal(names(s), c("drug", "task", "model", "metric", "n_folds", "mean", "sd",
                           "ci_low", "ci_high", "ci_method"))
  expect_equal(nrow(s), 3)
  r2 <- s[s$metric == "r2", ]
  expect_equal(r2$n_folds, 3)
  expect_equal(r2$mean, 0.2)
  expect_equal(r2$sd, 0.1)
  half <- stats::qt(0.975, 2) * 0.1 / sqrt(3)
  expect_equal(r2$ci_low, 0.2 - half)
  expect_equal(r2$ci_high, 0.2 + half)
  expect_equal(r2$ci_method, "t")
  rm <- s[s$metric == "rmse", ]
  expect_equal(rm$n_folds, 2)
  expect_equal(rm$mean, 0.6)
  expect_equal(rm$sd, stats::sd(c(0.5, 0.7)))
  ma <- s[s$metric == "mae", ]
  expect_equal(ma$mean, 0.4)
  expect_equal(ma$sd, 0)
  expect_equal(ma$ci_low, 0.4)
  expect_equal(ma$ci_high, 0.4)
})
