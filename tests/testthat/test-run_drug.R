# run_drug() end to end on the synthetic matrix with a tiny configuration:
# every file of section 12, consistent row counts, reproducibility and the
# leakage check.

m <- make_synthetic_matrix()
n <- length(m$auc)
n_models <- length(model_names())

test_that("default_config declares the spec defaults and the --fast column", {
  full <- default_config()
  expect_equal(c(full$k, full$repeats, full$k_inner, full$rf_budget), c(5L, 3L, 5L, 12L))
  expect_equal(full$n_genes, 200L)
  expect_equal(full$num_trees, 500L)
  expect_equal(full$seed, 2019L)
  expect_equal(full$models, model_names())
  expect_true(full$leakage_check)
  expect_false(full$fast)
  fast <- default_config(fast = TRUE)
  expect_equal(c(fast$k, fast$repeats, fast$k_inner, fast$rf_budget), c(3L, 1L, 3L, 6L))
  expect_false(fast$leakage_check)
  expect_true(fast$fast)
})

outdir <- tempfile("v1run")
dir.create(outdir)
results <- run_synthetic_example(outdir, seed = 1)

test_that("run_drug writes every per task file of section 12", {
  for (task in c("regression", "classification")) {
    task_dir <- file.path(outdir, "synthetic", task)
    for (f in task_files()) {
      expect_true(file.exists(file.path(task_dir, f)), info = paste(task, f))
    }
  }
})

test_that("metrics_by_fold and summary have the expected shape", {
  expected_metrics <- list(
    regression = c("r2", "rmse", "mae"),
    classification = c("accuracy", "balanced_accuracy", "auc", "sensitivity", "specificity")
  )
  for (task in c("regression", "classification")) {
    task_dir <- file.path(outdir, "synthetic", task)
    mbf <- read_csv_keep_names(file.path(task_dir, "metrics_by_fold.csv"))
    expect_equal(names(mbf)[1:6], c("drug", "task", "model", "repeat", "fold", "n_test"))
    expect_setequal(setdiff(names(mbf), names(mbf)[1:6]), expected_metrics[[task]])
    expect_equal(nrow(mbf), n_models * 2)
    expect_setequal(unique(mbf$model), model_names())
    expect_equal(sum(mbf$n_test[mbf$model == "ridge"]), n)
    expect_true(all(mbf$task == task))
    expect_true(all(mbf$drug == "synthetic"))

    s <- read_csv_keep_names(file.path(task_dir, "summary.csv"))
    expect_equal(names(s), c("drug", "task", "model", "metric", "n_folds", "mean", "sd",
                             "ci_low", "ci_high", "ci_method"))
    expect_equal(nrow(s), n_models * length(expected_metrics[[task]]))
    expect_equal(nrow(unique(s[, c("model", "metric")])), nrow(s))
    expect_true(all(s$n_folds <= 2))
    expect_equal(results[[task]]$summary$mean, s$mean)
  }
})

test_that("predictions, tuning, selected genes and timing are consistent", {
  for (task in c("regression", "classification")) {
    task_dir <- file.path(outdir, "synthetic", task)
    pred <- read_csv_keep_names(file.path(task_dir, "predictions.csv"))
    expect_equal(names(pred), c("cell_line", "repeat", "fold", "model", "observed", "predicted", "score"))
    expect_equal(nrow(pred), n * n_models)
    expect_setequal(pred$cell_line[pred$model == "rf"], m$cell_lines)
    if (task == "regression") {
      expect_true(is.numeric(pred$observed))
      expect_true(all(is.na(pred$score)))
    } else {
      expect_setequal(unique(pred$observed), c("resistant", "sensitive"))
      expect_true(all(pred$predicted %in% c("resistant", "sensitive")))
      expect_true(all(is.finite(pred$score)))
    }

    tuning <- read_csv_keep_names(file.path(task_dir, "tuning.csv"))
    expect_equal(names(tuning), c("drug", "task", "model", "repeat", "fold", "parameter", "value", "inner_score"))
    expect_setequal(unique(tuning$model), model_names())
    rf_params <- tuning$parameter[tuning$model == "rf" & tuning$fold == 1]
    expect_true(all(c("mtry", "min.node.size", "sample.fraction", "max.depth", "num.trees") %in% rf_params))
    expect_true(all(tuning$parameter[tuning$model == "enet"] != "" ))

    candidates <- read_csv_keep_names(file.path(task_dir, "tuning_candidates.csv"))
    expect_equal(names(candidates), c("drug", "task", "model", "repeat", "fold", "candidate",
                                      "parameter", "value", "inner_score"))
    rf_candidates <- candidates[candidates$model == "rf" & candidates$fold == 1, ]
    expect_equal(length(unique(rf_candidates$candidate)), 2)
    enet_candidates <- candidates[candidates$model == "enet" & candidates$fold == 1, ]
    expect_equal(length(unique(enet_candidates$candidate)), 2)

    genes <- read_csv_keep_names(file.path(task_dir, "selected_genes.csv"))
    expect_equal(names(genes), c("gene", "times_selected", "n_splits"))
    expect_true(all(genes$n_splits == 2))
    expect_true(all(genes$times_selected >= 1 & genes$times_selected <= 2))
    expect_equal(sum(genes$times_selected), 2 * 200)
    expect_true(all(colnames(m$x)[1:10] %in% genes$gene))

    timing <- read_csv_keep_names(file.path(task_dir, "timing.csv"))
    expect_equal(names(timing), c("model", "repeat", "fold", "seconds"))
    expect_equal(nrow(timing), n_models * 2)
    expect_true(all(timing$seconds >= 0))
  }
})

test_that("a second run with the same seed gives identical metrics_by_fold.csv", {
  outdir2 <- tempfile("v1run2")
  dir.create(outdir2)
  run_synthetic_example(outdir2, seed = 1)
  for (task in c("regression", "classification")) {
    a <- read_csv_keep_names(file.path(outdir, "synthetic", task, "metrics_by_fold.csv"))
    b <- read_csv_keep_names(file.path(outdir2, "synthetic", task, "metrics_by_fold.csv"))
    expect_identical(a, b)
    ta <- read_csv_keep_names(file.path(outdir, "synthetic", task, "tuning.csv"))
    tb <- read_csv_keep_names(file.path(outdir2, "synthetic", task, "tuning.csv"))
    expect_identical(ta, tb)
  }
})

test_that("the leakage check writes the comparison of the glmnet family", {
  outdir3 <- tempfile("v1leak")
  dir.create(outdir3)
  config <- tiny_config("classification", leakage_check = TRUE)
  res <- quiet_small_classes(run_drug(m$drug, "classification", m, config, outdir3))
  path <- file.path(outdir3, "synthetic", "classification", "leakage_check.csv")
  expect_true(file.exists(path))
  leak <- read_csv_keep_names(path)
  expect_equal(names(leak), c("model", "metric", "mean_fold_wise", "mean_global", "difference", "n_folds"))
  expect_setequal(unique(leak$model), c("ridge", "lasso", "enet"))
  expect_setequal(unique(leak$metric), c("accuracy", "balanced_accuracy", "auc", "sensitivity", "specificity"))
  expect_equal(leak$difference, leak$mean_global - leak$mean_fold_wise)
  expect_equal(res$leakage_check$mean_global, leak$mean_global)
  mbf <- res$metrics_by_fold
  ridge_auc <- mean(mbf$auc[mbf$model == "ridge"], na.rm = TRUE)
  expect_equal(leak$mean_fold_wise[leak$model == "ridge" & leak$metric == "auc"], ridge_auc)
})

test_that("write_run_log records the configuration, the packages and sessionInfo", {
  path <- file.path(tempfile("log"), "run_log.txt")
  config <- tiny_config("regression")
  write_run_log(path, config, c("--fast", "--drug", "synthetic"),
                list(`synthetic regression` = 1.5, total = 2), packages = c("glmnet", "nosuchpackage"))
  expect_true(file.exists(path))
  text <- readLines(path)
  expect_true(any(grepl("--fast --drug synthetic", text, fixed = TRUE)))
  expect_true(any(grepl("seed = 2019", text, fixed = TRUE)))
  expect_true(any(grepl("^  glmnet [0-9.]+", text)))
  expect_true(any(grepl("nosuchpackage not installed", text, fixed = TRUE)))
  expect_true(any(grepl("sessionInfo()", text, fixed = TRUE)))
  expect_true(any(grepl("R version", text, fixed = TRUE)))
})
