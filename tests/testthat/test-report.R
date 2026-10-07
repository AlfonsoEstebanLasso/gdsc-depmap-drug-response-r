# make_report() on the files written by run_drug() for the synthetic drug:
# figures, cross-drug summary tables and Markdown tables.

outdir <- tempfile("v1report")
dir.create(outdir)
run_synthetic_example(outdir, seed = 1)

test_that("make_report creates the figures and the summary files without error", {
  expect_error(make_report(outdir, drugs = "synthetic",
                           tasks = c("regression", "classification")), NA)
  reg_dir <- file.path(outdir, "synthetic", "regression")
  cls_dir <- file.path(outdir, "synthetic", "classification")
  expect_true(file.exists(file.path(reg_dir, "fig_metrics.png")))
  expect_true(file.exists(file.path(reg_dir, "fig_observed_vs_predicted.png")))
  expect_false(file.exists(file.path(reg_dir, "fig_roc.png")))
  expect_true(file.exists(file.path(cls_dir, "fig_metrics.png")))
  expect_true(file.exists(file.path(cls_dir, "fig_roc.png")))
  expect_false(file.exists(file.path(cls_dir, "fig_observed_vs_predicted.png")))
  for (f in c("summary.csv", "summary_regression.md", "summary_classification.md",
              "fig_summary_regression.png", "fig_summary_classification.png")) {
    expect_true(file.exists(file.path(outdir, f)), info = f)
  }
  expect_gt(file.size(file.path(reg_dir, "fig_metrics.png")), 1000)
})

test_that("summary.csv aggregates both tasks and the Markdown tables list every model", {
  s <- read_csv_keep_names(file.path(outdir, "summary.csv"))
  expect_equal(names(s), c("drug", "task", "model", "metric", "n_folds", "mean", "sd",
                           "ci_low", "ci_high", "ci_method"))
  expect_setequal(unique(s$task), c("regression", "classification"))
  expect_equal(nrow(s), length(model_names()) * (3 + 5))
  reg <- read_csv_keep_names(file.path(outdir, "synthetic", "regression", "summary.csv"))
  expect_equal(s$mean[s$task == "regression"], reg$mean)

  md <- readLines(file.path(outdir, "summary_regression.md"))
  expect_true(any(grepl("^## R squared", md)))
  expect_true(any(grepl("^## RMSE", md)))
  expect_true(any(grepl("^\\| Model \\| Synthetic \\|", md)))
  for (model in model_names()) {
    expect_true(any(grepl(paste0("^\\| ", model, " \\| "), md)), info = model)
  }
  expect_false(any(grepl("\u2013|\u2014", md)))
  md_cls <- readLines(file.path(outdir, "summary_classification.md"))
  expect_true(any(grepl("^## ROC AUC", md_cls)))
  expect_true(any(grepl("^## Balanced accuracy", md_cls)))
})

test_that("the plot functions return ggplot objects", {
  mbf <- read_csv_keep_names(file.path(outdir, "synthetic", "regression", "metrics_by_fold.csv"))
  expect_s3_class(plot_metric_boxes(mbf, "regression"), "ggplot")
  pred_reg <- read_csv_keep_names(file.path(outdir, "synthetic", "regression", "predictions.csv"))
  expect_s3_class(plot_observed_vs_predicted(pred_reg), "ggplot")
  pred_cls <- read_csv_keep_names(file.path(outdir, "synthetic", "classification", "predictions.csv"))
  expect_s3_class(plot_roc(pred_cls), "ggplot")
  s <- read_csv_keep_names(file.path(outdir, "summary.csv"))
  expect_s3_class(plot_summary(s, "classification"), "ggplot")
})

test_that("make_report discovers the drugs when none is given", {
  expect_error(make_report(outdir), NA)
  expect_equal(discover_drugs(outdir), "synthetic")
})
