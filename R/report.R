# R/report.R
#
# Purpose: figures and tables of the v1.0 outputs (docs/v1_spec.md, section
# 12). Everything is built from the CSV files written by run_drug(); no raw
# data file is read here.
#
# Inputs: <outdir>/<drug>/<task>/metrics_by_fold.csv, summary.csv and
# predictions.csv (the last one is local only and is used for the observed
# versus predicted and ROC figures).
#
# Outputs: per drug and task, fig_metrics.png (box plots of the fold metrics
# by model), fig_observed_vs_predicted.png (regression) or fig_roc.png
# (classification); across drugs, outputs/summary.csv, summary_<task>.md and
# fig_summary_<task>.png (mean and 95 percent interval per drug and model).
#
# Figures: ggplot2, PNG, 8 by 5 inches at 150 dpi, viridis palette (colour
# blind safe), sentence-case titles with the drug and the task.

THESIS_DRUGS <- c("erlotinib", "rapamycin", "sunitinib", "paclitaxel")
ALL_TASKS <- c("regression", "classification")

metric_labels <- function(metric) {
  labels <- c(
    r2 = "R squared", rmse = "RMSE", mae = "MAE",
    accuracy = "Accuracy", balanced_accuracy = "Balanced accuracy",
    auc = "ROC AUC", sensitivity = "Sensitivity", specificity = "Specificity"
  )
  out <- unname(labels[metric])
  out[is.na(out)] <- metric[is.na(out)]
  out
}

metric_order <- function(metrics) {
  known <- c("r2", "rmse", "mae", "accuracy", "balanced_accuracy", "auc",
             "sensitivity", "specificity")
  c(intersect(known, metrics), setdiff(metrics, known))
}

model_order <- function(models) {
  known <- if (exists("model_names", mode = "function")) model_names() else character(0)
  c(intersect(known, models), setdiff(models, known))
}

drug_order <- function(drugs) {
  drugs <- unique(tolower(drugs))
  c(intersect(THESIS_DRUGS, drugs), sort(setdiff(drugs, THESIS_DRUGS)))
}

drug_label <- function(drug) {
  drug <- as.character(drug)
  paste0(toupper(substring(drug, 1, 1)), substring(drug, 2))
}

read_report_csv <- function(path) {
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

discover_drugs <- function(outdir, tasks = ALL_TASKS) {
  dirs <- list.dirs(outdir, full.names = FALSE, recursive = FALSE)
  keep <- vapply(dirs, function(d) {
    any(file.exists(file.path(outdir, d, tasks, "summary.csv")))
  }, logical(1))
  drug_order(dirs[keep])
}

metrics_long <- function(metrics_by_fold) {
  id_cols <- c("drug", "task", "model", "repeat", "fold", "n_test")
  metrics <- metric_order(setdiff(names(metrics_by_fold), id_cols))
  rows <- lapply(metrics, function(m) {
    data.frame(
      model = metrics_by_fold$model,
      metric = m,
      value = as.numeric(metrics_by_fold[[m]]),
      stringsAsFactors = FALSE
    )
  })
  long <- do.call(rbind, rows)
  long$model <- factor(long$model, levels = model_order(unique(long$model)))
  long$metric <- factor(metric_labels(long$metric), levels = metric_labels(metrics))
  long
}

plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title.position = "plot"
    )
}

plot_metric_boxes <- function(metrics_by_fold, task) {
  long <- metrics_long(metrics_by_fold)
  long <- long[!is.na(long$value), , drop = FALSE]
  drugs <- unique(metrics_by_fold$drug)
  title <- sprintf("%s, %s: metrics per held-out fold",
                   paste(drug_label(drugs), collapse = ", "), task)
  ggplot2::ggplot(long, ggplot2::aes(x = .data$model, y = .data$value, fill = .data$model)) +
    ggplot2::geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.6) +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = 0.9, alpha = 0.6) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::scale_fill_viridis_d(end = 0.9) +
    ggplot2::labs(title = title, x = NULL, y = "Value on the held-out fold") +
    plot_theme() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

plot_observed_vs_predicted <- function(predictions, drug = NULL) {
  pred <- predictions[predictions[["repeat"]] == min(predictions[["repeat"]]), , drop = FALSE]
  pred$observed <- as.numeric(pred$observed)
  pred$predicted <- as.numeric(pred$predicted)
  pred$model <- factor(pred$model, levels = model_order(unique(pred$model)))
  lims <- range(c(pred$observed, pred$predicted), na.rm = TRUE)
  title <- sprintf("%sregression: observed versus predicted AUC (out of fold, repeat %d)",
                   if (is.null(drug)) "" else paste0(drug_label(drug), ", "),
                   min(predictions[["repeat"]]))
  ggplot2::ggplot(pred, ggplot2::aes(x = .data$observed, y = .data$predicted, colour = .data$model)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(size = 1.2, alpha = 0.7) +
    ggplot2::facet_wrap(~ model) +
    ggplot2::coord_equal(xlim = lims, ylim = lims) +
    ggplot2::scale_colour_viridis_d(end = 0.9) +
    ggplot2::labs(title = title, x = "Observed AUC", y = "Predicted AUC") +
    plot_theme() +
    ggplot2::theme(legend.position = "none")
}

roc_curve_points <- function(observed, score) {
  observed <- factor(as.character(observed), levels = c("resistant", "sensitive"))
  ok <- !is.na(observed) & !is.na(score)
  roc <- pROC::roc(observed[ok], as.numeric(score[ok]),
                   levels = c("resistant", "sensitive"), direction = "<", quiet = TRUE)
  pts <- data.frame(fpr = 1 - roc$specificities, tpr = roc$sensitivities)
  pts <- pts[order(pts$fpr, pts$tpr), , drop = FALSE]
  list(points = pts, auc = as.numeric(pROC::auc(roc)))
}

plot_roc <- function(predictions, drug = NULL) {
  pred <- predictions[predictions[["repeat"]] == min(predictions[["repeat"]]), , drop = FALSE]
  models <- model_order(unique(pred$model))
  curves <- list()
  labels <- character(0)
  for (m in models) {
    sub <- pred[pred$model == m, , drop = FALSE]
    rc <- roc_curve_points(sub$observed, sub$score)
    label <- sprintf("%s (AUC %.3f)", m, rc$auc)
    labels <- c(labels, label)
    curves[[m]] <- cbind(rc$points, model = label, stringsAsFactors = FALSE)
  }
  curves <- do.call(rbind, curves)
  curves$model <- factor(curves$model, levels = labels)
  title <- sprintf("%sclassification: ROC curves (pooled out of fold scores, repeat %d)",
                   if (is.null(drug)) "" else paste0(drug_label(drug), ", "),
                   min(predictions[["repeat"]]))
  ggplot2::ggplot(curves, ggplot2::aes(x = .data$fpr, y = .data$tpr, colour = .data$model)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_path(linewidth = 0.8) +
    ggplot2::coord_equal() +
    ggplot2::scale_colour_viridis_d(end = 0.9, name = "Model") +
    ggplot2::labs(title = title, x = "1 - specificity", y = "Sensitivity") +
    plot_theme()
}

plot_summary <- function(summary, task) {
  summary <- summary[summary$task == task, , drop = FALSE]
  summary$drug <- factor(drug_label(summary$drug), levels = drug_label(drug_order(summary$drug)))
  summary$model <- factor(summary$model, levels = model_order(unique(summary$model)))
  metrics <- metric_order(unique(summary$metric))
  summary$metric <- factor(metric_labels(summary$metric), levels = metric_labels(metrics))
  ggplot2::ggplot(summary, ggplot2::aes(x = .data$drug, y = .data$mean, colour = .data$model)) +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = .data$ci_low, ymax = .data$ci_high),
      position = ggplot2::position_dodge(width = 0.6), size = 0.3
    ) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::scale_colour_viridis_d(end = 0.9, name = "Model") +
    ggplot2::labs(
      title = sprintf("%s: mean and 95 percent interval over the folds, per drug and model",
                      drug_label(task)),
      x = NULL, y = "Mean over the folds"
    ) +
    plot_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

save_figure <- function(plot, path) {
  ggplot2::ggsave(path, plot, width = 8, height = 5, dpi = 150, bg = "white")
  invisible(path)
}

format_cell <- function(mean, low, high) {
  ifelse(is.na(mean), "NA",
         sprintf("%.3f (%.3f to %.3f)", mean, low, high))
}

markdown_table <- function(summary, metric, drugs) {
  sub <- summary[summary$metric == metric, , drop = FALSE]
  models <- model_order(unique(sub$model))
  header <- paste0("| Model | ", paste(drug_label(drugs), collapse = " | "), " |")
  rule <- paste0("|---|", paste(rep("---", length(drugs)), collapse = "|"), "|")
  rows <- vapply(models, function(m) {
    cells <- vapply(drugs, function(d) {
      row <- sub[sub$model == m & sub$drug == d, , drop = FALSE]
      if (nrow(row) == 0L) return("NA")
      format_cell(row$mean[1], row$ci_low[1], row$ci_high[1])
    }, character(1))
    paste0("| ", m, " | ", paste(cells, collapse = " | "), " |")
  }, character(1))
  c(header, rule, rows)
}

write_summary_tables <- function(outdir) {
  drugs <- discover_drugs(outdir)
  if (length(drugs) == 0L) {
    stop("no <drug>/<task>/summary.csv found under ", outdir)
  }
  parts <- list()
  for (d in drugs) {
    for (task in ALL_TASKS) {
      path <- file.path(outdir, d, task, "summary.csv")
      if (file.exists(path)) parts[[length(parts) + 1L]] <- read_report_csv(path)
    }
  }
  summary <- do.call(rbind, parts)
  summary$drug <- tolower(summary$drug)
  utils::write.csv(summary, file.path(outdir, "summary.csv"), row.names = FALSE)

  written <- file.path(outdir, "summary.csv")
  for (task in ALL_TASKS) {
    sub <- summary[summary$task == task, , drop = FALSE]
    if (nrow(sub) == 0L) next
    task_drugs <- drug_order(unique(sub$drug))
    n_folds <- sort(unique(sub$n_folds))
    lines <- c(
      sprintf("# %s: mean and 95 percent interval over the folds", drug_label(task)),
      "",
      sprintf("Cells show the mean over the %s held-out folds and a two-sided 95 percent %s interval (mean plus or minus qt(0.975, n - 1) times sd over sqrt(n)). Folds of different repeats share cell lines, so the interval is descriptive, not an exact inferential statement.",
              paste(n_folds, collapse = " or "),
              if (all(sub$ci_method == "t")) "t-based" else paste(unique(sub$ci_method), collapse = " or ")),
      ""
    )
    for (metric in metric_order(unique(sub$metric))) {
      lines <- c(lines, sprintf("## %s", metric_labels(metric)), "",
                 markdown_table(sub, metric, task_drugs), "")
    }
    md_path <- file.path(outdir, sprintf("summary_%s.md", task))
    writeLines(lines, md_path)
    written <- c(written, md_path)
  }
  invisible(written)
}

make_report <- function(outdir, drugs = NULL, tasks = ALL_TASKS) {
  tasks <- match.arg(tasks, ALL_TASKS, several.ok = TRUE)
  if (is.null(drugs)) drugs <- discover_drugs(outdir, tasks)
  drugs <- tolower(drugs)
  written <- character(0)

  for (d in drugs) {
    for (task in tasks) {
      task_dir <- file.path(outdir, d, task)
      metrics_path <- file.path(task_dir, "metrics_by_fold.csv")
      if (!file.exists(metrics_path)) next
      metrics_by_fold <- read_report_csv(metrics_path)
      written <- c(written, save_figure(plot_metric_boxes(metrics_by_fold, task),
                                        file.path(task_dir, "fig_metrics.png")))
      pred_path <- file.path(task_dir, "predictions.csv")
      if (file.exists(pred_path)) {
        predictions <- read_report_csv(pred_path)
        if (task == "regression") {
          written <- c(written, save_figure(plot_observed_vs_predicted(predictions, drug = d),
                                            file.path(task_dir, "fig_observed_vs_predicted.png")))
        } else {
          written <- c(written, save_figure(plot_roc(predictions, drug = d),
                                            file.path(task_dir, "fig_roc.png")))
        }
      }
    }
  }

  written <- c(written, write_summary_tables(outdir))
  summary <- read_report_csv(file.path(outdir, "summary.csv"))
  for (task in tasks) {
    if (!any(summary$task == task)) next
    written <- c(written, save_figure(plot_summary(summary, task),
                                      file.path(outdir, sprintf("fig_summary_%s.png", task))))
  }
  invisible(written)
}
