# R/run_drug.R
#
# Purpose: run the nested cross-validation of one drug and one task
# (regression on the AUC or classification of sensitive lines) and write the
# per drug and task outputs of docs/v1_spec.md, section 12.
#
# Inputs: the per-drug matrix of build_matrix() or load_drug_matrix()
# (list with drug, cell_lines, x, auc), a configuration from default_config()
# and an output directory.
#
# Outputs, under <outdir>/<drug>/<task>/: metrics_by_fold.csv, summary.csv,
# tuning.csv, tuning_candidates.csv, selected_genes.csv, timing.csv,
# predictions.csv (local only, not committed) and, with the leakage check,
# leakage_check.csv. Fitted objects are written only with save_models = TRUE.
#
# Transformations: outer folds from make_folds() (seed S + t), fold-wise
# preprocessing with preprocess_fold() (gene selection, class threshold and
# standardisation on the training rows only), inner folds from
# make_inner_folds() shared by every model family, fit_model() and
# predict_model() per family, evaluate() on the held-out rows and
# summarise_metrics() over the folds. Seeds follow section 4.3 of the spec:
# the base seed of a split is S + 1000 * t + j and fit_model() derives the
# family-specific offsets from it.
#
# Depends on the functions of R/build_matrix.R, R/cv.R, R/models.R and
# R/metrics.R, sourced beforehand by run_all.R or by the tests.

default_config <- function(fast = FALSE) {
  list(
    k = if (fast) 3L else 5L,
    repeats = if (fast) 1L else 3L,
    k_inner = if (fast) 3L else 5L,
    n_genes = 200L,
    rf_budget = if (fast) 6L else 12L,
    num_trees = 500L,
    seed = 2019L,
    threads = max(1L, parallel::detectCores() - 1L),
    models = model_names(),
    leakage_check = !fast,
    fast = fast,
    save_models = FALSE,
    grids = NULL
  )
}

task_index <- function(task) {
  c(regression = 1L, classification = 2L)[[task]]
}

split_seed <- function(config, task, j) {
  as.integer(config$seed + 1000L * task_index(task) + j)
}

task_dir_of <- function(outdir, drug, task) {
  file.path(outdir, tolower(drug), task)
}

read_output_csv <- function(path) {
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_output_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
}

grid_for_model <- function(model, task, config, seed) {
  if (!is.null(config$grids) && !is.null(config$grids[[model]])) {
    return(config$grids[[model]])
  }
  model_grid(model, task, rf_budget = config$rf_budget, seed = seed,
             num_trees = config$num_trees)
}

best_inner_score <- function(fitted, task) {
  if (!is.null(fitted$best$inner_score)) {
    return(as.numeric(fitted$best$inner_score))
  }
  scores <- fitted$tuning$inner_score
  if (is.null(scores) || all(is.na(scores))) {
    return(NA_real_)
  }
  if (task == "regression") min(scores, na.rm = TRUE) else max(scores, na.rm = TRUE)
}

value_as_character <- function(value) {
  value <- unlist(value, use.names = FALSE)
  if (is.numeric(value)) {
    return(format(value, digits = 10, trim = TRUE))
  }
  as.character(value)
}

selected_tuning_rows <- function(drug, task, model, rep, fold, fitted) {
  best <- fitted$best
  if (is.null(best) || length(best) == 0L) {
    return(NULL)
  }
  keep <- setdiff(names(best), "inner_score")
  data.frame(
    drug = drug, task = task, model = model, `repeat` = rep, fold = fold,
    parameter = keep,
    value = vapply(best[keep], value_as_character, character(1)),
    inner_score = best_inner_score(fitted, task),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

candidate_tuning_rows <- function(drug, task, model, rep, fold, fitted) {
  tuning <- fitted$tuning
  if (is.null(tuning) || nrow(tuning) == 0L) {
    return(NULL)
  }
  params <- setdiff(names(tuning), "inner_score")
  rows <- lapply(seq_len(nrow(tuning)), function(i) {
    data.frame(
      drug = drug, task = task, model = model, `repeat` = rep, fold = fold,
      candidate = i, parameter = params,
      value = vapply(params, function(p) value_as_character(tuning[[p]][i]), character(1)),
      inner_score = as.numeric(tuning$inner_score[i]),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

extract_predictions <- function(task, pred) {
  if (task == "regression") {
    list(predicted = as.numeric(pred), score = rep(NA_real_, length(pred)))
  } else {
    list(predicted = pred$class, score = as.numeric(pred$score))
  }
}

count_selected_genes <- function(genes_list) {
  counts <- table(unlist(genes_list, use.names = FALSE))
  out <- data.frame(
    gene = names(counts),
    times_selected = as.integer(counts),
    n_splits = length(genes_list),
    stringsAsFactors = FALSE
  )
  out[order(-out$times_selected, out$gene), , drop = FALSE]
}

run_drug <- function(drug, task, matrix, config, outdir) {
  task <- match.arg(task, c("regression", "classification"))
  drug <- tolower(drug)
  task_dir <- task_dir_of(outdir, drug, task)
  dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)

  x <- matrix$x
  cell_lines <- if (!is.null(matrix$cell_lines)) matrix$cell_lines else rownames(x)
  auc <- stats::setNames(as.numeric(matrix$auc), cell_lines)
  t_idx <- task_index(task)

  folds <- make_folds(auc, k = config$k, repeats = config$repeats, task = task,
                      seed = config$seed + t_idx)

  metrics_rows <- list()
  prediction_rows <- list()
  tuning_rows <- list()
  candidate_rows <- list()
  timing_rows <- list()
  genes_list <- list()

  j <- 0L
  for (rep in seq_len(config$repeats)) {
    for (fold in seq_len(config$k)) {
      j <- j + 1L
      idx <- split_rows(folds, rep, fold)
      pp <- preprocess_fold(x, auc, idx$train, idx$test, task, n_genes = config$n_genes)
      genes_list[[j]] <- pp$genes
      seed_j <- split_seed(config, task, j)
      foldid <- make_inner_folds(pp$y_train, config$k_inner, task, seed = seed_j)

      for (model in config$models) {
        grid <- grid_for_model(model, task, config, seed = seed_j + 100L)
        t0 <- proc.time()[["elapsed"]]
        fitted <- fit_model(model, pp$x_train, pp$y_train, task, foldid, grid,
                            seed = seed_j, threads = config$threads)
        pred <- predict_model(fitted, pp$x_test)
        seconds <- proc.time()[["elapsed"]] - t0

        out <- extract_predictions(task, pred)
        metrics <- evaluate(task, pp$y_test, out$predicted, out$score)

        metrics_rows[[length(metrics_rows) + 1L]] <- cbind(
          data.frame(drug = drug, task = task, model = model, `repeat` = rep,
                     fold = fold, n_test = length(idx$test),
                     stringsAsFactors = FALSE, check.names = FALSE),
          as.data.frame(as.list(metrics), check.names = FALSE)
        )
        prediction_rows[[length(prediction_rows) + 1L]] <- data.frame(
          cell_line = cell_lines[idx$test], `repeat` = rep, fold = fold, model = model,
          observed = if (task == "regression") as.numeric(pp$y_test) else as.character(pp$y_test),
          predicted = if (task == "regression") as.numeric(out$predicted) else as.character(out$predicted),
          score = out$score,
          stringsAsFactors = FALSE, check.names = FALSE
        )
        tuning_rows[[length(tuning_rows) + 1L]] <-
          selected_tuning_rows(drug, task, model, rep, fold, fitted)
        candidate_rows[[length(candidate_rows) + 1L]] <-
          candidate_tuning_rows(drug, task, model, rep, fold, fitted)
        timing_rows[[length(timing_rows) + 1L]] <- data.frame(
          model = model, `repeat` = rep, fold = fold, seconds = round(seconds, 3),
          stringsAsFactors = FALSE, check.names = FALSE
        )
        if (isTRUE(config$save_models)) {
          model_dir <- file.path(task_dir, "models")
          dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
          saveRDS(fitted, file.path(model_dir, sprintf("%s_r%d_f%d.rds", model, rep, fold)))
        }
        message(sprintf("[%s] %s %s repeat %d fold %d %s: %.1f s",
                        format(Sys.time(), "%H:%M:%S"), drug, task, rep, fold, model, seconds))
      }
    }
  }

  metrics_by_fold <- do.call(rbind, metrics_rows)
  predictions <- do.call(rbind, prediction_rows)
  tuning <- do.call(rbind, Filter(Negate(is.null), tuning_rows))
  tuning_candidates <- do.call(rbind, Filter(Negate(is.null), candidate_rows))
  timing <- do.call(rbind, timing_rows)
  selected_genes <- count_selected_genes(genes_list)
  summary <- summarise_metrics(metrics_by_fold)

  write_output_csv(metrics_by_fold, file.path(task_dir, "metrics_by_fold.csv"))
  write_output_csv(summary, file.path(task_dir, "summary.csv"))
  write_output_csv(tuning, file.path(task_dir, "tuning.csv"))
  write_output_csv(tuning_candidates, file.path(task_dir, "tuning_candidates.csv"))
  write_output_csv(selected_genes, file.path(task_dir, "selected_genes.csv"))
  write_output_csv(timing, file.path(task_dir, "timing.csv"))
  write_output_csv(predictions, file.path(task_dir, "predictions.csv"))

  leakage_check <- NULL
  if (isTRUE(config$leakage_check)) {
    leakage_check <- run_leakage_check(drug, task, matrix, folds, config, outdir)
  }

  invisible(list(
    metrics_by_fold = metrics_by_fold,
    summary = summary,
    tuning = tuning,
    tuning_candidates = tuning_candidates,
    predictions = predictions,
    timing = timing,
    selected_genes = selected_genes,
    leakage_check = leakage_check,
    folds = folds
  ))
}

# Leakage check of section 8: the glmnet family refitted on the same outer
# folds with the 2019 procedure (genes selected on all cell lines and, for
# classification, the class threshold computed on all cell lines). The
# difference column is mean_global minus mean_fold_wise.
run_leakage_check <- function(drug, task, matrix, folds, config, outdir) {
  task <- match.arg(task, c("regression", "classification"))
  drug <- tolower(drug)
  task_dir <- task_dir_of(outdir, drug, task)
  dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
  models <- intersect(c("ridge", "lasso", "enet"), config$models)
  if (length(models) == 0L) {
    return(NULL)
  }

  x <- matrix$x
  cell_lines <- if (!is.null(matrix$cell_lines)) matrix$cell_lines else rownames(x)
  auc <- stats::setNames(as.numeric(matrix$auc), cell_lines)
  global_threshold <- if (task == "classification") class_threshold(auc) else NULL

  global_rows <- list()
  j <- 0L
  for (rep in seq_len(config$repeats)) {
    for (fold in seq_len(config$k)) {
      j <- j + 1L
      idx <- split_rows(folds, rep, fold)
      pp <- preprocess_fold(x, auc, idx$train, idx$test, task, n_genes = config$n_genes,
                            selection = "global", threshold = global_threshold)
      seed_j <- split_seed(config, task, j)
      foldid <- make_inner_folds(pp$y_train, config$k_inner, task, seed = seed_j)
      for (model in models) {
        grid <- grid_for_model(model, task, config, seed = seed_j + 100L)
        fitted <- fit_model(model, pp$x_train, pp$y_train, task, foldid, grid,
                            seed = seed_j, threads = config$threads)
        out <- extract_predictions(task, predict_model(fitted, pp$x_test))
        metrics <- evaluate(task, pp$y_test, out$predicted, out$score)
        global_rows[[length(global_rows) + 1L]] <- cbind(
          data.frame(model = model, `repeat` = rep, fold = fold,
                     stringsAsFactors = FALSE, check.names = FALSE),
          as.data.frame(as.list(metrics), check.names = FALSE)
        )
      }
    }
  }
  global <- do.call(rbind, global_rows)

  fold_wise_path <- file.path(task_dir, "metrics_by_fold.csv")
  if (!file.exists(fold_wise_path)) {
    stop("metrics_by_fold.csv not found in ", task_dir,
         "; run_drug() must run before the leakage check")
  }
  fold_wise <- read_output_csv(fold_wise_path)
  metric_names <- setdiff(names(global), c("model", "repeat", "fold"))

  rows <- list()
  for (model in models) {
    for (metric in metric_names) {
      fw <- fold_wise[[metric]][fold_wise$model == model]
      gl <- global[[metric]][global$model == model]
      mean_fw <- mean(fw, na.rm = TRUE)
      mean_gl <- mean(gl, na.rm = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(
        model = model, metric = metric,
        mean_fold_wise = mean_fw, mean_global = mean_gl,
        difference = mean_gl - mean_fw,
        n_folds = sum(!is.na(gl)),
        stringsAsFactors = FALSE
      )
    }
  }
  leakage <- do.call(rbind, rows)
  write_output_csv(leakage, file.path(task_dir, "leakage_check.csv"))
  invisible(leakage)
}

format_config <- function(config) {
  vapply(names(config), function(name) {
    value <- config[[name]]
    text <- if (is.null(value)) {
      "NULL"
    } else if (is.list(value)) {
      paste0("list of ", length(value))
    } else {
      paste(format(value), collapse = ", ")
    }
    sprintf("  %s = %s", name, text)
  }, character(1))
}

package_version_line <- function(package) {
  version <- tryCatch(as.character(utils::packageVersion(package)),
                      error = function(e) "not installed")
  sprintf("  %s %s", package, version)
}

write_run_log <- function(path, config, args, timings, packages) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lines <- c(
    "v1.0 run log",
    paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Command line: Rscript run_all.R", paste(args, collapse = " ")),
    "",
    "Configuration:",
    format_config(config),
    "",
    "Seeds (docs/v1_spec.md, section 4.3):",
    sprintf("  base seed S = %s", config$seed),
    "  outer folds of task t: S + t (regression t = 1, classification t = 2)",
    "  split j of task t: inner folds S + 1000 * t + j, random forest grid draw + 100,",
    "  ranger seed + 1, svm and glmnet set.seed + 5",
    sprintf("  threads: %s", config$threads),
    "",
    paste("R version:", R.version.string),
    paste("Platform:", R.version$platform),
    "",
    "Package versions:",
    vapply(packages, package_version_line, character(1)),
    ""
  )
  if (!is.null(timings) && length(timings) > 0L) {
    timing_lines <- if (is.data.frame(timings)) {
      utils::capture.output(print(timings, row.names = FALSE))
    } else {
      vapply(names(timings), function(name) {
        sprintf("  %s: %.1f s", name, as.numeric(timings[[name]]))
      }, character(1))
    }
    lines <- c(lines, "Timings:", timing_lines, "")
  }
  lines <- c(lines, "sessionInfo():", utils::capture.output(print(utils::sessionInfo())))
  writeLines(lines, path)
  invisible(path)
}
