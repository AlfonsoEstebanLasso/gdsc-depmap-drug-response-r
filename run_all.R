# run_all.R
#
# Purpose: command line entry point of v1.0 (docs/v1_spec.md, section 11).
# Builds or loads the per-drug matrices, builds the outer folds once per drug
# and task (written to folds.csv and passed to run_drug()), runs the nested
# cross-validation of every requested drug and task, writes the per-fold
# metrics, the summaries, the chosen hyperparameters, the figures and a run
# log with the R and package versions (sessionInfo() without the operating
# system build, locale and time zone).
#
# Usage (from the repository root, so that .Rprofile activates renv):
#   Rscript run_all.R [--fast] [--drug erlotinib,paclitaxel]
#                     [--task regression|classification|both] [--outdir outputs]
#                     [--seed 2019] [--folds 5] [--repeats 3] [--inner-folds 5]
#                     [--rf-budget 12] [--threads N] [--no-leakage-check]
#                     [--refresh-cache] [--save-models]
#
# Inputs: data/raw/ (the three files of data/README.md, verified by
# scripts/00_download_data.R) or the cached matrices in data/derived/.
# Outputs: the layout of docs/v1_spec.md, section 12, under --outdir.
# Argument parsing uses base commandArgs() only.

source_v1 <- function(root = ".") {
  files <- c("build_matrix.R", "cv.R", "metrics.R", "models.R", "run_drug.R", "report.R")
  for (f in files) {
    source(file.path(root, "R", f), local = FALSE)
  }
  invisible(files)
}

usage_text <- function() {
  paste(
    "Usage: Rscript run_all.R [options]",
    "",
    "  --fast                 smoke configuration: 3 folds, 1 repeat, inner 3-fold,",
    "                         random forest budget 6, erlotinib and paclitaxel,",
    "                         no leakage check (explicit options override it)",
    "  --drug <names>         comma-separated GDSC drug names, case-insensitive",
    "                         (default: erlotinib,rapamycin,sunitinib,paclitaxel)",
    "  --task <task>          regression, classification or both (default: both)",
    "  --outdir <dir>         output directory (default: outputs)",
    "  --seed <int>           base seed (default: 2019)",
    "  --folds <int>          outer folds (default: 5)",
    "  --repeats <int>        outer repeats (default: 3)",
    "  --inner-folds <int>    inner folds for tuning (default: 5)",
    "  --rf-budget <int>      random forest configurations drawn (default: 12)",
    "  --threads <int>        ranger threads (default: all cores minus one)",
    "  --no-leakage-check     skip the leakage check of the glmnet family",
    "  --refresh-cache        rebuild data/derived/<drug>_matrix.rds",
    "  --save-models          write the fitted objects per fold (.rds, not committed)",
    "  --help                 show this text",
    sep = "\n"
  )
}

as_positive_integer <- function(value, name) {
  number <- suppressWarnings(as.numeric(value))
  if (is.na(number) || number < 1 || number != floor(number)) {
    stop(sprintf("%s must be a positive integer, got '%s'", name, value), call. = FALSE)
  }
  as.integer(number)
}

parse_args <- function(args) {
  opts <- list(
    fast = FALSE, drug = NULL, task = "both", outdir = "outputs", seed = NULL,
    folds = NULL, repeats = NULL, inner_folds = NULL, rf_budget = NULL,
    threads = NULL, leakage_check = NULL, refresh_cache = FALSE,
    save_models = FALSE, help = FALSE
  )
  flags <- c(
    "--fast" = "fast", "--no-leakage-check" = "no_leakage_check",
    "--refresh-cache" = "refresh_cache", "--save-models" = "save_models",
    "--help" = "help", "-h" = "help"
  )
  valued <- c(
    "--drug" = "drug", "--task" = "task", "--outdir" = "outdir", "--seed" = "seed",
    "--folds" = "folds", "--repeats" = "repeats", "--inner-folds" = "inner_folds",
    "--rf-budget" = "rf_budget", "--threads" = "threads"
  )
  i <- 1L
  while (i <= length(args)) {
    arg <- args[i]
    value <- NULL
    if (grepl("^--[a-z-]+=", arg)) {
      value <- sub("^[^=]*=", "", arg)
      arg <- sub("=.*$", "", arg)
    }
    if (arg %in% names(flags)) {
      if (!is.null(value)) stop(sprintf("%s takes no value", arg), call. = FALSE)
      key <- flags[[arg]]
      if (key == "no_leakage_check") opts$leakage_check <- FALSE else opts[[key]] <- TRUE
      i <- i + 1L
    } else if (arg %in% names(valued)) {
      if (is.null(value)) {
        if (i == length(args)) stop(sprintf("%s needs a value", arg), call. = FALSE)
        value <- args[i + 1L]
        i <- i + 2L
      } else {
        i <- i + 1L
      }
      key <- valued[[arg]]
      opts[[key]] <- if (key %in% c("seed", "folds", "repeats", "inner_folds", "rf_budget", "threads")) {
        as_positive_integer(value, arg)
      } else {
        value
      }
    } else {
      stop(sprintf("unknown argument '%s'\n%s", arg, usage_text()), call. = FALSE)
    }
  }
  opts$task <- match.arg(tolower(opts$task), c("both", "regression", "classification"))
  opts
}

build_config <- function(opts) {
  config <- default_config(fast = opts$fast)
  if (!is.null(opts$seed)) config$seed <- opts$seed
  if (!is.null(opts$folds)) config$k <- opts$folds
  if (!is.null(opts$repeats)) config$repeats <- opts$repeats
  if (!is.null(opts$inner_folds)) config$k_inner <- opts$inner_folds
  if (!is.null(opts$rf_budget)) config$rf_budget <- opts$rf_budget
  if (!is.null(opts$threads)) config$threads <- opts$threads
  if (!is.null(opts$leakage_check)) config$leakage_check <- opts$leakage_check
  config$save_models <- isTRUE(opts$save_models)
  config
}

requested_drugs <- function(opts) {
  if (!is.null(opts$drug)) {
    drugs <- tolower(trimws(strsplit(opts$drug, ",", fixed = TRUE)[[1]]))
    drugs <- unique(drugs[nzchar(drugs)])
    if (length(drugs) == 0L) stop("--drug needs at least one name", call. = FALSE)
    return(drugs)
  }
  if (opts$fast) c("erlotinib", "paclitaxel") else c("erlotinib", "rapamycin", "sunitinib", "paclitaxel")
}

requested_tasks <- function(opts) {
  if (opts$task == "both") c("regression", "classification") else opts$task
}

log_line <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...)))
}

raw_file_names <- function() {
  c(gdsc = "v17.3_fitted_dose_response.xlsx",
    metadata = "DepMap-2019q1-celllines.csv",
    tpm = "CCLE_depMap_19Q1_TPM.csv")
}

check_raw_files <- function(raw_dir) {
  paths <- file.path(raw_dir, raw_file_names())
  missing <- raw_file_names()[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing raw files in ", raw_dir, ": ", paste(missing, collapse = ", "),
         ". Run 'Rscript scripts/00_download_data.R' first (see data/README.md).",
         call. = FALSE)
  }
  invisible(paths)
}

# DepMap identifiers of every requested drug, so that the expression file is
# read once per run and only the needed rows are kept (spec, section 2.2).
depmap_ids_for <- function(drugs, gdsc, metadata) {
  ids <- character(0)
  for (drug in drugs) {
    name <- resolve_drug(drug, gdsc)
    rows <- gdsc[gdsc$DRUG_NAME == name, c("CELL_LINE_NAME", "COSMIC_ID"), drop = FALSE]
    joined <- merge(rows, metadata[, c("DepMap_ID", "COSMIC_ID"), drop = FALSE], by = "COSMIC_ID")
    ids <- c(ids, as.character(joined$DepMap_ID))
  }
  unique(ids)
}

load_matrices <- function(drugs, raw_dir, cache_dir, refresh) {
  cache_paths <- file.path(cache_dir, paste0(drugs, "_matrix.rds"))
  needs_build <- refresh | !file.exists(cache_paths)
  names(needs_build) <- drugs
  tpm <- NULL
  if (any(needs_build)) {
    paths <- check_raw_files(raw_dir)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    to_build <- drugs[needs_build]
    if (length(to_build) > 1L) {
      log_line("reading GDSC and DepMap metadata to resolve the cell lines of %d drugs",
               length(to_build))
      gdsc <- read_gdsc(paths[["gdsc"]])
      metadata <- read_depmap_metadata(paths[["metadata"]])
      keep_ids <- depmap_ids_for(to_build, gdsc, metadata)
      log_line("reading the expression file once for %d cell lines", length(keep_ids))
      tpm <- read_depmap_tpm(paths[["tpm"]], keep_ids = keep_ids)
    }
  }
  matrices <- list()
  for (drug in drugs) {
    log_line("%s: %s matrix", drug, if (needs_build[[drug]]) "building" else "loading cached")
    matrices[[drug]] <- load_drug_matrix(drug, raw_dir = raw_dir, cache_dir = cache_dir,
                                         refresh = needs_build[[drug]], tpm = tpm)
    log_line("%s: %d cell lines, %d genes", drug, nrow(matrices[[drug]]$x), ncol(matrices[[drug]]$x))
  }
  matrices
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opts <- parse_args(args)
  if (opts$help) {
    cat(usage_text(), "\n")
    return(invisible(0L))
  }
  source_v1(".")
  started <- Sys.time()
  config <- build_config(opts)
  drugs <- requested_drugs(opts)
  tasks <- requested_tasks(opts)
  outdir <- opts$outdir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  log_line("v1.0 run: drugs %s; tasks %s; k = %d, repeats = %d, inner folds = %d, rf budget = %d, seed = %d, threads = %d, leakage check = %s",
           paste(drugs, collapse = ", "), paste(tasks, collapse = ", "), config$k,
           config$repeats, config$k_inner, config$rf_budget, config$seed, config$threads,
           config$leakage_check)

  matrices <- load_matrices(drugs, raw_dir = "data/raw", cache_dir = "data/derived",
                            refresh = isTRUE(opts$refresh_cache))

  timings <- list()
  for (drug in drugs) {
    m <- matrices[[drug]]
    auc <- stats::setNames(as.numeric(m$auc), m$cell_lines)
    folds_regression <- make_folds(auc, k = config$k, repeats = config$repeats,
                                   task = "regression", seed = config$seed + task_index("regression"))
    folds_classification <- make_folds(auc, k = config$k, repeats = config$repeats,
                                       task = "classification", seed = config$seed + task_index("classification"))
    drug_dir <- file.path(outdir, drug)
    dir.create(drug_dir, recursive = TRUE, showWarnings = FALSE)
    write_folds(folds_regression, folds_classification, file.path(drug_dir, "folds.csv"))
    folds_by_task <- list(regression = folds_regression, classification = folds_classification)
    for (task in tasks) {
      t0 <- Sys.time()
      log_line("%s %s: starting %d outer splits", drug, task, config$k * config$repeats)
      # The same fold table that was just written: the saved indices and the
      # evaluated folds are one object (erratum 8).
      run_drug(drug, task, m, config, outdir, folds = folds_by_task[[task]])
      seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      timings[[paste(drug, task)]] <- seconds
      log_line("%s %s: done in %.1f s", drug, task, seconds)
    }
  }

  log_line("writing figures and summary tables")
  make_report(outdir, drugs, tasks)
  timings[["total"]] <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  write_run_log(file.path(outdir, "run_log.txt"), config, args, timings,
                packages = c("ranger", "glmnet", "e1071", "pROC", "caret", "ggplot2",
                             "viridisLite", "data.table", "readxl"))
  log_line("finished in %.1f minutes; outputs under %s", timings[["total"]] / 60, outdir)
  invisible(0L)
}

if (sys.nframe() == 0L) {
  main()
}
