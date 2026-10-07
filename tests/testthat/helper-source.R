# tests/testthat/helper-source.R
#
# Sources R/*.R from the repository root (no package structure) and defines
# the synthetic data used by every test. No raw data file is read by the
# tests.

find_project_root <- function(start = getwd()) {
  dir <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(dir, "run_all.R")) && dir.exists(file.path(dir, "R"))) {
      return(dir)
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) stop("repository root (run_all.R and R/) not found above ", start)
    dir <- parent
  }
}

project_root <- find_project_root()

for (f in c("build_matrix.R", "cv.R", "metrics.R", "models.R", "run_drug.R", "report.R")) {
  source(file.path(project_root, "R", f))
}

# Synthetic cell line by gene matrix: n lines, p genes, a planted linear
# signal between the first n_signal genes (given a larger variance so that
# fold-wise selection keeps them) and the AUC, which lies in (0, 1).
make_synthetic_matrix <- function(n = 60, p = 500, n_signal = 10, seed = 1,
                                  drug = "synthetic") {
  set.seed(seed)
  x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  x[, seq_len(n_signal)] <- 3 * x[, seq_len(n_signal)]
  beta <- rep(c(0.25, -0.25), length.out = n_signal)
  linear <- as.numeric(x[, seq_len(n_signal), drop = FALSE] %*% beta)
  auc <- stats::plogis(linear + stats::rnorm(n, sd = 0.5))
  cell_lines <- sprintf("CL%03d", seq_len(n))
  genes <- sprintf("GENE%d (ENSG%08d)", seq_len(p), seq_len(p))
  dimnames(x) <- list(cell_lines, genes)
  names(auc) <- cell_lines
  list(drug = drug, cell_lines = cell_lines, x = x, auc = auc,
       depmap_ids = sprintf("ACH-%06d", seq_len(n)))
}

# Small grids so that the model tests run in seconds (the inner resampling
# keeps three folds because cv.glmnet refuses fewer): two random forest
# candidates with 50 trees, the single ridge and lasso candidates, two
# elastic net alphas and two SVM candidates, all taken from the declared
# grids of model_grid().
tiny_grids <- function(task, num_trees = 50) {
  list(
    rf = model_grid("rf", task, rf_budget = 2, seed = 1, num_trees = num_trees),
    ridge = model_grid("ridge", task),
    lasso = model_grid("lasso", task),
    enet = utils::head(model_grid("enet", task), 2),
    svm_linear = utils::head(model_grid("svm_linear", task), 2)
  )
}

tiny_config <- function(task, leakage_check = FALSE) {
  config <- default_config()
  config$k <- 2L
  config$repeats <- 1L
  config$k_inner <- 3L
  config$rf_budget <- 2L
  config$num_trees <- 50L
  config$threads <- 1L
  config$leakage_check <- leakage_check
  config$grids <- tiny_grids(task, num_trees = 50)
  config
}

# glmnet warns when a binomial class has fewer than eight observations or
# an inner fold has fewer than ten, which is expected on the small synthetic
# folds; those two warnings alone are muffled so that the suite stays
# readable.
quiet_small_classes <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    msg <- conditionMessage(w)
    if (grepl("fewer than 8", msg, fixed = TRUE) || grepl("Too few", msg, fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

# Runs both tasks of the synthetic drug into outdir with the tiny
# configuration; used by the run_drug and report tests.
run_synthetic_example <- function(outdir, seed = 1, leakage_check = FALSE) {
  m <- make_synthetic_matrix(seed = seed)
  results <- list()
  for (task in c("regression", "classification")) {
    config <- tiny_config(task, leakage_check = leakage_check)
    results[[task]] <- quiet_small_classes(run_drug(m$drug, task, m, config, outdir))
  }
  invisible(results)
}

task_files <- function(leakage_check = FALSE) {
  files <- c("metrics_by_fold.csv", "summary.csv", "tuning.csv", "tuning_candidates.csv",
             "selected_genes.csv", "timing.csv", "predictions.csv")
  if (leakage_check) files <- c(files, "leakage_check.csv")
  files
}

read_csv_keep_names <- function(path) {
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}
