# ----------------------------------------------------------------------------
# R/cv.R: cross-validation partitions and fold-wise preprocessing (v1.0,
# data builder).
#
# Purpose: repeated stratified k-fold partitions with saved indices (errata 4
# and 8), inner folds shared by every model family (section 4.2 of
# docs/v1_spec.md), and the preprocessing that is fitted on the training rows
# of each outer split only: selection of the 200 most variable genes, class
# threshold at the first quartile of the training AUC (sensitive = low AUC,
# erratum 5) and standardisation with training mean and standard deviation
# (errata 6 and 7).
#
# Inputs: the list returned by build_matrix() (x: cell lines by genes, auc:
# named numeric vector aligned with the rows).
# Outputs: fold tables (data.frame row, cell_line, repeat, fold), integer
# inner foldid vectors, and the preprocessed split (x_train, x_test, y_train,
# y_test, auc_train, auc_test, genes, center, scale, threshold).
# Side effect: write_folds() writes outputs/<drug>/folds.csv (names and fold
# numbers only, no AUC or expression values).
#
# Seeds: set by the caller through the seed arguments (section 4.3 of the
# specification); nothing here calls set.seed() without an explicit seed.
# ----------------------------------------------------------------------------

.tasks <- c("regression", "classification")
.class_levels <- c("resistant", "sensitive")

.check_task <- function(task) {
  task <- match.arg(task, .tasks)
  task
}

# Sample variance of every column, base R, same definition as
# resample::colVars() used in 2019: colSums((x - colMeans(x))^2) / (n - 1).
col_vars <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (n < 2) {
    stop("col_vars() needs at least two rows.", call. = FALSE)
  }
  centred <- x - matrix(colMeans(x), nrow = n, ncol = ncol(x), byrow = TRUE)
  v <- colSums(centred^2) / (n - 1)
  names(v) <- colnames(x)
  v
}

# Factor used to stratify the folds (section 4.1):
#   regression: quartile group of the AUC (cut at the three quartiles,
#               include.lowest = TRUE);
#   classification: the indicator auc < first quartile. This provisional
#               indicator only balances the folds; the labels used for
#               training and evaluation are recomputed on each training fold.
stratum_for <- function(auc, task) {
  task <- .check_task(task)
  auc <- as.numeric(auc)
  if (anyNA(auc)) stop("'auc' has missing values.", call. = FALSE)
  if (task == "regression") {
    breaks <- unique(quantile(auc, probs = c(0, 0.25, 0.5, 0.75, 1),
                              names = FALSE))
    if (length(breaks) < 2) {
      return(factor(rep("q1", length(auc))))
    }
    g <- cut(auc, breaks = breaks, include.lowest = TRUE)
    factor(paste0("q", as.integer(g)))
  } else {
    class_rule(auc, class_threshold(auc))
  }
}

# Repeated stratified k-fold partition built with caret::createMultiFolds().
# auc: numeric vector named by cell line (names(auc) are the cell_line
# values; the row index is used when names are missing).
# Returns data.frame(row, cell_line, repeat, fold): fold is the fold in which
# the row is held out, one row per cell line and repeat, ordered by repeat,
# fold and row.
make_folds <- function(auc, k, repeats, task, seed) {
  task <- .check_task(task)
  n <- length(auc)
  if (n < k) stop("Fewer rows than folds.", call. = FALSE)
  cell_line <- names(auc)
  if (is.null(cell_line)) cell_line <- as.character(seq_len(n))
  if (anyDuplicated(cell_line)) {
    stop("Cell line names must be unique.", call. = FALSE)
  }
  stratum <- stratum_for(auc, task)

  set.seed(seed)
  train_sets <- caret::createMultiFolds(stratum, k = k, times = repeats)

  # Names are "Fold<f>.Rep<r>"; parse them instead of relying on the order.
  nm <- names(train_sets)
  f <- as.integer(sub("^Fold0*([0-9]+)\\.Rep0*([0-9]+)$", "\\1", nm))
  r <- as.integer(sub("^Fold0*([0-9]+)\\.Rep0*([0-9]+)$", "\\2", nm))
  if (anyNA(f) || anyNA(r)) {
    stop("Unexpected fold names from caret::createMultiFolds().", call. = FALSE)
  }

  pieces <- lapply(seq_along(train_sets), function(i) {
    held_out <- setdiff(seq_len(n), train_sets[[i]])
    data.frame(row = held_out, cell_line = cell_line[held_out],
               "repeat" = r[i], fold = f[i],
               stringsAsFactors = FALSE, check.names = FALSE)
  })
  folds <- do.call(rbind, pieces)
  folds <- folds[order(folds[["repeat"]], folds$fold, folds$row), ,
                 drop = FALSE]
  rownames(folds) <- NULL

  # Every repeat must partition all rows exactly once.
  for (rr in unique(folds[["repeat"]])) {
    rows_r <- folds$row[folds[["repeat"]] == rr]
    if (length(rows_r) != n || anyDuplicated(rows_r) ||
        !setequal(rows_r, seq_len(n))) {
      stop("Repeat ", rr, " does not partition the rows.", call. = FALSE)
    }
  }
  attr(folds, "task") <- task
  attr(folds, "seed") <- seed
  folds
}

# Training and held-out row indices of one outer split.
split_rows <- function(folds, rep, fold) {
  sel <- folds[["repeat"]] == rep
  if (!any(sel)) stop("Repeat ", rep, " not found.", call. = FALSE)
  rows_r <- folds$row[sel]
  test <- sort(folds$row[sel & folds$fold == fold])
  if (length(test) == 0) {
    stop("Fold ", fold, " of repeat ", rep, " not found.", call. = FALSE)
  }
  train <- sort(setdiff(rows_r, test))
  list(train = train, test = test)
}

# Inner stratified k-fold assignment of the training rows, shared by every
# model family of the split (section 4.2). y_train is the training AUC
# (regression) or the fold-wise class factor (classification); a numeric
# y_train for classification is converted with the fold-wise rule first.
# Returns an integer vector in 1..k_inner of length n_train.
make_inner_folds <- function(y_train, k_inner, task, seed) {
  task <- .check_task(task)
  n <- length(y_train)
  if (n < k_inner) stop("Fewer training rows than inner folds.", call. = FALSE)
  if (task == "classification") {
    stratum <- if (is.factor(y_train)) {
      factor(as.character(y_train), levels = .class_levels)
    } else {
      stratum_for(y_train, "classification")
    }
  } else {
    stratum <- stratum_for(y_train, "regression")
  }
  set.seed(seed)
  foldid <- caret::createFolds(stratum, k = k_inner, list = FALSE,
                               returnTrain = FALSE)
  foldid <- as.integer(foldid)
  if (length(foldid) != n || length(unique(foldid)) != k_inner) {
    stop("Inner folds could not be built with k_inner = ", k_inner, ".",
         call. = FALSE)
  }
  foldid
}

# The n_genes genes with the largest sample variance over the rows of
# x_train (ties broken by column order). Genes with zero variance are never
# selected.
select_genes <- function(x_train, n_genes = 200) {
  v <- col_vars(x_train)
  candidates <- which(is.finite(v) & v > 0)
  if (length(candidates) < n_genes) {
    stop("Only ", length(candidates), " genes with positive variance; ",
         n_genes, " requested.", call. = FALSE)
  }
  # order() is stable: ties keep the column order.
  ord <- candidates[order(-v[candidates])]
  colnames(x_train)[ord[seq_len(n_genes)]]
}

# First quartile of the training AUC (R default quantile type 7).
class_threshold <- function(auc_train) {
  auc_train <- as.numeric(auc_train)
  if (anyNA(auc_train)) stop("'auc_train' has missing values.", call. = FALSE)
  quantile(auc_train, probs = 0.25, names = FALSE, type = 7)
}

# Sensitive when AUC is strictly below the threshold, resistant otherwise.
# Levels: resistant, sensitive (sensitive is the positive class).
class_rule <- function(auc, threshold) {
  auc <- as.numeric(auc)
  lab <- ifelse(auc < threshold, "sensitive", "resistant")
  factor(lab, levels = .class_levels)
}

# Preprocess one outer split using the training rows only, then apply to
# the held-out rows:
#   1. gene selection on x[train, ] (selection = "fold") or, for the leakage
#      check only, on all rows of x (selection = "global");
#   2. class label (classification): threshold on auc[train], or the explicit
#      threshold given (leakage check only);
#   3. standardisation of the selected genes with training mean and sd.
preprocess_fold <- function(x, auc, train, test, task, n_genes = 200,
                            selection = c("fold", "global"),
                            threshold = NULL) {
  task <- .check_task(task)
  selection <- match.arg(selection)
  if (length(intersect(train, test)) > 0) {
    stop("'train' and 'test' overlap.", call. = FALSE)
  }
  if (length(auc) != nrow(x)) {
    stop("'auc' must have one value per row of 'x'.", call. = FALSE)
  }

  genes <- if (selection == "fold") {
    select_genes(x[train, , drop = FALSE], n_genes)
  } else {
    select_genes(x, n_genes)
  }

  x_tr <- x[train, genes, drop = FALSE]
  center <- colMeans(x_tr)
  scale <- sqrt(col_vars(x_tr))
  names(center) <- genes
  names(scale) <- genes
  if (any(!is.finite(scale) | scale <= 0)) {
    stop("A selected gene has zero variance on the training rows.",
         call. = FALSE)
  }
  x_train <- .standardise(x_tr, center, scale)
  x_test <- .standardise(x[test, genes, drop = FALSE], center, scale)

  auc_train <- as.numeric(auc[train])
  auc_test <- as.numeric(auc[test])
  names(auc_train) <- rownames(x)[train]
  names(auc_test) <- rownames(x)[test]

  if (task == "classification") {
    thr <- if (is.null(threshold)) class_threshold(auc_train) else {
      as.numeric(threshold)
    }
    y_train <- class_rule(auc_train, thr)
    y_test <- class_rule(auc_test, thr)
  } else {
    thr <- NA_real_
    y_train <- auc_train
    y_test <- auc_test
  }

  list(
    x_train   = x_train,
    x_test    = x_test,
    y_train   = y_train,
    y_test    = y_test,
    auc_train = auc_train,
    auc_test  = auc_test,
    genes     = genes,
    center    = center,
    scale     = scale,
    threshold = thr,
    selection = selection
  )
}

# Centre and scale with given parameters; returns a plain matrix with the
# same dimnames (no "scaled:" attributes).
.standardise <- function(m, center, scale) {
  out <- scale(m, center = center, scale = scale)
  attr(out, "scaled:center") <- NULL
  attr(out, "scaled:scale") <- NULL
  dimnames(out) <- dimnames(m)
  out
}

# outputs/<drug>/folds.csv: cell_line, task, repeat, fold. Either table may
# be NULL when only one task is run. Names and fold numbers only.
write_folds <- function(folds_regression, folds_classification, path) {
  pieces <- list()
  if (!is.null(folds_regression)) {
    pieces[[length(pieces) + 1]] <- data.frame(
      cell_line = folds_regression$cell_line, task = "regression",
      "repeat" = folds_regression[["repeat"]], fold = folds_regression$fold,
      stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!is.null(folds_classification)) {
    pieces[[length(pieces) + 1]] <- data.frame(
      cell_line = folds_classification$cell_line, task = "classification",
      "repeat" = folds_classification[["repeat"]],
      fold = folds_classification$fold,
      stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (length(pieces) == 0) stop("Nothing to write.", call. = FALSE)
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  dir <- dirname(path)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  utils::write.csv(out, path, row.names = FALSE, quote = TRUE)
  invisible(out)
}

# Read a folds.csv written by write_folds() and rebuild the per-task tables
# (list(regression = , classification = ), with the row index recovered from
# cell_lines, the row names of the matrix).
read_folds <- function(path, cell_lines) {
  tab <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  out <- list()
  for (task in unique(tab$task)) {
    t <- tab[tab$task == task, , drop = FALSE]
    row <- match(t$cell_line, cell_lines)
    if (anyNA(row)) {
      stop("folds.csv names a cell line absent from the matrix.", call. = FALSE)
    }
    f <- data.frame(row = row, cell_line = t$cell_line,
                    "repeat" = t[["repeat"]], fold = t$fold,
                    stringsAsFactors = FALSE, check.names = FALSE)
    f <- f[order(f[["repeat"]], f$fold, f$row), , drop = FALSE]
    rownames(f) <- NULL
    attr(f, "task") <- task
    out[[task]] <- f
  }
  out
}
