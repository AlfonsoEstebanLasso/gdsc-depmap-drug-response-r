# Outer and inner folds, stratification and the saved partition file.

m <- make_synthetic_matrix()

check_partition <- function(folds, n, k, repeats, strata) {
  expect_true(all(c("row", "cell_line", "repeat", "fold") %in% names(folds)))
  expect_equal(nrow(folds), n * repeats)
  expect_equal(sort(unique(folds[["repeat"]])), seq_len(repeats))
  for (r in seq_len(repeats)) {
    sub <- folds[folds[["repeat"]] == r, , drop = FALSE]
    expect_equal(sort(sub$row), seq_len(n))
    expect_equal(sort(unique(sub$fold)), seq_len(k))
    counts <- table(strata[sub$row], sub$fold)
    expect_true(all(counts > 0))
  }
}

test_that("make_folds partitions every repeat into k stratified disjoint folds", {
  for (task in c("regression", "classification")) {
    folds <- make_folds(m$auc, k = 3, repeats = 2, task = task, seed = 5)
    strata <- stratum_for(m$auc, task)
    expect_true(is.factor(strata))
    expect_equal(length(strata), length(m$auc))
    check_partition(folds, n = length(m$auc), k = 3, repeats = 2, strata = strata)
    expect_equal(folds$cell_line, m$cell_lines[folds$row])
  }
})

test_that("make_folds is reproducible for the same seed and changes with another", {
  a <- make_folds(m$auc, k = 3, repeats = 2, task = "regression", seed = 5)
  b <- make_folds(m$auc, k = 3, repeats = 2, task = "regression", seed = 5)
  c <- make_folds(m$auc, k = 3, repeats = 2, task = "regression", seed = 6)
  expect_identical(a, b)
  # The assignment of rows to folds changes with the seed (the data frame
  # may be sorted by fold, so the fold column alone is not informative).
  key <- function(f) sort(paste(f[["repeat"]], f$row, f$fold))
  expect_false(identical(key(a), key(c)))
})

test_that("split_rows returns disjoint and exhaustive train and test rows", {
  folds <- make_folds(m$auc, k = 3, repeats = 1, task = "regression", seed = 5)
  for (f in 1:3) {
    idx <- split_rows(folds, 1, f)
    expect_length(intersect(idx$train, idx$test), 0)
    expect_equal(sort(c(idx$train, idx$test)), seq_along(m$auc))
    expect_equal(sort(idx$test), sort(folds$row[folds$fold == f]))
  }
})

test_that("write_folds round trip keeps both tasks", {
  reg <- make_folds(m$auc, k = 3, repeats = 2, task = "regression", seed = 5)
  cls <- make_folds(m$auc, k = 3, repeats = 2, task = "classification", seed = 6)
  path <- file.path(tempfile("folds"), "folds.csv")
  dir.create(dirname(path), recursive = TRUE)
  write_folds(reg, cls, path)
  expect_true(file.exists(path))
  saved <- read_csv_keep_names(path)
  expect_equal(names(saved), c("cell_line", "task", "repeat", "fold"))
  expect_equal(nrow(saved), 2 * nrow(reg))
  expect_setequal(unique(saved$task), c("regression", "classification"))
  saved_reg <- saved[saved$task == "regression", ]
  key_saved <- paste(saved_reg$cell_line, saved_reg[["repeat"]], saved_reg$fold)
  key_reg <- paste(reg$cell_line, reg[["repeat"]], reg$fold)
  expect_setequal(key_saved, key_reg)
  saved_cls <- saved[saved$task == "classification", ]
  key_saved <- paste(saved_cls$cell_line, saved_cls[["repeat"]], saved_cls$fold)
  key_cls <- paste(cls$cell_line, cls[["repeat"]], cls$fold)
  expect_setequal(key_saved, key_cls)
  # No AUC or expression values in the partition file.
  expect_false(any(c("auc", "AUC") %in% names(saved)))
})

test_that("make_inner_folds has the right length and keeps every stratum in each fold", {
  train <- 1:45
  y_reg <- m$auc[train]
  foldid <- make_inner_folds(y_reg, k_inner = 3, task = "regression", seed = 9)
  expect_length(foldid, length(train))
  expect_equal(sort(unique(foldid)), 1:3)
  expect_true(all(table(stratum_for(y_reg, "regression"), foldid) > 0))

  y_cls <- class_rule(m$auc[train], class_threshold(m$auc[train]))
  foldid_cls <- make_inner_folds(y_cls, k_inner = 3, task = "classification", seed = 9)
  expect_length(foldid_cls, length(train))
  expect_true(all(table(y_cls, foldid_cls) > 0))

  expect_identical(foldid, make_inner_folds(y_reg, k_inner = 3, task = "regression", seed = 9))
})
