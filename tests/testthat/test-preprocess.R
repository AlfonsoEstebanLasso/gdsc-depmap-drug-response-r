# Fold-wise preprocessing: gene selection, class threshold and
# standardisation computed on the training rows only.

m <- make_synthetic_matrix()
train <- 1:45
test <- 46:60

test_that("col_vars equals the sample variance per column", {
  expect_equal(unname(col_vars(m$x)), unname(apply(m$x, 2, stats::var)))
  expect_equal(names(col_vars(m$x)), colnames(m$x))
})

test_that("a gene with huge variance only in the held-out rows is not selected", {
  x <- m$x
  train_vars <- col_vars(x[train, , drop = FALSE])
  leaky <- names(which.min(train_vars))
  x[test, leaky] <- x[test, leaky] + 1000
  expect_true(leaky %in% select_genes(x, 200))
  expect_false(leaky %in% select_genes(x[train, , drop = FALSE], 200))

  pp <- preprocess_fold(x, m$auc, train, test, "regression", n_genes = 200)
  expect_false(leaky %in% pp$genes)
  expect_setequal(pp$genes, select_genes(x[train, , drop = FALSE], 200))

  global <- preprocess_fold(x, m$auc, train, test, "regression", n_genes = 200,
                            selection = "global")
  expect_true(leaky %in% global$genes)
  expect_setequal(global$genes, select_genes(x, 200))
})

test_that("select_genes keeps the n most variable genes in column order", {
  genes <- select_genes(m$x[train, , drop = FALSE], 200)
  expect_length(genes, 200)
  expect_true(all(genes %in% colnames(m$x)))
  vars <- col_vars(m$x[train, , drop = FALSE])
  expect_true(min(vars[genes]) >= max(vars[setdiff(names(vars), genes)]))
  # The ten signal genes were given a larger variance and must be selected.
  expect_true(all(colnames(m$x)[1:10] %in% genes))
})

test_that("the class threshold is the training first quartile and labels follow it", {
  expect_equal(class_threshold(m$auc[train]), unname(stats::quantile(m$auc[train], 0.25)))
  pp <- preprocess_fold(m$x, m$auc, train, test, "classification", n_genes = 200)
  expect_equal(pp$threshold, unname(stats::quantile(m$auc[train], 0.25)))
  expect_true(is.factor(pp$y_train))
  expect_equal(levels(pp$y_train), c("resistant", "sensitive"))
  expect_equal(levels(pp$y_test), c("resistant", "sensitive"))
  expect_equal(as.character(pp$y_train),
               unname(ifelse(m$auc[train] < pp$threshold, "sensitive", "resistant")))
  expect_equal(as.character(pp$y_test),
               unname(ifelse(m$auc[test] < pp$threshold, "sensitive", "resistant")))
  labels <- class_rule(c(0.1, 0.5, 0.9), 0.5)
  expect_equal(as.character(labels), c("sensitive", "resistant", "resistant"))
  expect_equal(levels(labels), c("resistant", "sensitive"))
  explicit <- preprocess_fold(m$x, m$auc, train, test, "classification", n_genes = 200,
                              threshold = 0.5)
  expect_equal(explicit$threshold, 0.5)
})

test_that("the held-out matrix uses the training centre and scale with the same columns", {
  pp <- preprocess_fold(m$x, m$auc, train, test, "regression", n_genes = 200)
  expect_equal(ncol(pp$x_train), 200)
  expect_equal(ncol(pp$x_test), 200)
  expect_equal(colnames(pp$x_train), pp$genes)
  expect_equal(colnames(pp$x_test), pp$genes)
  expect_equal(nrow(pp$x_train), length(train))
  expect_equal(nrow(pp$x_test), length(test))
  expect_equal(unname(colMeans(pp$x_train)), rep(0, 200), tolerance = 1e-10)
  expect_equal(unname(apply(pp$x_train, 2, stats::sd)), rep(1, 200), tolerance = 1e-10)
  centre <- colMeans(m$x[train, pp$genes, drop = FALSE])
  scale <- apply(m$x[train, pp$genes, drop = FALSE], 2, stats::sd)
  expect_equal(unname(pp$center[pp$genes]), unname(centre))
  expect_equal(unname(pp$scale[pp$genes]), unname(scale))
  expected_test <- sweep(sweep(m$x[test, pp$genes, drop = FALSE], 2, centre), 2, scale, "/")
  expect_equal(unname(pp$x_test), unname(expected_test))
  expect_equal(unname(pp$auc_train), unname(m$auc[train]))
  expect_equal(unname(pp$auc_test), unname(m$auc[test]))
  expect_equal(unname(pp$y_train), unname(m$auc[train]))
})
