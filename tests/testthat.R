# tests/testthat.R
#
# Runs the v1.0 test suite without a package structure:
#   Rscript tests/testthat.R
# or, equivalently, from the repository root:
#   Rscript -e "testthat::test_dir('tests/testthat')"
# Every test uses synthetic data generated with a fixed seed; no raw data
# file is read.

script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1L) dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/")) else getwd()
})

test_dir_path <- if (dir.exists(file.path(script_dir, "testthat"))) {
  file.path(script_dir, "testthat")
} else {
  file.path(script_dir, "tests", "testthat")
}

library(testthat)
results <- test_dir(test_dir_path, reporter = "summary", stop_on_failure = TRUE)
invisible(results)
