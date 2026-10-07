# build_matrix() on tiny synthetic GDSC, metadata and TPM inputs.

synthetic_inputs <- function() {
  gdsc <- data.frame(
    DRUG_NAME = c("DrugA", "DrugA", "DrugA", "DrugA", "DrugB"),
    CELL_LINE_NAME = c("B-LINE", "A-LINE", "C-LINE", "D-LINE", "A-LINE"),
    COSMIC_ID = c(1, 2, 3, 4, 2),
    DRUG_ID = c(1L, 1L, 1L, 1L, 2L),
    AUC = c(0.9, 0.5, 0.7, 0.6, 0.8),
    stringsAsFactors = FALSE
  )
  # COSMIC 3 has no DepMap identifier and COSMIC 4 has no expression row.
  metadata <- data.frame(
    DepMap_ID = c("ACH-000001", "ACH-000002", "ACH-000004", "ACH-000009"),
    COSMIC_ID = c(1, 2, 4, 9),
    stringsAsFactors = FALSE
  )
  genes <- c("G5 (ENSG5)", "G1 (ENSG1)", "G3 (ENSG3)", "G2 (ENSG2)", "G4 (ENSG4)")
  set.seed(7)
  x <- matrix(round(stats::runif(15, 0, 10), 2), nrow = 3, ncol = 5, dimnames = list(NULL, genes))
  tpm <- list(ids = c("ACH-000002", "ACH-000001", "ACH-000009"), x = x)
  list(gdsc = gdsc, metadata = metadata, tpm = tpm)
}

test_that("build_matrix keeps only matched lines, sorted by name, with aligned AUC", {
  inp <- synthetic_inputs()
  out <- build_matrix("druga", inp$gdsc, inp$metadata, inp$tpm)
  expect_equal(tolower(out$drug), "druga")
  expect_equal(out$cell_lines, c("A-LINE", "B-LINE"))
  expect_equal(rownames(out$x), c("A-LINE", "B-LINE"))
  expect_equal(nrow(out$x), 2L)
  expect_equal(colnames(out$x), colnames(inp$tpm$x))
  expect_equal(unname(out$auc), c(0.5, 0.9))
  expect_equal(length(out$auc), nrow(out$x))
  # A-LINE is COSMIC 2, DepMap ACH-000002, first TPM row; B-LINE is the second row.
  expect_equal(unname(out$x["A-LINE", ]), unname(inp$tpm$x[1, ]))
  expect_equal(unname(out$x["B-LINE", ]), unname(inp$tpm$x[2, ]))
  expect_true(is.numeric(out$x))
  expect_setequal(out$depmap_ids, c("ACH-000002", "ACH-000001"))
})

test_that("build_matrix stops on an ambiguous COSMIC_ID", {
  inp <- synthetic_inputs()
  metadata <- rbind(inp$metadata, data.frame(DepMap_ID = "ACH-000022", COSMIC_ID = 2))
  tpm <- inp$tpm
  tpm$ids <- c(tpm$ids, "ACH-000022")
  tpm$x <- rbind(tpm$x, tpm$x[1, ])
  expect_error(build_matrix("DrugA", inp$gdsc, metadata, tpm), "2")
})

test_that("build_matrix stops on an unknown drug and on missing values", {
  inp <- synthetic_inputs()
  expect_error(build_matrix("nosuchdrug", inp$gdsc, inp$metadata, inp$tpm))
  tpm <- inp$tpm
  tpm$x[1, 2] <- NA
  expect_error(build_matrix("DrugA", inp$gdsc, inp$metadata, tpm))
  gdsc <- inp$gdsc
  gdsc$AUC[2] <- NA
  expect_error(build_matrix("DrugA", gdsc, inp$metadata, inp$tpm))
})

test_that("resolve_drug is case-insensitive and rejects absent or ambiguous names", {
  inp <- synthetic_inputs()
  expect_equal(resolve_drug("DRUGA", inp$gdsc), "DrugA")
  expect_equal(resolve_drug("drugb", inp$gdsc), "DrugB")
  expect_error(resolve_drug("drugz", inp$gdsc))
  ambiguous <- rbind(inp$gdsc, data.frame(DRUG_NAME = "DrugA", CELL_LINE_NAME = "E-LINE",
                                          COSMIC_ID = 5, DRUG_ID = 99L, AUC = 0.4))
  expect_error(resolve_drug("DrugA", ambiguous))
})
