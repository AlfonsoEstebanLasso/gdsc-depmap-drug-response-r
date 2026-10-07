# ----------------------------------------------------------------------------
# Purpose:  Rebuild the two DepMap Public 19Q1 input files expected by the
#           legacy scripts from Bioconductor's ExperimentHub (package depmap),
#           because the DepMap portal sits behind an interactive check and
#           figshare does not host the 19Q1 expression release.
# Inputs:   network access to ExperimentHub (records EH2264, TPM_19Q1, and
#           EH2266, metadata_19Q1); packages ExperimentHub and depmap
#           (install with --install, which uses BiocManager).
# Outputs:  data/raw/CCLE_depMap_19Q1_TPM.csv and
#           data/raw/DepMap-2019q1-celllines.csv, in the layout of the portal
#           files (same rows, columns and order), plus
#           data/raw/depmap_19Q1_rebuild.txt with the checks performed.
# Checks:   dimensions, order of lines and genes, a digest of the expression
#           values rounded to six decimals and a digest of the COSMIC_ID key,
#           all compared with the values recorded from the portal files used
#           in 2019. The byte-level SHA-256 of the portal files cannot be
#           reproduced (number formatting differs); the content is identical.
#           The ExperimentHub metadata lacks one line of the portal file
#           (ACH-001825, without COSMIC_ID), with no effect on the join.
# Usage:    Rscript scripts/01_depmap_from_bioconductor.R [--install]
# Exit:     0 when both files are written and every check passes; 1 otherwise.
# Depends:  data.table, digest (renv.lock); ExperimentHub, depmap (Bioconductor,
#           not recorded in renv.lock, installed on demand).
# ----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
raw_dir <- file.path("data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

if ("--install" %in% args) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(c("ExperimentHub", "depmap"), update = FALSE, ask = FALSE)
}
for (pkg in c("ExperimentHub", "depmap", "data.table", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package %s is not installed. Run this script with --install, or install ExperimentHub and depmap with BiocManager.", pkg), call. = FALSE)
  }
}
suppressPackageStartupMessages({ library(ExperimentHub); library(data.table); library(digest) })

expected <- list(
  n_lines = 1165L, n_genes = 57820L, n_cosmic = 982L,
  matrix = "5c04ea6dced9723dc67ddd596f0881968df0382e8a896fc3caf29e9d15e18a07",
  genes  = "0e41eda854e8d8ec12841c37f153c7c9c69f7c990006f9e07b98013745094ab6",
  lines  = "244af96058a619b15a51ece2d65641a71f5915cd72067a4145c5807670034254",
  key    = "505740f46f016bd7174896521fcaacb98a995dd831f103210e377a0d5310fdcb"
)
string_digest <- function(x) digest(paste(x, collapse = "
"), algo = "sha256", serialize = FALSE)

options(timeout = 3600)
eh <- ExperimentHub()
message("Fetching metadata_19Q1 (EH2266) and TPM_19Q1 (EH2264) from ExperimentHub; the expression record is about 500 MB and is cached locally.")
meta <- eh[["EH2266"]]
tpm <- eh[["EH2264"]]

# 1. Cell line metadata in the layout of DepMap-2019q1-celllines.csv
celllines <- data.table(
  DepMap_ID = meta$depmap_id, CCLE_Name = meta$cell_line, Aliases = meta$aliases,
  COSMIC_ID = meta$cosmic_id, "Sanger ID" = meta$sanger_id,
  "Primary Disease" = meta$primary_disease, "Subtype Disease" = meta$subtype_disease,
  Gender = meta$gender, Source = meta$source)
fwrite(celllines, file.path(raw_dir, "DepMap-2019q1-celllines.csv"))

# 2. Expression matrix in the layout of CCLE_depMap_19Q1_TPM.csv: one row per
#    cell line (unnamed first column with the DepMap_ID), one column per gene
#    named "SYMBOL (ENSG...)", lines and genes in their order of appearance,
#    which is the order of the portal file.
long <- as.data.table(tpm)[, .(depmap_id, gene, expression)]
line_order <- unique(long$depmap_id)
gene_order <- unique(long$gene)
wide <- dcast(long, depmap_id ~ gene, value.var = "expression")
setcolorder(wide, c("depmap_id", gene_order))
wide <- wide[match(line_order, depmap_id)]
setnames(wide, "depmap_id", "")
fwrite(wide, file.path(raw_dir, "CCLE_depMap_19Q1_TPM.csv"))

# 3. Checks against the values recorded from the portal files
m <- as.matrix(wide[, -1]); storage.mode(m) <- "double"
key <- celllines[!is.na(COSMIC_ID) & COSMIC_ID != "", .(DepMap_ID, COSMIC_ID = as.character(COSMIC_ID))][order(DepMap_ID)]
checks <- c(
  lines  = nrow(m) == expected$n_lines,
  genes  = ncol(m) == expected$n_genes,
  matrix = identical(digest(round(m, 6), algo = "sha256"), expected$matrix),
  gene_order = identical(string_digest(colnames(m)), expected$genes),
  line_order = identical(string_digest(wide[[1]]), expected$lines),
  cosmic = nrow(key) == expected$n_cosmic && identical(digest(key, algo = "sha256"), expected$key)
)
note <- c(
  sprintf("DepMap Public 19Q1 rebuilt from ExperimentHub on %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "Records: EH2266 (metadata_19Q1), EH2264 (TPM_19Q1); package depmap.",
  sprintf("Expression matrix: %d cell lines x %d genes.", nrow(m), ncol(m)),
  sprintf("Cell line metadata: %d rows, %d with COSMIC_ID.", nrow(celllines), nrow(key)),
  paste0("Check ", names(checks), ": ", ifelse(checks, "OK", "FAILED")),
  "Byte-level SHA-256 of the portal files is not reproduced (number formatting); content checks above compare values, order and the COSMIC_ID key."
)
writeLines(note, file.path(raw_dir, "depmap_19Q1_rebuild.txt"))
cat(note, sep = "
")
if (!all(checks)) {
  stop("At least one content check failed; the rebuilt files do not match the 2019 inputs.", call. = FALSE)
}
cat("Both DepMap 19Q1 files rebuilt and verified.
")
