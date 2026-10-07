#!/usr/bin/env Rscript
# 00_download_data.R
#
# Purpose:  download the public input data used by the legacy scripts into
#           data/raw/ and verify file size and SHA-256. The script aborts with a
#           clear message if a downloaded (or pre-existing) file does not match
#           the expected size or hash.
# Inputs:   none (network access to the GDSC download server and, optionally,
#           the figshare API).
# Outputs:  data/raw/v17.3_fitted_dose_response.xlsx
#           data/raw/CCLE_depMap_19Q1_TPM.csv
#           data/raw/DepMap-2019q1-celllines.csv
# Usage:    Rscript scripts/00_download_data.R
#           Rscript scripts/00_download_data.R --verify-only <directory>
# Exit:     0 all files present and verified; 1 size or SHA-256 mismatch
#           (error); 2 one or more files missing.
# Depends:  digest, jsonlite (recorded in renv.lock).
#
# DepMap source: the two DepMap files belong to the "DepMap Public 19Q1" release.
# At the time of writing, the figshare record 7655150 ("DepMap Achilles 19Q1
# Public") does NOT contain them, so no figshare download is attempted unless a
# record id is given in DEPMAP_FIGSHARE_ARTICLE (constant below or environment
# variable of the same name). See data/README.md for the alternative routes.
# Files placed manually in data/raw/ are verified in any case.

suppressPackageStartupMessages({
  library(digest)
  library(jsonlite)
})

DEPMAP_FIGSHARE_ARTICLE <- Sys.getenv("DEPMAP_FIGSHARE_ARTICLE", unset = "")
GDSC_BASE_URL <- "https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-7.0/"
FIGSHARE_API  <- "https://api.figshare.com/v2/articles/"

files <- data.frame(
  name   = c("v17.3_fitted_dose_response.xlsx",
             "CCLE_depMap_19Q1_TPM.csv",
             "DepMap-2019q1-celllines.csv"),
  size   = c(16064522, 568839199, 154145),
  sha256 = c("3c48e4dae5d55aaa7ad84018172308cbc1d43c4bd0fba83db0040bc7fb672bae",
             "8385d0928cfbd8524398509f89b3b6ce23e41e29251db1f7281b9743284bc45c",
             "1cd583848ca9470ae2d233f8e7d97030d5b9fad78980f2263af537ec4ab16dee"),
  source = c("gdsc", "depmap", "depmap"),
  stringsAsFactors = FALSE
)

# ---- helpers ----------------------------------------------------------------

verify_file <- function(path, size, sha256) {
  if (!file.exists(path)) return(list(status = "missing", detail = "file not found"))
  actual_size <- file.info(path)$size
  if (actual_size != size) {
    return(list(status = "mismatch",
                detail = sprintf("size %d bytes, expected %d bytes", actual_size, size)))
  }
  actual_sha <- digest(file = path, algo = "sha256")
  if (!identical(tolower(actual_sha), tolower(sha256))) {
    return(list(status = "mismatch",
                detail = sprintf("SHA-256 %s, expected %s", actual_sha, sha256)))
  }
  list(status = "ok", detail = sprintf("%d bytes, SHA-256 verified", actual_size))
}

download_to <- function(url, path) {
  old <- options(timeout = max(7200, getOption("timeout")))
  on.exit(options(old), add = TRUE)
  tmp <- paste0(path, ".part")
  message(sprintf("Downloading %s", url))
  status <- tryCatch(download.file(url, tmp, mode = "wb", quiet = FALSE),
                     error = function(e) { message(conditionMessage(e)); 1L })
  if (!identical(as.integer(status), 0L) || !file.exists(tmp)) {
    if (file.exists(tmp)) unlink(tmp)
    return(FALSE)
  }
  file.rename(tmp, path)
}

figshare_download_urls <- function(article_id, wanted) {
  record <- tryCatch(fromJSON(paste0(FIGSHARE_API, article_id)),
                     error = function(e) NULL)
  if (is.null(record) || is.null(record$files)) {
    message(sprintf("Could not read figshare record %s.", article_id))
    return(NULL)
  }
  available <- record$files[record$files$name %in% wanted, c("name", "download_url"), drop = FALSE]
  missing <- setdiff(wanted, available$name)
  if (length(missing) > 0) {
    message(sprintf("figshare record %s (%s) does not contain: %s",
                    article_id, record$title, paste(missing, collapse = ", ")))
  }
  available
}

# ---- arguments --------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
verify_only <- length(args) >= 1 && args[1] == "--verify-only"
if (verify_only && (length(args) < 2 || !dir.exists(args[2]))) {
  stop("Usage: Rscript scripts/00_download_data.R --verify-only <existing directory>", call. = FALSE)
}
raw_dir <- if (verify_only) args[2] else file.path("data", "raw")
if (!verify_only) dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

# ---- download ---------------------------------------------------------------

if (!verify_only) {
  for (i in seq_len(nrow(files))) {
    f <- files[i, ]
    path <- file.path(raw_dir, f$name)
    if (verify_file(path, f$size, f$sha256)$status == "ok") {
      message(sprintf("%s already present and verified; skipping download.", f$name))
      next
    }
    if (f$source == "gdsc") {
      if (!isTRUE(download_to(paste0(GDSC_BASE_URL, f$name), path))) {
        message(sprintf("Download of %s failed; check the network and the GDSC server.", f$name))
      }
    }
  }
  depmap <- files[files$source == "depmap", ]
  pending <- depmap$name[vapply(seq_len(nrow(depmap)), function(i) {
    verify_file(file.path(raw_dir, depmap$name[i]), depmap$size[i], depmap$sha256[i])$status != "ok"
  }, logical(1))]
  if (length(pending) > 0) {
    if (nzchar(DEPMAP_FIGSHARE_ARTICLE)) {
      urls <- figshare_download_urls(DEPMAP_FIGSHARE_ARTICLE, pending)
      if (!is.null(urls)) {
        for (j in seq_len(nrow(urls))) {
          if (!isTRUE(download_to(urls$download_url[j], file.path(raw_dir, urls$name[j])))) {
            message(sprintf("Download of %s failed; check the network and the figshare record.", urls$name[j]))
          }
        }
      }
    } else {
      message("")
      message("DepMap Public 19Q1 files are not downloaded automatically: no figshare")
      message("record id is configured (DEPMAP_FIGSHARE_ARTICLE). Obtain them from the")
      message("DepMap portal (release DepMap Public 19Q1) or through the Bioconductor")
      message("package depmap (ExperimentHub), place them in data/raw/ and re-run this")
      message("script to verify their size and SHA-256. See data/README.md.")
      message("")
    }
  }
}

# ---- verify -----------------------------------------------------------------

message(sprintf("Verifying files in %s", normalizePath(raw_dir, mustWork = FALSE)))
results <- lapply(seq_len(nrow(files)), function(i) {
  f <- files[i, ]
  verify_file(file.path(raw_dir, f$name), f$size, f$sha256)
})
for (i in seq_len(nrow(files))) {
  message(sprintf("  [%s] %s: %s", toupper(results[[i]]$status), files$name[i], results[[i]]$detail))
}
statuses <- vapply(results, function(r) r$status, character(1))
if (any(statuses == "mismatch")) {
  stop("One or more files do not match the expected size or SHA-256. ",
       "Do not use them: delete the mismatching files and re-run the script. ",
       "If the mismatch persists, the file on the server has changed and the ",
       "expected values in this script need to be reviewed.", call. = FALSE)
}
if (any(statuses == "missing")) {
  message("Some files are missing. The legacy scripts cannot run until all three files are present. The two DepMap files can be rebuilt from Bioconductor with: Rscript scripts/01_depmap_from_bioconductor.R --install (see data/README.md).")
  quit(status = 2)
}
message("All input files are present and verified.")
