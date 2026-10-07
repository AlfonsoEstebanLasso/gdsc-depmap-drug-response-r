# ----------------------------------------------------------------------------
# R/build_matrix.R: per-drug cell line by gene matrix (v1.0, data builder).
#
# Purpose: reproduce the join of legacy/creacionmatriz.R (GDSC v17.3 AUC and
# DepMap 19Q1 TPM joined through COSMIC_ID and DepMap_ID) as plain functions,
# one drug at a time, WITHOUT any gene selection or class threshold on the
# full data. Gene selection and thresholds move inside each training fold
# (see R/cv.R and docs/v1_spec.md, sections 2.2 and 5).
#
# Inputs (data/raw, verified by scripts/00_download_data.R):
#   v17.3_fitted_dose_response.xlsx  DRUG_NAME, CELL_LINE_NAME, COSMIC_ID, DRUG_ID, AUC
#   DepMap-2019q1-celllines.csv      DepMap_ID, COSMIC_ID
#   CCLE_depMap_19Q1_TPM.csv         DepMap_ID (unnamed first column) and 57820 genes
#
# Output of build_matrix() / load_drug_matrix(): a list with
#   drug        lower-case key used in paths and cache names
#   drug_name   DRUG_NAME as spelled in the GDSC file
#   drug_id     the single DRUG_ID of that drug
#   cell_lines  CELL_LINE_NAME, sorted (C locale), one per row of x
#   x           numeric matrix, cell lines by 57820 genes, log2(TPM + 1), file order
#   auc         numeric vector aligned with the rows of x, named by cell line
#   depmap_ids  DepMap_ID of every row, aligned
#   cosmic_ids  COSMIC_ID of every row, aligned
#
# Transformations: none on the values. The AUC is used as it is (a fraction
# in (0, 1], higher means more resistant) and the expression values are the
# log2(TPM + 1) of the file. Rows are sorted by cell line name. Any missing
# value in x or auc stops with an error.
#
# Cache: data/derived/<drug>_matrix.rds (ignored by git through data/*).
# ----------------------------------------------------------------------------

# Expected cell line counts after the join, for the four thesis drugs
# (docs/v1_spec.md, section 2.2). Used only for an informative message.
.expected_cell_lines <- c(erlotinib = 206, rapamycin = 208,
                          sunitinib = 229, paclitaxel = 230)

.raw_file_names <- c(
  gdsc     = "v17.3_fitted_dose_response.xlsx",
  metadata = "DepMap-2019q1-celllines.csv",
  tpm      = "CCLE_depMap_19Q1_TPM.csv"
)

# Stop with a helpful message when a raw file is missing.
.check_raw_file <- function(path) {
  if (!file.exists(path)) {
    stop("Raw file not found: '", basename(path), "'. Place the three input ",
         "files in data/raw and verify them with ",
         "'Rscript scripts/00_download_data.R' (see data/README.md).",
         call. = FALSE)
  }
  invisible(path)
}

# GDSC fitted dose-response table: the columns used by v1.0 only.
read_gdsc <- function(path) {
  .check_raw_file(path)
  g <- readxl::read_excel(path, sheet = 1)
  needed <- c("DRUG_NAME", "CELL_LINE_NAME", "COSMIC_ID", "DRUG_ID", "AUC")
  missing_cols <- setdiff(needed, names(g))
  if (length(missing_cols) > 0) {
    stop("GDSC file lacks the columns: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }
  data.frame(
    DRUG_NAME      = as.character(g$DRUG_NAME),
    CELL_LINE_NAME = as.character(g$CELL_LINE_NAME),
    COSMIC_ID      = as.integer(g$COSMIC_ID),
    DRUG_ID        = as.integer(g$DRUG_ID),
    AUC            = as.numeric(g$AUC),
    stringsAsFactors = FALSE
  )
}

# DepMap cell line metadata: DepMap_ID and COSMIC_ID, rows without a
# COSMIC_ID dropped (they can never match a GDSC line).
read_depmap_metadata <- function(path) {
  .check_raw_file(path)
  m <- data.table::fread(path, select = c("DepMap_ID", "COSMIC_ID"),
                         data.table = FALSE, showProgress = FALSE)
  m <- data.frame(
    DepMap_ID = as.character(m$DepMap_ID),
    COSMIC_ID = suppressWarnings(as.integer(m$COSMIC_ID)),
    stringsAsFactors = FALSE
  )
  m[!is.na(m$COSMIC_ID) & !is.na(m$DepMap_ID) & nzchar(m$DepMap_ID), ,
    drop = FALSE]
}

# DepMap expression file: read once per run with data.table::fread(), name
# the first column DepMap_ID, keep the rows whose DepMap_ID is in keep_ids
# (all rows when NULL) and convert the gene columns to a numeric matrix.
read_depmap_tpm <- function(path, keep_ids = NULL) {
  .check_raw_file(path)
  dt <- data.table::fread(path, header = TRUE, showProgress = FALSE)
  data.table::setnames(dt, 1L, "DepMap_ID")
  ids <- as.character(dt[[1L]])
  keep <- if (is.null(keep_ids)) seq_along(ids) else which(ids %in% keep_ids)
  if (length(keep) == 0) {
    stop("No row of the expression file matches the requested DepMap_ID values.",
         call. = FALSE)
  }
  dt <- dt[keep]
  ids <- ids[keep]
  genes <- names(dt)[-1L]
  x <- as.matrix(dt[, -1L, with = FALSE])
  storage.mode(x) <- "double"
  dimnames(x) <- list(ids, genes)
  rm(dt)
  list(ids = ids, x = x)
}

# DRUG_NAME as spelled in GDSC; case-insensitive; error if absent or if the
# name maps to more than one DRUG_ID.
resolve_drug <- function(name, gdsc) {
  if (length(name) != 1 || is.na(name) || !nzchar(trimws(name))) {
    stop("'name' must be one non-empty drug name.", call. = FALSE)
  }
  key <- tolower(trimws(name))
  hit <- gdsc[tolower(gdsc$DRUG_NAME) == key, c("DRUG_NAME", "DRUG_ID"),
              drop = FALSE]
  if (nrow(hit) == 0) {
    stop("Drug '", name, "' is not present in the GDSC file.", call. = FALSE)
  }
  ids <- sort(unique(hit$DRUG_ID))
  if (length(ids) > 1) {
    stop("Drug '", name, "' maps to more than one DRUG_ID (",
         paste(ids, collapse = ", "), ").", call. = FALSE)
  }
  spellings <- unique(hit$DRUG_NAME)
  if (length(spellings) > 1) {
    # Identical names up to case: keep the most frequent spelling.
    spellings <- names(sort(table(hit$DRUG_NAME), decreasing = TRUE))[1L]
  }
  spellings[1L]
}

# Steps 1 to 3 of docs/v1_spec.md section 2.2: the rows of one drug joined
# with the metadata on COSMIC_ID. Returns a data.frame with CELL_LINE_NAME,
# COSMIC_ID, DepMap_ID, AUC (one row per cell line) and the attributes
# drug_name and drug_id.
match_drug_lines <- function(drug, gdsc, metadata) {
  drug_name <- resolve_drug(drug, gdsc)
  rows <- gdsc[gdsc$DRUG_NAME == drug_name, , drop = FALSE]
  drug_id <- unique(rows$DRUG_ID)[1L]
  rows <- rows[!is.na(rows$COSMIC_ID), , drop = FALSE]
  dup <- unique(rows$COSMIC_ID[duplicated(rows$COSMIC_ID)])
  if (length(dup) > 0) {
    stop("Drug '", drug_name, "' has more than one row per cell line in GDSC ",
         "(COSMIC_ID: ", paste(dup, collapse = ", "), ").", call. = FALSE)
  }
  joined <- merge(
    rows[, c("CELL_LINE_NAME", "COSMIC_ID", "AUC")],
    metadata[, c("DepMap_ID", "COSMIC_ID")],
    by = "COSMIC_ID", all = FALSE
  )
  dup <- unique(joined$COSMIC_ID[duplicated(joined$COSMIC_ID)])
  if (length(dup) > 0) {
    stop("COSMIC_ID maps to several DepMap_ID for drug '", drug_name, "': ",
         paste(dup, collapse = ", "), ".", call. = FALSE)
  }
  dup <- unique(joined$DepMap_ID[duplicated(joined$DepMap_ID)])
  if (length(dup) > 0) {
    stop("DepMap_ID maps to several cell lines for drug '", drug_name, "': ",
         paste(dup, collapse = ", "), ".", call. = FALSE)
  }
  out <- data.frame(
    CELL_LINE_NAME = as.character(joined$CELL_LINE_NAME),
    COSMIC_ID      = as.integer(joined$COSMIC_ID),
    DepMap_ID      = as.character(joined$DepMap_ID),
    AUC            = as.numeric(joined$AUC),
    stringsAsFactors = FALSE
  )
  attr(out, "drug_name") <- drug_name
  attr(out, "drug_id") <- drug_id
  out
}

# DepMap_ID values needed by one or several drugs: the keep_ids argument of
# read_depmap_tpm() when a run covers several drugs.
drug_depmap_ids <- function(drugs, gdsc, metadata) {
  ids <- unlist(lapply(drugs, function(d) {
    match_drug_lines(d, gdsc, metadata)$DepMap_ID
  }), use.names = FALSE)
  sort(unique(ids))
}

# Steps 4 to 6: inner join with the expression rows, numeric matrix with one
# row per cell line (row names CELL_LINE_NAME, sorted by name), AUC aligned.
# tpm is the list returned by read_depmap_tpm() (or a numeric matrix whose
# row names are DepMap_ID values).
build_matrix <- function(drug, gdsc, metadata, tpm) {
  if (is.matrix(tpm)) {
    if (is.null(rownames(tpm))) {
      stop("'tpm' given as a matrix must have DepMap_ID row names.",
           call. = FALSE)
    }
    tpm <- list(ids = rownames(tpm), x = tpm)
  }
  if (!is.list(tpm) || is.null(tpm$ids) || is.null(tpm$x)) {
    stop("'tpm' must be the list returned by read_depmap_tpm().", call. = FALSE)
  }
  if (anyDuplicated(tpm$ids)) {
    stop("The expression table has duplicated DepMap_ID values.", call. = FALSE)
  }
  lines <- match_drug_lines(drug, gdsc, metadata)
  drug_name <- attr(lines, "drug_name")
  drug_id <- attr(lines, "drug_id")

  pos <- match(lines$DepMap_ID, tpm$ids)
  lines <- lines[!is.na(pos), , drop = FALSE]
  pos <- pos[!is.na(pos)]
  if (nrow(lines) == 0) {
    stop("No cell line of drug '", drug_name, "' has expression data.",
         call. = FALSE)
  }
  if (anyDuplicated(lines$CELL_LINE_NAME)) {
    stop("Duplicated CELL_LINE_NAME after the join for drug '", drug_name, "'.",
         call. = FALSE)
  }

  # Sort by cell line name in the C locale so the row order does not depend
  # on the session locale.
  ord <- order(lines$CELL_LINE_NAME, method = "radix")
  lines <- lines[ord, , drop = FALSE]
  pos <- pos[ord]

  x <- tpm$x[pos, , drop = FALSE]
  storage.mode(x) <- "double"
  rownames(x) <- lines$CELL_LINE_NAME
  auc <- lines$AUC
  names(auc) <- lines$CELL_LINE_NAME

  if (anyNA(x)) {
    stop("Missing values in the expression matrix of drug '", drug_name, "'.",
         call. = FALSE)
  }
  if (anyNA(auc)) {
    stop("Missing AUC values for drug '", drug_name, "'.", call. = FALSE)
  }

  key <- tolower(drug_name)
  if (key %in% names(.expected_cell_lines) &&
      nrow(x) != .expected_cell_lines[[key]]) {
    message("Note: ", nrow(x), " cell lines for ", key, ", expected ",
            .expected_cell_lines[[key]], " with the verified raw files.")
  }

  list(
    drug       = key,
    drug_name  = drug_name,
    drug_id    = drug_id,
    cell_lines = lines$CELL_LINE_NAME,
    x          = x,
    auc        = auc,
    depmap_ids = lines$DepMap_ID,
    cosmic_ids = lines$COSMIC_ID
  )
}

# Build the matrix of one drug from data/raw, or read it from the cache.
# tpm: the list of read_depmap_tpm() when the expression file was already
# read for this run (several drugs); when NULL, only the rows needed by this
# drug are kept after reading the file.
load_drug_matrix <- function(drug, raw_dir = "data/raw",
                             cache_dir = "data/derived", refresh = FALSE,
                             tpm = NULL) {
  key <- tolower(trimws(drug))
  cache_path <- file.path(cache_dir, paste0(key, "_matrix.rds"))
  if (!refresh && file.exists(cache_path)) {
    out <- readRDS(cache_path)
    if (is.list(out) && is.matrix(out$x) && length(out$auc) == nrow(out$x)) {
      return(out)
    }
    warning("Cache '", cache_path, "' is not a valid matrix object; rebuilding.",
            call. = FALSE)
  }

  paths <- file.path(raw_dir, .raw_file_names)
  names(paths) <- names(.raw_file_names)
  for (p in paths) .check_raw_file(p)

  gdsc <- read_gdsc(paths[["gdsc"]])
  metadata <- read_depmap_metadata(paths[["metadata"]])
  if (is.null(tpm)) {
    keep_ids <- drug_depmap_ids(drug, gdsc, metadata)
    tpm <- read_depmap_tpm(paths[["tpm"]], keep_ids = keep_ids)
  }
  out <- build_matrix(drug, gdsc, metadata, tpm)

  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  saveRDS(out, cache_path)
  out
}
