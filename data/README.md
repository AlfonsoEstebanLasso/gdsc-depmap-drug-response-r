# Input data

The legacy scripts need three public files. **This repository does not
redistribute any data.** `scripts/00_download_data.R` downloads what it can,
verifies size and SHA-256, and aborts if a file does not match. Files placed
in `data/raw/` by hand are verified in the same way. The whole `data/`
directory is ignored by git.

| File | Expected size (bytes) | SHA-256 | Source |
|---|---|---|---|
| `v17.3_fitted_dose_response.xlsx` | 16064522 | `3c48e4dae5d55aaa7ad84018172308cbc1d43c4bd0fba83db0040bc7fb672bae` | GDSC release 7.0 |
| `CCLE_depMap_19Q1_TPM.csv` | 568839199 | `8385d0928cfbd8524398509f89b3b6ce23e41e29251db1f7281b9743284bc45c` | DepMap Public 19Q1 |
| `DepMap-2019q1-celllines.csv` | 154145 | `1cd583848ca9470ae2d233f8e7d97030d5b9fad78980f2263af537ec4ab16dee` | DepMap Public 19Q1 |

## GDSC (Genomics of Drug Sensitivity in Cancer)

- File: `v17.3_fitted_dose_response.xlsx`, fitted dose-response parameters
  (including the AUC used by the thesis) of release 7.0 (March 2018).
- Download server: `https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-7.0/`
  (downloaded automatically by the script; checked on 2026-10-06).
- Terms of use, as published by the GDSC project: a non-exclusive,
  non-transferable licence for internal research and educational purposes
  only; the data must not be sold or redistributed. Check the current terms
  on the GDSC website before use and cite the GDSC resource in any work that
  uses the data.

## DepMap Public 19Q1 (Broad Institute)

- Files: `CCLE_depMap_19Q1_TPM.csv` (RNA-seq expression, log2(TPM + 1)) and
  `DepMap-2019q1-celllines.csv` (cell line metadata, including the mapping
  between DepMap_ID and COSMIC_ID used by `legacy/creacionmatriz.R`).
- Licence: DepMap releases are published under the Creative Commons
  Attribution 4.0 (CC BY 4.0) licence. Cite DepMap and the release (DepMap
  Public 19Q1) in any work that uses the data.
- Where to get them. On 2026-10-06 the figshare record 7655150 ("DepMap
  Achilles 19Q1 Public") did **not** contain these two files (it holds the
  Achilles CRISPR screens of the same release), and no other 2019 DepMap
  record returned by the figshare API contained them either. Two routes
  remain:
  1. The DepMap portal (`depmap.org`), release "DepMap Public 19Q1",
     downloads section. Place the two files in `data/raw/` and run
     `Rscript scripts/00_download_data.R` to verify them.
  2. The Bioconductor package `depmap`, which serves the 19Q1 TPM table
     through ExperimentHub (see `depmap::depmap_TPM()` and the package
     documentation for the 19Q1 release). The object delivered this way is an
     R data frame in long format, not the original CSV, so the SHA-256 above
     cannot be checked against it and `legacy/creacionmatriz.R` would need
     the data written back to the expected CSV layout.
  If a figshare record containing the original files is identified later,
  set its id in `DEPMAP_FIGSHARE_ARTICLE` (constant in the script or
  environment variable) and the script will download and verify them.

## Derived files

`legacy/creacionmatriz.R` writes `farmaco.rds`, `expresion.rds` and
`finalmatrix.rds` to the working directory, and the model scripts write
fitted models as `.RData`. All of them are regenerated from the inputs and
are ignored by git.
