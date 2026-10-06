# Drug response prediction in cancer cell lines (GDSC v17.3 and DepMap 19Q1)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](LICENSE)
![R 4.3](https://img.shields.io/badge/R-4.3-276DC3?style=flat-square&logo=r&logoColor=white)
![renv](https://img.shields.io/badge/renv-lockfile-1F6FB2?style=flat-square)
![h2o](https://img.shields.io/badge/h2o-random%20forest-6F42C1?style=flat-square)
![caret](https://img.shields.io/badge/caret-model%20tuning-2C7FB8?style=flat-square)
![glmnet](https://img.shields.io/badge/glmnet-penalized%20regression-2E8B57?style=flat-square)

**Predicting the drug response of cancer cell lines from gene expression with random forests, penalized regression and linear SVMs in R: the 2019 thesis code as deposited, plus a verified data download, a renv environment and an errata list.**

R code of a Master's thesis that predicts the response of cancer cell lines
to a drug from their gene expression: the area under the dose-response curve
(AUC) of GDSC release 7.0 (v17.3) as response, the RNA-seq expression of
DepMap Public 19Q1 as predictors. The nine scripts are published unchanged
except for text encoding and a provenance header; this release (v0.1.0) adds
a reproducible data download, a `renv` environment and a list of errata of
the 2019 code.

**Disclaimer:** this code is intended for research and teaching only. It is
not a clinical tool and must not be used to guide the treatment of any patient.

Master's thesis: MSc in Bioinformatics and Biostatistics, Universitat
Oberta de Catalunya and Universitat de Barcelona (UOC-UB), June 2019. Author:
Alfonso Esteban Lasso.

> The scripts in `legacy/` are published as deposited in June 2019 and keep
> their original Spanish comments; everything written for this release (this
> README, `legacy/README.md`, `data/README.md`, `docs/errata.md` and
> `scripts/00_download_data.R`) is in English.

## Objective

The code predicts the response of cancer cell lines to a drug, measured as the
area under the dose-response curve (AUC) of GDSC release 7.0 (v17.3), from the
RNA-seq expression of DepMap Public 19Q1. For each drug it builds a cell line
by gene matrix (joined by COSMIC_ID), keeps the 200 most variable genes, and
fits random forests (h2o and caret), ridge, lasso and elastic net (glmnet) and
linear SVM (e1071), both as regression on the AUC and as classification after
discretising the AUC at its third quartile, with an 80/20 split and a 60/20/20
split (cross-validation inside the training set where the package supports
it: 10-fold repeated five times in caret, 10-fold in glmnet and in the SVM
tuning; the h2o searches use training-frame metrics). The thesis analysed four
drugs: erlotinib, rapamycin, sunitinib and paclitaxel.

## Data & methods

This repository does **not** redistribute any data. `scripts/00_download_data.R`
downloads the GDSC file and verifies all three inputs by size and SHA-256;
files placed in `data/raw/` by hand are verified in the same way, and
everything under `data/` except its README is ignored by git.

- **GDSC release 7.0** (Wellcome Sanger Institute and Massachusetts General
  Hospital): `v17.3_fitted_dose_response.xlsx`, downloaded from the
  [GDSC download server](https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-7.0/);
  non-exclusive licence for internal research and education, no resale.
- **DepMap Public 19Q1** (Broad Institute): `CCLE_depMap_19Q1_TPM.csv` and
  `DepMap-2019q1-celllines.csv`, CC BY 4.0. These two files are not in the
  figshare record of the Achilles 19Q1 release; `data/README.md` explains how
  to obtain them (the [DepMap portal](https://depmap.org) or the Bioconductor
  `depmap` package) and the script verifies them once placed in `data/raw/`.

Exit codes of the script: 0 when all three files are present and verified, 1
when a file does not match its expected size or SHA-256 (the script stops with
an error), 2 when a file is missing. The mode
`Rscript scripts/00_download_data.R --verify-only <directory>` only checks
files already present in a directory, without downloading anything.

Methods: the nine legacy scripts implement the pipeline summarised under
Objective, one drug at a time. `creacionmatriz.R` builds the cell line by gene
matrix with the AUC column of the chosen drug, and each of the eight model
scripts fits its models twice, with the 80/20 split and with the 60/20/20
split (hyperparameters tuned on the 20 percent validation subset, final fit
on the remaining 80 percent, test on the last 20 percent). `legacy/README.md`
maps each script to the thesis sections.

## Known limitations and errata

The 2019 code has seven methodological defects and one reproducibility
caveat, documented with file and line in the table of
[`docs/errata.md`](docs/errata.md) (titled "Errata of the 2019 code") and
left uncorrected in v0.1.0:

1. The R-squared is computed as the squared correlation between the observed
   value and the residual, which is not the coefficient of determination; the
   R-squared equal to 1 highlighted in the thesis is an artefact of that
   calculation.
2. `randomForest()` receives tuned arguments that it does not understand and
   silently ignores.
3. The 80/20 "elastic net" calls `cv.glmnet()` without `alpha` and is a lasso.
4. `caret::train()` performs a grid search, not the random search described.
5. The class label probably assigns "sensitive" to the most resistant quartile.
6. Gene selection and the class threshold are computed before the train/test
   split (information leakage).
7. The SVM regression uses 2000 genes instead of 200 and saves the wrong
   object; in the penalized logistic script three `predict()` calls have
   misplaced parentheses (without effect on the ROC curves).
8. The 2019 partitions and resampling are not reproduced exactly by a current
   R: the 60/20/20 split and the internal folds depend on `sample()` (changed
   in R 3.6.0) and on the random numbers consumed by earlier fits, and the
   h2o searches have no seed and a time budget.

The thesis results should be read with these limitations in mind.

`docs/errata.md` also lists further observations of the 2026 review (section
"Additional observations"), for example the sorting of the h2o grid in the
classification script and the sign of the SVM decision values used for its
ROC curves.

## Tech stack

R 4.3 with a `renv` lockfile (`renv.lock`, `.Rprofile`), caret, glmnet,
randomForest, e1071, h2o (needs a Java runtime), caTools, ROCR, miscTools,
readxl, data.table, plyr, resample, ggplot2; digest and jsonlite for the
download script.

## How to run

1. Install R 4.3 or later, Rtools on Windows, and a Java runtime supported by
   the installed h2o version (needed only by the two h2o scripts; the
   environment builds without it).

2. Restore the package environment:

   ```bash
   Rscript -e "renv::restore()"
   ```

3. Download and verify the data:

   ```bash
   Rscript scripts/00_download_data.R
   ```

   The DepMap files must be placed in `data/raw/` by hand (the script says so
   and verifies them). See the exit codes under Data & methods.

4. Run the legacy scripts from `data/raw/`, which holds the three input
   files, in an R session started at the project root so that `renv` is active
   (open `gdsc-depmap-drug-response-r.Rproj` or run R there). Source
   `creacionmatriz.R` first, which builds the matrix, and then any model
   script in the same way:

   ```r
   setwd("data/raw")
   source("../../legacy/creacionmatriz.R")
   source("../../legacy/h2ocontinuosRF.R")
   ```

Note: the scripts read and write by relative file name in the working
directory (`.rds` matrices and `.RData` fitted models, ignored by git).
Figures and console summaries are produced interactively; the scripts do not
export any figure or table. The 2019 random partitions are not reproduced
exactly by a current R: see item 8 of `docs/errata.md` and `legacy/README.md`.

## Repository structure

```
gdsc-depmap-drug-response-r/
├── legacy/                                             # the nine scripts deposited with the thesis (UTF-8, provenance header)
│   ├── README.md                                       # script-to-thesis-section map, how to run them, original environment
│   ├── creacionmatriz.R                                # cell line by gene matrix joined by COSMIC_ID (run first)
│   ├── h2ocontinuosRF.R                                # random forest regression, h2o search (needs Java)
│   ├── h2odiscretizadosRF.R                            # random forest classification, h2o search (needs Java)
│   ├── randomsearchcontinuosRF.R                       # random forest regression tuned with caret
│   ├── randomsearchdiscretizadosRF.R                   # random forest classification tuned with caret
│   ├── ridgelassoelasticcontinuosrandomsearch.R        # ridge, lasso and elastic net regression (glmnet)
│   ├── ridgelassoelasticdiscretizadosrandomsearch.R    # penalized logistic regression (glmnet)
│   ├── SVMLinealcontinuos.R                            # linear SVM regression (e1071)
│   └── SVMLinealdiscretizados.R                        # linear SVM classification (e1071)
├── scripts/
│   └── 00_download_data.R                              # downloads the GDSC file, verifies size and SHA-256 of the three inputs
├── data/
│   ├── README.md                                       # sources, licences, expected sizes and hashes
│   └── raw/                                            # input files (created locally, ignored by git)
├── docs/
│   └── errata.md                                       # eight items of the 2019 code and additional observations
├── renv/                                               # activate.R, settings.json and its own .gitignore (the project library is not versioned)
│   ├── activate.R
│   └── settings.json
├── renv.lock                                           # package versions of the reproducible environment (R 4.3)
├── .Rprofile                                           # activates renv when R starts at the project root
├── .gitignore                                          # ignores data/ (except its README), *.rds, *.RData, *.csv, *.xlsx and renv/library/
├── gdsc-depmap-drug-response-r.Rproj                   # RStudio project file (open it to start R at the root)
├── CITATION.cff                                        # citation metadata for the code and the thesis
├── LICENSE                                             # MIT, the author's code only
└── README.md
```

## Provenance note

Supervised at the Bioinformatics Unit of the Spanish National Cancer Research
Centre (CNIO) by Dr. Héctor Tejero Franco, with Jeroni Luna (UOC) as academic
tutor. The thesis (in Spanish) is deposited at
<https://hdl.handle.net/10609/97486> under a CC BY-NC-ND 3.0 ES licence; it is
linked here, not copied.

The scripts in `legacy/` are published for transparency and traceability, as
deposited in June 2019, with Spanish comments and the hyperparameters chosen
at the time written by hand. No figures or tables are exported; the only
files they write are the `.rds` matrices and the `.RData` fitted models.
`legacy/README.md` maps each script to the thesis sections and describes the
original environment (R 3.5.1 on Ubuntu 18.04 and the package versions as
listed in the thesis).

What this release is, and what it is **not**:

- **v0.1.0 (this release):** the code deposited with the thesis, unchanged
  except for text encoding and a provenance header, plus a reproducible data
  download, a `renv` environment and a list of errata of the 2019 code. It is
  **not** a corrected version: the defects listed in the errata are left as
  they are.
- **Roadmap, v1.0:** a corrected reimplementation (functions instead of
  repeated blocks, correct metrics, gene selection and thresholds inside the
  training folds, exported figures and tables), kept separate from the legacy
  code.

## License

The author's code is released under the MIT License (see `LICENSE`). The data
remain under the terms of their providers and are not part of this
repository. The thesis is linked under its own CC BY-NC-ND 3.0 ES licence and
is not redistributed here.

## How to cite

Use the metadata in `CITATION.cff` (GitHub shows a "Cite this repository"
button) and cite the thesis:

Esteban Lasso, A. (2019). *Predicción de respuesta a fármacos
quimioterapéuticos a partir de datos genómicos*. Master's thesis, Universitat
Oberta de Catalunya and Universitat de Barcelona.
<https://hdl.handle.net/10609/97486>

## Acknowledgements

Model fitting follows the public documentation and tutorials of the time for
caret, glmnet, randomForest, e1071 and h2o. Thanks to the GDSC and DepMap
projects for making their data publicly available.
