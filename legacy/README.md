# Legacy scripts (MSc thesis annex, June 2019)

These nine scripts are the code deposited with the thesis "Predicción de
respuesta a fármacos quimioterapéuticos a partir de datos genómicos" (MSc in
Bioinformatics and Biostatistics, UOC-UB, June 2019), available at
<https://hdl.handle.net/10609/97486>. They are published **unchanged except for
text encoding and a provenance header**: each file was converted from its
original encoding (ISO-8859-1, or plain ASCII for `creacionmatriz.R`) to
UTF-8, keeps its CRLF line endings, and starts with a ten-line comment header
that records the original file name, date, encoding and SHA-256. Removing the
header and converting back to the original encoding gives a byte-identical
copy of the deposited file (verified for all nine files).

Line numbers quoted in `docs/errata.md` refer to these files; subtract 10 to
obtain the line in the deposited original.

## Scripts and thesis sections

| Script | What it does | Thesis section (pages of the deposited PDF) |
|---|---|---|
| `creacionmatriz.R` | Builds the cell line by gene matrix with the AUC column of the chosen drug, joining GDSC and DepMap by COSMIC_ID; writes `farmaco.rds`, `expresion.rds` and `finalmatrix.rds` | Objective 1; Methods, pages 9-10 |
| `h2ocontinuosRF.R` | Random forest regression with hyperparameter search in h2o | Methods, pages 13-16; Results, pages 22-24 and 40-42 |
| `h2odiscretizadosRF.R` | Random forest classification with h2o | Results, pages 21-22 and 38-40 |
| `randomsearchcontinuosRF.R` | Random forest regression tuned with caret | Results, pages 26-27 and 44-46 |
| `randomsearchdiscretizadosRF.R` | Random forest classification tuned with caret | Results, pages 24-26 and 42-44 |
| `ridgelassoelasticcontinuosrandomsearch.R` | Ridge, lasso and elastic net regression | Methods, pages 16-18; Results, pages 28-34 and 48-58 |
| `ridgelassoelasticdiscretizadosrandomsearch.R` | Penalized logistic regression (ridge, lasso, elastic net) | Results, pages 27-33 and 46-56 |
| `SVMLinealcontinuos.R` | Linear SVM regression | Methods, pages 18-19; Results, pages 36-38 and 60-62 |
| `SVMLinealdiscretizados.R` | Linear SVM classification | Results, pages 35-36 and 58-60 |

Every model script works on one drug at a time: the drug name is set near the
top of `creacionmatriz.R`, and the hyperparameters chosen in 2019 are written
by hand in the model scripts. Regression scripts use the AUC as a continuous
response; classification scripts discretise it at the third quartile. Each
script fits the models twice, with an 80/20 split and with a 60/20/20 split
(hyperparameters tuned on the 20 percent validation subset, final fit on the
remaining 80 percent, test on the last 20 percent). Cross-validation is used
inside the training set where the package supports it: 10-fold repeated five
times in caret, 10-fold in glmnet and in the SVM tuning; the h2o searches use
training-frame metrics.

## How the scripts expect to be run

1. Obtain the input files (see `../data/README.md`).
2. Start R at the project root, so that `renv` is activated by `.Rprofile`
   (open the `.Rproj` file or run R from that directory), then move to the
   directory that holds the three input files and source the scripts:
   `setwd("data/raw"); source("../../legacy/creacionmatriz.R")` first, then
   any of the model scripts. They read and write by relative file name in the
   working directory (`.rds` matrices and `.RData` fitted models), which git
   ignores.
3. Figures and console summaries are produced interactively; the scripts do
   not export any figure or table, only the `.rds` and `.RData` files above.

## Original environment

The thesis reports R 3.5.1 on Ubuntu 18.04 with RStudio and the following
package versions: ROCR 1.0.7, e1071 1.7.1, glmnet 2.0.16, caret 6.0,
miscTools 0.6.22, randomForest 4.6, h2o 3.22.1, readxl 1.3.0 and
data.table 1.12.0. The `renv.lock` of this repository records a current
environment (R 4.3) with newer versions of the same packages, so numeric
results may differ from the thesis.

**Random partitions.** The 80/20 split uses `caTools::sample.split()`, which
is not affected by the R 3.6.0 change of the default `sample()` algorithm and
is reproduced with the same seed. The 60/20/20 split, `createDataPartition()`
and the internal folds of `cv.glmnet`, caret and `tune` use `sample()` and
depend both on `sample.kind` and on the random numbers consumed by the fits
run before them; the h2o searches have no seed and a time budget. Calling
`RNGkind(sample.kind = "Rounding")` before the scripts is therefore
necessary but not sufficient to reproduce the 2019 numbers (item 8 of
`../docs/errata.md`).

**h2o needs Java.** `h2ocontinuosRF.R` and `h2odiscretizadosRF.R` start a
local h2o server, which requires a Java runtime supported by the installed
h2o version. Java is a system dependency and is not managed by renv.
