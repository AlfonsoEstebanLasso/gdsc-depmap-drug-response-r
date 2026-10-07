# v1.0 specification: corrected reimplementation

This document is the binding specification of v1.0 of this repository. v1.0
is a corrected reimplementation of the 2019 thesis pipeline, kept separate
from `legacy/`, which fixes the eight items of `docs/errata.md`. Every file
written for v1.0 must follow this specification; where the specification
leaves a choice open, the implementer records the choice in the section
"Implementation record" at the end of this file.

Scope reminder: research and teaching only. The models predict the drug
response of cancer cell lines, not of patients, and nothing here is a clinical
tool.

Conventions for everything written for v1.0: English, sentence-case headings,
no em or en dashes, no emojis, no absolute paths, no e-mail addresses, no
third-party names other than the ones already cited in `README.md`. Seeds are
fixed everywhere. R 4.3 with the project `renv`; `Rscript` is run from the
repository root so that `.Rprofile` activates `renv`.

## 1. Mapping of the errata to v1.0

| Erratum | Fix in v1.0 | Where |
|---|---|---|
| 1. R squared computed as cor(y, residual)^2 | R squared is 1 minus SSres over SStot on the held-out fold, with RMSE and MAE | `R/metrics.R` |
| 2. `randomForest()` ignores the tuned h2o arguments | The random forest is fitted with `ranger` using the hyperparameters that were tuned (`mtry`, `min.node.size`, `sample.fraction`, `max.depth`, `num.trees`); h2o is not used | `R/models.R` |
| 3. The 80/20 "elastic net" is a lasso; two lambda rules | `alpha` tuned on a declared grid for both tasks; `lambda.min` everywhere | `R/models.R` |
| 4. Grid search labelled random search; mixed resampling | Explicitly declared grids; the random forest grid is sampled at random with a seed and a fixed budget; the same inner folds for every model family | `R/cv.R`, `R/models.R` |
| 5. "Sensitive" assigned to the highest AUC quartile | Sensitive is AUC below the first quartile of the training fold; documented | `R/cv.R` |
| 6. Gene selection and threshold computed on the whole dataset | Both computed inside each outer training fold; a leakage check compares fold-wise and global selection | `R/cv.R`, `R/run_drug.R` |
| 7. SVM regression on 2000 genes; wrong saved object | The same 200-gene matrix for every model; the fitted objects and tuning tables are saved per fold | `R/models.R`, `R/run_drug.R` |
| 8. Partitions not reproducible | Seeds fixed, partition indices saved to `outputs/<drug>/folds.csv`, no h2o, package versions in `renv.lock` | `run_all.R`, `renv.lock` |

The additional observations of `docs/errata.md` are covered as follows: the
ROC curves use scores oriented towards the sensitive class (section 7.4), the
60/20/20 scheme is replaced by cross-validation (section 4), the ridge model
with a degenerate lambda cannot occur because lambda is chosen by inner
cross-validation on the glmnet path, and hyperparameters are never copied by
hand.

## 2. Data flow

### 2.1 Inputs

The three raw files described in `data/README.md`, placed in `data/raw/` and
verified by `scripts/00_download_data.R` (sizes and SHA-256 unchanged):

| File | Content used by v1.0 |
|---|---|
| `v17.3_fitted_dose_response.xlsx` | Single sheet (`Sheet1`), 224202 rows, 13 columns. Columns used: `DRUG_NAME`, `CELL_LINE_NAME`, `COSMIC_ID`, `AUC`. `DRUG_ID` is read to check that the drug name maps to one identifier |
| `DepMap-2019q1-celllines.csv` | 1677 rows. Columns used: `DepMap_ID`, `COSMIC_ID` (982 non-missing, one duplicated value) |
| `CCLE_depMap_19Q1_TPM.csv` | 1165 rows (cell lines) by 57821 columns: an unnamed first column with the DepMap identifier and 57820 genes named `SYMBOL (ENSG...)`, values log2(TPM + 1) |

### 2.2 Per-drug matrix (reusing the logic of `legacy/creacionmatriz.R`)

The legacy script joins the GDSC table and the DepMap metadata by `COSMIC_ID`,
keeps the unique pairs (`CELL_LINE_NAME`, `DepMap_ID`), joins the expression
rows by `DepMap_ID`, and finally joins the drug's (`CELL_LINE_NAME`, `AUC`)
rows by `CELL_LINE_NAME` with an inner join. v1.0 reproduces the same result
in this order:

1. Read the GDSC sheet with `readxl::read_excel()`; keep the four columns
   above. Match the requested drug against `DRUG_NAME` ignoring case (`--drug
   erlotinib` matches `Erlotinib`). Stop with an error if the name is absent
   or maps to more than one `DRUG_ID` (each of the four thesis drugs maps to
   exactly one: 1, 3, 5 and 11).
2. Keep the rows of that drug: one row per cell line (no duplicated
   `COSMIC_ID` within a drug in this release; stop with an error otherwise).
3. Inner join with the metadata on `COSMIC_ID` (missing `COSMIC_ID` rows are
   dropped by the join). If one `COSMIC_ID` maps to several `DepMap_ID`, stop
   with an error that names the identifier; this does not happen for the four
   thesis drugs.
4. Read the expression file once per run with `data.table::fread()`, name the
   first column `DepMap_ID`, keep only the rows whose `DepMap_ID` appears in
   any requested drug, and convert the gene columns to a numeric matrix.
5. Inner join on `DepMap_ID`. The result is a numeric matrix `x` with one row
   per cell line (row names `CELL_LINE_NAME`, rows sorted by name) and 57820
   gene columns in the file order, plus a numeric vector `auc` aligned with
   the rows. The AUC is used as it is: a fraction in (0, 1], higher means more
   resistant. No log or other transform of the response, as in 2019.
6. Stop with an error if any value of `x` or `auc` is missing (none is
   expected; the 2019 script only counted missing values).

Expected sizes with the verified files (cell lines matched with expression):
erlotinib 206, rapamycin 208, sunitinib 229, paclitaxel 230.

The 200 most variable genes are **not** selected here. Selection moves inside
each outer training fold (section 5.2). The full matrix of each drug is
cached as `data/derived/<drug>_matrix.rds` (ignored by git through the
`data/*` rule); `--refresh-cache` rebuilds it.

## 3. Tasks and class rule

- **Regression**: the response is the AUC.
- **Classification**: the response is `sensitive` when AUC is strictly below
  the first quartile (`quantile(auc_train, 0.25)`, R default type 7) of the
  **outer training fold**, and `resistant` otherwise. The factor levels are
  `c("resistant", "sensitive")`; `sensitive` is the positive class for
  sensitivity, specificity and ROC AUC. The threshold is computed once per
  outer split on its training rows, applied to the training rows and to the
  held-out rows of that split, and reused unchanged by the inner resampling of
  that split (inner folds are subsets of the training rows, so no held-out
  information enters). Expected first quartiles on the full matrices, for
  orientation only: erlotinib 0.955, rapamycin 0.77, sunitinib 0.80,
  paclitaxel 0.56.

This replaces the 2019 rule (AUC at or above the third quartile labelled
sensitive), which is biologically inverted; see erratum 5.

## 4. Cross-validation design

### 4.1 Outer resampling

Repeated stratified k-fold cross-validation, built with
`caret::createMultiFolds()` on a stratification factor:

| Setting | Default | `--fast` |
|---|---|---|
| Outer folds `k` | 5 | 3 |
| Repeats `r` | 3 | 1 |
| Drugs | erlotinib, rapamycin, sunitinib, paclitaxel | erlotinib, paclitaxel |
| Inner folds `k_inner` | 5 | 3 |
| Random forest search budget | 12 configurations | 6 |
| Leakage check (section 8) | on | off |

Stratification factor:

- Regression: the AUC quartile group of each cell line, computed on all cell
  lines of the drug (`cut()` at the three quartiles, `include.lowest = TRUE`).
- Classification: the indicator `auc < quantile(auc, 0.25)` computed on all
  cell lines. This provisional indicator is used **only** to balance the
  folds; the labels used for training and evaluation are always recomputed
  from the training fold (section 3).

Stratifying on the outcome is standard practice and does not leak outcome
values into the predictors. The two tasks therefore have different fold
assignments; both are saved.

### 4.2 Inner resampling for tuning

Inside every outer training fold, one round of stratified `k_inner`-fold
cross-validation produces an integer vector `foldid` (length equal to the
training rows). The **same** `foldid` is passed to every model family of that
outer fold: to `cv.glmnet(foldid = ...)` and to the hand-rolled loops of the
random forest and the SVM. Stratification uses the same rule as the outer
folds, applied to the training rows (quartile groups of the training AUC for
regression; the fold-wise class for classification).

### 4.3 Seeds

One base seed `S` (`--seed`, default 2019). Task index `t` is 1 for
regression and 2 for classification; outer splits are numbered `j = 1..k*r`
in the order returned by `createMultiFolds()` (repeat 1 folds 1..k, then
repeat 2, and so on).

| Use | Seed |
|---|---|
| Outer folds of task `t` | `set.seed(S + t)` before `createMultiFolds()` |
| Inner folds of split `j`, task `t` | `set.seed(S + 1000 * t + j)` |
| Random forest grid draw | `set.seed(S + 1000 * t + j + 100)` |
| `ranger(seed = )` | `S + 1000 * t + j + 1` |
| SVM and glmnet fits | `set.seed(S + 1000 * t + j + 5)` before the first fit of each family (glmnet is deterministic given `foldid`) |
| Bootstrap, if any | `set.seed(S + 7)` |

Thread count does not change the results of `ranger` when `seed` is set; the
thread count is still recorded in the run log.

### 4.4 Saved partitions

`outputs/<drug>/folds.csv` with columns `cell_line`, `task`, `repeat`,
`fold` (the fold in which the cell line is held out). One row per cell line,
task and repeat. This file contains names and fold numbers only, no AUC or
expression values, and is committed.

## 5. Fold-wise preprocessing

Performed by `preprocess_fold()` on each outer split, using training rows
only, then applied to the held-out rows:

1. **Gene selection**: sample variance of every gene over the training rows
   (same definition as `resample::colVars()`, implemented in base R as
   `colSums((x - colMeans(x))^2) / (n - 1)`); keep the 200 genes with the
   largest variance (ties broken by column order). Every model family of
   the split, the SVM included, receives this 200-gene matrix.
2. **Class label** (classification only): rule of section 3.
3. **Standardisation**: centre and scale each selected gene with the training
   mean and standard deviation; genes with zero training variance cannot be
   selected, so no division by zero occurs. The held-out matrix is transformed
   with the training parameters and has the same columns in the same order.
   `glmnet(standardize = FALSE)` and `svm(scale = FALSE)` are used so that the
   data are scaled exactly once; `ranger` is scale invariant.
4. **Response scaling for the SVM regression**: the training AUC is centred
   and scaled (training mean and standard deviation) before `svm()`, and the
   predictions are transformed back; the epsilon grid refers to the scaled
   response. No other model scales the response.

The 200 selected genes of every split are counted in
`outputs/<drug>/<task>/selected_genes.csv` (`gene`, `times_selected`,
`n_splits`): an aggregate, committed.

## 6. Implementation choice: hand-rolled nested loops

v1.0 does **not** use `caret::train()`. Tuning is hand-rolled around the
packages themselves, with the same `foldid` for every family (section 4.2):

- `caret` is used only for `createMultiFolds()`.
- `ranger` for the random forest, `glmnet::cv.glmnet()` for ridge, lasso and
  elastic net, `e1071::svm()` for the linear SVM, `pROC` for ROC AUC,
  `ggplot2` for figures, `data.table` and `readxl` for reading.
- `h2o` is not called by v1.0. The substitution is documented in the README
  and in `docs/comparison_2019_vs_v1.md`: the 2019 h2o grid searches had no
  seed and a time budget, and their tuned values were then dropped by
  `randomForest()`; `ranger` fits the tuned values directly and is
  deterministic given a seed, without a Java dependency.

Reasons: `caret`'s `rf` and `ranger` methods do not expose
`sample.fraction` or `max.depth`, `cv.glmnet` already tunes the lambda path
with a user-supplied `foldid`, and a hand-rolled loop makes the seeds, the
fold reuse and the saved tuning tables explicit and testable.

## 7. Model families and tuning

All grids are declared in `model_grid()` and written to the tuning tables, so
every reported hyperparameter is traceable to a declared candidate. The inner
selection metric is the same for every family: **RMSE** for regression
(lower is better) and **ROC AUC** for classification (higher is better),
computed on the inner validation folds and averaged over the `k_inner` folds.
Ties are broken by the first candidate in grid order.

### 7.1 Random forest (`ranger`), model name `rf`

Declared grid (mapping of the 2019 h2o names in parentheses):

| Parameter | Values | 2019 h2o grid |
|---|---|---|
| `num.trees` (ntrees) | 500 (fixed) | 200, 350, 500 |
| `mtry` (mtries) | 5, 10, 14, 20, 35, 50 | 15, 25, 35 |
| `min.node.size` (min_rows) | 1, 3, 5, 10 | 1, 3, 5 |
| `sample.fraction` (sample_rate) | 0.55, 0.632, 0.75 | 0.55, 0.632, 0.75 |
| `max.depth` (max_depth) | 0 (unlimited), 20, 30, 40 | 20 to 40 by 5 |
| nbins | no equivalent in `ranger`, dropped | 10 to 30 by 5 |

The full grid has 288 configurations. A random search draws `--rf-budget`
of them without replacement (default 12, `--fast` 6) with the seed of section
4.3; the drawn candidates are evaluated on the inner folds and the best one is
refitted on the whole outer training fold with `num.trees = 500`. The refit
uses exactly the tuned values (erratum 2), which the tests verify on the
`ranger` object. Regression uses `splitrule = "variance"`; classification uses
`probability = TRUE` (probability forest, `splitrule = "gini"`), the predicted
class is `sensitive` when its probability is at least 0.5, and the ROC score is
that probability. `replace = FALSE` with `sample.fraction` reproduces the
h2o meaning of `sample_rate`. `num.threads = --threads` (default: all cores
minus one). `importance = "none"`.

### 7.2 Ridge, lasso and elastic net (`glmnet`), model names `ridge`, `lasso`, `enet`

`cv.glmnet(x, y, alpha, foldid, nlambda = 100, standardize = FALSE,
type.measure = "mse")` for regression and `family = "binomial",
type.measure = "auc"` for classification. The selected lambda is
**`lambda.min`** for every model and task (erratum 3).

| Model | alpha |
|---|---|
| `ridge` | 0 |
| `lasso` | 1 |
| `enet` | grid 0.1, 0.2, ..., 0.9; the alpha with the best inner `cvm` at its `lambda.min` is selected (erratum 3), using the same `foldid` for every alpha |

The refit is the `cv.glmnet` object itself (its `glmnet.fit` on the whole
training fold), predictions use `s = "lambda.min"`; classification
probabilities use `type = "response"` (probability of the second level,
`sensitive`), the predicted class is `sensitive` at probability 0.5 or
above. The three `predict()` calls with misplaced parentheses of 2019 have no
counterpart here.

### 7.3 Linear SVM (`e1071`), model name `svm_linear`

`svm(x, y, kernel = "linear", scale = FALSE)` on the standardised 200-gene
matrix (erratum 7).

| Task | Grid |
|---|---|
| Regression (`type = "eps-regression"`) | `cost` in 0.01, 0.1, 1, 10 crossed with `epsilon` in 0.05, 0.1, 0.2 on the scaled response (12 configurations) |
| Classification (`type = "C-classification"`) | `cost` in 0.01, 0.1, 1, 10, 100 (5 configurations) |

The 2019 grid was `cost` in 0.1, 1, 10, 100, 1000 with the default epsilon
on the raw AUC; the upper costs are dropped because they are slow to converge
on standardised data and were never selected in 2019 (0.1 was chosen every
time). No class weights are used in any family; the class imbalance (about
one sensitive in four) is reported through balanced accuracy.

### 7.4 Classification scores and their orientation

Every ROC AUC is computed with
`pROC::roc(response, score, levels = c("resistant", "sensitive"), direction =
"<", quiet = TRUE)`, so the score must increase with the probability of
`sensitive`. For `rf` and the `glmnet` models the score is that probability.
For the SVM the score is the decision value oriented by its column name: the
decision value returned by `e1071` is positive for the class named first in
`colnames(attr(pred, "decision.values"))`; when that name starts with
`resistant` the score is the negated decision value. A unit test checks the
orientation on synthetic data with a planted signal (AUC well above 0.5).

## 8. Leakage check (erratum 6, quantification)

With `--leakage-check` (default on in the full run, off in `--fast`), the
three `glmnet` models of both tasks are refitted on the **same outer folds**
with the 2019 procedure for selection and threshold: the 200 genes chosen on
all cell lines of the drug and, for classification, the class threshold
computed on all cell lines. Standardisation stays fold-wise. The result is
`outputs/<drug>/<task>/leakage_check.csv` with columns `model`, `metric`,
`mean_fold_wise`, `mean_global`, `difference`, `n_folds`. Only the `glmnet`
family is used to keep the budget small; the difference measures the effect of
global selection on the held-out estimates. The random forest and the SVM
are not refitted.

## 9. Metrics and confidence intervals

Computed by `evaluate()` on each held-out fold:

| Task | Metrics |
|---|---|
| Regression | `r2` = 1 minus SSres over SStot, SStot around the mean of the held-out fold (can be negative); `rmse`; `mae` |
| Classification | `accuracy`; `balanced_accuracy` = (sensitivity + specificity) / 2; `auc` (ROC, section 7.4); `sensitivity` = recall of `sensitive`; `specificity` = recall of `resistant` |

If a held-out fold has no cell line of one class (not expected with
stratification), `sensitivity` or `specificity` and `auc` are `NA` for that
fold and the summary reports the number of folds used.

`summarise_metrics()` aggregates over the `k * r` folds per drug, task, model
and metric: `mean`, `sd`, `n_folds`, and a two-sided 95 percent **t-based**
interval `mean plus or minus qt(0.975, n_folds - 1) * sd / sqrt(n_folds)`
(`ci_method = "t"`). The folds of different repeats share cell lines, so the
interval is descriptive, not an exact inferential statement; the caveat is
printed in the README results section.

## 10. Function signatures

All functions are plain R functions in `R/`, sourced by `run_all.R` and by the
tests (no package structure). Names and arguments are binding; additional
internal helpers are allowed.

### `R/build_matrix.R` (data builder)

```r
read_gdsc(path)                       # data.frame: DRUG_NAME, CELL_LINE_NAME, COSMIC_ID, DRUG_ID, AUC
read_depmap_metadata(path)            # data.frame: DepMap_ID, COSMIC_ID (rows with missing COSMIC_ID dropped)
read_depmap_tpm(path, keep_ids = NULL) # list(ids = character, x = numeric matrix rows by genes)
resolve_drug(name, gdsc)              # the DRUG_NAME as spelled in GDSC; error if absent or ambiguous
build_matrix(drug, gdsc, metadata, tpm)
#   -> list(drug, cell_lines, x = matrix [n x 57820, dimnames], auc = numeric n, depmap_ids)
load_drug_matrix(drug, raw_dir = "data/raw", cache_dir = "data/derived",
                 refresh = FALSE, tpm = NULL)
#   -> same list, read from the cache when present
```

`build_matrix()` takes data frames and a matrix, not paths, so the tests feed
it synthetic inputs.

### `R/cv.R` (data builder)

```r
col_vars(x)                                   # sample variance per column, base R
stratum_for(auc, task)                        # factor used to stratify (section 4.1)
make_folds(auc, k, repeats, task, seed)       # data.frame(row, cell_line, repeat, fold)
split_rows(folds, rep, fold)                  # list(train = integer, test = integer)
make_inner_folds(y_train, k_inner, task, seed) # integer foldid of length n_train
select_genes(x_train, n_genes = 200)          # character vector of gene names
class_threshold(auc_train)                    # quantile(auc_train, 0.25)
class_rule(auc, threshold)                    # factor, levels resistant, sensitive
preprocess_fold(x, auc, train, test, task, n_genes = 200,
                selection = c("fold", "global"), threshold = NULL)
#   -> list(x_train, x_test, y_train, y_test, auc_train, auc_test,
#           genes, center, scale, threshold)
write_folds(folds_regression, folds_classification, path)   # outputs/<drug>/folds.csv
```

`selection = "global"` and an explicit `threshold` exist only for the leakage
check of section 8.

### `R/models.R` (models builder)

```r
model_names()                                  # c("rf", "ridge", "lasso", "enet", "svm_linear")
model_grid(model, task, rf_budget = 12, seed = NULL, num_trees = 500)
#   -> data.frame of candidate configurations (rf: the random draw of the declared grid)
fit_model(model, x_train, y_train, task, foldid, grid, seed, threads = 1)
#   -> list(model, task, fit, best = named list, tuning = data.frame(candidate columns, inner_score), elapsed)
predict_model(fitted, x_test)
#   -> regression: numeric vector; classification: list(prob = numeric, class = factor, score = numeric)
inner_score(task, observed, predicted)         # RMSE or ROC AUC on inner validation rows
```

`fit_model()` runs the inner loop, refits the best candidate on all training
rows and returns both. For the SVM regression the response scaling of section
5 happens inside `fit_model()` and `predict_model()`.

### `R/metrics.R` (models builder)

```r
r_squared(observed, predicted)                 # 1 - SSres / SStot
rmse(observed, predicted)
mae(observed, predicted)
roc_auc(observed, score)                       # pROC, levels and direction fixed (section 7.4)
classification_metrics(observed, predicted_class, score)
evaluate(task, observed, predicted, score = NULL)   # named numeric vector of the task metrics
summarise_metrics(metrics_by_fold, conf_level = 0.95)
#   -> data.frame(drug, task, model, metric, n_folds, mean, sd, ci_low, ci_high, ci_method)
```

### `R/run_drug.R` (runner and tests builder)

```r
default_config(fast = FALSE)
#   -> list(k, repeats, k_inner, n_genes, rf_budget, num_trees, seed, threads,
#           models, leakage_check, fast)
run_drug(drug, task, matrix, config, outdir)
#   -> list(metrics_by_fold, summary, tuning, predictions, timing, selected_genes, leakage_check)
#      and writes the files of section 12 under outdir/<drug>/<task>/
run_leakage_check(drug, task, matrix, folds, config, outdir)
write_run_log(path, config, args, timings, packages)
```

### `R/report.R` (runner and tests builder)

```r
make_report(outdir, drugs, tasks)              # per drug and task figures, cross-drug tables and figures
plot_metric_boxes(metrics_by_fold, task)       # ggplot object
plot_observed_vs_predicted(predictions)        # regression, out-of-fold predictions of repeat 1
plot_roc(predictions)                          # classification, pooled scores of repeat 1
write_summary_tables(outdir)                   # outputs/summary.csv and outputs/summary_<task>.md
```

### `run_all.R` (runner and tests builder)

Command line interface only; it sources `R/*.R` in the order
`build_matrix.R`, `cv.R`, `metrics.R`, `models.R`, `run_drug.R`,
`report.R`, parses the arguments with base `commandArgs()` (no extra
package), and calls `load_drug_matrix()`, `run_drug()` and `make_report()`.

## 11. Command line

```bash
Rscript run_all.R [--fast] [--drug erlotinib,paclitaxel] [--task regression|classification|both]
                  [--outdir outputs] [--seed 2019] [--folds 5] [--repeats 3] [--inner-folds 5]
                  [--rf-budget 12] [--threads N] [--no-leakage-check] [--refresh-cache]
```

- `--fast`: the `--fast` column of section 4.1; explicit arguments override
  it. The two drugs are erlotinib (narrow AUC distribution) and paclitaxel
  (wide).
- `--drug`: comma-separated, case-insensitive, any `DRUG_NAME` of the GDSC
  file that maps to one `DRUG_ID`; default the four thesis drugs.
- `--task`: default `both`.
- `--outdir`: default `outputs`; created if missing; existing files of the
  same drug and task are overwritten.
- `--threads`: default `max(1, parallel::detectCores() - 1)`.
- Exit code 0 on success, 1 on any error (R's default for `stop()`); missing
  raw files produce an error that points to `scripts/00_download_data.R`.

Full run: `Rscript run_all.R`. Tests: `Rscript -e
"testthat::test_dir('tests/testthat')"`.

## 12. Output naming scheme

Drug names in paths are lower case; task names are `regression` and
`classification`; model names are those of `model_names()`.

```
outputs/
  run_log.txt                         command line, config, seeds, R and package versions, timings
  summary.csv                         all drugs, tasks, models and metrics (section 9 columns)
  summary_regression.md               the same, as Markdown tables for the README
  summary_classification.md
  fig_summary_regression.png          mean and 95 percent interval per drug and model
  fig_summary_classification.png
  <drug>/
    folds.csv                         cell_line, task, repeat, fold
    <task>/
      metrics_by_fold.csv             drug, task, model, repeat, fold, n_test, one column per metric
      summary.csv                     section 9 columns for this drug and task
      tuning.csv                      drug, task, model, repeat, fold, parameter, value, inner_score (one row per selected parameter)
      tuning_candidates.csv           every evaluated candidate with its inner score
      selected_genes.csv              gene, times_selected, n_splits
      timing.csv                      model, repeat, fold, seconds
      predictions.csv                 cell_line, repeat, fold, model, observed, predicted, score (NOT committed)
      leakage_check.csv               section 8 (full run only)
      fig_metrics.png                 box plots of the fold metrics by model
      fig_observed_vs_predicted.png   regression only
      fig_roc.png                     classification only
```

Figures: `ggplot2`, PNG, 8 by 5 inches at 150 dpi, a colour-blind safe
palette (`viridisLite`, already in `renv.lock`), sentence-case titles with
the drug and task. Fitted model objects are not written to disk by default
(`--save-models` writes `outputs/<drug>/<task>/models/<model>_r<repeat>_f<fold>.rds`,
ignored by git through `*.rds`).

`predictions.csv` holds observed AUC values per cell line, which are GDSC
data, so it stays local. Committed files are aggregates, fold assignments and
figures.

## 13. `.gitignore` exceptions

`.gitignore` currently ignores `*.csv`, `*.rds` and `*.RData` globally and
everything under `data/` except its README. The runner and tests builder
appends:

```
# v1.0 outputs: aggregate tables, fold assignments and figures are committed;
# per cell line predictions and fitted models are not
!outputs/**/*.csv
outputs/**/predictions.csv
outputs/**/models/
```

Negations must come after the global `*.csv` rule, and the re-ignore of
`predictions.csv` after the negation. `*.png` is not ignored, so figures need
no exception. `data/derived/` is covered by `data/*`.

## 14. Runtime budget and `--fast`

Budget: the full run (four drugs, both tasks, 5 folds by 3 repeats, inner
5-fold, random forest budget 12, leakage check) must finish in under about
45 minutes on the reference desktop (12 logical cores). Estimated work per
outer split and task with about 165 to 185 training rows and 200 genes: 61
`ranger` fits of 500 trees (12 candidates by 5 inner folds plus the refit),
11 `cv.glmnet` paths, and 61 (regression) or 26 (classification) `svm` fits;
about 10 to 25 seconds, hence 15 splits by 2 tasks by 4 drugs is roughly 20
to 50 minutes.

Procedure, binding: the runner and tests builder times `--fast` first,
extrapolates (factor 5 for splits, 2 for drugs, 2 for the random forest
budget, plus the leakage check), then runs the full configuration. If the
full run exceeds 45 minutes, the levers are applied in this order, one at a
time, and the final values are recorded in section 17 and in the README:
(1) `--repeats 2`; (2) `--rf-budget 8`; (3) `--inner-folds 3`. The default
values written in `default_config()` are the ones finally used, so that
`Rscript run_all.R` reproduces the committed outputs.

`--fast` exists for smoke runs and for the first timing; its outputs are not
committed.

## 15. Tests (`tests/testthat/`)

Run with `testthat::test_dir()`, no package structure; a
`helper-source.R` sources `R/*.R`. All tests use synthetic data generated in
the tests with a fixed seed (a matrix of 60 cell lines by 500 genes with a
planted linear signal between ten genes and the AUC, AUC values in (0, 1)),
small grids and `num_trees = 50`, so the whole suite runs in well under a
minute. No raw file is read by the tests.

| File | What it checks |
|---|---|
| `test-build_matrix.R` | `build_matrix()` on tiny synthetic GDSC, metadata and TPM inputs: only matched lines kept, one row per cell line, rows sorted by name, gene order preserved, AUC aligned; error on an ambiguous `COSMIC_ID`, on an unknown drug and on missing values; `resolve_drug()` is case-insensitive |
| `test-cv.R` | `make_folds()`: every repeat partitions all rows into `k` disjoint folds; each fold has at least one row of every stratum; identical output for the same seed and different for another; `write_folds()` round trip; `make_inner_folds()` lengths and stratification |
| `test-preprocess.R` | `col_vars()` equals `apply(x, 2, var)`; a gene with huge variance only in the held-out rows is not selected; the threshold equals the training first quartile; the held-out matrix uses training centre and scale and has the same 200 columns in the same order; `selection = "global"` reproduces the 2019 selection |
| `test-models.R` | For each of the five models and both tasks: `fit_model()` returns a fit, a `best` list whose values belong to the declared grid, and a `tuning` table with one inner score per candidate; `predict_model()` returns the right length and type; classification probabilities lie in [0, 1]; the `ranger` object carries exactly the tuned `mtry`, `min.node.size`, `sample.fraction`, `max.depth` and `num.trees`; the elastic net alpha is in the grid and the lambda is `lambda.min`; the SVM uses 200 columns; the same `foldid` is used by every family |
| `test-metrics.R` | `r_squared()` is 1 for a perfect prediction, 0 for the mean, negative for a worse one, and equals `1 - SSres / SStot` on a hand-computed example; `rmse()` and `mae()` on known values; classification metrics on a hand-built confusion matrix; `roc_auc()` is above 0.5 for a score that increases with `sensitive` and the SVM orientation of section 7.4 holds; `summarise_metrics()` reproduces a hand-computed t interval and handles `NA` |
| `test-run_drug.R` | `run_drug()` on the synthetic matrix with `k = 2`, `repeats = 1`, `k_inner = 3` (`cv.glmnet()` refuses fewer than three inner folds, so the value 2 of the first draft of this table was not usable), tiny grids, both tasks: writes every file of section 12 to a temporary directory, the summary has one row per model and metric, `predictions.csv` has `n * models` rows, and a second run with the same seed gives identical `metrics_by_fold.csv` |
| `test-report.R` | `make_report()` on the files written by the previous test creates the PNG figures and `summary.csv` without error |

## 16. Builders and files

| Builder | Writes | Must not touch |
|---|---|---|
| data | `R/build_matrix.R`, `R/cv.R` | anything else |
| models | `R/models.R`, `R/metrics.R`; installs `ranger` with `renv::install("ranger")` and runs `renv::snapshot()` (`renv.lock` updated; `ranger` must be referenced with `library()` or `ranger::` in `R/models.R` so the implicit snapshot records it) | `R/cv.R` |
| runner and tests | `run_all.R`, `R/run_drug.R`, `R/report.R`, `tests/testthat.R`, `tests/testthat/helper-source.R`, `tests/testthat/test-*.R`, `.gitignore` (section 13), `renv.lock` after `renv::install("testthat")` and `renv::snapshot()`; runs the test suite, the `--fast` run and the full run; commits `outputs/` as in section 12; fills section 17 of this file | `legacy/`, `scripts/`, `docs/errata.md` |
| docs | `README.md` (v1.0 sections: what v1.0 is, how to run, results tables from `outputs/summary_*.md`, substitution of h2o, repository structure, provenance roadmap turned into "v1.0 (this release)"), `docs/comparison_2019_vs_v1.md`, a v1.0 note appended to `data/README.md` (the `data/derived/` cache and the fact that no new input is needed), `CITATION.cff` version and date | `docs/errata.md`, `legacy/`, existing statements of `data/README.md` |

Never modified by anyone: `legacy/*.R`, `LICENSE`, `docs/errata.md`,
`scripts/00_download_data.R`, the existing text of `data/README.md`. No git
command that rewrites history. No data file committed; `data/raw/` and
`data/derived/` stay ignored.

Dependencies to add to `renv.lock`: `ranger` and `testthat` (from CRAN).
Everything else used by v1.0 (`caret`, `glmnet`, `e1071`, `pROC`,
`ggplot2`, `viridisLite`, `data.table`, `readxl`) is already recorded.

## 17. Implementation record

Filled after the full run.

| Item | Value |
|---|---|
| Date of the full run | 2026-10-07 |
| Final `k`, `repeats`, `k_inner`, `rf_budget` | 5, 3, 5, 12 (the defaults of `default_config()`) |
| Levers applied and why | none: the full run finished under the 45 minute budget |
| Wall time of `--fast` and of the full run | `--fast` 1.2 minutes (71.7 s); full run 39.6 minutes (2377.9 s, of which 253 to 412 s per drug and task, inner loop and leakage check included, plus the report) |
| Threads | 11 (`parallel::detectCores() - 1` on a 12 logical core desktop) |
| `ranger`, `glmnet`, `e1071`, `pROC`, `caret` versions | 0.18.0, 4.1.8, 1.7.14, 1.18.5, 6.0.94 (from `renv.lock`) |
| Tests passed | 555 expectations in 7 files, 0 failures (`testthat` 3.3.2) |

Choices recorded by the builders where the specification left room or was
adjusted:

- `build_matrix()` and `load_drug_matrix()` return the list of section 10
  plus `drug_name`, `drug_id` and `cosmic_ids`; `drug` is the lower-case key
  used in paths and in the cache name `data/derived/<drug>_matrix.rds`.
  `preprocess_fold()` also returns `selection`; `threshold` is `NA` for
  regression.
- `tuning_candidates.csv` is written in long format (`drug`, `task`, `model`,
  `repeat`, `fold`, `candidate`, `parameter`, `value`, `inner_score`),
  because the families have different parameter columns.
- `fit_model()` receives the split seed `S + 1000 t + j` and applies the
  family offsets of section 4.3 itself (`ranger` seed + 1, `set.seed` + 5
  before the first `svm` and `glmnet` fit); the random forest draw uses
  + 100 in the runner. `ranger` does not store `sample.fraction` as a field
  of the fitted object, so the test of the tuned values reads it from the
  call.
- `predict_model()` returns `prob = NA` for the SVM classification (the
  decision value is the score, section 7.4); `predictions.csv` carries the
  score column.
- The `run_drug` test uses `k_inner = 3` (section 15), because `cv.glmnet()`
  stops with fewer than three folds.
- `renv::install("testthat")` upgraded `R6`, `cli`, `jsonlite`, `rlang` and
  `withr` in the project library; `renv.lock` records the new versions.

## 18. Material for `docs/comparison_2019_vs_v1.md`

The docs builder writes the comparison from this section and from
`outputs/summary.csv`. A local copy of the deposited thesis PDF was found
next to the raw files and read for the numbers below; it is not part of the
repository. The thesis reports, per drug, model, split and task, a single
test-set value: accuracy for classification and an "R2" for regression.
Nothing else is available from 2019: no ROC AUC values (only curves), no
sensitivity or specificity, no RMSE or MAE, no fold-level variability and no
confidence intervals. The 2019 "R2" is `cor(y, y - prediction)^2` (erratum
1), not a coefficient of determination, so it is **not comparable** with the
v1.0 R squared; the 2019 accuracies are comparable only loosely, because the
2019 class label is inverted (erratum 5), the splits differ (section 4) and
the 2019 selection and threshold saw the test lines (erratum 6).

### 18.1 Classification, accuracy on the test set (2019, pages 39 to 60 of the deposited PDF)

| Model (2019 label) | Split | Erlotinib | Rapamycin | Sunitinib | Paclitaxel |
|---|---|---|---|---|---|
| Random forest, h2o search | 80/20 | 0.7561 | 0.7073 | 0.7778 | 0.7609 |
| Random forest, h2o search | 60/20/20 | 0.7826 | 0.7073 | 0.8043 | 0.8095 |
| Random forest, caret search | 80/20 | 0.7561 | 0.7073 | 0.7556 | 0.7609 |
| Random forest, caret search | 60/20/20 | 0.8333 | 0.7381 | 0.7826 | 0.7174 |
| Ridge logistic | 80/20 | 0.7561 | 0.7561 | 0.7556 | 0.7556 |
| Ridge logistic | 60/20/20 | 0.7561 | 0.7561 | 0.7556 | 0.7556 |
| Lasso logistic | 80/20 | 0.7073 | 0.7778 | 0.7556 | 0.7778 |
| Lasso logistic | 60/20/20 | 0.8764 | 0.7561 | 0.8222 | 0.6889 |
| Elastic net logistic (a lasso in the 80/20 split, erratum 3) | 80/20 | 0.7561 | 0.7317 | 0.7556 | 0.7778 |
| Elastic net logistic | 60/20/20 | 0.7561 | 0.7561 | 0.7556 | 0.6889 |
| Linear SVM | 80/20 | 0.6829 | 0.7317 | 0.7111 | 0.6304 |
| Linear SVM | 60/20/20 | 0.6429 | 0.6190 | 0.6087 | 0.7174 |

Note: 0.7561 (31 of 41), 0.7556 (34 of 45) and similar values equal the
proportion of the majority class in the test sets, that is, a model that
predicts "resistant" for every line (see the ridge rows and the additional
observations of `docs/errata.md`). The abstract's headline value, 87.64
percent for erlotinib, is the 60/20/20 lasso logistic model.

### 18.2 Regression, "R2" as computed in 2019 (not a coefficient of determination)

| Model (2019 label) | Split | Erlotinib | Rapamycin | Sunitinib | Paclitaxel |
|---|---|---|---|---|---|
| Random forest, h2o search | 80/20 | 0.9646 | 0.7004 | 0.7144 | 0.8850 |
| Random forest, h2o search | 60/20/20 | 0.2785 | 0.7374 | 0.9131 | 0.9306 |
| Random forest, caret search | 80/20 | 0.9870 | 0.6235 | 0.7892 | 0.9113 |
| Random forest, caret search | 60/20/20 | 0.9827 | 0.9017 | 0.9490 | 0.9730 |
| Ridge | 80/20 | 1 | 1 | 1 | 0.9993 |
| Ridge | 60/20/20 | 0.9801 | 0.9924 | 0.9207 | 1.0000 |
| Lasso | 80/20 | 1 | 0.9974 | 0.9996 | 1 |
| Lasso | 60/20/20 | 0.9689 | 0.9810 | 0.8940 | 1 |
| Elastic net (a lasso in the 80/20 split, erratum 3) | 80/20 | 1 | 0.9661 | 0.9920 | 1 |
| Elastic net | 60/20/20 | 0.9618 | 0.9852 | 0.9972 | 1 |
| Linear SVM (2000 genes, erratum 7) | 80/20 | 0.8942 | 0.3801 | 0.2404 | 0.5558 |
| Linear SVM (2000 genes, erratum 7) | 60/20/20 | 0.0374 | 0.1500 | 0.5116 | 0.7149 |

Values of 1 correspond to nearly constant predictions (the most penalised
models), which is the artefact described in erratum 1. The thesis also gives
the `mtry` selected by caret per drug (classification 80/20: 2, 33, 2, 43;
classification 60/20/20: 2, 33, 2, 2; regression 80/20: 2, 74, 2, 22, for
erlotinib, rapamycin, sunitinib and paclitaxel respectively); the other
per-drug hyperparameters appear only in figures and are not transcribed.

### 18.3 Hyperparameters hard-coded in the legacy scripts (erlotinib)

| Script | Final fit |
|---|---|
| `h2ocontinuosRF.R` | 80/20: `ntree = 200`, `mtry = 25` (the other arguments ignored); 60/20/20: `ntree = 350`, `mtry = 35` |
| `h2odiscretizadosRF.R` | 80/20 and 60/20/20: `ntree = 500`, `mtry = 15` (the other arguments ignored) |
| `randomsearchcontinuosRF.R` | `mtry = 2` in both splits |
| `randomsearchdiscretizadosRF.R` | 80/20: `mtry = 43`; 60/20/20: `mtry = 2` |
| `ridgelassoelasticcontinuosrandomsearch.R` | 80/20: `cv.glmnet` with `lambda.1se`; 60/20/20: ridge `lambda = 4.328761`, lasso `lambda = 0.02477076`, elastic net `alpha = 0.8105263`, `lambda = 0.07328962` |
| `ridgelassoelasticdiscretizadosrandomsearch.R` | 80/20: `cv.glmnet` with `lambda.min`; 60/20/20: ridge `lambda = 1000`, lasso `lambda = 0.04977024`, elastic net `alpha = 0.1473684`, `lambda = 0.2952464` |
| `SVMLinealcontinuos.R`, `SVMLinealdiscretizados.R` | `cost = 0.1` in every split (grid 0.1 to 1000) |

The 2019 grids: h2o `ntrees` 200, 350, 500; `mtries` 15, 25, 35;
`max_depth` 20 to 40 by 5; `min_rows` 1, 3, 5; `nbins` 10 to 30 by 5;
`sample_rate` 0.55, 0.632, 0.75, `RandomDiscrete` with a 30 minute budget
and no seed; caret `tuneLength = 20` (a grid, erratum 4); glmnet lambda
`10^seq(-3, 3, length = 100)` with caret in the 60/20/20 split.

### 18.4 Structure of the comparison document

1. Purpose and caveats (this section's first paragraph).
2. Method differences, one row per erratum (section 1 of this file).
3. The 2019 tables of 18.1 and 18.2, as they are.
4. The v1.0 tables from `outputs/summary_classification.md` and
   `outputs/summary_regression.md` (mean and 95 percent interval over the
   folds), with the thesis drugs in the same order.
5. The leakage check table (section 8) per drug and task.
6. What cannot be compared and why (the list in the first paragraph), and a
   short reading of the differences without over-interpretation.
