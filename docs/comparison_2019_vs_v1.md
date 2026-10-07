# Comparison of the 2019 thesis results and v1.0

This document places the numbers reported in the 2019 thesis next to the
numbers produced by the v1.0 reimplementation, and explains why most of them
cannot be compared directly. The 2019 figures were transcribed from a local
copy of the deposited thesis PDF (the thesis is linked from `README.md`, not
redistributed here) and checked against the text of its results section. The
v1.0 figures come from `outputs/summary.csv`, written by `run_all.R` in the
full run of 2026-10-07 (sections 4, 5 and 7 were completed from those files
after that run).

Scope reminder: research and teaching only. The models predict the drug
response of cancer cell lines, not of patients, and nothing here is a clinical
tool.

## 1. Purpose and caveats

The thesis reports, per drug, model, split and task, a single test-set value:
accuracy for classification and an "R2" for regression. Nothing else is
available from 2019: no ROC AUC values (only curves), no sensitivity or
specificity, no RMSE or MAE, no fold-level variability and no confidence
intervals. Two further facts limit the comparison:

- The 2019 "R2" is `cor(y, y - prediction)^2` (erratum 1 of
  `docs/errata.md`), not a coefficient of determination, so it is **not
  comparable** with the v1.0 R squared (1 minus SSres over SStot on the
  held-out fold). It approaches 1 whenever the prediction is nearly constant.
- The 2019 accuracies are comparable only loosely: the 2019 class label is
  inverted (erratum 5), the splits differ (one 80/20 or 60/20/20 test set in
  2019 against repeated stratified cross-validation in v1.0, section 4 of
  `docs/v1_spec.md`), and the 2019 gene selection and class threshold saw the
  test lines (erratum 6).

The 2019 numbers are therefore reproduced below as historical values, and the
v1.0 numbers are read on their own terms. Section 6 lists what cannot be
compared and section 7 gives a short reading of the differences.

## 2. Method differences, one row per erratum

| Erratum (`docs/errata.md`) | 2019 code | v1.0 | Effect on the comparison |
|---|---|---|---|
| 1. R squared | `cor(y, y - prediction)^2` on the test set | 1 minus SSres over SStot on the held-out fold, plus RMSE and MAE (`R/metrics.R`) | The regression numbers of the two versions measure different things; values of 1 in 2019 are an artefact |
| 2. Random forest refit | h2o grid search, then `randomForest()` with the tuned arguments silently dropped | `ranger` fitted with the tuned `mtry`, `min.node.size`, `sample.fraction`, `max.depth` and `num.trees` (`R/models.R`); h2o not used | The 2019 random forests are default `randomForest` fits with the tuned `ntree` and `mtry` only |
| 3. Elastic net and lambda rule | The 80/20 "elastic net" is a lasso (`alpha` not passed); `lambda.1se` in regression and `lambda.min` in classification | `alpha` tuned on the grid 0.1 to 0.9 for both tasks; `lambda.min` everywhere | The 2019 80/20 elastic net rows are a second lasso |
| 4. Search procedure | `caret::train()` grid search labelled random search; two calls with bootstrap resampling instead of cross-validation | Explicitly declared grids; the random forest grid sampled at random with a seed and a fixed budget; the same inner folds for every model family (`R/cv.R`, `R/models.R`) | The 2019 hyperparameters cannot all be traced to the declared search (for example `mtry = 43`) |
| 5. Class label | "Sensitive" is AUC at or above the third quartile (the most resistant quartile) | Sensitive is AUC strictly below the first quartile of the training fold (`R/cv.R`) | Accuracies are unaffected by the inversion, but sensitivity, specificity and the ROC reading are |
| 6. Leakage | The 200 most variable genes and the class threshold computed on all cell lines before the split | Both computed inside each outer training fold; a leakage check quantifies the difference on the `glmnet` family (section 5) | The 2019 test estimates are not strictly out of sample |
| 7. Feature set and saved objects | SVM regression on 2000 genes; the wrong 60/20/20 SVM object saved; misplaced parentheses in three `predict()` calls | The same 200-gene matrix for every model including the SVM; fitted objects and tuning tables saved per fold (`R/models.R`, `R/run_drug.R`) | The 2019 SVM regression rows use a different feature set from the other models |
| 8. Reproducibility | Partitions depend on `sample()` and on the random numbers consumed by earlier fits; h2o searches without seed and with a time budget | Seeds fixed, fold assignments written to `outputs/<drug>/folds.csv`, no h2o, package versions in `renv.lock` (`run_all.R`) | The 2019 numbers cannot be reproduced exactly; the v1.0 numbers can |

Other differences that are not errata: v1.0 replaces the two 2019 splits by
repeated stratified 5-fold cross-validation with 3 repeats (3 folds and 1
repeat with `--fast`), tunes every model with the same inner 5-fold
resampling, standardises the selected genes with training-fold parameters,
drops the upper SVM costs (100 and 1000, never selected in 2019) and reports
per-fold metrics with a mean and a t-based 95 percent interval.

## 3. The 2019 numbers as reported in the thesis

The model labels are those of the thesis. "h2o search" and "caret search"
are the two random forest scripts (`h2o*.R` and `randomsearch*.R` in
`legacy/`). The page ranges refer to the deposited PDF; `legacy/README.md`
maps each script to its thesis sections.

### 3.1 Classification, accuracy on the test set (pages 39 to 60)

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
proportion of the majority class in the test sets, that is, the accuracy of a
model that predicts "resistant" for every cell line (see the ridge rows and the
additional observations of `docs/errata.md`). The headline value of the
abstract, 87.64 percent for erlotinib, is the 60/20/20 lasso logistic model.
The thesis gives no ROC AUC values, only curves.

### 3.2 Regression, "R2" as computed in 2019 (pages 40 to 61; not a coefficient of determination)

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
models), which is the artefact described in erratum 1. The thesis reports no
RMSE or MAE on the test sets.

### 3.3 Hyperparameters reported or hard-coded in 2019

The thesis gives the `mtry` selected by caret per drug (erlotinib, rapamycin,
sunitinib and paclitaxel, in this order): classification 80/20: 2, 33, 2, 43;
classification 60/20/20: 2, 33, 2, 2; regression 80/20: 2, 74, 2, 22. The
other per-drug hyperparameters appear only in figures and are not
transcribed. The legacy scripts, which are tied to erlotinib, hard-code the
following final fits:

| Script | Final fit |
|---|---|
| `h2ocontinuosRF.R` | 80/20: `ntree = 200`, `mtry = 25` (the other arguments ignored, erratum 2); 60/20/20: `ntree = 350`, `mtry = 35` |
| `h2odiscretizadosRF.R` | 80/20 and 60/20/20: `ntree = 500`, `mtry = 15` (the other arguments ignored) |
| `randomsearchcontinuosRF.R` | `mtry = 2` in both splits |
| `randomsearchdiscretizadosRF.R` | 80/20: `mtry = 43`; 60/20/20: `mtry = 2` |
| `ridgelassoelasticcontinuosrandomsearch.R` | 80/20: `cv.glmnet` with `lambda.1se`; 60/20/20: ridge `lambda = 4.328761`, lasso `lambda = 0.02477076`, elastic net `alpha = 0.8105263`, `lambda = 0.07328962` |
| `ridgelassoelasticdiscretizadosrandomsearch.R` | 80/20: `cv.glmnet` with `lambda.min`; 60/20/20: ridge `lambda = 1000`, lasso `lambda = 0.04977024`, elastic net `alpha = 0.1473684`, `lambda = 0.2952464` |
| `SVMLinealcontinuos.R`, `SVMLinealdiscretizados.R` | `cost = 0.1` in every split (grid 0.1 to 1000) |

The 2019 grids: h2o `ntrees` 200, 350, 500; `mtries` 15, 25, 35; `max_depth`
20 to 40 by 5; `min_rows` 1, 3, 5; `nbins` 10 to 30 by 5; `sample_rate` 0.55,
0.632, 0.75, strategy `RandomDiscrete` with a 30 minute budget and no seed;
caret `tuneLength = 20` (a grid, erratum 4); glmnet lambda
`10^seq(-3, 3, length = 100)` with caret in the 60/20/20 split.

### 3.4 Substitution of h2o by ranger

The 2019 random forest searches ran in h2o, which needs a Java runtime, had
no seed and a time budget, so the set of models explored depended on the
clock; their tuned values were then dropped by `randomForest()` (erratum 2).
v1.0 does not call h2o. The random forest is fitted with `ranger`, which
accepts the equivalents of the h2o parameters directly (`mtry` for `mtries`,
`min.node.size` for `min_rows`, `sample.fraction` with `replace = FALSE` for
`sample_rate`, `max.depth` for `max_depth`, `num.trees` for `ntrees`; `nbins`
has no equivalent and is dropped), is deterministic given a seed and has no
Java dependency. The v1.0 grid, its random draw and the mapping of every
parameter are declared in section 7.1 of `docs/v1_spec.md`.

## 4. v1.0 results

Mean over the 15 outer folds (5 folds by 3 repeats) and two-sided 95 percent
t-based interval, from `outputs/summary.csv` (the same numbers, one table per
metric, are in `outputs/summary_classification.md` and
`outputs/summary_regression.md`), with the thesis drugs in the same order
(erlotinib, rapamycin, sunitinib, paclitaxel) and the model names of
`model_names()`. The interval is descriptive, because the folds of different
repeats share cell lines. In every held-out fold about one cell line in four
is `sensitive` (share on the full matrices: 0.252, 0.250, 0.249 and 0.252),
so a model that predicts `resistant` for every line scores an accuracy of
about 0.75, a balanced accuracy of 0.5 and a sensitivity of 0.

### 4.1 Classification (sensitive = AUC below the first quartile of the training fold)

| Drug | Model | Accuracy | Balanced accuracy | ROC AUC | Sensitivity | Specificity | n folds |
|---|---|---|---|---|---|---|---|
| erlotinib | rf | 0.791 (0.772 to 0.811) | 0.617 (0.585 to 0.648) | 0.706 (0.630 to 0.781) | 0.264 (0.205 to 0.323) | 0.970 (0.954 to 0.985) | 15 |
| erlotinib | ridge | 0.752 (0.739 to 0.766) | 0.512 (0.494 to 0.530) | 0.703 (0.618 to 0.788) | 0.024 (-0.011 to 0.060) | 1.000 (1.000 to 1.000) | 15 |
| erlotinib | lasso | 0.746 (0.723 to 0.769) | 0.542 (0.512 to 0.571) | 0.624 (0.539 to 0.708) | 0.124 (0.044 to 0.205) | 0.959 (0.919 to 0.998) | 15 |
| erlotinib | enet | 0.744 (0.726 to 0.762) | 0.536 (0.504 to 0.568) | 0.642 (0.557 to 0.726) | 0.111 (0.024 to 0.198) | 0.961 (0.928 to 0.993) | 15 |
| erlotinib | svm_linear | 0.757 (0.720 to 0.794) | 0.599 (0.571 to 0.627) | 0.650 (0.587 to 0.714) | 0.279 (0.215 to 0.343) | 0.919 (0.863 to 0.976) | 15 |
| rapamycin | rf | 0.724 (0.703 to 0.746) | 0.506 (0.485 to 0.528) | 0.640 (0.602 to 0.677) | 0.063 (0.021 to 0.104) | 0.950 (0.923 to 0.978) | 15 |
| rapamycin | ridge | 0.745 (0.738 to 0.753) | 0.500 (0.500 to 0.500) | 0.591 (0.549 to 0.633) | 0.000 (0.000 to 0.000) | 1.000 (1.000 to 1.000) | 15 |
| rapamycin | lasso | 0.705 (0.677 to 0.733) | 0.496 (0.474 to 0.518) | 0.562 (0.498 to 0.626) | 0.070 (0.023 to 0.116) | 0.923 (0.877 to 0.969) | 15 |
| rapamycin | enet | 0.704 (0.684 to 0.724) | 0.504 (0.485 to 0.524) | 0.596 (0.539 to 0.653) | 0.096 (0.047 to 0.145) | 0.913 (0.876 to 0.949) | 15 |
| rapamycin | svm_linear | 0.694 (0.662 to 0.727) | 0.519 (0.483 to 0.555) | 0.605 (0.545 to 0.665) | 0.159 (0.061 to 0.258) | 0.879 (0.819 to 0.938) | 15 |
| sunitinib | rf | 0.738 (0.717 to 0.759) | 0.543 (0.513 to 0.572) | 0.694 (0.654 to 0.735) | 0.154 (0.097 to 0.210) | 0.932 (0.906 to 0.958) | 15 |
| sunitinib | ridge | 0.750 (0.736 to 0.764) | 0.500 (0.500 to 0.500) | 0.696 (0.657 to 0.736) | 0.000 (0.000 to 0.000) | 1.000 (1.000 to 1.000) | 15 |
| sunitinib | lasso | 0.745 (0.725 to 0.765) | 0.516 (0.500 to 0.533) | 0.705 (0.672 to 0.738) | 0.056 (0.013 to 0.100) | 0.977 (0.954 to 0.999) | 15 |
| sunitinib | enet | 0.753 (0.735 to 0.770) | 0.513 (0.492 to 0.534) | 0.706 (0.671 to 0.742) | 0.032 (-0.016 to 0.081) | 0.994 (0.985 to 1.003) | 15 |
| sunitinib | svm_linear | 0.710 (0.666 to 0.754) | 0.540 (0.499 to 0.581) | 0.626 (0.577 to 0.674) | 0.199 (0.130 to 0.267) | 0.881 (0.825 to 0.938) | 15 |
| paclitaxel | rf | 0.754 (0.742 to 0.765) | 0.515 (0.503 to 0.527) | 0.647 (0.602 to 0.692) | 0.039 (0.015 to 0.064) | 0.990 (0.983 to 0.998) | 15 |
| paclitaxel | ridge | 0.751 (0.741 to 0.761) | 0.500 (0.500 to 0.500) | 0.611 (0.570 to 0.652) | 0.000 (0.000 to 0.000) | 1.000 (1.000 to 1.000) | 15 |
| paclitaxel | lasso | 0.693 (0.656 to 0.730) | 0.481 (0.455 to 0.508) | 0.542 (0.495 to 0.589) | 0.058 (0.015 to 0.102) | 0.904 (0.850 to 0.958) | 15 |
| paclitaxel | enet | 0.703 (0.674 to 0.733) | 0.496 (0.478 to 0.513) | 0.570 (0.528 to 0.611) | 0.085 (0.036 to 0.133) | 0.907 (0.864 to 0.950) | 15 |
| paclitaxel | svm_linear | 0.665 (0.631 to 0.699) | 0.488 (0.466 to 0.511) | 0.563 (0.516 to 0.610) | 0.137 (0.048 to 0.226) | 0.840 (0.773 to 0.907) | 15 |

### 4.2 Regression (R squared = 1 minus SSres over SStot on the held-out fold)

| Drug | Model | R squared | RMSE | MAE | n folds |
|---|---|---|---|---|---|
| erlotinib | rf | 0.087 (0.013 to 0.161) | 0.074 (0.061 to 0.088) | 0.041 (0.038 to 0.045) | 15 |
| erlotinib | ridge | 0.120 (0.045 to 0.196) | 0.073 (0.059 to 0.087) | 0.040 (0.037 to 0.044) | 15 |
| erlotinib | lasso | 0.095 (-0.009 to 0.199) | 0.074 (0.060 to 0.087) | 0.041 (0.038 to 0.044) | 15 |
| erlotinib | enet | 0.115 (0.024 to 0.206) | 0.073 (0.060 to 0.087) | 0.041 (0.037 to 0.044) | 15 |
| erlotinib | svm_linear | 0.137 (0.075 to 0.199) | 0.073 (0.059 to 0.087) | 0.038 (0.034 to 0.041) | 15 |
| rapamycin | rf | 0.111 (0.071 to 0.152) | 0.207 (0.195 to 0.219) | 0.152 (0.143 to 0.162) | 15 |
| rapamycin | ridge | 0.087 (0.049 to 0.125) | 0.210 (0.198 to 0.222) | 0.152 (0.143 to 0.160) | 15 |
| rapamycin | lasso | 0.093 (0.051 to 0.135) | 0.209 (0.196 to 0.223) | 0.152 (0.143 to 0.161) | 15 |
| rapamycin | enet | 0.088 (0.047 to 0.128) | 0.210 (0.197 to 0.223) | 0.153 (0.144 to 0.161) | 15 |
| rapamycin | svm_linear | 0.009 (-0.048 to 0.067) | 0.218 (0.206 to 0.231) | 0.144 (0.135 to 0.153) | 15 |
| sunitinib | rf | 0.059 (-0.093 to 0.211) | 0.143 (0.132 to 0.154) | 0.097 (0.091 to 0.102) | 15 |
| sunitinib | ridge | 0.095 (-0.003 to 0.193) | 0.142 (0.130 to 0.154) | 0.094 (0.090 to 0.099) | 15 |
| sunitinib | lasso | 0.023 (-0.124 to 0.170) | 0.147 (0.135 to 0.159) | 0.099 (0.093 to 0.104) | 15 |
| sunitinib | enet | 0.063 (-0.062 to 0.189) | 0.144 (0.132 to 0.156) | 0.096 (0.091 to 0.101) | 15 |
| sunitinib | svm_linear | -0.012 (-0.119 to 0.095) | 0.151 (0.137 to 0.164) | 0.097 (0.092 to 0.103) | 15 |
| paclitaxel | rf | 0.053 (0.016 to 0.090) | 0.193 (0.189 to 0.198) | 0.165 (0.160 to 0.170) | 15 |
| paclitaxel | ridge | 0.043 (0.013 to 0.073) | 0.195 (0.190 to 0.199) | 0.167 (0.162 to 0.171) | 15 |
| paclitaxel | lasso | 0.029 (0.005 to 0.053) | 0.196 (0.191 to 0.201) | 0.168 (0.164 to 0.171) | 15 |
| paclitaxel | enet | 0.042 (0.014 to 0.070) | 0.195 (0.190 to 0.199) | 0.166 (0.162 to 0.170) | 15 |
| paclitaxel | svm_linear | -0.103 (-0.153 to -0.054) | 0.209 (0.205 to 0.213) | 0.172 (0.168 to 0.176) | 15 |

Run record (from section 17 of `docs/v1_spec.md` and `outputs/run_log.txt`):
full run of 2026-10-07 with the defaults `k = 5`, `repeats = 3`, `k_inner =
5` and `rf_budget = 12`, no lever of section 14 of the specification applied,
39.6 minutes of wall time with 11 `ranger` threads (the `--fast` smoke run
took 1.2 minutes); `ranger` 0.18.0, `glmnet` 4.1.8, `e1071` 1.7.14, `pROC`
1.18.5 and `caret` 6.0.94 as recorded in `renv.lock`; 555 unit tests passed.

## 5. Leakage check (erratum 6, quantification)

With the leakage check on (default in the full run), the three `glmnet`
models of both tasks are refitted on the same outer folds with the 2019
procedure for selection and threshold: the 200 genes chosen on all cell lines
of the drug and, for classification, the class threshold computed on all cell
lines. Standardisation stays fold-wise. The difference between the fold-wise
and the global procedure measures the effect of global selection on the
held-out estimates (`outputs/<drug>/<task>/leakage_check.csv`).

The table lists the main metric of each task (ROC AUC and R squared) for the
three `glmnet` models; `outputs/<drug>/<task>/leakage_check.csv` holds every
metric. The difference is the global mean minus the fold-wise mean.

| Drug | Task | Model | Metric | Mean, fold-wise | Mean, global (2019 procedure) | Difference | n folds |
|---|---|---|---|---|---|---|---|
| erlotinib | classification | ridge | auc | 0.703 | 0.699 | -0.004 | 15 |
| erlotinib | classification | lasso | auc | 0.624 | 0.619 | -0.005 | 15 |
| erlotinib | classification | enet | auc | 0.642 | 0.649 | +0.008 | 15 |
| erlotinib | regression | ridge | r2 | 0.120 | 0.129 | +0.009 | 15 |
| erlotinib | regression | lasso | r2 | 0.095 | 0.090 | -0.005 | 15 |
| erlotinib | regression | enet | r2 | 0.115 | 0.108 | -0.007 | 15 |
| rapamycin | classification | ridge | auc | 0.591 | 0.583 | -0.008 | 15 |
| rapamycin | classification | lasso | auc | 0.562 | 0.558 | -0.004 | 15 |
| rapamycin | classification | enet | auc | 0.596 | 0.611 | +0.015 | 15 |
| rapamycin | regression | ridge | r2 | 0.087 | 0.089 | +0.002 | 15 |
| rapamycin | regression | lasso | r2 | 0.093 | 0.100 | +0.007 | 15 |
| rapamycin | regression | enet | r2 | 0.088 | 0.093 | +0.005 | 15 |
| sunitinib | classification | ridge | auc | 0.696 | 0.696 | -0.000 | 15 |
| sunitinib | classification | lasso | auc | 0.705 | 0.704 | -0.001 | 15 |
| sunitinib | classification | enet | auc | 0.706 | 0.711 | +0.004 | 15 |
| sunitinib | regression | ridge | r2 | 0.095 | 0.088 | -0.007 | 15 |
| sunitinib | regression | lasso | r2 | 0.023 | 0.026 | +0.003 | 15 |
| sunitinib | regression | enet | r2 | 0.063 | 0.043 | -0.020 | 15 |
| paclitaxel | classification | ridge | auc | 0.611 | 0.605 | -0.006 | 15 |
| paclitaxel | classification | lasso | auc | 0.542 | 0.557 | +0.015 | 15 |
| paclitaxel | classification | enet | auc | 0.570 | 0.577 | +0.007 | 15 |
| paclitaxel | regression | ridge | r2 | 0.043 | 0.046 | +0.003 | 15 |
| paclitaxel | regression | lasso | r2 | 0.029 | 0.035 | +0.006 | 15 |
| paclitaxel | regression | enet | r2 | 0.042 | 0.054 | +0.011 | 15 |

The largest absolute differences over all metrics of the CSV files are 0.033
(rapamycin, elastic net, sensitivity) and 0.032 (sunitinib, lasso,
sensitivity); for ROC AUC they lie between -0.008 and +0.015 and for
R squared between -0.020 and +0.011, with no consistent sign.

The random forest and the SVM are not refitted, to keep the budget small.

## 6. What cannot be compared and why

| Quantity | 2019 | v1.0 | Comparable |
|---|---|---|---|
| Regression "R2" | `cor(y, y - prediction)^2` on one test set | 1 minus SSres over SStot per held-out fold | No: different definitions; the 2019 values near 1 are artefacts of erratum 1 |
| RMSE, MAE | Not reported on the test sets | Per fold, mean and interval | No 2019 value |
| Accuracy | One test set, class label inverted, selection and threshold with leakage | Per fold, label "sensitive = low AUC", fold-wise selection and threshold | Loosely: both are proportions of correct labels, but the class definition, the splits and the leakage differ |
| Balanced accuracy, sensitivity, specificity | Not reported (confusion matrices and ROC curves appear only as figures) | Per fold | No 2019 value; the 2019 reading of sensitivity would in any case be inverted (erratum 5) |
| ROC AUC | Curves only, no values; the SVM curves may be inverted | Per fold, with the score oriented towards the sensitive class | No 2019 value |
| Random forest | h2o search then `randomForest()` refit that drops the tuned values; caret grid search | `ranger` with the tuned values; declared grid sampled with a seed | Model family only; not the same estimator |
| Elastic net | A lasso in the 80/20 split; `lambda.1se` or `lambda.min` depending on the task | `alpha` tuned, `lambda.min` everywhere | 80/20 rows are lasso fits |
| Linear SVM regression | 2000 genes | 200 genes | Different feature set |
| Variability | A single number per cell | 15 folds (or the final `k` times `repeats`) with mean, sd and interval | The 2019 numbers carry no uncertainty |
| Partitions | Not reproducible (erratum 8) | Saved in `outputs/<drug>/folds.csv` | The 2019 test sets cannot be rebuilt |

## 7. Reading of the differences

What can be said before looking at the v1.0 numbers:

- The 2019 regression results do not show that the penalised models fit the
  AUC almost perfectly; they show that those models predicted an almost
  constant value, which the 2019 "R2" rewards with values near 1. A v1.0
  R squared near zero or negative for the same models on the same drug is the
  expected correction, not a regression in quality.
- Most 2019 classification accuracies sit at the majority-class proportion
  (about 0.75), which a constant "resistant" prediction attains. The v1.0
  balanced accuracy, sensitivity and ROC AUC make that case visible: a model
  that predicts the majority class scores 0.5 in balanced accuracy and ROC AUC
  whatever its accuracy.
- Erlotinib has the narrowest AUC distribution of the four drugs (first
  quartile about 0.955 on the full matrix), so its regression task has little
  variance to explain and its classification task separates cell lines that
  differ by a few hundredths of AUC. The 2019 headline of 87.64 percent for
  erlotinib should be read with that in mind.

Reading of the v1.0 tables (section 4) and of the leakage check (section 5):

- **Regression.** The out-of-fold R squared lies between -0.10 and 0.14 for
  every drug and model. The best means per drug are 0.137 (erlotinib, linear
  SVM) and 0.120 (erlotinib, ridge), 0.111 (rapamycin, random forest), 0.095
  (sunitinib, ridge) and 0.053 (paclitaxel, random forest); for sunitinib
  every interval includes zero. The 200 most variable genes therefore explain
  at most about a tenth of the AUC variance of unseen cell lines. This is the
  expected correction of the 2019 values near 1 (erratum 1), not a loss of
  quality: the 2019 "R2" rewarded nearly constant predictions. RMSE and MAE
  track the width of each AUC distribution (erlotinib 0.07, paclitaxel 0.19
  to 0.21) and separate the models very little.
- **Classification.** Accuracy sits at about 0.75 for most cells, the
  majority-class proportion, exactly as in the 2019 tables; ridge predicts
  `resistant` for every line of rapamycin, sunitinib and paclitaxel (balanced
  accuracy 0.500, sensitivity 0) while its ROC AUC is 0.59 to 0.70, so its
  scores rank the lines better than chance but never cross the 0.5
  probability threshold. The models whose balanced accuracy interval excludes
  0.5 are, for erlotinib, the random forest (0.617), the linear SVM (0.599),
  the lasso (0.542) and the elastic net (0.536), and, for sunitinib, the
  random forest (0.543); no model beats the trivial classifier in balanced
  accuracy for rapamycin or paclitaxel. In ROC AUC most models are above 0.5
  with the interval excluding it: the best per drug are the random forest
  for erlotinib (0.706), rapamycin (0.640) and paclitaxel (0.647) and the
  elastic net, lasso, ridge and random forest for sunitinib (0.69 to 0.71).
  The 2019 headline of 87.64 percent accuracy for erlotinib has no
  counterpart: the v1.0 erlotinib accuracies lie between 0.74 and 0.79, the
  highest being the random forest at 0.791 (0.772 to 0.811).
- **Leakage check.** Selecting the 200 genes and the class threshold on all
  cell lines, as in 2019, changes the fold-wise means of the `glmnet` family
  by at most 0.015 in ROC AUC and 0.020 in R squared, in both directions.
  Variance-based selection does not use the response, so a small effect was
  expected; the 2019 procedure was wrong in principle (erratum 6) but, on
  these data and for this family, it did not inflate the estimates in a
  measurable way. The random forest and the SVM were not refitted.
- **Linear SVM on 200 genes.** With the same feature set as the other models
  (erratum 7), the SVM behaves like the penalised models for erlotinib (best
  R squared, second balanced accuracy) but worse than them for the drugs
  with a wide AUC distribution in regression (R squared about zero for
  rapamycin and sunitinib, negative for paclitaxel). The lowest cost of the
  grid (0.01) was chosen in nearly every split, with epsilon 0.2 in
  regression, so the grid edge may be limiting it. In classification it
  gives the highest sensitivity of every family (0.14 to 0.28) at the cost of
  specificity.
- **Hyperparameters.** The random forest mostly selected small `mtry` values
  (5 or 10 in a majority of regression splits) and `num.trees = 500`; the
  elastic net alpha was 0.1 (close to ridge) in most regression splits and
  spread over the grid in classification; `lambda.min` was used everywhere.
  Every selected value is in `outputs/<drug>/<task>/tuning.csv` and every
  evaluated candidate in `tuning_candidates.csv`.

None of these numbers supports a clinical reading; they describe how well
the 2019 design predicts the GDSC AUC of unseen cell lines once the errata are
fixed.
