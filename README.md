# Bayesian Network Workflow for qAOP

This repository provides illustrative R code for constructing and evaluating
Bayesian network (BN) models in quantitative Adverse Outcome Pathways (qAOPs).
It is a modular workflow template from preprocessing through probabilistic
inference and external-group validation.

## Scope

- No raw data, derived data, manuscript results, or fitted models are included.
- Scripts require a user-provided data set and are not directly executable.
- Example names and thresholds are placeholders, not reported results.
- The focus is methodology and reusable pipeline structure.

## Repository structure

```text
R_code/
|-- 01_preprocessing.R
|-- 02_exploration.R
|-- 03_bn_structure.R
|-- 04_parameter_learning.R
|-- 05_validation.R
|-- 06_inference.R
|-- 07_in_vivo_analysis.R
|-- 08_resampling_validation.R
|-- 09_grouped_validation.R
|-- packages.R
|-- packages_list.txt
`-- utils.R
```

## Workflow

1. Preprocessing and distribution diagnostics.
2. Exploratory association and exposure-effect analyses.
3. Biologically informed DAG definition.
4. Maximum-likelihood learning of Gaussian BN parameters.
5. Basic k-fold validation with R-squared, NMSE, AUC, and calibration plots.
6. Likelihood-weighted inference and threshold-based probability queries.
7. In vivo analysis template using rank tests, logistic models, and ROC curves.
8. Reusable k-fold prediction and bootstrap uncertainty summaries.
9. Leave-one-study/chemical/group-out validation and label-permutation testing.

## Expected data interface

Users supply a complete-case data frame with one column per DAG node. Grouped
validation also requires a grouping column such as `Study` or `Chemical`.
Optional metadata columns can be retained in prediction outputs.

## Dependencies

R (>= 4.0), `bnlearn`, and `pROC` are central to the validation pipeline.
Additional packages are listed in `R_code/packages_list.txt`.

## Illustrative usage

```r
source("R_code/08_resampling_validation.R")
source("R_code/09_grouped_validation.R")

dag <- bnlearn::model2network(
  "[Exposure][KeyEvent1|Exposure][KeyEvent2|KeyEvent1]"
)
targets <- c("KeyEvent1", "KeyEvent2")
thresholds <- c(KeyEvent1 = 1.5, KeyEvent2 = 1.5)

cv_results <- bn_kfold_cv(
  my_data, dag, targets, thresholds,
  metadata = c("Study", "Chemical"), seed = 123
)
calculate_bn_metrics(cv_results, targets)

grouped_results <- leave_one_group_out_cv(
  my_data, dag, targets, "Study", thresholds
)
```

The permutation pipeline may be computationally expensive because each
permutation refits the BN and uses posterior sampling. Select repetitions and
sample counts appropriate to the intended analysis.

If you use or adapt this workflow, please cite the associated publication.
