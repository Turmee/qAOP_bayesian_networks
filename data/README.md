# Data

This directory contains the two derived datasets used in the liver-fibrosis
analyses. They are shared to make the analytical inputs transparent. The
repository's R scripts remain illustrative workflow examples rather than a
complete reproduction package.

## `liver_fibrosis_bn_analysis_data.rds`

This R data frame is the final complete-case dataset used for the Gaussian
Bayesian network, k-fold cross-validation, bootstrap evaluation, and
leave-one-study/chemical-out analyses.

- **Dimensions:** 26 rows and 21 columns
- **BN variables:** `ATP_perc`, `aSMA`, and `Col1a1`
- **Model structure:** `ATP_perc -> aSMA -> Col1a1`
- **Grouping metadata:** `Chemical` and `Study_ID`
- **Coverage:** six chemicals and four studies contribute complete cases

### Preprocessing

The starting literature-curated table contained 99 rows. Outliers were
identified separately for `ATP_perc`, `aSMA`, and `Col1a1` using Tukey's rule:
values below Q1 - 1.5 x IQR or above Q3 + 1.5 x IQR were flagged. A row was
removed if any available BN endpoint was flagged as an outlier; missing values
were retained during this screening. This left 87 rows.

`ATP_perc` was converted to numeric. Binary high-response indicators were
defined as values strictly greater than a fold-change threshold of 1.5.
`Study_ID` was generated from the
source DOI as anonymized factor labels (`Study_1`, `Study_2`, and so on).
Finally, rows missing any of `ATP_perc`, `aSMA`, or `Col1a1` were removed,
leaving 26 rows with complete observations for these three BN variables. Other
supporting columns retained in the dataset may contain missing values.

The RDS also retains provenance and supporting columns from preprocessing,
including DOI, chemical, exposure, concentration, gene labels, extracellular
matrix measurements, and endpoint-specific outlier flags. See the file itself
for the complete column set.

## `tg_gates.csv`

This table was used for the complementary in vivo analysis linking ACTA2 and
COL1A1 responses to liver pathology annotations.

- **Dimensions:** 3,528 observations and 22 columns (including the exported row
  index)
- **Exposure information:** experiment, concentration, duration, and single or
  repeated exposure
- **Molecular endpoints:** ACTA2 and COL1A1 log2 fold changes, fold changes,
  p-values, and adjusted p-values
- **Pathology indicators:** single-cell necrosis, necrosis, fibrinoid necrosis,
  and fibrosis

### Preprocessing

Blank entries in the four pathology columns were interpreted as missing values,
and the pathology variables were treated as categorical variables. Fold-change
columns were checked against the corresponding log2 fold changes using
`fold change = 2^(log2 fold change)`. Statistical-significance indicators used
in downstream analyses were derived from the adjusted p-value columns and are
not stored as additional columns in this file.

The project files available with this repository do not document the complete
upstream extraction and aggregation procedure that produced this TG-GATES
table. The description above is therefore limited to preprocessing that could
be verified from the analysis scripts.

## Missing values and file formats

Missing values in `tg_gates.csv` are represented as `NA`. The Bayesian-network
dataset is stored as an RDS file to preserve R column classes and factor levels;
load it with `readRDS()`.

## File integrity

SHA-256 checksums of the published files:

```text
506b00ddedf8fe6409453d3d15f9fccf84290b3f59792920c095ec90d08f9bd9  liver_fibrosis_bn_analysis_data.rds
f8e70ffe03451332c23df35cfa30db7b3adfb4c73a30809b9aa25e0f5e157b13  tg_gates.csv
```

## Citation

When reusing these data, please cite the original studies identified in the
data and the associated publication:

Durnik R, Juchelkova T, Hecht H, Winkelman LMT, Beltman JB, Coumoul X,
Jornod F, Audouze K, Blaha L, Bajard L. (2026). Towards Bayesian-based
quantitative adverse outcome pathways using in vitro data from open literature
and continuous variables: a case example for liver fibrosis. *Toxicological
Sciences*, 209(8). https://doi.org/10.1093/toxsci/kfag090
