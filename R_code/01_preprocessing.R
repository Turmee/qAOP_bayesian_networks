# ============================================================
# 01_preprocessing.R
# Data preprocessing template for Bayesian network workflow
# ============================================================

# -------------------------------
# Setup
# -------------------------------

source("packages.R")
source("utils.R")

# ============================================================
# Example workflow (template)
# ============================================================

# The following is an illustrative workflow.
# Replace `data` with your own dataset.

# ------------------------------------------------------------
# Example expected structure:
# data <- data.frame(
#   ATP_perc = numeric,
#   aSMA     = numeric,
#   Col1a1   = numeric
# )
# ------------------------------------------------------------

# -------------------------------
# Load data
# -------------------------------

# Example (replace with your path)
data <- read.csv2("data/my_data.csv")

# Variable names used in the workflow
vars <- c("ATP_perc", "aSMA", "Col1a1")

# --- Outlier detection
data$aSMA_outlier <- detect_outliers(data$aSMA)
data$ATP_outlier <- detect_outliers(data$ATP_perc)
data$Col1a1_outlier <- detect_outliers(data$Col1a1)


# dataset without outliers 
data_no_outliers <- data[(data$ATP_outlier == FALSE | is.na(data$ATP_outlier) == TRUE) & (data$aSMA_outlier == FALSE | is.na(data$aSMA_outlier) == TRUE) & (data$Col1a1_outlier == FALSE | is.na(data$Col1a1_outlier) == TRUE),]

# -------------------------------
# Shuffle data
# -------------------------------

set.seed(123)  # for reproducibility
data_out_shuffle <- data_out[sample(nrow(data_out)), ]

# ------------------------------------------------------------
# Output datasets (example)
# ------------------------------------------------------------
# datasets <- list(
#   raw           = data,
#   no_outliers   = data_no_outliers,
#   shuffled      = data_out_shuffle
# )


# ============================================================
# Notes
# ============================================================

# - Outlier detection uses Tukey’s rule (1.5 × IQR)
# - Threshold variables are used for probabilistic classification
# - Multiple dataset variants can be used to assess robustness
# - This script is data-agnostic and requires user-provided input