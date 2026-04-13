# ============================================================
# 04_parameter_learning.R
# Parameter learning for Gaussian Bayesian Network (GBN)
# ============================================================

# -------------------------------
# Setup
# -------------------------------

source("packages.R")
source("utils.R")
source("03_bn_structure.R")

# ============================================================
# Example workflow (template)
# ============================================================

# The following is an illustrative workflow.
# Replace `data` with your processed dataset.

# ------------------------------------------------------------
# Expected structure:
# data <- data.frame(
#   ATP_perc = numeric,
#   aSMA     = numeric,
#   Col1a1   = numeric,
# )
# ------------------------------------------------------------

# -------------------------------
# Load data
# -------------------------------

# Example:
data <- read.csv2("data/my_processed_data.csv")

# Extract only the variables used in the DAG
dag_vars <- nodes(bn_structure)
data_fitted <- data[, dag_vars]

# Ensure numeric and remove NAs (MLE fitting requires complete cases)
data_fitted <- as.data.frame(lapply(data_fitted, as.numeric))
data_complete <- na.omit(data_fitted)


# ============================================================
# 1. Fit Gaussian Bayesian Network (MLE)
# ============================================================

# Structure is defined in 03_bn_structure.R as 'bn_structure'

fitted_bn <- bn.fit(
  bn_structure,
  data = data_complete,
  method = "mle"
)

# ============================================================
# 2. Inspect learned parameters (optional)
# ============================================================

# Regression: aSMA ~ ATP_perc
# regression coefficients represent linear Gaussian relationships

# Example:
 fitted_bn$aSMA

# Regression: Col1a1 ~ aSMA
 fitted_bn$Col1a1

# ============================================================
# 3. Extract coefficients (optional utility view)
# ============================================================

 extract_params(fitted_bn, "aSMA")
 extract_params(fitted_bn, "Col1a1")

# ============================================================
# 4. Store model object
# ============================================================

bn_fitted <- fitted_bn

 
# ============================================================
# 5. Model assumption checks
# ============================================================
 # Save old par settings
 old_par <- par(no.readonly = TRUE)
 par(mfrow = c(2,2))
 
 # Use the fitted object directly to extract the "lm" flavor
 # This ensures diagnostics match the EXACT model in the BN
 plot(as.lm(fitted_bn$aSMA), main = "aSMA Diagnostics")
 plot(as.lm(fitted_bn$Col1a1), main = "Col1a1 Diagnostics")
 
 # Reset plotting window
 par(old_par)

# ============================================================
# Notes
# ============================================================

# - Parameters estimated using Maximum Likelihood Estimation (MLE)
# - Gaussian BN assumes linear relationships with Gaussian residuals
# - Each node is modeled conditionally on its parents
# - No validation or prediction performed in this step
# - Output is used in:
#     05_validation.R
#     06_inference.R