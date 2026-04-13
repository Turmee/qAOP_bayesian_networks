# ============================================================
# 02_exploration.R
# Data exploration template for Bayesian network workflow
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
# Replace `data` with your processed dataset.

# ------------------------------------------------------------
# Expected structure:
# data <- data.frame(
#   ATP_perc = numeric,
#   aSMA     = numeric,
#   Col1a1   = numeric,
#   time_exposure = numeric
# )
# ------------------------------------------------------------

# -------------------------------
# Load processed data
# -------------------------------

# Example:
data <- read.csv2("data/my_processed_data.csv")

vars <- c("ATP_perc", "aSMA", "Col1a1")

# ============================================================
# 1. Distribution diagnostics
# ============================================================

# ATP
normality_score(data$ATP_perc, name = "ATP_perc")

# aSMA
normality_score(data$aSMA, name = "aSMA")

# COL1A1
normality_score(data$Col1a1, name = "Col1a1")

# Optional transformation view (biological interpretation)
normality_score(max(data$ATP_perc, na.rm = TRUE) - data$ATP_perc, name = "Cell injury (max(ATP) - ATP)")

# ============================================================
# 2. Correlation analysis (Pearson + bootstrap CI)
# ============================================================

# ATP vs aSMA
corr_ATP_aSMA <- get_corr_boot(data, "ATP_perc", "aSMA", nboot = 10000)

# ATP vs Col1a1
corr_ATP_COL  <- get_corr_boot(data, "ATP_perc", "Col1a1", nboot = 10000)

# aSMA vs Col1a1
corr_aSMA_COL <- get_corr_boot(data, "aSMA", "Col1a1", nboot = 10000)

# ============================================================
# 3. Optional: MANOVA-based threshold exploration
# ============================================================

# Exposure time column must exist in dataset (if applicable)

# Example:
 best_split <- find_best_multivar_split(
   data = data,
   time_col = "time_exposure",
   param_cols = vars)

 print(best_split)

 # Add group label to data
 data$time_exposure_cat <- ifelse(data$time_exposure <= best_split$best_threshold, "short", "long")
 
 ggplot(data, aes(x = ATP_perc, y = aSMA)) + 
   geom_point(size = 2, aes(color = time_exposure_cat)) +
   geom_smooth(aes(color = time_exposure_cat), method = "lm", se = TRUE)
   
 
# ============================================================
# Notes
# ============================================================

# - Distribution diagnostics assess Gaussian assumptions for BN modeling
# - Correlation structure informs network plausibility
# - Bootstrap CI quantifies uncertainty in correlations
# - MANOVA step is optional and depends on available metadata
# - This script is exploratory and does not define the BN structure