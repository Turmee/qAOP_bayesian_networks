# ============================================================
# 07_in_vivo_analysis.R
# In vivo analysis: Fold change analysis, ROC threshold (Youden’s J), Logistic + Wilcoxon analysis, Visualization# ============================================================

# -------------------------------
# Setup
# -------------------------------

source("packages.R")
source("utils.R")

# ============================================================
# 1. Load in vivo data
# ============================================================

data_in_vivo <- read.csv2("data/in_vivo_data.csv")

# Expected columns (adapt if needed):
# - acta2_fc
# - col1a1_fc
# - Fibrosis (0/1 or no/yes)
# - exposure_type

# Ensure correct types
data_in_vivo$acta2_fc  <- as.numeric(data_in_vivo$acta2_fc)
data_in_vivo$col1a1_fc <- as.numeric(data_in_vivo$col1a1_fc)

# Convert fibrosis to binary if needed
data_in_vivo$Fibrosis <- ifelse(data_in_vivo$Fibrosis %in% c("yes", 1), 1, 0)

# ============================================================
# 2. ROC analysis (Youden threshold)
# ============================================================

roc_obj <- roc(data_in_vivo$Fibrosis, data_in_vivo$col1a1_fc)

youden_coords <- coords(
  roc_obj, 
  x = "best", 
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity")
)

# Store threshold
col1a1_threshold <- as.numeric(youden_coords[1, "threshold"])

message(paste("Optimal In Vivo Threshold (Youden):", round(col1a1_threshold, 2)))

# ============================================================
# 3. Run fold-change analysis
# ============================================================

results_fc <- analyze_fc_combined(
  data = data_in_vivo,
  genes = c("acta2_fc", "col1a1_fc"),
  necrosis_vars = c("Fibrosis"),
  pval_threshold = 0.05,
  plot_results = TRUE
)

# ============================================================
# 4. Plot FC vs damage
# ============================================================

plot_fc_vs_damages(
  data_in_vivo,
  damage_vars = c("Fibrosis")
)

# ============================================================
# Notes
# ============================================================

# - ROC threshold is derived independently from BN
# - Logistic regression evaluates association with damage endpoints
# - Wilcoxon tests provide non-parametric confirmation
# - Plots show distribution + threshold alignment
# - Results can be compared to in vitro BN predictions