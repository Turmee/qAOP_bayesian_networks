# ============================================================
# packages.R: Dependency Management
# ============================================================

required_packages <- c(
  "dplyr",      # Data manipulation
  "tidyr",      # Data reshaping
  "ggplot2",    # Visualization
  "bnlearn",    # Bayesian Network modeling
  "e1071",      # Statistical moments (skewness/kurtosis)
  "nortest",    # Normality tests
  "pROC",       # ROC and AUC analysis
  "broom",      # Tidy model outputs
  "gridExtra",  # Plot arrangement
  "cowplot"     # Legend handling
)

# Identify missing packages
missing <- required_packages[!required_packages %in% rownames(installed.packages())]

# Prompt for installation if missing
if (length(missing) > 0) {
  message("The following packages are missing: ", paste(missing, collapse = ", "))
  # Optional: Uncomment the next line to allow auto-installation
  # install.packages(missing)
}

# Load all packages quietly
suppressPackageStartupMessages(
  invisible(lapply(required_packages, library, character.only = TRUE))
)

message("All dependencies loaded successfully.")