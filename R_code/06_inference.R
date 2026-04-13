# ============================================================
# 06_inference.R
# Probabilistic inference and qKER analysis (GBN)
# ============================================================

# -------------------------------
# Setup
# -------------------------------

source("packages.R")
source("utils.R")
source("03_bn_structure.R")
source("04_parameter_learning.R")
source("05_validation.R")

# -------------------------------
#  Load CV results from validation step
# -------------------------------
# Expected structure:
# cv_results <- data.frame(
#   ATP_perc = numeric,
#   prob_Col1a1_cpdist = numeric,
#   fold   = numeric
# )
# 
data <- read.csv2("data/my_processed_data.csv")

# Ensure numeric
data <- as.data.frame(lapply(data, as.numeric))


# ============================================================
# 1. Construct biological predictor (cell injury)
# ============================================================

cv_results <- cv_results %>%
  mutate(
    cell_injury = max(ATP_perc, na.rm = TRUE) - ATP_perc
  )


# ============================================================
# 2. qKER visualization (cross-validated probability curve)
# ============================================================

p1 <- ggplot(cv_results,
             aes(x = cell_injury,
                 y = prob_Col1a1_cpdist)) +
  
  # CV scatter (uncertainty across folds)
  geom_point(aes(color = factor(fold)),
             alpha = 0.5,
             size = 2) +
  
  # Smooth biological response curve
  geom_smooth(method = "loess",
              se = TRUE,
              color = "black",
              fill = "gray80") +
  
  labs(
    x = "Cell injury (max ATP - ATP)",
    y = "P(Col1a1 ≥ 1.5)",
    color = "Fold",
    title = "Cross-validated probabilistic response (qKER)",
    subtitle = "Probabilistic response curve across 5-fold CV"
  ) +
  
  theme_minimal(base_size = 14)

print(p1)

# ============================================================
# 4. Fold-level robustness view (optional but useful)
# ============================================================

p2 <- ggplot(cv_results,
             aes(x = cell_injury,
                 y = prob_Col1a1_cpdist)) +
  
  geom_point(alpha = 0.3, size = 1.5) +
  
  geom_smooth(aes(group = fold),
              method = "loess",
              se = FALSE,
              alpha = 0.4) +
  
  labs(
    x = "Cell injury",
    y = "P(Col1a1 ≥ 1.5)",
    title = "Fold-wise variability of qKER"
  ) +
  
  theme_minimal(base_size = 14)

print(p2)

# ============================================================
# 5. Optional: aggregated response curve
# ============================================================

cv_summary <- cv_results %>%
  group_by(cell_injury) %>%
  summarise(
    mean_prob = mean(prob_Col1a1_cpdist, na.rm = TRUE),
    sd_prob   = sd(prob_Col1a1_cpdist, na.rm = TRUE),
    .groups = "drop"
  )

p3 <- ggplot(cv_summary,
             aes(x = cell_injury,
                 y = mean_prob)) +
  
  geom_line(linewidth = 1.2) +
  
  geom_ribbon(aes(ymin = mean_prob - sd_prob,
                  ymax = mean_prob + sd_prob),
              alpha = 0.2) +
  
  labs(
    x = "Cell injury",
    y = "Mean P(Col1a1 ≥ 1.5)",
    title = "Aggregated cross-validated qKER (mean ± SD)"
  ) +
  
  theme_minimal(base_size = 14)

print(p3)

# ============================================================
# 6. Analytical understanding (Mechanistic Insight)
# ============================================================

# 1. Extract Path Coefficients (Effect Sizes)
# In a GBN, these are the 'slopes' of the biological relationships
slopes <- list(
  MIE_to_KE = coef(fitted_bn$aSMA)["ATP_perc"],
  KE_to_AO  = coef(fitted_bn$Col1a1)["aSMA"]
)

# 2. Calculate Total Causal Effect (MIE -> AO)
# For linear chains, this is the product of the path coefficients
total_path_effect <- slopes$MIE_to_KE * slopes$KE_to_AO

# 3. Variance Propagation Analysis
# This tells us how much noise comes from the first step vs the second
var_from_aSMA   <- (slopes$KE_to_AO^2 * fitted_bn$aSMA$sd^2)
var_from_Col1a1 <- fitted_bn$Col1a1$sd^2
total_variance  <- var_from_aSMA + var_from_Col1a1

# 4. Sensitivity (Percentage of AO variance explained by KE)
prop_variance_from_upstream <- var_from_aSMA / total_variance

# --- Output Summary ---
message("--- Mechanistic Summary ---")
message(paste("For every 1 unit drop in ATP, Col1a1 changes by:", round(total_path_effect, 3)))
message(paste("Uncertainty in AO driven by upstream KE:", round(prop_variance_from_upstream * 100, 1), "%"))

# ============================================================
# Notes
# ============================================================

# - Uses cross-validated CPD probabilities (no training bias)
# - Cell injury defined as inverse ATP signal
# - LOESS provides smooth mechanistic response curve
# - Fold coloring shows model stability
# - Suitable for qAOP / qKER visualization
