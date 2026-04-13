# ============================================================
# utils.R: Shared functions for qAOP BN Workflow
# ============================================================

# ------------------------------------------------------------
# 1. Statistical Diagnostics & Preprocessing
# ------------------------------------------------------------

# Detect outliers using Tukey's rule
detect_outliers <- function(x) {
  if (!is.numeric(x)) return(rep(FALSE, length(x)))
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_val
  upper_bound <- Q3 + 1.5 * IQR_val
  return(x < lower_bound | x > upper_bound)
}

# NMSE Calculation
nmse <- function(obs, pred) {
  # Handles potential NAs in validation
  valid <- !is.na(obs) & !is.na(pred)
  mean((obs[valid] - pred[valid])^2) / var(obs[valid])
}

# ------------------------------------------------------------
# 2. Distributional Assessment (Normality)
# ------------------------------------------------------------

normality_score <- function(x, n_bins = 100, plot = TRUE, name = NULL) {
  library(e1071)
  library(nortest)
  
  var_name <- if (is.null(name)) deparse(substitute(x)) else name
  x <- na.omit(as.numeric(x))
  n <- length(x)
  
  # Divergence Calculations
  m <- mean(x); s <- sd(x)
  hist_data <- hist(x, breaks = n_bins, plot = FALSE)
  p <- hist_data$density
  q <- dnorm(hist_data$mids, mean = m, sd = s)
  bin_width <- diff(hist_data$breaks[1:2])
  
  valid <- p > 0 & q > 0
  p <- p[valid]; q <- q[valid]
  kl_div <- sum(p * log(p / q)) * bin_width
  m_mix  <- 0.5 * (p + q)
  js_div <- 0.5 * sum(p * log(p / m_mix)) * bin_width + 0.5 * sum(q * log(q / m_mix)) * bin_width
  
  # Visualization
  if (plot) {
    par(mfrow = c(1, 2))
    h <- hist(x, breaks = n_bins, freq = FALSE, col = "lightgrey", border = "white",
              main = paste0("JS Div = ", round(js_div, 3)), xlab = var_name)
    lines(density(x), col = "red", lwd = 2)
    curve(dnorm(x, m, s), add = TRUE, col = "black", lwd = 2)
    
    qqnorm(x, pch = 21, bg = "lightgrey", col = "darkgrey")
    qqline(x, col = "red", lwd = 2)
    par(mfrow = c(1, 1))
  }
  
  return(data.frame(variable = var_name, n = n, js_divergence = js_div, 
                    shapiro_p = if(n <= 5000) shapiro.test(x)$p.value else NA))
}

# ------------------------------------------------------------
# 3. BN Specific Utilities
# ------------------------------------------------------------

# Extract coefficients from a specific node in a fitted BN
extract_params <- function(fitted_bn, node) {
  return(coefficients(fitted_bn[[node]]))
}

# Bootstrap correlation with confidence intervals
get_corr_boot <- function(df, xvar, yvar, nboot = 1000) {
  df <- df %>% filter(!is.na(.data[[xvar]]), !is.na(.data[[yvar]]))
  if (nrow(df) < 5) return(c(NA, NA, NA))
  
  boot_r <- replicate(nboot, {
    s <- df[sample(nrow(df), replace = TRUE), ]
    cor(s[[xvar]], s[[yvar]], use = "complete.obs")
  })
  return(c(mean = mean(boot_r, na.rm=T), quantile(boot_r, c(0.025, 0.975), na.rm=T)))
}

# ------------------------------------------------------------
# 4. In Vivo Validation Utilities
# ------------------------------------------------------------

# Arranging plots with one legend
grid_arrange_shared_legend <- function(plot_list, ncol = 2) {
  plots <- lapply(plot_list, function(x) x + theme(legend.position = "none"))
  legend <- cowplot::get_legend(plot_list[[1]] + theme(legend.position = "right"))
  gridExtra::grid.arrange(
    do.call(gridExtra::arrangeGrob, c(plots, ncol = ncol)),
    legend,
    ncol = 2,
    widths = c(5, 1)
  )
}

# Unified plotting for Fold Change vs Binary Damage
plot_fc_vs_damages <- function(data, damage_vars) {
  plot_list <- list()
  col_map <- c("no" = "grey40", "yes" = "firebrick")
  
  for (damage in damage_vars) {
    for (gene in c("acta2_fc", "col1a1_fc")) {
      df <- data %>% 
        filter(!is.na(.data[[gene]]), !is.na(.data[[damage]])) %>%
        mutate(Damage = factor(.data[[damage]], levels = c("no", "yes")))
      
      p <- ggplot(df, aes(x = .data[[gene]], y = as.numeric(Damage) - 1, color = Damage)) +
        geom_jitter(height = 0.05, width = 0, alpha = 0.7) +
        geom_vline(xintercept = 1.5, linetype = "dashed", color = "firebrick") +
        scale_color_manual(values = col_map) +
        scale_y_continuous(breaks = c(0, 1), labels = c("no", "yes")) +
        labs(x = gene, y = damage, title = paste(gene, "vs", damage)) +
        theme_minimal()
      
      plot_list[[paste(damage, gene, sep = "_")]] <- p
    }
  }
  grid_arrange_shared_legend(plot_list, ncol = 2)
}