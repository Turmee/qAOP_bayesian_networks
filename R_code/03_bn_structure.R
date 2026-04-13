# ============================================================
# 03_bn_structure.R
# Bayesian network structure definition (AOP-based DAG)
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
#   time_exposure = numeric,
#   time_exposure_cat = factor (optional grouping variable)
# )
# ------------------------------------------------------------

# -------------------------------
# Load data
# -------------------------------

# Example:

data <- read.csv2("data/my_processed_data.csv")

# ============================================================
# 1. Define Bayesian network structure (DAG)
# ============================================================

# AOP-based causal structure:
# ATP depletion (cell injury) -> HSC activation -> ECM production

dag <- model2network("[ATP_perc][aSMA|ATP_perc][Col1a1|aSMA]")

# Visual inspection (optional)
graphviz.plot(dag, shape = "ellipse", main = "qAOP Structural Model")

# ============================================================
# 2. Biological justification (explicit mapping)
# ============================================================

# Node interpretation:
# ATP_perc  -> hepatocellular energy status (cell injury)
# aSMA      -> hepatic stellate cell activation marker
# Col1a1    -> extracellular matrix production (fibrosis marker)

# Edge justification:
# ATP_perc -> aSMA:
#   injury (ATP depletion) induces stellate cell activation

# aSMA -> Col1a1:
#   activation state drives collagen expression

# This structure reflects:
# - AOP ID: liver fibrosis pathway (AOP38 context)
# - known mechanistic sequence of fibrotic progression

# ============================================================
# 3. Structural object output
# ============================================================

# Store DAG for downstream steps
bn_structure <- dag

# Optional: list form for reproducibility
structure_list <- list(
  nodes = nodes(dag),
  arcs  = arcs(dag)
)

# ============================================================
# Notes
# ============================================================

# - Structure is defined a priori (no data-driven learning)
# - DAG reflects mechanistic AOP assumptions
# - This structure is fixed and across all downstream analysis:
#     04_parameter_learning.R
#     05_validation.R
#     06_inference.R
# - Gaussian assumptions are applied in parameter learning stage