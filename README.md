# Bayesian Network Workflow for qAOP

This repository provides **illustrative R code** for constructing and evaluating Bayesian network (BN) models in the context of **quantitative Adverse Outcome Pathways (qAOPs)**

The code is designed as a **modular workflow template**, covering key steps from preprocessing to probabilistic inference.

**Important**  
- No raw data or results are included  
- Scripts are **not directly executable** without user-provided data  
- The focus is on **methodology and structure**, not reproducibility of published results  

## 📁 Repository Structure
R/

├── 01_preprocessing.R

├── 02_exploration.R

├── 03_bn_structure.R

├── 04_parameter_learning.R

├── 05_validation.R

├── 06_inference.R

├── 07_in_vivo_analysis.R

└── utils.R


Each script represents one step of the analytical workflow.

---

## 🧭 Workflow Overview

### 1. Data preprocessing  
- Outlier detection using Tukey’s rule (1.5 × IQR)  
- Distribution diagnostics:
  - Kullback–Leibler divergence  
  - Q–Q plot correlation  
  - Shapiro–Wilk test  

---

### 2. Exploratory analysis  
- Distributional assessment  
- Correlation analysis with bootstrap confidence intervals  
- Exploration of exposure effects (e.g., MANOVA, regression)

---

### 3. Bayesian network structure  
- A priori structure based on biological knowledge (AOP framework)  
- Example structure used in the workflow: ATPperc → aSMA → COL1A1


---

### 4. Parameter learning  
- Gaussian Bayesian Network (GBN)  
- Parameters estimated via maximum likelihood (`bnlearn`)  
- Each node modeled as a linear function of its parents with Gaussian residuals  

---

### 5. Model validation  
- k-fold cross-validation  
- Performance metrics:
  - R² (posterior predictive correlation)  
  - Normalized Mean Squared Error (NMSE)  
- Calibration via predicted vs observed plots  
- Posterior sampling using `cpdist()`  

---

### 6. Probabilistic inference  
- Likelihood-weighted sampling  
- Conditional probability queries  
- Threshold-based predictions (e.g., low vs high expression)  

---

### 7. In vivo data analysis (template)  
- Wilcoxon rank-sum tests  
- Logistic regression (univariate and multivariate)  
- ROC analysis and threshold selection (Youden’s J)  

---

## 🧠 Key Concepts

**Performance vs Uncertainty**  
- Performance: agreement between predictions and observations  
- Uncertainty: dispersion of posterior predictive distributions  

**Gaussian Bayesian Networks**  
- Linear relationships between variables  
- Normally distributed residuals  

---

## 📦 Dependencies

- R (≥ 4.0)
- `bnlearn`
- `ggplot2`
- `dplyr`
- `boot`
- `pROC`

---

## 📌 Usage

This repository is intended as a **template**.

Typical workflow:

```r
source("R/01_preprocessing.R")
source("R/03_bn_structure.R")
source("R/04_parameter_learning.R")
source("R/05_validation.R")
source("R/06_inference.R")
```

If you use or adapt this workflow, please cite.
