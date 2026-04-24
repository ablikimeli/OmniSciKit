# OmniSciKit

[![R-CMD-check](https://github.com/yourusername/OmniSciKit/workflows/R-CMD-check/badge.svg)](https://github.com/yourusername/OmniSciKit/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

OmniSciKit is a comprehensive R package for scientific research, providing a unified framework for:

- **Bioinformatics**: Differential expression, pathway enrichment, omics data analysis
- **Mendelian Randomization**: Causal inference methods with genetic instruments
- **Machine Learning**: Prediction, classification, and model evaluation
- **Epidemiology**: Infectious disease modeling and outbreak analysis
- **Statistics**: Comprehensive statistical tests and regression analysis
- **Sample Size**: Power and sample size calculations for study design

## Installation

### From GitHub

```r
# install.packages("devtools")
devtools::install_github("yourusername/OmniSciKit")
```

### From Source

```r
# Download OmniSciKit_0.1.0.tar.gz
install.packages("path/to/OmniSciKit_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Quick Start

```r
library(OmniSciKit)

# Check dependencies
check_deps()

# Load demo data
data(omni_demo)
head(omni_demo)

# Statistical analysis
result <- stats_ttest(omni_demo$baseline_score, omni_demo$follow_up_score)
print(result)

# Sample size calculation
ss_ttest(delta = 5, sd = 10, power = 0.8)

# Machine learning
model <- ml_logistic(omni_ml, outcome ~ .)
summary(model)
```

## Documentation

- [Installation Guide](vignettes/installation.Rmd)
- [Bioinformatics Analysis](vignettes/omics.Rmd)
- [Mendelian Randomization](vignettes/mr.Rmd)
- [Machine Learning](vignettes/ml.Rmd)
- [Epidemiological Modeling](vignettes/epi.Rmd)
- [Statistical Analysis](vignettes/stats.Rmd)
- [Sample Size Calculation](vignettes/ss.Rmd)

## Features

### 1. Bioinformatics Module (`omics`)
- Differential expression analysis (DESeq2, edgeR, limma)
- Pathway enrichment (GO, KEGG, GSEA)
- Quality control and normalization
- GEO and TCGA data download
- Visualization (volcano plots, heatmaps, PCA)

### 2. Mendelian Randomization Module (`mr`)
- IVW, Egger, Weighted Median methods
- Sensitivity analysis (leave-one-out, radial)
- Pleiotropy and heterogeneity tests
- Power calculation
- Forest and funnel plots

### 3. Machine Learning Module (`ml`)
- Logistic regression, Random Forest, XGBoost
- SVM, Naive Bayes, KNN
- Cross-validation and hyperparameter tuning
- Feature selection and importance
- Model evaluation (ROC, calibration)

### 4. Epidemiology Module (`epi`)
- SIR, SEIR, SIRS, SEIRS models
- R0 and Rt estimation
- Growth rate and doubling time
- Vaccine impact simulation
- Outbreak detection

### 5. Statistics Module (`stats`)
- T-tests, ANOVA, non-parametric tests
- Regression analysis (linear, logistic, Cox)
- Survival analysis (KM curves, log-rank)
- Diagnostic test evaluation
- Bootstrap and permutation

### 6. Sample Size Module (`ss`)
- T-test, ANOVA, Chi-square
- Correlation and regression
- Survival analysis
- Non-inferiority and equivalence
- Adaptive designs

## Dependencies

### Required
- R (>= 4.1.0)
- ggplot2, dplyr, tidyr, tibble
- stats, utils, graphics, grDevices, methods

### Suggested
- DESeq2, edgeR, limma
- survival, survminer
- randomForest, xgboost, glmnet
- clusterProfiler, org.Hs.eg.db
- TwoSampleMR, MendelianRandomization

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

## License

This package is licensed under the MIT License. See [LICENSE](LICENSE) file for details.

## Citation

```r
citation("OmniSciKit")
```

## Contact

- Author: Ablikim Ali
- Email: ablikim.ali@example.com
- GitHub: https://github.com/yourusername/OmniSciKit

## Acknowledgments

This package was developed with support from the scientific community and builds upon many excellent R packages.
