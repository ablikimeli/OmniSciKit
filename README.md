# OmniSciKit

Omni-Science Toolkit for Biomedical Research

## Description

A comprehensive R package for biomedical research covering:
- Statistical Analysis (stats_*)
- Sample Size Calculation (ss_*)
- Mendelian Randomization (mr_*)
- Machine Learning (ml_*)
- Epidemiology (epi_*)
- Omics Analysis (omics_*)

## Installation

`
install.packages("OmniSciKit_0.1.0.tar.gz", repos = NULL, type = "source")
`

## Usage

`
library(OmniSciKit)

# Sample size calculation
result <- ss_ttest(delta = 0.5, sd = 1, power = 0.8)
print(result\)

# Statistical test
x <- rnorm(100)
result <- stats_ttest(x, mu = 0)
print(result\)
`

## Author

Ablikim Ali <medical_stat@163.com>

## License

MIT