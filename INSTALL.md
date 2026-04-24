# Installation Guide

## System Requirements

- R (>= 4.1.0)
- Rtools (Windows) or Xcode Command Line Tools (macOS)
- At least 2GB free disk space

## Installation Methods

### Method 1: Install from Pre-built Tarball (Recommended)

1. Download `OmniSciKit_0.1.0.tar.gz`
2. In R, run:

```r
install.packages("path/to/OmniSciKit_0.1.0.tar.gz", repos = NULL, type = "source")
```

### Method 2: Install from GitHub

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install from GitHub
devtools::install_github("yourusername/OmniSciKit")
```

### Method 3: Build from Source

1. Clone the repository:
```bash
git clone https://github.com/yourusername/OmniSciKit.git
```

2. In R:
```r
devtools::install("path/to/OmniSciKit")
```

## Dependencies

### Required Packages
The following packages will be installed automatically:

- ggplot2 (>= 3.3.0)
- dplyr (>= 1.0.0)
- tidyr (>= 1.1.0)
- tibble (>= 3.1.0)
- rlang (>= 0.4.0)
- cli (>= 3.0.0)

### Suggested Packages
For full functionality, you may want to install:

```r
# Bioinformatics
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "limma", "clusterProfiler"))

# Machine Learning
install.packages(c("randomForest", "xgboost", "glmnet", "e1071", "class"))

# Survival Analysis
install.packages(c("survival", "survminer"))

# Mendelian Randomization
install.packages(c("TwoSampleMR", "MendelianRandomization"))
```

## Verification

After installation, verify the package works:

```r
library(OmniSciKit)

# Check version
packageVersion("OmniSciKit")

# Run examples
check_deps()
data(omni_demo)
head(omni_demo)

# Run tests
devtools::test()
```

## Troubleshooting

### Windows Users

If you encounter compilation errors:
1. Install Rtools from https://cran.r-project.org/bin/windows/Rtools/
2. Ensure Rtools is in your PATH

### macOS Users

If you encounter compilation errors:
1. Install Xcode Command Line Tools:
   ```bash
   xcode-select --install
   ```

### Linux Users

Install required system libraries:
```bash
# Ubuntu/Debian
sudo apt-get install build-essential libcurl4-openssl-dev libssl-dev libxml2-dev

# Fedora/RHEL
sudo dnf install R-devel curl-devel openssl-devel libxml2-devel
```

## Getting Help

If you encounter issues:

1. Check the [FAQ](https://github.com/yourusername/OmniSciKit/wiki/FAQ)
2. Search existing [issues](https://github.com/yourusername/OmniSciKit/issues)
3. Open a new issue with:
   - Your operating system
   - R version
   - Error message
   - Minimal reproducible example

## Uninstallation

```r
remove.packages("OmniSciKit")
```
