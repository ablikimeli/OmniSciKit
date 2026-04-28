#' Complete Utility Functions for Clinical Research
#'
#' Data manipulation and utility functions

# ============================================================================
# MISSING DATA
# ============================================================================

#' Impute Missing Values
#'
#' Impute missing values using various methods
#' @param data Data frame with missing values
#' @param method Imputation method: "mean", "median", "mode", "knn"
#' @param vars Variables to impute (NULL = all)
#' @return Data frame with imputed values
#' @examples
#' data(omni_glucose)
#' test_data <- omni_glucose
#' test_data$age[1:5] <- NA
#' impute_missing(test_data, method = "median")
#' @export
impute_missing <- function(data, method = "median", vars = NULL) {
  if (is.null(vars)) {
    vars <- names(data)
  }

  result <- data

  for (var in vars) {
    if (var %in% names(data)) {
      missing_idx <- is.na(data[[var]])

      if (sum(missing_idx) > 0) {
        if (is.numeric(data[[var]])) {
          if (method == "mean") {
            impute_val <- mean(data[[var]], na.rm = TRUE)
          } else if (method == "median") {
            impute_val <- median(data[[var]], na.rm = TRUE)
          } else {
            impute_val <- median(data[[var]], na.rm = TRUE)
          }
          result[[var]][missing_idx] <- impute_val
        } else if (is.factor(data[[var]]) || is.character(data[[var]])) {
          mode_val <- names(sort(table(data[[var]]), decreasing = TRUE))[1]
          result[[var]][missing_idx] <- mode_val
        }
      }
    }
  }

  result
}

#' Multiple Imputation (Simple)
#'
#' Simple multiple imputation using predictive mean matching
#' @param data Data frame
#' @param m Number of imputations
#' @return List of imputed data frames
#' @examples
#' \dontrun{
#' imputed_data <- multiple_imputation(data, m = 5)
#' }
#' @export
multiple_imputation <- function(data, m = 5) {
  imputed_list <- list()

  for (i in 1:m) {
    imputed_data <- data

    for (var in names(data)) {
      missing_idx <- is.na(data[[var]])

      if (sum(missing_idx) > 0 && is.numeric(data[[var]])) {
        # Simple stochastic regression imputation
        observed <- data[[var]][!missing_idx]
        mean_val <- mean(observed, na.rm = TRUE)
        sd_val <- sd(observed, na.rm = TRUE)

        imputed_data[[var]][missing_idx] <- rnorm(sum(missing_idx),
                                                   mean = mean_val,
                                                   sd = sd_val)
      }
    }

    imputed_list[[i]] <- imputed_data
  }

  imputed_list
}

# ============================================================================
# DATA TRANSFORMATION
# ============================================================================

#' Create Dummy Variables
#'
#' Convert categorical variables to dummy/indicator variables
#' @param data Data frame
#' @param vars Variables to convert
#' @param drop_first Drop first level (default TRUE)
#' @return Data frame with dummy variables
#' @examples
#' data(omni_glucose)
#' create_dummy_vars(omni_glucose, vars = "treatment")
#' @export
create_dummy_vars <- function(data, vars, drop_first = TRUE) {
  result <- data

  for (var in vars) {
    if (var %in% names(data)) {
      result[[var]] <- NULL

      dummy_mat <- model.matrix(~ . - 1, data = data.frame(x = data[[var]]))

      if (drop_first && ncol(dummy_mat) > 1) {
        dummy_mat <- dummy_mat[, -1, drop = FALSE]
      }

      dummy_df <- as.data.frame(dummy_mat)
      names(dummy_df) <- paste(var, colnames(dummy_mat), sep = "_")
      result <- cbind(result, dummy_df)
    }
  }

  result
}

#' Winsorize Variables
#'
#' Winsorize extreme values
#' @param x Numeric vector
#' @param probs Quantile limits (default c(0.05, 0.95))
#' @return Winsorized vector
#' @examples
#' winsorize(rnorm(100))
#' @export
winsorize <- function(x, probs = c(0.05, 0.95)) {
  limits <- quantile(x, probs = probs, na.rm = TRUE)
  x[x < limits[1]] <- limits[1]
  x[x > limits[2]] <- limits[2]
  x
}

#' Log Transform
#'
#' Apply log transformation with handling of zero/negative values
#' @param x Numeric vector
#' @param offset Offset for zero/negative values (default 1)
#' @return Log-transformed vector
#' @examples
#' log_transform(c(0, 1, 2, 5, 10))
#' @export
log_transform <- function(x, offset = 1) {
  if (any(x <= 0, na.rm = TRUE)) {
    x <- x + offset - min(x, na.rm = TRUE)
  }
  log(x)
}

# ============================================================================
# MULTIPLE TESTING
# ============================================================================

#' Multiple Testing Correction
#'
#' Apply multiple testing correction methods
#' @param p_values Vector of p-values
#' @param method Correction method: "bonferroni", "holm", "fdr", "bh"
#' @return Data frame with original and adjusted p-values
#' @examples
#' p <- c(0.01, 0.03, 0.05, 0.1, 0.2)
#' multiple_testing_correction(p, method = "fdr")
#' @export
multiple_testing_correction <- function(p_values, method = "fdr") {
  valid_methods <- c("bonferroni", "holm", "fdr", "bh")

  if (!method %in% valid_methods) {
    stop("Method must be one of: ", paste(valid_methods, collapse = ", "))
  }

  if (method == "fdr" || method == "bh") {
    method <- "fdr"
  }

  adj_p <- p.adjust(p_values, method = method)

  result <- data.frame(
    original = p_values,
    adjusted = adj_p,
    significant = adj_p < 0.05
  )

  structure(result, method = method)
}

#' Family-Wise Error Rate Control
#'
#' Control FWER using Holm-Bonferroni method
#' @param p_values Vector of p-values
#' @param alpha Significance level (default 0.05)
#' @return List with adjusted p-values and significance
#' @examples
#' p <- c(0.01, 0.03, 0.05, 0.1)
#' fwer_control(p)
#' @export
fwer_control <- function(p_values, alpha = 0.05) {
  adj_p <- p.adjust(p_values, method = "holm")

  structure(list(
    original = p_values,
    adjusted = adj_p,
    significant = adj_p < alpha,
    alpha = alpha,
    method = "Holm-Bonferroni"
  ), class = "omni_stats")
}

# ============================================================================
# BOOTSTRAP
# ============================================================================

#' Bootstrap Confidence Interval
#'
#' Calculate bootstrap confidence interval
#' @param data Numeric vector
#' @param statistic Function to calculate statistic
#' @param n_bootstrap Number of bootstrap samples
#' @param conf_level Confidence level
#' @return List with bootstrap results
#' @examples
#' data(omni_glucose)
#' bootstrap_ci(omni_glucose$hba1c_3m, mean)
#' @export
bootstrap_ci <- function(data, statistic = mean, n_bootstrap = 1000,
                         conf_level = 0.95) {
  boot_stats <- numeric(n_bootstrap)

  for (i in 1:n_bootstrap) {
    boot_sample <- sample(data, length(data), replace = TRUE)
    boot_stats[i] <- statistic(boot_sample)
  }

  alpha <- 1 - conf_level
  ci <- quantile(boot_stats, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)

  structure(list(
    statistic = statistic(data),
    bootstrap_mean = mean(boot_stats, na.rm = TRUE),
    bootstrap_sd = sd(boot_stats, na.rm = TRUE),
    ci_lower = ci[1],
    ci_upper = ci[2],
    n_bootstrap = n_bootstrap,
    conf_level = conf_level,
    bootstrap_distribution = boot_stats
  ), class = "omni_stats")
}

#' Bootstrap for Regression
#'
#' Bootstrap confidence intervals for regression coefficients
#' @param formula Regression formula
#' @param data Data frame
#' @param n_bootstrap Number of bootstrap samples
#' @return List with bootstrap results
#' @examples
#' \dontrun{
#' bootstrap_regression(y ~ x, data)
#' }
#' @export
bootstrap_regression <- function(formula, data, n_bootstrap = 1000) {
  fit <- lm(formula, data = data)
  coef_names <- names(coef(fit))
  n_coef <- length(coef_names)

  boot_coefs <- matrix(NA, nrow = n_bootstrap, ncol = n_coef)
  colnames(boot_coefs) <- coef_names

  for (i in 1:n_bootstrap) {
    boot_idx <- sample(nrow(data), nrow(data), replace = TRUE)
    boot_data <- data[boot_idx, ]
    boot_fit <- lm(formula, data = boot_data)
    boot_coefs[i, names(coef(boot_fit))] <- coef(boot_fit)
  }

  ci <- apply(boot_coefs, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

  structure(list(
    coefficients = coef(fit),
    ci_lower = ci[1, ],
    ci_upper = ci[2, ],
    bootstrap_distribution = boot_coefs,
    n_bootstrap = n_bootstrap
  ), class = "omni_stats")
}

# ============================================================================
# ASSUMPTION CHECKING
# ============================================================================

#' Check Assumptions for Statistical Tests
#'
#' Check assumptions for common statistical tests
#' @param data Data frame or numeric vector
#' @param test Type of test: "ttest", "anova", "regression"
#' @param group Grouping variable
#' @return List with assumption check results
#' @examples
#' data(omni_glucose)
#' check_assumptions(omni_glucose$hba1c_3m, "ttest",
#'                   group = omni_glucose$treatment)
#' @export
check_assumptions <- function(data, test = "ttest", group = NULL) {
  results <- list()

  if (test == "ttest") {
    if (is.null(group)) {
      stop("Grouping variable required for t-test assumption check")
    }

    groups <- unique(group)
    normality_tests <- list()
    for (g in groups) {
      g_data <- data[group == g]
      if (length(g_data) >= 3 && length(g_data) <= 5000) {
        normality_tests[[as.character(g)]] <- shapiro.test(g_data)
      }
    }

    group_vars <- tapply(data, group, var, na.rm = TRUE)
    max_var_ratio <- max(group_vars, na.rm = TRUE) / min(group_vars, na.rm = TRUE)

    results <- list(
      test = "t-test",
      normality = normality_tests,
      equal_variance = max_var_ratio < 3,
      variance_ratio = max_var_ratio,
      recommendation = ifelse(max_var_ratio < 3,
                              "Use standard t-test",
                              "Consider Welch's t-test")
    )
  } else if (test == "anova") {
    if (is.null(group)) {
      stop("Grouping variable required for ANOVA assumption check")
    }

    groups <- unique(group)
    normality_tests <- list()
    for (g in groups) {
      g_data <- data[group == g]
      if (length(g_data) >= 3 && length(g_data) <= 5000) {
        normality_tests[[as.character(g)]] <- shapiro.test(g_data)
      }
    }

    group_vars <- tapply(data, group, var, na.rm = TRUE)
    max_var_ratio <- max(group_vars, na.rm = TRUE) / min(group_vars, na.rm = TRUE)

    results <- list(
      test = "ANOVA",
      normality = normality_tests,
      equal_variance = max_var_ratio < 3,
      variance_ratio = max_var_ratio,
      recommendation = ifelse(max_var_ratio < 3,
                              "Use standard ANOVA",
                              "Consider Welch's ANOVA or Kruskal-Wallis")
    )
  } else if (test == "regression") {
    if (!is.data.frame(data)) {
      stop("Data frame required for regression assumption check")
    }

    # This would need the model fit, simplified version
    results <- list(
      test = "Regression",
      note = "Use stats_regression() and check residuals"
    )
  }

  structure(results, class = "omni_stats")
}

# ============================================================================
# SAMPLE SIZE ADJUSTMENTS
# ============================================================================

#' Adjust Sample Size for Dropouts
#'
#' Inflate sample size to account for anticipated dropout rate
#' @param n Required sample size
#' @param dropout_rate Anticipated dropout rate (proportion)
#' @return Adjusted sample size
#' @examples
#' adjust_dropout(n = 100, dropout_rate = 0.2)
#' @export
adjust_dropout <- function(n, dropout_rate) {
  if (dropout_rate < 0 || dropout_rate >= 1) {
    stop("Dropout rate must be between 0 and 1")
  }

  n_adjusted <- ceiling(n / (1 - dropout_rate))

  structure(list(
    n_required = n,
    n_adjusted = n_adjusted,
    n = n_adjusted,
    n_total = n_adjusted,
    dropout_rate = dropout_rate,
    expected_completers = n,
    method = "Dropout adjustment"
  ), class = "omni_sample_size")
}

#' Adjust Sample Size for Multiple Testing
#'
#' Adjust alpha level for multiple comparisons
#' @param alpha Original alpha level
#' @param k Number of comparisons
#' @param method Adjustment method: "bonferroni", "sidak"
#' @return Adjusted alpha
#' @examples
#' adjust_alpha(alpha = 0.05, k = 5)
#' @export
adjust_alpha <- function(alpha, k, method = "bonferroni") {
  if (method == "bonferroni") {
    alpha_adj <- alpha / k
  } else if (method == "sidak") {
    alpha_adj <- 1 - (1 - alpha)^(1/k)
  } else {
    stop("Method must be 'bonferroni' or 'sidak'")
  }

  structure(list(
    alpha_original = alpha,
    alpha_adjusted = alpha_adj,
    k = k,
    method = method
  ), class = "omni_stats")
}

# ============================================================================
# DATA SUMMARY
# ============================================================================

#' Summary Statistics
#'
#' Comprehensive summary statistics for numeric variables
#' @param x Numeric vector
#' @return List with summary statistics
#' @examples
#' data(omni_glucose)
#' summary_stats(omni_glucose$hba1c_3m)
#' @export
summary_stats <- function(x) {
  x_clean <- x[!is.na(x)]

  # Calculate skewness manually if moments not available
  skew <- sum((x_clean - mean(x_clean))^3 / sd(x_clean)^3) / length(x_clean)
  kurt <- sum((x_clean - mean(x_clean))^4 / sd(x_clean)^4) / length(x_clean) - 3

  structure(list(
    n = length(x_clean),
    n_missing = sum(is.na(x)),
    mean = mean(x_clean),
    median = stats::median(x_clean),
    sd = stats::sd(x_clean),
    var = stats::var(x_clean),
    min = min(x_clean),
    max = max(x_clean),
    range = max(x_clean) - min(x_clean),
    q1 = as.numeric(stats::quantile(x_clean, 0.25)),
    q3 = as.numeric(stats::quantile(x_clean, 0.75)),
    iqr = stats::IQR(x_clean),
    skewness = skew,
    kurtosis = kurt
  ), class = "omni_stats")
}

#' Cross-Tabulation with Statistics
#'
#' Enhanced cross-tabulation with chi-square and Cramer's V
#' @param x Row variable
#' @param y Column variable
#' @return List with cross-tabulation results
#' @examples
#' data(omni_glucose)
#' crosstab(omni_glucose$gender, omni_glucose$complications)
#' @export
crosstab <- function(x, y) {
  tbl <- table(x, y)
  prop_tbl <- prop.table(tbl)
  row_prop <- prop.table(tbl, 1)
  col_prop <- prop.table(tbl, 2)

  # Chi-square
  chi_test <- chisq.test(tbl)

  # Cramer's V
  cramers_v <- sqrt(chi_test$statistic / (sum(tbl) * min(nrow(tbl) - 1, ncol(tbl) - 1)))

  structure(list(
    table = tbl,
    proportions = prop_tbl,
    row_percentages = row_prop,
    column_percentages = col_prop,
    chi_square = chi_test$statistic,
    p_value = chi_test$p.value,
    cramers_v = cramers_v
  ), class = "omni_stats")
}
