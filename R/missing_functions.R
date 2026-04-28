#' Missing Functions for OmniSciKit
#'
#' Supplementary functions completing the OmniSciKit package interface.
#' Includes basic statistical tests, sample size wrappers, Table One,
#' and additional visualizations.

# ============================================================================
# BASIC STATISTICAL TESTS
# ============================================================================

#' Two-Sample T-Test
#'
#' Perform independent two-sample t-test
#' @param x Numeric vector for group 1
#' @param y Numeric vector for group 2
#' @param alternative "two.sided", "greater", or "less"
#' @param conf.level Confidence level (default 0.95)
#' @return List with t-test results
#' @examples
#' data(omni_glucose)
#' g1 <- omni_glucose$hba1c_3m[omni_glucose$treatment == "胰岛素"]
#' g2 <- omni_glucose$hba1c_3m[omni_glucose$treatment == "口服药"]
#' stats_ttest(g1, g2)
#' @export
stats_ttest <- function(x, y, alternative = "two.sided", conf.level = 0.95) {
  x <- na.omit(x)
  y <- na.omit(y)

  res <- t.test(x, y, alternative = alternative, conf.level = conf.level)

  structure(list(
    method = "Two-sample t-test",
    statistic = res$statistic,
    parameter = res$parameter,
    p_value = res$p.value,
    estimate = c(mean_x = mean(x), mean_y = mean(y),
                 diff = mean(x) - mean(y)),
    ci_lower = res$conf.int[1],
    ci_upper = res$conf.int[2],
    conf_level = conf.level,
    alternative = alternative
  ), class = "omni_stats")
}

#' Wilcoxon Rank-Sum Test
#'
#' Perform Wilcoxon/Mann-Whitney rank-sum test for two independent samples
#' @param x Numeric vector for group 1
#' @param y Numeric vector for group 2
#' @param alternative "two.sided", "greater", or "less"
#' @return List with test results
#' @examples
#' data(omni_glucose)
#' g1 <- omni_glucose$hba1c_3m[omni_glucose$treatment == "胰岛素"]
#' g2 <- omni_glucose$hba1c_3m[omni_glucose$treatment == "口服药"]
#' stats_wilcoxon(g1, g2)
#' @export
stats_wilcoxon <- function(x, y, alternative = "two.sided") {
  x <- na.omit(x)
  y <- na.omit(y)

  res <- wilcox.test(x, y, alternative = alternative, exact = FALSE)

  structure(list(
    method = "Wilcoxon rank-sum test",
    statistic = res$statistic,
    p_value = res$p.value,
    alternative = alternative,
    n1 = length(x),
    n2 = length(y)
  ), class = "omni_stats")
}

#' One-Way ANOVA
#'
#' Perform one-way analysis of variance
#' @param response Numeric response variable
#' @param group Factor grouping variable
#' @return List with ANOVA results
#' @examples
#' data(omni_glucose)
#' stats_anova(omni_glucose$hba1c_3m, omni_glucose$treatment)
#' @export
stats_anova <- function(response, group) {
  valid <- complete.cases(response, group)
  response <- response[valid]
  group <- as.factor(group[valid])

  fit <- aov(response ~ group)
  s <- summary(fit)

  # Extract ANOVA table
  anova_table <- as.data.frame(s[[1]])

  # Effect size (eta-squared)
  ss_between <- anova_table["group", "Sum Sq"]
  ss_total <- ss_between + anova_table["Residuals", "Sum Sq"]
  eta_sq <- ss_between / ss_total

  # Pairwise comparisons
  tukey <- TukeyHSD(fit)

  structure(list(
    method = "One-way ANOVA",
    statistic = anova_table["group", "F value"],
    df1 = anova_table["group", "Df"],
    df2 = anova_table["Residuals", "Df"],
    p_value = anova_table["group", "Pr(>F)"],
    eta_squared = eta_sq,
    anova_table = anova_table,
    posthoc = tukey,
    n = length(response),
    groups = levels(group)
  ), class = "omni_stats")
}

#' Correlation Test
#'
#' Calculate Pearson or Spearman correlation
#' @param x Numeric vector
#' @param y Numeric vector
#' @param method "pearson" or "spearman"
#' @return List with correlation results
#' @examples
#' data(omni_glucose)
#' stats_correlation(omni_glucose$age, omni_glucose$bmi)
#' @export
stats_correlation <- function(x, y, method = "pearson") {
  valid <- complete.cases(x, y)
  x <- x[valid]
  y <- y[valid]

  res <- cor.test(x, y, method = method)

  structure(list(
    method = paste0(method, " correlation"),
    estimate = res$estimate,
    statistic = res$statistic,
    parameter = res$parameter,
    p_value = res$p.value,
    ci_lower = if (method == "pearson") res$conf.int[1] else NA,
    ci_upper = if (method == "pearson") res$conf.int[2] else NA,
    n = length(x)
  ), class = "omni_stats")
}

#' Linear Regression
#'
#' Perform linear regression analysis
#' @param formula Model formula
#' @param data Data frame
#' @return List with regression results
#' @examples
#' data(omni_glucose)
#' stats_regression(hba1c_3m ~ age + bmi, data = omni_glucose)
#' @export
stats_regression <- function(formula, data) {
  fit <- lm(formula, data = data)
  s <- summary(fit)

  coef_table <- as.data.frame(s$coefficients)

  structure(list(
    method = "Linear Regression",
    coefficients = coef_table,
    r_squared = s$r.squared,
    adj_r_squared = s$adj.r.squared,
    f_statistic = s$fstatistic[1],
    f_df1 = s$fstatistic[2],
    f_df2 = s$fstatistic[3],
    f_p_value = pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3],
                   lower.tail = FALSE),
    sigma = s$sigma,
    n = length(fit$residuals),
    formula = deparse(formula),
    fit = fit
  ), class = "omni_stats")
}

#' Chi-Square Test
#'
#' Perform chi-square test of independence
#' @param x Factor or table
#' @param y Factor (optional, for 2-variable case)
#' @return List with test results
#' @examples
#' data(omni_glucose)
#' stats_chisq(omni_glucose$treatment, omni_glucose$complications)
#' @export
stats_chisq <- function(x, y = NULL) {
  if (!is.null(y)) {
    tbl <- table(x, y)
  } else if (is.table(x) || is.matrix(x)) {
    tbl <- x
  } else {
    stop("Provide two factors or a contingency table")
  }

  res <- chisq.test(tbl)

  # Cramer's V
  n <- sum(tbl)
  cramers_v <- sqrt(res$statistic / (n * min(nrow(tbl) - 1, ncol(tbl) - 1)))

  structure(list(
    method = "Chi-square test",
    statistic = res$statistic,
    parameter = res$parameter,
    p_value = res$p.value,
    observed = tbl,
    expected = res$expected,
    cramers_v = as.numeric(cramers_v),
    n = n
  ), class = "omni_stats")
}

#' Normality Test
#'
#' Perform Shapiro-Wilk normality test with Q-Q plot option
#' @param x Numeric vector
#' @param plot Logical, show Q-Q plot (default FALSE)
#' @return List with normality test results
#' @examples
#' data(omni_glucose)
#' stats_normality(omni_glucose$hba1c_3m)
#' @export
stats_normality <- function(x, plot = FALSE) {
  x_clean <- na.omit(x)

  if (length(x_clean) < 3) {
    stop("At least 3 non-missing values required")
  }

  if (length(x_clean) > 5000) {
    x_clean <- sample(x_clean, 5000)
  }

  shapiro <- shapiro.test(x_clean)

  # Skewness and kurtosis
  n <- length(x_clean)
  skew <- sum((x_clean - mean(x_clean))^3 / sd(x_clean)^3) / n
  kurt <- sum((x_clean - mean(x_clean))^4 / sd(x_clean)^4) / n - 3

  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 required for Q-Q plot")
    } else {
      qq_data <- data.frame(
        theoretical = qnorm(ppoints(n))[rank(x_clean)],
        sample = sort(x_clean)
      )
      p <- ggplot2::ggplot(qq_data, ggplot2::aes(x = theoretical, y = sample)) +
        ggplot2::geom_point(color = "steelblue") +
        ggplot2::geom_abline(intercept = mean(x_clean),
                            slope = sd(x_clean), color = "red") +
        ggplot2::labs(title = "Q-Q Plot",
                      x = "Theoretical Quantiles",
                      y = "Sample Quantiles") +
        ggplot2::theme_minimal()
      print(p)
    }
  }

  structure(list(
    method = "Shapiro-Wilk normality test",
    statistic = shapiro$statistic,
    p_value = shapiro$p.value,
    is_normal = shapiro$p.value > 0.05,
    skewness = skew,
    kurtosis = kurt,
    n = length(x_clean)
  ), class = "omni_stats")
}

# ============================================================================
# TABLE ONE - BASELINE CHARACTERISTICS
# ============================================================================

#' Table One: Baseline Characteristics
#'
#' Create a baseline characteristics table (Table 1) for clinical studies
#' @param data Data frame
#' @param vars Character vector of variable names (NULL = all)
#' @param group Grouping variable name (e.g., "treatment")
#' @param continuous_fun Function for continuous summaries (default mean/sd)
#' @param categorical_fun Function for categorical summaries (default n/pct)
#' @param digits Number of decimal places (default 2)
#' @param pvalue Logical, include p-values (default TRUE)
#' @return Data frame formatted as Table 1
#' @examples
#' data(omni_glucose)
#' table_one(omni_glucose, vars = c("age", "bmi", "hba1c_3m", "gender", "complications"),
#'           group = "treatment")
#' @export
table_one <- function(data, vars = NULL, group = NULL,
                      continuous_fun = NULL, categorical_fun = NULL,
                      digits = 2, pvalue = TRUE) {

  if (is.null(vars)) {
    vars <- setdiff(names(data), group)
  }

  if (is.null(continuous_fun)) {
    continuous_fun <- function(x) {
      paste0(round(mean(x, na.rm = TRUE), digits), " (",
             round(sd(x, na.rm = TRUE), digits), ")")
    }
  }

  if (is.null(categorical_fun)) {
    categorical_fun <- function(x) {
      tbl <- table(x, useNA = "ifany")
      freq <- paste0(names(tbl), ": ", tbl, " (",
                     round(100 * tbl / sum(tbl), 1), "%)")
      paste(freq, collapse = "; ")
    }
  }

  result_list <- list()

  for (var in vars) {
    if (!var %in% names(data)) next

    is_num <- is.numeric(data[[var]])
    is_cat <- is.factor(data[[var]]) || is.character(data[[var]])

    if (is.null(group) || !group %in% names(data)) {
      # Overall only
      if (is_num) {
        result_list[[var]] <- data.frame(
          Variable = var,
          Overall = continuous_fun(data[[var]]),
          stringsAsFactors = FALSE
        )
      } else if (is_cat) {
        result_list[[var]] <- data.frame(
          Variable = var,
          Overall = categorical_fun(data[[var]]),
          stringsAsFactors = FALSE
        )
      }
    } else {
      # By group
      group_levels <- unique(data[[group]])
      group_levels <- sort(as.character(group_levels))

      row_data <- data.frame(Variable = var, stringsAsFactors = FALSE)

      for (gl in group_levels) {
        gl_data <- data[[var]][data[[group]] == gl]
        if (is_num) {
          row_data[[gl]] <- continuous_fun(gl_data)
        } else if (is_cat) {
          row_data[[gl]] <- categorical_fun(gl_data)
        }
      }

      # P-value
      if (pvalue) {
        if (is_num && length(unique(data[[group]])) == 2) {
          pv <- tryCatch(
            t.test(data[[var]] ~ data[[group]])$p.value,
            error = function(e) NA
          )
        } else if (is_cat) {
          pv <- tryCatch(
            suppressWarnings(chisq.test(data[[var]], data[[group]]))$p.value,
            error = function(e) NA
          )
        } else {
          pv <- NA
        }

        row_data[["P-value"]] <- ifelse(is.na(pv), "",
                                        ifelse(pv < 0.001, "<0.001",
                                               round(pv, 4)))
      }

      result_list[[var]] <- row_data
    }
  }

  result <- do.call(rbind, result_list)
  row.names(result) <- NULL

  structure(result, class = c("omni_table_one", "data.frame"))
}

# ============================================================================
# BASIC VISUALIZATIONS
# ============================================================================

#' Box Plot
#'
#' Create a box plot for comparing groups
#' @param x Grouping variable
#' @param y Numeric response variable
#' @param title Plot title
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @return ggplot object
#' @examples
#' data(omni_glucose)
#' plot_boxplot(omni_glucose$treatment, omni_glucose$hba1c_3m,
#'              title = "HbA1c by Treatment Group")
#' @export
plot_boxplot <- function(x, y, title = NULL, xlab = NULL, ylab = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  valid <- complete.cases(x, y)
  plot_data <- data.frame(
    group = as.factor(x[valid]),
    value = y[valid]
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = value, fill = group)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.color = "red", outlier.alpha = 0.5) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.3, color = "gray30") +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 18,
                          size = 3, color = "darkred") +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(title = title %||% "Box Plot",
                  x = xlab %||% deparse(substitute(x)),
                  y = ylab %||% deparse(substitute(y))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  p
}

#' Scatter Plot with Regression
#'
#' Create a scatter plot with optional regression line
#' @param x Numeric predictor
#' @param y Numeric response
#' @param group Optional grouping variable (for colors)
#' @param reg_line Logical, add regression line (default TRUE)
#' @param conf_int Logical, show confidence band (default TRUE)
#' @param title Plot title
#' @return ggplot object
#' @examples
#' data(omni_glucose)
#' plot_scatter(omni_glucose$age, omni_glucose$hba1c_3m,
#'              title = "Age vs HbA1c")
#' @export
plot_scatter <- function(x, y, group = NULL, reg_line = TRUE,
                          conf_int = TRUE, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  valid <- complete.cases(x, y)
  plot_data <- data.frame(
    x = x[valid],
    y = y[valid]
  )

  if (!is.null(group)) {
    plot_data$group <- as.factor(group[valid])
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y))

  if (!is.null(group)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = group), alpha = 0.6)
    if (reg_line) {
      p <- p + ggplot2::geom_smooth(ggplot2::aes(color = group),
                                     method = "lm", se = conf_int,
                                     formula = y ~ x, alpha = 0.2)
    }
  } else {
    p <- p + ggplot2::geom_point(color = "steelblue", alpha = 0.6)
    if (reg_line) {
      p <- p + ggplot2::geom_smooth(method = "lm", se = conf_int,
                                     color = "red", formula = y ~ x,
                                     alpha = 0.2)
    }
  }

  p + ggplot2::labs(title = title %||% "Scatter Plot",
                    x = deparse(substitute(x)),
                    y = deparse(substitute(y))) +
    ggplot2::theme_minimal()
}

#' Correlation Matrix Heatmap
#'
#' Plot correlation matrix as a heatmap
#' @param data Data frame with numeric variables
#' @param method "pearson" or "spearman"
#' @param title Plot title
#' @param show_values Logical, display correlation values (default TRUE)
#' @param palette Color palette: "red", "blue", or "green"
#' @return ggplot object
#' @examples
#' data(omni_glucose)
#' plot_correlation_matrix(omni_glucose[, c("age", "bmi", "hba1c_3m", "hba1c_baseline")])
#' @export
plot_correlation_matrix <- function(data, method = "pearson", title = NULL,
                                     show_values = TRUE, palette = "blue") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  # Select only numeric columns
  num_data <- data[, sapply(data, is.numeric)]
  if (ncol(num_data) < 2) {
    stop("At least 2 numeric variables required")
  }

  # Correlation matrix
  cormat <- cor(num_data, use = "pairwise.complete.obs", method = method)

  # Reshape for ggplot
  var_names <- colnames(cormat)
  plot_data <- expand.grid(Var1 = var_names, Var2 = var_names)
  plot_data$value <- as.vector(cormat)

  # Color scale
  if (palette == "red") {
    colors <- c("#B2182B", "#EF8A62", "#FDDBC7", "#F7F7F7",
                "#D1E5F0", "#67A9CF", "#2166AC")
  } else if (palette == "green") {
    colors <- c("#A6611A", "#DFC27D", "#F5F5F5", "#80CDC1", "#018571")
  } else {
    colors <- c("#0571B0", "#92C5DE", "#F7F7F7", "#F4A582", "#CA0020")
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradientn(colors = colors, limits = c(-1, 1),
                                   name = "Correlation")

  if (show_values) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = round(value, 2)),
                                 size = 3.5, color = "gray20")
  }

  p + ggplot2::labs(title = title %||% "Correlation Matrix",
                    x = "", y = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                        hjust = 1, vjust = 1),
                   panel.grid = ggplot2::element_blank())
}

# ============================================================================
# HELPER: NULL coalescing operator
# ============================================================================

#' NULL coalescing operator
#'
#' Returns left-hand side if not NULL, otherwise right-hand side
#' @param x Left value
#' @param y Right value (fallback)
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}

# ============================================================================
# SAMPLE SIZE CONVENIENCE WRAPPERS
# ============================================================================

#' Sample Size for Two-Sample T-Test (Equal Variances)
#'
#' Convenience wrapper for two-sample t-test sample size
#' @param delta Expected difference between groups
#' @param sd Standard deviation (pooled)
#' @param ratio Allocation ratio n2/n1 (default 1)
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_ttest(delta = 0.5, sd = 1)
#' @export
ss_ttest <- function(delta, sd, ratio = 1, power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  n1 <- ceiling((z_alpha + z_beta)^2 * sd^2 * (1 + 1/ratio) / delta^2)
  n2 <- ceiling(n1 * ratio)
  n_total <- n1 + n2

  structure(list(
    n_per_group = n1,
    n1 = n1,
    n2 = n2,
    n_total = n_total,
    delta = delta,
    sd = sd,
    ratio = ratio,
    power = power,
    alpha = alpha,
    method = "Two-sample t-test (equal variances)",
    formula = "n₁ = (Z_{α/2} + Z_β)² × σ² × (1+1/ρ) / δ²",
    reference = "Chow et al. (2017) Sample Size Calculations in Clinical Research, 3rd ed."
  ), class = "omni_sample_size")
}

#' Sample Size for One-Way ANOVA
#'
#' Convenience wrapper for ANOVA sample size
#' @param k Number of groups
#' @param f Cohen's f effect size
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_anova(k = 3, f = 0.25)
#' @export
ss_anova <- function(k, f, power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  n_per_group <- ceiling((z_alpha + z_beta)^2 / f^2)
  n_total <- n_per_group * k

  structure(list(
    n_per_group = n_per_group,
    n_total = n_total,
    k = k,
    f = f,
    power = power,
    alpha = alpha,
    method = "One-way ANOVA",
    formula = "n = (Z_{α/2} + Z_β)² / f²",
    reference = "Cohen (1988) Statistical Power Analysis, 2nd ed."
  ), class = "omni_sample_size")
}

#' Sample Size for Correlation Test
#'
#' Convenience wrapper for Pearson correlation sample size
#' @param r Expected correlation coefficient
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @param alternative "two.sided" or "one.sided"
#' @return List with sample size results
#' @examples
#' ss_correlation(r = 0.3)
#' @export
ss_correlation <- function(r, power = 0.8, alpha = 0.05,
                            alternative = "two.sided") {
  if (alternative == "two.sided") {
    z_alpha <- qnorm(1 - alpha/2)
  } else {
    z_alpha <- qnorm(1 - alpha)
  }
  z_beta <- qnorm(power)

  z_r <- 0.5 * log((1 + r) / (1 - r))
  n <- ceiling((z_alpha + z_beta)^2 / z_r^2 + 3)

  structure(list(
    n = n,
    r = r,
    power = power,
    alpha = alpha,
    alternative = alternative,
    method = "Correlation test",
    formula = "n = (Z_{α/2} + Z_β)² / z(r)² + 3, z(r) = 0.5×log((1+r)/(1-r))",
    reference = "Cohen (1988) Statistical Power Analysis, 2nd ed."
  ), class = "omni_sample_size")
}

#' Sample Size for Two Proportions
#'
#' Convenience wrapper for comparing two proportions sample size
#' @param p1 Proportion in group 1
#' @param p2 Proportion in group 2
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @param alternative "two.sided" or "one.sided"
#' @return List with sample size results
#' @examples
#' ss_proportion(p1 = 0.6, p2 = 0.4)
#' @export
ss_proportion <- function(p1, p2, power = 0.8, alpha = 0.05,
                           alternative = "two.sided") {
  if (alternative == "two.sided") {
    z_alpha <- qnorm(1 - alpha/2)
  } else {
    z_alpha <- qnorm(1 - alpha)
  }
  z_beta <- qnorm(power)

  p_bar <- (p1 + p2) / 2
  n_per_group <- ceiling(
    (z_alpha * sqrt(2 * p_bar * (1 - p_bar)) +
       z_beta * sqrt(p1 * (1 - p1) + p2 * (1 - p2)))^2 / (p1 - p2)^2
  )
  n_total <- n_per_group * 2

  structure(list(
    n_per_group = n_per_group,
    n_total = n_total,
    p1 = p1,
    p2 = p2,
    power = power,
    alpha = alpha,
    method = "Two proportions test",
    formula = "n = (Z_{α/2}√(2p̄(1-p̄)) + Z_β√(p₁(1-p₁)+p₂(1-p₂)))² / (p₁-p₂)²",
    reference = "Chow et al. (2017) Sample Size Calculations in Clinical Research, 3rd ed."
  ), class = "omni_sample_size")
}

#' Sample Size for Chi-Square Test
#'
#' Convenience wrapper for chi-square test sample size
#' @param w Cohen's w effect size
#' @param df Degrees of freedom
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_chisq(w = 0.3, df = 1)
#' @export
ss_chisq <- function(w, df = 1, power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  n <- ceiling((z_alpha + z_beta)^2 / w^2)

  structure(list(
    n = n,
    w = w,
    df = df,
    power = power,
    alpha = alpha,
    method = "Chi-square test",
    formula = "n = (Z_{α/2} + Z_β)² / w²",
    reference = "Cohen (1988) Statistical Power Analysis, 2nd ed."
  ), class = "omni_sample_size")
}

#' Sample Size for Diagnostic Test (Sensitivity/Specificity)
#'
#' Convenience wrapper for diagnostic test sample size
#' @param sens Expected sensitivity
#' @param spec Expected specificity
#' @param prev Disease prevalence
#' @param w Margin of error (default 0.05)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_diagnostic(sens = 0.85, spec = 0.9, prev = 0.3)
#' @export
ss_diagnostic <- function(sens, spec, prev, w = 0.05, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)

  n_sens <- ceiling(z_alpha^2 * sens * (1 - sens) / w^2)
  n_spec <- ceiling(z_alpha^2 * spec * (1 - spec) / w^2)

  n_disease <- ceiling(n_sens)
  n_healthy <- ceiling(n_spec)

  # Adjust for prevalence
  n_total <- ceiling(max(n_disease / prev, n_healthy / (1 - prev)))

  structure(list(
    n_disease = n_disease,
    n_healthy = n_healthy,
    n_total = n_total,
    n_sens = n_sens,
    n_spec = n_spec,
    prevalence = prev,
    sens = sens,
    spec = spec,
    w = w,
    alpha = alpha,
    method = "Diagnostic test",
    formula = "n_disease = Z² × sens × (1-sens) / w², n_healthy = Z² × spec × (1-spec) / w²",
    reference = "Buderer (1996) Academic Emergency Medicine, 3(9), 895-900"
  ), class = "omni_sample_size")
}

# ============================================================================
# ADDITIONAL PLOTS
# ============================================================================

#' ROC Curve (Alias)
#'
#' Alias for plot_roc_enhanced. See \code{\link{plot_roc_enhanced}} for details.
#' @param predictor Predictor variable
#' @param response Response variable (binary)
#' @param title Plot title
#' @param ci Show confidence interval (default TRUE)
#' @return List with AUC and plot
#' @examples
#' data(omni_diagnostic)
#' disease_binary <- ifelse(omni_diagnostic$true_disease == "阳性", 1, 0)
#' plot_roc(omni_diagnostic$test1_result, disease_binary)
#' @export
plot_roc <- function(predictor, response, title = "ROC Curve", ci = TRUE) {
  plot_roc_enhanced(predictor, response, title, ci)
}

#' Sample Size Curve
#'
#' Plot required sample size across effect sizes
#' @param effect_sizes Vector of effect sizes
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @param test_type Type of test: "ttest", "proportion", "correlation"
#' @param ... Additional arguments passed to internal calculators
#' @return ggplot object
#' @examples
#' plot_sample_size_curve(effect_sizes = seq(0.2, 1.0, by = 0.1), power = 0.8)
#' @export
plot_sample_size_curve <- function(effect_sizes = seq(0.2, 1.0, by = 0.1),
                                    power = 0.8, alpha = 0.05,
                                    test_type = "ttest", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  dots <- list(...)
  sample_sizes <- numeric(length(effect_sizes))

  for (i in seq_along(effect_sizes)) {
    es <- effect_sizes[i]
    tryCatch({
      if (test_type == "ttest") {
        res <- ss_ttest(delta = es, sd = 1, power = power, alpha = alpha)
        sample_sizes[i] <- res$n_total
      } else if (test_type == "proportion") {
        p2 <- if (!is.null(dots$p2)) dots$p2 else 0.5
        res <- ss_proportion(p1 = es, p2 = p2, power = power, alpha = alpha)
        sample_sizes[i] <- res$n_total
      } else if (test_type == "correlation") {
        res <- ss_correlation(r = es, power = power, alpha = alpha)
        sample_sizes[i] <- res$n
      }
    }, error = function(e) {
      sample_sizes[i] <- NA
    })
  }

  df <- data.frame(effect_size = effect_sizes, n = sample_sizes)

  ggplot2::ggplot(df, ggplot2::aes(x = effect_size, y = n)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(color = "steelblue", size = 2) +
    ggplot2::labs(
      title = "Sample Size Curve",
      subtitle = paste("Power =", power, ", Alpha =", alpha),
      x = "Effect Size",
      y = "Required Sample Size"
    ) +
    ggplot2::theme_minimal()
}
