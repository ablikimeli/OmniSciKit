#' Complete Power Analysis Functions
#'
#' Power calculations for various study designs

# ============================================================================
# BASIC POWER FUNCTIONS
# ============================================================================

#' Power for Two-Sample T-Test
#'
#' Calculate statistical power for independent two-sample t-test
#' @param n Sample size per group
#' @param delta Expected difference between groups
#' @param sd Standard deviation (pooled)
#' @param alpha Significance level (default 0.05)
#' @param alternative "two.sided" or "one.sided"
#' @return List with power analysis results
#' @examples
#' power_ttest(n = 64, delta = 0.5, sd = 1)
#' @export
power_ttest <- function(n, delta, sd, alpha = 0.05, alternative = "two.sided") {
  if (alternative == "two.sided") {
    z_alpha <- qnorm(1 - alpha/2)
  } else {
    z_alpha <- qnorm(1 - alpha)
  }

  se <- sd * sqrt(2/n)
  z_beta <- (delta / se) - z_alpha
  power <- pnorm(z_beta)

  structure(list(
    power = power,
    n = n,
    delta = delta,
    sd = sd,
    alpha = alpha,
    alternative = alternative,
    method = "Two-sample t-test"
  ), class = "omni_power")
}

#' Power for One-Way ANOVA
#'
#' Calculate statistical power for one-way ANOVA
#' @param n Sample size per group
#' @param k Number of groups
#' @param f Cohen's f effect size
#' @param alpha Significance level (default 0.05)
#' @return List with power analysis results
#' @examples
#' power_anova(n = 50, k = 3, f = 0.25)
#' @export
power_anova <- function(n, k, f, alpha = 0.05) {
  df1 <- k - 1
  df2 <- k * (n - 1)
  ncp <- f^2 * (df1 + df2 + 1)

  crit <- qf(1 - alpha, df1, df2)
  power <- 1 - pf(crit, df1, df2, ncp = ncp)

  structure(list(
    power = power,
    n = n,
    k = k,
    f = f,
    alpha = alpha,
    method = "One-way ANOVA"
  ), class = "omni_power")
}

#' Power for Correlation
#'
#' Calculate statistical power for Pearson correlation
#' @param n Sample size
#' @param r Expected correlation coefficient
#' @param alpha Significance level (default 0.05)
#' @param alternative "two.sided" or "one.sided"
#' @return List with power analysis results
#' @examples
#' power_correlation(n = 100, r = 0.3)
#' @export
power_correlation <- function(n, r, alpha = 0.05, alternative = "two.sided") {
  if (alternative == "two.sided") {
    z_alpha <- qnorm(1 - alpha/2)
  } else {
    z_alpha <- qnorm(1 - alpha)
  }

  z_r <- 0.5 * log((1 + r) / (1 - r))
  se <- 1 / sqrt(n - 3)
  z_beta <- z_r / se - z_alpha
  power <- pnorm(z_beta)

  structure(list(
    power = power,
    n = n,
    r = r,
    alpha = alpha,
    alternative = alternative,
    method = "Pearson correlation"
  ), class = "omni_power")
}

#' Power for Two Proportions
#'
#' Calculate statistical power for comparing two proportions
#' @param n Sample size per group
#' @param p1 Proportion in group 1
#' @param p2 Proportion in group 2
#' @param alpha Significance level (default 0.05)
#' @param alternative "two.sided" or "one.sided"
#' @return List with power analysis results
#' @examples
#' power_proportion(n = 100, p1 = 0.3, p2 = 0.5)
#' @export
power_proportion <- function(n, p1, p2, alpha = 0.05, alternative = "two.sided") {
  if (alternative == "two.sided") {
    z_alpha <- qnorm(1 - alpha/2)
  } else {
    z_alpha <- qnorm(1 - alpha)
  }

  p_bar <- (p1 + p2) / 2
  se_null <- sqrt(2 * p_bar * (1 - p_bar) / n)
  se_alt <- sqrt(p1 * (1 - p1) / n + p2 * (1 - p2) / n)

  z_beta <- (abs(p1 - p2) - z_alpha * se_null) / se_alt
  power <- pnorm(z_beta)

  structure(list(
    power = power,
    n = n,
    p1 = p1,
    p2 = p2,
    alpha = alpha,
    alternative = alternative,
    method = "Two proportions"
  ), class = "omni_power")
}

# ============================================================================
# ADVANCED POWER FUNCTIONS
# ============================================================================

#' Power for Survival Analysis (Log-rank Test)
#'
#' Calculate power for survival analysis
#' @param n Total sample size
#' @param hr Hazard ratio
#' @param p_event Probability of event
#' @param p_t Proportion in treatment group (default 0.5)
#' @param alpha Significance level (default 0.05)
#' @return List with power results
#' @examples
#' power_survival(n = 200, hr = 0.7, p_event = 0.6)
#' @export
power_survival <- function(n, hr, p_event, p_t = 0.5, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)

  d <- n * p_event

  z_beta <- sqrt(d * p_t * (1 - p_t)) * abs(log(hr)) - z_alpha
  power <- pnorm(z_beta)

  structure(list(
    method = "Survival analysis power",
    n = n,
    hr = hr,
    p_event = p_event,
    expected_events = d,
    power = power,
    alpha = alpha
  ), class = "omni_power")
}

#' Power for Cox Regression
#'
#' Calculate power for Cox proportional hazards regression
#' @param n Total sample size
#' @param hr Hazard ratio
#' @param p_event Event probability
#' @param r2 R-squared with other covariates (default 0)
#' @param alpha Significance level (default 0.05)
#' @return List with power results
#' @examples
#' power_cox_regression(n = 300, hr = 1.5, p_event = 0.4)
#' @export
power_cox_regression <- function(n, hr, p_event, r2 = 0, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)

  d <- n * p_event
  log_hr <- log(hr)

  z_beta <- sqrt(d * (1 - r2)) * abs(log_hr) - z_alpha
  power <- pnorm(z_beta)

  structure(list(
    method = "Cox regression power",
    n = n,
    hr = hr,
    p_event = p_event,
    events = d,
    power = power,
    alpha = alpha
  ), class = "omni_power")
}

#' Power for Logistic Regression
#'
#' Calculate power for logistic regression
#' @param n Sample size
#' @param p Event probability
#' @param or Odds ratio
#' @param r2 R-squared with other covariates (default 0)
#' @param alpha Significance level (default 0.05)
#' @return List with power results
#' @examples
#' power_logistic_regression(n = 200, p = 0.2, or = 2.0)
#' @export
power_logistic_regression <- function(n, p, or, r2 = 0, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)

  log_or <- log(or)
  b <- p * (1 - p) * log_or^2

  z_beta <- sqrt(n * (1 - r2) * b) - z_alpha
  power <- pnorm(z_beta)

  structure(list(
    method = "Logistic regression power",
    n = n,
    p = p,
    or = or,
    power = power,
    alpha = alpha
  ), class = "omni_power")
}

#' Power for Equivalence Test (Means)
#'
#' Calculate power for equivalence trial
#' @param n Sample size per group
#' @param delta Equivalence margin
#' @param sd Standard deviation
#' @param true_diff True difference (default 0)
#' @param alpha Significance level (default 0.05)
#' @return List with power results
#' @examples
#' power_equivalence_mean(n = 100, delta = 5, sd = 10)
#' @export
power_equivalence_mean <- function(n, delta, sd, true_diff = 0, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha)

  se <- sd * sqrt(2/n)
  margin <- delta - abs(true_diff)

  power <- 2 * pnorm(margin / se - z_alpha) - 1

  structure(list(
    method = "Equivalence test power (means)",
    n = n,
    delta = delta,
    sd = sd,
    true_diff = true_diff,
    power = max(0, power),
    alpha = alpha
  ), class = "omni_power")
}

#' Power for Non-inferiority Test (Means)
#'
#' Calculate power for non-inferiority trial
#' @param n Sample size per group
#' @param margin Non-inferiority margin
#' @param sd Standard deviation
#' @param true_diff True difference (default 0)
#' @param alpha Significance level (default 0.05)
#' @return List with power results
#' @examples
#' power_noninferiority_mean(n = 100, margin = 3, sd = 8)
#' @export
power_noninferiority_mean <- function(n, margin, sd, true_diff = 0, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha)

  se <- sd * sqrt(2/n)
  eff_margin <- margin - abs(true_diff)

  z_beta <- eff_margin / se - z_alpha
  power <- pnorm(z_beta)

  structure(list(
    method = "Non-inferiority test power (means)",
    n = n,
    margin = margin,
    sd = sd,
    true_diff = true_diff,
    power = power,
    alpha = alpha
  ), class = "omni_power")
}

#' Power for Cluster Randomized Trial
#'
#' Calculate power for cluster randomized design
#' @param n_clusters Number of clusters per group
#' @param cluster_size Average cluster size
#' @param delta Effect size
#' @param sd Standard deviation
#' @param icc Intraclass correlation
#' @param alpha Significance level (default 0.05)
#' @return List with power results
#' @examples
#' power_cluster_randomized(n_clusters = 20, cluster_size = 30,
#'                          delta = 0.5, sd = 1, icc = 0.05)
#' @export
power_cluster_randomized <- function(n_clusters, cluster_size, delta, sd,
                                      icc, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)

  # Design effect
  de <- 1 + (cluster_size - 1) * icc

  # Effective sample size
  n_effective <- n_clusters * cluster_size / de

  se <- sd * sqrt(2 / n_effective)
  z_beta <- delta / se - z_alpha
  power <- pnorm(z_beta)

  structure(list(
    method = "Cluster randomized trial power",
    n_clusters = n_clusters,
    cluster_size = cluster_size,
    effective_n = n_effective,
    design_effect = de,
    delta = delta,
    power = power,
    alpha = alpha
  ), class = "omni_power")
}

# ============================================================================
# POWER CURVES
# ============================================================================

#' Plot Power Curve
#'
#' Visualize power across different sample sizes
#' @param n_range Range of sample sizes
#' @param effect_size Effect size
#' @param alpha Significance level (default 0.05)
#' @param test_type Type of test
#' @param ... Additional parameters
#' @return ggplot object
#' @examples
#' plot_power_curve(n_range = seq(10, 200, by = 10), effect_size = 0.5)
#' @export
plot_power_curve <- function(n_range = seq(10, 200, by = 10),
                              effect_size = 0.5,
                              alpha = 0.05,
                              test_type = "ttest", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  powers <- numeric(length(n_range))

  for (i in seq_along(n_range)) {
    n <- n_range[i]
    tryCatch({
      dots <- list(...)
      if (test_type == "ttest") {
        result <- power_ttest(n = n, delta = effect_size, sd = 1, alpha = alpha)
      } else if (test_type == "anova") {
        k <- if (!is.null(dots$k)) dots$k else 3
        result <- power_anova(n = n, k = k, f = effect_size, alpha = alpha)
      } else if (test_type == "correlation") {
        result <- power_correlation(n = n, r = effect_size, alpha = alpha)
      } else if (test_type == "proportion") {
        p1 <- if (!is.null(dots$p1)) dots$p1 else 0.3
        p2 <- if (!is.null(dots$p2)) dots$p2 else 0.5
        result <- power_proportion(n = n, p1 = p1, p2 = p2, alpha = alpha)
      } else if (test_type == "survival") {
        hr <- if (!is.null(dots$hr)) dots$hr else 0.7
        p_event <- if (!is.null(dots$p_event)) dots$p_event else 0.6
        result <- power_survival(n = n, hr = hr, p_event = p_event, alpha = alpha)
      }
      powers[i] <- result$power
    }, error = function(e) {
      powers[i] <- NA
    })
  }

  df <- data.frame(n = n_range, power = powers)

  ggplot2::ggplot(df, ggplot2::aes(x = n, y = power)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(color = "steelblue", size = 2) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed",
                        color = "red", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0.9, linetype = "dashed",
                        color = "orange", alpha = 0.7) +
    ggplot2::labs(
      title = "Power Curve",
      subtitle = paste("Effect size =", effect_size, ", Alpha =", alpha),
      x = "Sample Size per Group",
      y = "Statistical Power"
    ) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_minimal()
}

#' Plot Effect Size Curve
#'
#' Visualize power across different effect sizes
#' @param effect_range Range of effect sizes
#' @param n Sample size
#' @param alpha Significance level
#' @param test_type Type of test
#' @return ggplot object
#' @examples
#' plot_effect_size_curve(effect_range = seq(0.1, 1, by = 0.1), n = 100)
#' @export
plot_effect_size_curve <- function(effect_range = seq(0.1, 1, by = 0.1),
                                    n = 100, alpha = 0.05,
                                    test_type = "ttest") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  powers <- numeric(length(effect_range))

  for (i in seq_along(effect_range)) {
    es <- effect_range[i]
    tryCatch({
      if (test_type == "ttest") {
        result <- power_ttest(n = n, delta = es, sd = 1, alpha = alpha)
      } else if (test_type == "correlation") {
        result <- power_correlation(n = n, r = es, alpha = alpha)
      }
      powers[i] <- result$power
    }, error = function(e) {
      powers[i] <- NA
    })
  }

  df <- data.frame(effect_size = effect_range, power = powers)

  ggplot2::ggplot(df, ggplot2::aes(x = effect_size, y = power)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(color = "steelblue", size = 2) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed",
                        color = "red", alpha = 0.7) +
    ggplot2::labs(
      title = "Power vs Effect Size",
      subtitle = paste("N =", n, ", Alpha =", alpha),
      x = "Effect Size",
      y = "Statistical Power"
    ) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_minimal()
}
