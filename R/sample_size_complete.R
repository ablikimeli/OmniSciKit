#' Complete Sample Size Functions for Clinical Research
#'
#' Comprehensive sample size calculations covering all major clinical trial designs
#' including superiority, non-inferiority, equivalence, and various study types.

# ============================================================================
# SUPERIORITY TRIALS (差异性检验)
# ============================================================================

#' One-Sample T-Test Sample Size
#'
#' Calculate sample size for one-sample t-test
#' @param mu0 Null hypothesis mean
#' @param mu1 Alternative hypothesis mean
#' @param sd Standard deviation
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @param alternative "two.sided" or "one.sided"
#' @return List with sample size results
#' @examples
#' ss_ttest_one_sample(mu0 = 100, mu1 = 105, sd = 15)
#' @export
ss_ttest_one_sample <- function(mu0, mu1, sd, power = 0.8, alpha = 0.05,
                                 alternative = "two.sided") {
  delta <- abs(mu1 - mu0)
  if (alternative == "two.sided") {
    z_alpha <- qnorm(1 - alpha/2)
  } else {
    z_alpha <- qnorm(1 - alpha)
  }
  z_beta <- qnorm(power)

  n <- ceiling(((z_alpha + z_beta) * sd / delta)^2)

  structure(list(
    n = n,
    mu0 = mu0,
    mu1 = mu1,
    delta = delta,
    sd = sd,
    power = power,
    alpha = alpha,
    alternative = alternative,
    method = "One-sample t-test"
  ), class = "omni_sample_size")
}

#' Paired T-Test Sample Size
#'
#' Calculate sample size for paired t-test
#' @param delta Expected mean difference
#' @param sd_diff Standard deviation of differences
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @param alternative "two.sided" or "one.sided"
#' @return List with sample size results
#' @examples
#' ss_ttest_paired(delta = 5, sd_diff = 10)
#' @export
ss_ttest_paired <- function(delta, sd_diff, power = 0.8, alpha = 0.05,
                             alternative = "two.sided") {
  if (alternative == "two.sided") {
    z_alpha <- qnorm(1 - alpha/2)
  } else {
    z_alpha <- qnorm(1 - alpha)
  }
  z_beta <- qnorm(power)

  n <- ceiling(((z_alpha + z_beta) * sd_diff / delta)^2)

  structure(list(
    n = n,
    delta = delta,
    sd_diff = sd_diff,
    power = power,
    alpha = alpha,
    alternative = alternative,
    method = "Paired t-test"
  ), class = "omni_sample_size")
}

#' Two-Sample T-Test Sample Size (Welch's)
#'
#' Calculate sample size for Welch's t-test (unequal variances)
#' @param delta Expected difference between groups
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param ratio Allocation ratio n2/n1 (default 1)
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_ttest_welch(delta = 5, sd1 = 10, sd2 = 15, ratio = 1)
#' @export
ss_ttest_welch <- function(delta, sd1, sd2, ratio = 1, power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  # Variance term
  var_term <- sd1^2 + sd2^2 / ratio

  n1 <- ceiling((z_alpha + z_beta)^2 * var_term / delta^2)
  n2 <- ceiling(n1 * ratio)
  n_total <- n1 + n2

  structure(list(
    n1 = n1,
    n2 = n2,
    n_total = n_total,
    delta = delta,
    sd1 = sd1,
    sd2 = sd2,
    ratio = ratio,
    power = power,
    alpha = alpha,
    method = "Welch's t-test"
  ), class = "omni_sample_size")
}

#' Mann-Whitney U Test Sample Size
#'
#' Calculate sample size for Mann-Whitney U test
#' @param p Probability that observation from group 1 > observation from group 2
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_mann_whitney(p = 0.6)
#' @export
ss_mann_whitney <- function(p, power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  # Effect size
  omega <- p - 0.5

  n_per_group <- ceiling((z_alpha + z_beta)^2 / (12 * omega^2))
  n_total <- n_per_group * 2

  structure(list(
    n_per_group = n_per_group,
    n_total = n_total,
    p = p,
    omega = omega,
    power = power,
    alpha = alpha,
    method = "Mann-Whitney U test"
  ), class = "omni_sample_size")
}

#' McNemar's Test Sample Size
#'
#' Calculate sample size for paired proportions (McNemar's test)
#' @param p12 Proportion with response in treatment but not control
#' @param p21 Proportion with response in control but not treatment
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_mcnemar(p12 = 0.2, p21 = 0.1)
#' @export
ss_mcnemar <- function(p12, p21, power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  # Discordant pairs
  pd <- p12 + p21

  n <- ceiling((z_alpha * sqrt(pd) + z_beta * sqrt(pd - (p12 - p21)^2))^2 /
                 (p12 - p21)^2)

  structure(list(
    n = n,
    p12 = p12,
    p21 = p21,
    discordant = pd,
    power = power,
    alpha = alpha,
    method = "McNemar's test"
  ), class = "omni_sample_size")
}

# ============================================================================
# NON-INFERIORITY TRIALS (非劣效设计)
# ============================================================================

#' Non-inferiority Trial Sample Size (Means) - Full
#'
#' Calculate sample size for non-inferiority trial comparing two means
#' Based on Chow, Shao & Wang formula
#' @param margin Non-inferiority margin (positive value, clinically acceptable difference)
#' @param sd Standard deviation
#' @param true_diff Expected true difference (test - control, default 0)
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05, one-sided)
#' @param ratio Allocation ratio n2/n1 (default 1)
#' @return List with sample size results
#' @examples
#' # Non-inferiority margin = 3, SD=8
#' ss_noninferiority_mean(margin = 3, sd = 8)
#'
#' # Expected test is slightly worse by 1 unit
#' ss_noninferiority_mean(margin = 5, sd = 10, true_diff = 1)
#' @export
ss_noninferiority_mean <- function(margin, sd, true_diff = 0,
                                    power = 0.8, alpha = 0.05, ratio = 1) {
  z_alpha <- qnorm(1 - alpha)
  z_beta <- qnorm(power)

  # Effective margin
  eff_margin <- margin - abs(true_diff)
  if (eff_margin <= 0) {
    stop("Effective margin must be positive. Check margin and true_diff values.")
  }

  # Variance inflation for unequal allocation
  var_inflation <- (1 + 1/ratio)

  n1 <- ceiling((z_alpha + z_beta)^2 * sd^2 * var_inflation / eff_margin^2)
  n2 <- ceiling(n1 * ratio)
  n_total <- n1 + n2

  structure(list(
    n1 = n1,
    n2 = n2,
    n_total = n_total,
    margin = margin,
    effective_margin = eff_margin,
    sd = sd,
    true_diff = true_diff,
    ratio = ratio,
    power = power,
    alpha = alpha,
    method = "Non-inferiority test (means)"
  ), class = "omni_sample_size")
}

#' Non-inferiority Trial Sample Size (Proportions) - Full
#'
#' Calculate sample size for non-inferiority trial comparing two proportions
#' Based on Farrington-Manning score test
#' @param p_control Expected proportion in control group
#' @param margin Non-inferiority margin (positive)
#' @param p_test Expected proportion in test group (default = p_control)
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05, one-sided)
#' @param ratio Allocation ratio n2/n1 (default 1)
#' @return List with sample size results
#' @examples
#' # Control rate 80%, non-inferiority margin 10%
#' ss_noninferiority_prop(p_control = 0.8, margin = 0.1)
#' @export
ss_noninferiority_prop <- function(p_control, margin, p_test = NULL,
                                    power = 0.8, alpha = 0.05, ratio = 1) {
  if (is.null(p_test)) p_test <- p_control

  z_alpha <- qnorm(1 - alpha)
  z_beta <- qnorm(power)

  # Effective difference
  eff_diff <- margin - abs(p_test - p_control)
  if (eff_diff <= 0) {
    stop("Effective difference must be positive")
  }

  # Pooled proportion under null
  p0 <- (p_control + ratio * p_test) / (1 + ratio)

  # Variance terms
  var_null <- p0 * (1 - p0) * (1 + 1/ratio)
  var_alt <- p_control * (1 - p_control) + p_test * (1 - p_test) / ratio

  n1 <- ceiling((z_alpha * sqrt(var_null) + z_beta * sqrt(var_alt))^2 / eff_diff^2)
  n2 <- ceiling(n1 * ratio)
  n_total <- n1 + n2

  structure(list(
    n1 = n1,
    n2 = n2,
    n_total = n_total,
    p_control = p_control,
    p_test = p_test,
    margin = margin,
    effective_diff = eff_diff,
    ratio = ratio,
    power = power,
    alpha = alpha,
    method = "Non-inferiority test (proportions)"
  ), class = "omni_sample_size")
}

#' Non-inferiority Trial Sample Size (Survival)
#'
#' Calculate sample size for non-inferiority trial with time-to-event endpoint
#' @param hr_margin Non-inferiority margin for hazard ratio (e.g., 1.25)
#' @param hr_expected Expected hazard ratio (default 1.0)
#' @param p_event Probability of event
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_noninferiority_survival(hr_margin = 1.25, p_event = 0.6)
#' @export
ss_noninferiority_survival <- function(hr_margin, hr_expected = 1.0,
                                        p_event = 0.5, power = 0.8,
                                        alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha)
  z_beta <- qnorm(power)

  # Number of events needed
  log_margin <- log(hr_margin)
  log_expected <- log(hr_expected)

  d <- ceiling((z_alpha + z_beta)^2 / (log_margin - log_expected)^2)

  # Total sample size
  n_total <- ceiling(d / p_event)

  structure(list(
    n_total = n_total,
    events_needed = d,
    hr_margin = hr_margin,
    hr_expected = hr_expected,
    p_event = p_event,
    power = power,
    alpha = alpha,
    method = "Non-inferiority test (survival)"
  ), class = "omni_sample_size")
}

# ============================================================================
# EQUIVALENCE TRIALS (等效性设计)
# ============================================================================

#' Equivalence Trial Sample Size (Means) - Full
#'
#' Calculate sample size for equivalence trial comparing two means
#' Using two one-sided tests (TOST) approach
#' @param delta Equivalence margin (clinically acceptable difference)
#' @param sd Standard deviation
#' @param true_diff True difference (default 0)
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @param ratio Allocation ratio n2/n1 (default 1)
#' @return List with sample size results
#' @examples
#' # Equivalence margin = 5, SD=10
#' ss_equivalence_mean(delta = 5, sd = 10)
#' @export
ss_equivalence_mean <- function(delta, sd, true_diff = 0,
                                 power = 0.8, alpha = 0.05, ratio = 1) {
  # For equivalence, alpha is split between two one-sided tests
  z_alpha <- qnorm(1 - alpha)
  z_beta <- qnorm(power)

  # Effective margin
  eff_delta <- delta - abs(true_diff)
  if (eff_delta <= 0) {
    stop("Effective margin must be positive")
  }

  # Variance inflation
  var_inflation <- (1 + 1/ratio)

  n1 <- ceiling((z_alpha + z_beta)^2 * sd^2 * var_inflation / eff_delta^2)
  n2 <- ceiling(n1 * ratio)
  n_total <- n1 + n2

  structure(list(
    n1 = n1,
    n2 = n2,
    n_total = n_total,
    delta = delta,
    effective_delta = eff_delta,
    sd = sd,
    true_diff = true_diff,
    ratio = ratio,
    power = power,
    alpha = alpha,
    method = "Equivalence test (means)"
  ), class = "omni_sample_size")
}

#' Equivalence Trial Sample Size (Proportions)
#'
#' Calculate sample size for equivalence trial comparing two proportions
#' @param p1 Proportion in group 1
#' @param p2 Proportion in group 2
#' @param delta Equivalence margin
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_equivalence_prop(p1 = 0.8, p2 = 0.82, delta = 0.15)
#' @export
ss_equivalence_prop <- function(p1, p2, delta, power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha)
  z_beta <- qnorm(power)

  # Average proportion
  p_avg <- (p1 + p2) / 2

  # Effective difference
  eff_diff <- delta - abs(p1 - p2)
  if (eff_diff <= 0) {
    stop("Effective difference must be positive")
  }

  n_per_group <- ceiling(
    (z_alpha + z_beta)^2 * 2 * p_avg * (1 - p_avg) / eff_diff^2
  )
  n_total <- n_per_group * 2

  structure(list(
    n_per_group = n_per_group,
    n_total = n_total,
    p1 = p1,
    p2 = p2,
    delta = delta,
    effective_diff = eff_diff,
    power = power,
    alpha = alpha,
    method = "Equivalence test (proportions)"
  ), class = "omni_sample_size")
}

#' Equivalence Trial Sample Size (Ratio of Means)
#'
#' Calculate sample size for equivalence based on ratio of means (log-transformed)
#' Common for bioequivalence studies
#' @param theta True ratio (test/reference)
#' @param delta Bioequivalence margin (default log(1.25) ~ 0.2231, corresponding to 80-125 percent)
#' @param cv Coefficient of variation
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' # Bioequivalence: CV=20 percent, true ratio 95 percent
#' ss_equivalence_ratio(theta = 0.95, cv = 0.2)
#' @export
ss_equivalence_ratio <- function(theta, delta = log(1.25), cv,
                                  power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha)
  z_beta <- qnorm(power)

  # Log-scale parameters
  log_theta <- log(theta)
  sigma <- sqrt(log(1 + cv^2))

  # Sample size per sequence (for 2x2 crossover)
  n <- ceiling((z_alpha + z_beta)^2 * sigma^2 / (delta - abs(log_theta))^2)

  # Total for 2x2 crossover
  n_total <- 2 * n

  structure(list(
    n_per_sequence = n,
    n_total = n_total,
    theta = theta,
    delta = delta,
    cv = cv,
    sigma_log = sigma,
    power = power,
    alpha = alpha,
    method = "Equivalence test (ratio of means, 2x2 crossover)"
  ), class = "omni_sample_size")
}

# ============================================================================
# SUPERIORITY TRIALS (优效性设计)
# ============================================================================

#' Superiority Trial Sample Size (Means)
#'
#' Calculate sample size for superiority trial comparing two means
#' @param delta Superiority margin (minimum clinically important difference)
#' @param sd Standard deviation
#' @param true_diff Expected true difference
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @param ratio Allocation ratio n2/n1 (default 1)
#' @return List with sample size results
#' @examples
#' ss_superiority_mean(delta = 5, sd = 10, true_diff = 8)
#' @export
ss_superiority_mean <- function(delta, sd, true_diff,
                                 power = 0.8, alpha = 0.05, ratio = 1) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  # Effective difference
  eff_diff <- true_diff - delta
  if (eff_diff <= 0) {
    stop("True difference must be greater than superiority margin")
  }

  var_inflation <- (1 + 1/ratio)

  n1 <- ceiling((z_alpha + z_beta)^2 * sd^2 * var_inflation / eff_diff^2)
  n2 <- ceiling(n1 * ratio)
  n_total <- n1 + n2

  structure(list(
    n1 = n1,
    n2 = n2,
    n_total = n_total,
    delta = delta,
    true_diff = true_diff,
    effective_diff = eff_diff,
    sd = sd,
    ratio = ratio,
    power = power,
    alpha = alpha,
    method = "Superiority test (means)"
  ), class = "omni_sample_size")
}

#' Superiority Trial Sample Size (Proportions)
#'
#' Calculate sample size for superiority trial comparing two proportions
#' @param p_control Expected proportion in control
#' @param p_test Expected proportion in test
#' @param delta Superiority margin
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_superiority_prop(p_control = 0.6, p_test = 0.75, delta = 0.1)
#' @export
ss_superiority_prop <- function(p_control, p_test, delta,
                                 power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  # Effective difference
  eff_diff <- (p_test - p_control) - delta
  if (eff_diff <= 0) {
    stop("Test proportion must exceed control by more than margin")
  }

  p_avg <- (p_control + p_test) / 2

  n_per_group <- ceiling(
    (z_alpha * sqrt(2 * p_avg * (1 - p_avg)) +
       z_beta * sqrt(p_control * (1 - p_control) + p_test * (1 - p_test)))^2 /
      eff_diff^2
  )
  n_total <- n_per_group * 2

  structure(list(
    n_per_group = n_per_group,
    n_total = n_total,
    p_control = p_control,
    p_test = p_test,
    delta = delta,
    effective_diff = eff_diff,
    power = power,
    alpha = alpha,
    method = "Superiority test (proportions)"
  ), class = "omni_sample_size")
}

# ============================================================================
# SURVIVAL ANALYSIS
# ============================================================================

#' Survival Analysis Sample Size (Log-rank Test) - Enhanced
#'
#' Calculate sample size for comparing two survival curves using log-rank test
#' @param hr Hazard ratio (treatment/control)
#' @param p_t Proportion of patients in treatment group (default 0.5)
#' @param median_control Median survival time in control group
#' @param accrual Accrual period duration
#' @param follow_up Follow-up period duration
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @param alternative "two.sided" or "one.sided"
#' @return List with sample size results
#' @examples
#' ss_survival_logrank(hr = 0.7, median_control = 12, accrual = 12, follow_up = 18)
#' @export
ss_survival_logrank <- function(hr, p_t = 0.5, median_control = NULL,
                                 accrual = 12, follow_up = 12,
                                 power = 0.8, alpha = 0.05,
                                 alternative = "two.sided") {
  if (is.null(median_control)) {
    stop("median_control must be provided")
  }

  if (alternative == "two.sided") {
    z_alpha <- qnorm(1 - alpha/2)
  } else {
    z_alpha <- qnorm(1 - alpha)
  }
  z_beta <- qnorm(power)

  # Calculate event probabilities using exponential model
  lambda_c <- log(2) / median_control
  lambda_t <- lambda_c * hr

  # Overall event probability
  p_event <- 1 - (1 / (accrual * (lambda_c + lambda_t))) *
    (exp(-(accrual + follow_up) * lambda_c) - exp(-follow_up * lambda_c) +
       exp(-(accrual + follow_up) * lambda_t) - exp(-follow_up * lambda_t))

  # Number of events needed (Freedman formula)
  d <- ceiling((z_alpha + z_beta)^2 / (p_t * (1 - p_t) * (log(hr))^2))

  # Total sample size
  n_total <- ceiling(d / p_event)
  n_treatment <- ceiling(n_total * p_t)
  n_control <- n_total - n_treatment

  structure(list(
    n_total = n_total,
    n_treatment = n_treatment,
    n_control = n_control,
    events_needed = d,
    event_prob = p_event,
    hr = hr,
    median_control = median_control,
    accrual = accrual,
    follow_up = follow_up,
    power = power,
    alpha = alpha,
    method = "Log-rank test (survival)"
  ), class = "omni_sample_size")
}

#' Cox Regression Sample Size
#'
#' Calculate sample size for Cox proportional hazards regression
#' Based on Hsieh & Lavori (2000)
#' @param hr Hazard ratio for the predictor of interest
#' @param r2 R-squared from regression of predictor on other covariates
#' @param p_event Event probability (proportion experiencing event)
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_cox_regression(hr = 2.0, r2 = 0, p_event = 0.3)
#' @export
ss_cox_regression <- function(hr, r2 = 0, p_event = 0.5,
                               power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  # Events needed
  log_hr <- log(hr)
  d <- ceiling((z_alpha + z_beta)^2 / ((1 - r2) * p_event * (1 - p_event) * log_hr^2))

  # Total sample size
  n_total <- ceiling(d / p_event)

  structure(list(
    n_total = n_total,
    events_needed = d,
    p_event = p_event,
    hr = hr,
    r2 = r2,
    power = power,
    alpha = alpha,
    method = "Cox regression"
  ), class = "omni_sample_size")
}

# ============================================================================
# REGRESSION MODELS
# ============================================================================

#' Linear Regression Sample Size
#'
#' Calculate sample size for multiple linear regression
#' @param r2_squared Expected R-squared for full model
#' @param r2_reduced R-squared for reduced model (default 0)
#' @param n_predictors Number of predictors being tested
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_linear_regression(r2_squared = 0.15, n_predictors = 3)
#' @export
ss_linear_regression <- function(r2_squared, r2_reduced = 0, n_predictors = 1,
                                  power = 0.8, alpha = 0.05) {
  # Cohen's f2
  f2 <- (r2_squared - r2_reduced) / (1 - r2_squared)

  if (f2 <= 0) {
    stop("f2 must be positive. Check R-squared values.")
  }

  # Using power approximation
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  n <- ceiling((z_alpha + z_beta)^2 / f2) + n_predictors + 1

  structure(list(
    n = n,
    r2_squared = r2_squared,
    r2_reduced = r2_reduced,
    f2 = f2,
    n_predictors = n_predictors,
    power = power,
    alpha = alpha,
    method = "Linear regression"
  ), class = "omni_sample_size")
}

#' Logistic Regression Sample Size
#'
#' Calculate sample size for multiple logistic regression
#' Based on Hsieh, Bloch & Larsen (1998)
#' @param p Event probability
#' @param r2 R-squared from regression of predictor on other covariates
#' @param or Odds ratio for the predictor of interest
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_logistic_regression(p = 0.2, or = 2.0, r2 = 0)
#' @export
ss_logistic_regression <- function(p, or, r2 = 0,
                                    power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  log_or <- log(or)
  b <- p * (1 - p) * log_or^2

  if (b <= 0) {
    stop("Invalid parameters")
  }

  n <- ceiling((z_alpha + z_beta)^2 / ((1 - r2) * b))

  structure(list(
    n = n,
    p = p,
    or = or,
    r2 = r2,
    power = power,
    alpha = alpha,
    method = "Logistic regression"
  ), class = "omni_sample_size")
}

#' Poisson Regression Sample Size
#'
#' Calculate sample size for Poisson regression
#' @param lambda0 Baseline rate
#' @param rr Rate ratio of interest
#' @param exposure Mean exposure time
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_poisson_regression(lambda0 = 0.5, rr = 0.7, exposure = 1)
#' @export
ss_poisson_regression <- function(lambda0, rr, exposure = 1,
                                   power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  # Variance of log rate ratio
  var_log_rr <- (1 / (lambda0 * exposure) + 1 / (lambda0 * rr * exposure))

  n <- ceiling((z_alpha + z_beta)^2 * var_log_rr / log(rr)^2)

  structure(list(
    n = n,
    lambda0 = lambda0,
    rr = rr,
    exposure = exposure,
    power = power,
    alpha = alpha,
    method = "Poisson regression"
  ), class = "omni_sample_size")
}

# ============================================================================
# SPECIAL DESIGNS
# ============================================================================

#' Crossover Design Sample Size
#'
#' Calculate sample size for 2x2 crossover design
#' @param delta Expected treatment difference
#' @param sd_within Within-subject standard deviation
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @param alternative "two.sided" or "one.sided"
#' @return List with sample size results
#' @examples
#' ss_crossover(delta = 3, sd_within = 5)
#' @export
ss_crossover <- function(delta, sd_within, power = 0.8,
                          alpha = 0.05, alternative = "two.sided") {
  if (alternative == "two.sided") {
    z_alpha <- qnorm(1 - alpha/2)
  } else {
    z_alpha <- qnorm(1 - alpha)
  }
  z_beta <- qnorm(power)

  # Total sample size (both sequences)
  n_total <- ceiling(((z_alpha + z_beta)^2 * sd_within^2) / delta^2)

  # Ensure even number for equal allocation
  if (n_total %% 2 != 0) n_total <- n_total + 1

  structure(list(
    n_total = n_total,
    n_per_sequence = n_total / 2,
    delta = delta,
    sd_within = sd_within,
    power = power,
    alpha = alpha,
    alternative = alternative,
    method = "2x2 Crossover design"
  ), class = "omni_sample_size")
}

#' Group Sequential Design Sample Size
#'
#' Calculate sample size for group sequential trial
#' @param delta Effect size
#' @param sd Standard deviation
#' @param n_interim Number of planned interim analyses
#' @param alpha_spending Alpha spending function: "pocock" or "obf"
#' @param power Desired power (default 0.8)
#' @param alpha Overall significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_group_sequential(delta = 0.5, sd = 1, n_interim = 3)
#' @export
ss_group_sequential <- function(delta, sd, n_interim = 1,
                                 alpha_spending = "obf",
                                 power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  # Fixed sample size
  n_fixed <- ceiling(((z_alpha + z_beta)^2 * 2 * sd^2) / delta^2)

  # Inflation factor
  if (alpha_spending == "pocock") {
    inflation <- 1 + 0.3 * log(n_interim + 1)
  } else {
    inflation <- 1 + 0.1 * log(n_interim + 1)
  }

  n_total <- ceiling(n_fixed * inflation)
  n_per_group <- ceiling(n_total / 2)

  structure(list(
    n_per_group = n_per_group,
    n_total = n_total,
    n_fixed = n_fixed,
    inflation_factor = inflation,
    n_interim = n_interim,
    alpha_spending = alpha_spending,
    delta = delta,
    sd = sd,
    power = power,
    alpha = alpha,
    method = "Group sequential design"
  ), class = "omni_sample_size")
}

#' Repeated Measures ANOVA Sample Size
#'
#' Calculate sample size for repeated measures design
#' @param delta Expected difference between groups
#' @param sd_within Within-subject standard deviation
#' @param n_timepoints Number of time points
#' @param corr Within-subject correlation (default 0.5)
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_repeated_measures(delta = 3, sd_within = 5, n_timepoints = 4)
#' @export
ss_repeated_measures <- function(delta, sd_within, n_timepoints,
                                  corr = 0.5, power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  # Adjusted SD for repeated measures
  sd_adj <- sd_within * sqrt(1 - corr)

  n_per_group <- ceiling(((z_alpha + z_beta)^2 * 2 * sd_adj^2) / delta^2)
  n_total <- n_per_group * 2

  structure(list(
    n_per_group = n_per_group,
    n_total = n_total,
    delta = delta,
    sd_within = sd_within,
    sd_adjusted = sd_adj,
    n_timepoints = n_timepoints,
    corr = corr,
    power = power,
    alpha = alpha,
    method = "Repeated measures ANOVA"
  ), class = "omni_sample_size")
}

#' Cluster Randomized Trial Sample Size
#'
#' Calculate sample size for cluster randomized trial
#' @param delta Effect size
#' @param sd Standard deviation
#' @param icc Intraclass correlation coefficient
#' @param cluster_size Average cluster size
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_cluster_randomized(delta = 0.5, sd = 1, icc = 0.05, cluster_size = 20)
#' @export
ss_cluster_randomized <- function(delta, sd, icc, cluster_size,
                                   power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  # Design effect
  de <- 1 + (cluster_size - 1) * icc

  # Individual sample size
  n_individual <- ceiling(((z_alpha + z_beta)^2 * 2 * sd^2) / delta^2)

  # Adjusted for clustering
  n_total <- ceiling(n_individual * de)
  n_per_group <- ceiling(n_total / 2)

  # Number of clusters
  n_clusters <- ceiling(n_total / cluster_size)
  clusters_per_group <- ceiling(n_clusters / 2)

  structure(list(
    n_per_group = n_per_group,
    n_total = n_total,
    n_clusters = n_clusters,
    clusters_per_group = clusters_per_group,
    cluster_size = cluster_size,
    design_effect = de,
    icc = icc,
    delta = delta,
    power = power,
    alpha = alpha,
    method = "Cluster randomized trial"
  ), class = "omni_sample_size")
}

#' Multiple Comparisons Sample Size
#'
#' Calculate sample size with Bonferroni correction
#' @param k Number of comparisons
#' @param delta Effect size
#' @param sd Standard deviation
#' @param power Desired power (default 0.8)
#' @param alpha Overall significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_multiple_comparisons(k = 5, delta = 0.5, sd = 1)
#' @export
ss_multiple_comparisons <- function(k, delta, sd,
                                     power = 0.8, alpha = 0.05) {
  # Bonferroni corrected alpha
  alpha_adj <- alpha / k

  z_alpha <- qnorm(1 - alpha_adj/2)
  z_beta <- qnorm(power)

  n_per_group <- ceiling(((z_alpha + z_beta)^2 * 2 * sd^2) / delta^2)
  n_total <- n_per_group * 2

  structure(list(
    n_per_group = n_per_group,
    n_total = n_total,
    k = k,
    alpha_original = alpha,
    alpha_adjusted = alpha_adj,
    delta = delta,
    sd = sd,
    power = power,
    method = "Multiple comparisons (Bonferroni)"
  ), class = "omni_sample_size")
}

# ============================================================================
# CLINICAL PREDICTION MODELS
# ============================================================================

#' Clinical Prediction Model Sample Size (Binary)
#'
#' Calculate sample size for developing binary prediction model
#' Based on Riley et al. methodology
#' @param prevalence Event prevalence
#' @param cstatistic Expected C-statistic (AUC)
#' @param parameters Number of candidate parameters
#' @param shrinkage Desired shrinkage (default 0.9)
#' @param r2_pseudo Expected pseudo R-squared (optional)
#' @return List with sample size results
#' @examples
#' ss_prediction_binary(prevalence = 0.3, cstatistic = 0.75, parameters = 10)
#' @export
ss_prediction_binary <- function(prevalence, cstatistic, parameters,
                                  shrinkage = 0.9, r2_pseudo = NULL) {
  # Criterion 1: Shrinkage
  n_shrinkage <- ceiling(parameters / ((shrinkage - 1) * log(1 - prevalence / shrinkage)))

  # Criterion 2: Small absolute difference in predicted risk
  n_small_abs <- ceiling(parameters * 10 / prevalence)

  # Criterion 3: C-statistic (simplified)
  # Using rough approximation
  n_cstat <- ceiling(parameters * 15 / (prevalence * (1 - prevalence)))

  # Maximum of criteria
  n_total <- max(n_shrinkage, n_small_abs, n_cstat, parameters * 10)

  structure(list(
    n_total = n_total,
    n_events = ceiling(n_total * prevalence),
    n_per_parameter = ceiling(n_total / parameters),
    events_per_parameter = ceiling(n_total * prevalence / parameters),
    prevalence = prevalence,
    cstatistic = cstatistic,
    parameters = parameters,
    shrinkage = shrinkage,
    method = "Clinical prediction model (binary)"
  ), class = "omni_sample_size")
}

#' Clinical Prediction Model Sample Size (Continuous)
#'
#' Calculate sample size for developing continuous prediction model
#' @param r2 Expected R-squared
#' @param parameters Number of candidate parameters
#' @param shrinkage Desired shrinkage (default 0.9)
#' @return List with sample size results
#' @examples
#' ss_prediction_continuous(r2 = 0.3, parameters = 8)
#' @export
ss_prediction_continuous <- function(r2, parameters, shrinkage = 0.9) {
  # Riley et al. formula
  n_shrinkage <- ceiling(parameters / (1 - shrinkage))

  # For R-squared
  n_r2 <- ceiling(parameters * 10 / r2)

  n_total <- max(n_shrinkage, n_r2, parameters * 10)

  structure(list(
    n_total = n_total,
    n_per_parameter = ceiling(n_total / parameters),
    r2 = r2,
    parameters = parameters,
    shrinkage = shrinkage,
    method = "Clinical prediction model (continuous)"
  ), class = "omni_sample_size")
}

#' Clinical Prediction Model Sample Size (Survival)
#'
#' Calculate sample size for developing survival prediction model
#' @param cstatistic Expected C-index
#' @param parameters Number of candidate parameters
#' @param p_cens Censoring proportion
#' @param shrinkage Desired shrinkage (default 0.9)
#' @return List with sample size results
#' @examples
#' ss_prediction_survival(cstatistic = 0.7, parameters = 10, p_cens = 0.5)
#' @export
ss_prediction_survival <- function(cstatistic, parameters, p_cens = 0.5,
                                    shrinkage = 0.9) {
  # Events needed
  n_events <- ceiling(parameters * 10 / (1 - p_cens))

  # Total sample size
  n_total <- ceiling(n_events / (1 - p_cens))

  structure(list(
    n_total = n_total,
    n_events = n_events,
    n_per_parameter = ceiling(n_events / parameters),
    cstatistic = cstatistic,
    parameters = parameters,
    p_cens = p_cens,
    shrinkage = shrinkage,
    method = "Clinical prediction model (survival)"
  ), class = "omni_sample_size")
}

# ============================================================================
# DIAGNOSTIC TESTS
# ============================================================================

#' Diagnostic Test Sample Size (Sensitivity)
#'
#' Calculate sample size for estimating sensitivity
#' @param sens Expected sensitivity
#' @param w Margin of error (default 0.05)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_diagnostic_sens(sens = 0.9, w = 0.05)
#' @export
ss_diagnostic_sens <- function(sens, w = 0.05, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  n <- ceiling((z_alpha^2 * sens * (1 - sens)) / w^2)

  structure(list(
    n = n,
    sens = sens,
    w = w,
    alpha = alpha,
    method = "Diagnostic test (sensitivity)"
  ), class = "omni_sample_size")
}

#' Diagnostic Test Sample Size (Specificity)
#'
#' Calculate sample size for estimating specificity
#' @param spec Expected specificity
#' @param w Margin of error (default 0.05)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_diagnostic_spec(spec = 0.85, w = 0.05)
#' @export
ss_diagnostic_spec <- function(spec, w = 0.05, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  n <- ceiling((z_alpha^2 * spec * (1 - spec)) / w^2)

  structure(list(
    n = n,
    spec = spec,
    w = w,
    alpha = alpha,
    method = "Diagnostic test (specificity)"
  ), class = "omni_sample_size")
}

#' Diagnostic Test Sample Size (AUC)
#'
#' Calculate sample size for comparing AUC to null hypothesis (0.5)
#' @param auc Expected AUC
#' @param ratio Ratio of positive to negative cases
#' @param power Desired power (default 0.8)
#' @param alpha Significance level (default 0.05)
#' @return List with sample size results
#' @examples
#' ss_diagnostic_auc(auc = 0.8, ratio = 1)
#' @export
ss_diagnostic_auc <- function(auc, ratio = 1, power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)

  # Hanley-McNeil formula
  q1 <- auc / (2 - auc)
  q2 <- 2 * auc^2 / (1 + auc)

  n_pos <- ceiling((z_alpha * sqrt(q1) + z_beta * sqrt(q2))^2 / (auc - 0.5)^2)
  n_neg <- ceiling(n_pos * ratio)
  n_total <- n_pos + n_neg

  structure(list(
    n_total = n_total,
    n_positive = n_pos,
    n_negative = n_neg,
    auc = auc,
    ratio = ratio,
    power = power,
    alpha = alpha,
    method = "Diagnostic test (AUC)"
  ), class = "omni_sample_size")
}

# ============================================================================
# PRECISION ANALYSIS
# ============================================================================

#' Sample Size for Precision (Confidence Interval Width)
#'
#' Calculate sample size based on desired precision
#' @param sd Standard deviation
#' @param width Desired confidence interval width
#' @param conf_level Confidence level (default 0.95)
#' @return List with sample size results
#' @examples
#' ss_precision(sd = 10, width = 5)
#' @export
ss_precision <- function(sd, width, conf_level = 0.95) {
  alpha <- 1 - conf_level
  z <- qnorm(1 - alpha/2)

  n <- ceiling((2 * z * sd / width)^2)

  structure(list(
    n = n,
    sd = sd,
    width = width,
    conf_level = conf_level,
    method = "Precision-based sample size"
  ), class = "omni_sample_size")
}

#' Sample Size for Precision (Proportion)
#'
#' Calculate sample size for estimating proportion with desired precision
#' @param p Expected proportion
#' @param width Desired confidence interval width
#' @param conf_level Confidence level (default 0.95)
#' @return List with sample size results
#' @examples
#' ss_precision_prop(p = 0.5, width = 0.1)
#' @export
ss_precision_prop <- function(p, width, conf_level = 0.95) {
  alpha <- 1 - conf_level
  z <- qnorm(1 - alpha/2)

  n <- ceiling(4 * z^2 * p * (1 - p) / width^2)

  structure(list(
    n = n,
    p = p,
    width = width,
    conf_level = conf_level,
    method = "Precision-based sample size (proportion)"
  ), class = "omni_sample_size")
}
