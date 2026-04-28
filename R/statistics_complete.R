#' Complete Statistical Analysis Functions
#'
#' Advanced statistical functions for clinical research

# ============================================================================
# REGRESSION MODELS
# ============================================================================

#' Cox Proportional Hazards Regression
#'
#' Perform Cox regression analysis for survival data
#' @param formula Formula: Surv(time, status) ~ predictors
#' @param data Data frame
#' @return List with Cox regression results
#' @examples
#' \dontrun{
#' data(omni_cancer)
#' stats_cox_regression(Surv(time, status) ~ age + treatment, data = omni_cancer)
#' }
#' @export
stats_cox_regression <- function(formula, data) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("survival package required for Cox regression")
  }

  fit <- survival::coxph(formula, data = data)
  s <- summary(fit)

  # Extract coefficients
  coef_table <- as.data.frame(s$coefficients)
  confint <- as.data.frame(s$conf.int)

  structure(list(
    method = "Cox Proportional Hazards Regression",
    coefficients = coef_table,
    confint = confint,
    concordance = s$concordance["C"],
    loglik = fit$loglik,
    n = s$n,
    nevent = s$nevent,
    formula = deparse(formula),
    fit = fit
  ), class = "omni_stats")
}

#' Logistic Regression
#'
#' Perform logistic regression analysis
#' @param formula Model formula
#' @param data Data frame
#' @return List with logistic regression results
#' @examples
#' data(omni_glucose)
#' omni_glucose$poor_control <- ifelse(omni_glucose$hba1c_3m > 7, 1, 0)
#' stats_logistic_regression(poor_control ~ age + bmi, data = omni_glucose)
#' @export
stats_logistic_regression <- function(formula, data) {
  fit <- glm(formula, data = data, family = binomial(link = "logit"))
  s <- summary(fit)

  # Odds ratios
  coef_table <- as.data.frame(s$coefficients)
  or <- exp(coef_table[, "Estimate"])
  ci_lower <- exp(coef_table[, "Estimate"] - 1.96 * coef_table[, "Std. Error"])
  ci_upper <- exp(coef_table[, "Estimate"] + 1.96 * coef_table[, "Std. Error"])

  or_table <- data.frame(
    OR = or,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    p_value = coef_table[, "Pr(>|z|)"]
  )

  structure(list(
    method = "Logistic Regression",
    coefficients = coef_table,
    odds_ratios = or_table,
    aic = fit$aic,
    deviance = fit$deviance,
    null_deviance = fit$null.deviance,
    n = length(fit$y),
    converged = fit$converged,
    fit = fit
  ), class = "omni_stats")
}

#' Poisson Regression
#'
#' Perform Poisson regression for count data
#' @param formula Model formula
#' @param data Data frame
#' @param offset Optional offset variable
#' @return List with Poisson regression results
#' @examples
#' \dontrun{
#' stats_poisson_regression(counts ~ treatment + offset(log(time)), data = mydata)
#' }
#' @export
stats_poisson_regression <- function(formula, data, offset = NULL) {
  if (!is.null(offset)) {
    fit <- glm(formula, data = data, family = poisson(link = "log"), offset = offset)
  } else {
    fit <- glm(formula, data = data, family = poisson(link = "log"))
  }

  s <- summary(fit)

  structure(list(
    method = "Poisson Regression",
    coefficients = as.data.frame(s$coefficients),
    aic = fit$aic,
    deviance = fit$deviance,
    null_deviance = fit$null.deviance,
    n = length(fit$y),
    converged = fit$converged,
    fit = fit
  ), class = "omni_stats")
}

# ============================================================================
# SURVIVAL ANALYSIS
# ============================================================================

#' Kaplan-Meier Survival Analysis
#'
#' Calculate Kaplan-Meier survival estimates
#' @param time Survival time
#' @param status Event status (1 = event, 0 = censored)
#' @param group Optional grouping variable
#' @return List with KM estimates
#' @examples
#' data(omni_cancer)
#' stats_km_survival(omni_cancer$time, omni_cancer$status, omni_cancer$treatment)
#' @export
stats_km_survival <- function(time, status, group = NULL) {
  valid <- complete.cases(time, status)
  time <- time[valid]
  status <- status[valid]

  if (!is.null(group)) {
    group <- group[valid]
    groups <- unique(group)
    results <- list()

    for (g in groups) {
      g_idx <- group == g
      results[[as.character(g)]] <- calculate_km_internal(time[g_idx], status[g_idx])
    }

    structure(list(
      method = "Kaplan-Meier (by group)",
      groups = groups,
      estimates = results,
      n = length(time),
      events = sum(status)
    ), class = "omni_stats")
  } else {
    km <- calculate_km_internal(time, status)

    structure(list(
      method = "Kaplan-Meier",
      time = km$time,
      survival = km$survival,
      std_err = km$std_err,
      n = length(time),
      events = sum(status),
      median_survival = km$median
    ), class = "omni_stats")
  }
}

#' Internal KM calculation
calculate_km_internal <- function(time, status) {
  ord <- order(time)
  time <- time[ord]
  status <- status[ord]

  unique_times <- unique(time[status == 1])
  unique_times <- sort(unique_times)

  survival <- numeric(length(unique_times))
  std_err <- numeric(length(unique_times))
  n_risk <- numeric(length(unique_times))

  s <- 1
  se_sq <- 0

  for (i in seq_along(unique_times)) {
    t <- unique_times[i]
    at_risk <- sum(time >= t)
    events <- sum(time == t & status == 1)

    if (at_risk > 0 && events > 0) {
      s <- s * (1 - events / at_risk)
      se_sq <- se_sq + events / (at_risk * (at_risk - events))
    }

    survival[i] <- s
    std_err[i] <- s * sqrt(se_sq)
    n_risk[i] <- at_risk
  }

  # Median survival
  median_surv <- NA
  if (any(survival <= 0.5)) {
    idx <- which(survival <= 0.5)[1]
    if (idx == 1) {
      median_surv <- unique_times[1]
    } else {
      median_surv <- unique_times[idx - 1] +
        (0.5 - survival[idx - 1]) * (unique_times[idx] - unique_times[idx - 1]) /
        (survival[idx] - survival[idx - 1])
    }
  }

  list(
    time = unique_times,
    survival = survival,
    std_err = std_err,
    n_risk = n_risk,
    median = median_surv
  )
}

#' Log-rank Test
#'
#' Compare survival curves between groups
#' @param time Survival time
#' @param status Event status (1 = event, 0 = censored)
#' @param group Grouping variable
#' @return List with log-rank test results
#' @examples
#' data(omni_cancer)
#' stats_logrank_test(omni_cancer$time, omni_cancer$status, omni_cancer$treatment)
#' @export
stats_logrank_test <- function(time, status, group) {
  valid <- complete.cases(time, status, group)
  time <- time[valid]
  status <- status[valid]
  group <- as.factor(group[valid])

  if (requireNamespace("survival", quietly = TRUE)) {
    surv_obj <- survival::Surv(time, status)
    fit <- survival::survdiff(surv_obj ~ group)

    structure(list(
      method = "Log-rank test",
      statistic = fit$chisq,
      df = length(fit$n) - 1,
      p_value = pchisq(fit$chisq, df = length(fit$n) - 1, lower.tail = FALSE),
      observed = fit$obs,
      expected = fit$exp,
      groups = levels(group),
      n = length(time),
      events = sum(status)
    ), class = "omni_stats")
  } else {
    # Simple approximation for 2 groups
    groups <- levels(group)
    if (length(groups) != 2) {
      stop("Simple log-rank requires exactly 2 groups without survival package")
    }

    g1_idx <- group == groups[1]
    g2_idx <- group == groups[2]

    o1 <- sum(status[g1_idx])
    o2 <- sum(status[g2_idx])
    n1 <- sum(g1_idx)
    n2 <- sum(g2_idx)
    e1 <- sum(status) * n1 / length(status)
    e2 <- sum(status) * n2 / length(status)

    chisq <- ((o1 - e1)^2 / e1 + (o2 - e2)^2 / e2)
    p_value <- pchisq(chisq, df = 1, lower.tail = FALSE)

    structure(list(
      method = "Log-rank test (approximation)",
      statistic = chisq,
      df = 1,
      p_value = p_value,
      observed = c(o1, o2),
      expected = c(e1, e2),
      groups = groups,
      n = length(time),
      events = sum(status)
    ), class = "omni_stats")
  }
}

# ============================================================================
# EFFECT SIZE MEASURES
# ============================================================================

#' Odds Ratio
#'
#' Calculate odds ratio and confidence interval from 2x2 table
#' @param a Exposed cases
#' @param b Exposed non-cases
#' @param c Unexposed cases
#' @param d Unexposed non-cases
#' @return List with OR results
#' @examples
#' odds_ratio(a = 30, b = 20, c = 15, d = 35)
#' @export
odds_ratio <- function(a, b, c, d) {
  or <- (a * d) / (b * c)
  se_log_or <- sqrt(1/a + 1/b + 1/c + 1/d)
  ci_lower <- exp(log(or) - 1.96 * se_log_or)
  ci_upper <- exp(log(or) + 1.96 * se_log_or)

  # Chi-square test
  total <- a + b + c + d
  e_a <- (a + b) * (a + c) / total
  e_b <- (a + b) * (b + d) / total
  e_c <- (c + d) * (a + c) / total
  e_d <- (c + d) * (b + d) / total

  chisq <- sum((c(a, b, c, d) - c(e_a, e_b, e_c, e_d))^2 / c(e_a, e_b, e_c, e_d))
  p_value <- pchisq(chisq, df = 1, lower.tail = FALSE)

  structure(list(
    method = "Odds Ratio",
    or = or,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    log_or = log(or),
    se_log_or = se_log_or,
    chisq = chisq,
    p_value = p_value,
    table = matrix(c(a, b, c, d), nrow = 2,
                   dimnames = list(c("Case", "Control"),
                                   c("Exposed", "Unexposed")))
  ), class = "omni_stats")
}

#' Relative Risk
#'
#' Calculate relative risk and confidence interval
#' @param a Exposed cases
#' @param b Exposed total
#' @param c Unexposed cases
#' @param d Unexposed total
#' @return List with RR results
#' @examples
#' relative_risk(a = 30, b = 50, c = 15, d = 50)
#' @export
relative_risk <- function(a, b, c, d) {
  p1 <- a / b
  p2 <- c / d
  rr <- p1 / p2

  se_log_rr <- sqrt((1 - p1) / (b * p1) + (1 - p2) / (d * p2))
  ci_lower <- exp(log(rr) - 1.96 * se_log_rr)
  ci_upper <- exp(log(rr) + 1.96 * se_log_rr)

  structure(list(
    method = "Relative Risk",
    rr = rr,
    p1 = p1,
    p2 = p2,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    risk_difference = p1 - p2,
    nnt = ifelse(abs(p1 - p2) > 0, ceiling(1 / abs(p1 - p2)), Inf)
  ), class = "omni_stats")
}

#' Hazard Ratio
#'
#' Calculate hazard ratio from survival data
#' @param time Survival time
#' @param status Event status
#' @param group Binary grouping variable
#' @return List with HR results
#' @examples
#' data(omni_cancer)
#' hazard_ratio(omni_cancer$time, omni_cancer$status, omni_cancer$treatment)
#' @export
hazard_ratio <- function(time, status, group) {
  valid <- complete.cases(time, status, group)
  time <- time[valid]
  status <- status[valid]
  group <- as.factor(group[valid])

  groups <- levels(group)
  if (length(groups) != 2) {
    stop("Hazard ratio requires exactly 2 groups")
  }

  g1_idx <- group == groups[1]
  g2_idx <- group == groups[2]

  # Person-time
  pt1 <- sum(time[g1_idx])
  pt2 <- sum(time[g2_idx])

  # Events
  e1 <- sum(status[g1_idx])
  e2 <- sum(status[g2_idx])

  # Hazard rates
  h1 <- e1 / pt1
  h2 <- e2 / pt2

  hr <- h1 / h2
  se_log_hr <- sqrt(1/e1 + 1/e2)
  ci_lower <- exp(log(hr) - 1.96 * se_log_hr)
  ci_upper <- exp(log(hr) + 1.96 * se_log_hr)

  structure(list(
    method = "Hazard Ratio",
    hr = hr,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    hazard1 = h1,
    hazard2 = h2,
    events = c(e1, e2),
    person_time = c(pt1, pt2),
    groups = groups
  ), class = "omni_stats")
}

#' Number Needed to Treat
#'
#' Calculate NNT with confidence interval
#' @param p_treatment Proportion with outcome in treatment group
#' @param p_control Proportion with outcome in control group
#' @return List with NNT results
#' @examples
#' nnt(p_treatment = 0.15, p_control = 0.25)
#' @export
nnt <- function(p_treatment, p_control) {
  arr <- p_control - p_treatment  # Absolute risk reduction
  nnt_val <- 1 / arr

  # CI for ARR
  se_arr <- sqrt(p_treatment * (1 - p_treatment) / 100 +
                   p_control * (1 - p_control) / 100)
  ci_lower_arr <- arr - 1.96 * se_arr
  ci_upper_arr <- arr + 1.96 * se_arr

  structure(list(
    method = "Number Needed to Treat",
    nnt = nnt_val,
    arr = arr,
    ci_lower = ifelse(ci_upper_arr > 0, 1 / ci_upper_arr, Inf),
    ci_upper = ifelse(ci_lower_arr > 0, 1 / ci_lower_arr, Inf),
    p_treatment = p_treatment,
    p_control = p_control
  ), class = "omni_stats")
}

# ============================================================================
# DIAGNOSTIC TESTS
# ============================================================================

#' Diagnostic Test Performance
#'
#' Calculate sensitivity, specificity, PPV, NPV
#' @param test Test result (positive/negative)
#' @param reference Reference standard (positive/negative)
#' @return List with diagnostic performance metrics
#' @examples
#' data(omni_diagnostic)
#' diag_test_performance(omni_diagnostic$test1_result,
#'                       omni_diagnostic$true_disease)
#' @export
diag_test_performance <- function(test, reference) {
  valid <- complete.cases(test, reference)
  test <- test[valid]
  reference <- reference[valid]

  # Convert to binary if needed
  if (is.character(test) || is.factor(test)) {
    test_levels <- unique(test)
    if (length(test_levels) != 2) {
      stop("Test must have exactly 2 levels")
    }
    test_bin <- as.numeric(factor(test)) - 1
  } else {
    test_bin <- test
  }

  if (is.character(reference) || is.factor(reference)) {
    ref_levels <- unique(reference)
    if (length(ref_levels) != 2) {
      stop("Reference must have exactly 2 levels")
    }
    ref_bin <- as.numeric(factor(reference)) - 1
  } else {
    ref_bin <- reference
  }

  # 2x2 table
  tp <- sum(test_bin == 1 & ref_bin == 1)
  fp <- sum(test_bin == 1 & ref_bin == 0)
  fn <- sum(test_bin == 0 & ref_bin == 1)
  tn <- sum(test_bin == 0 & ref_bin == 0)

  # Metrics
  sens <- tp / (tp + fn)
  spec <- tn / (tn + fp)
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)
  prev <- (tp + fn) / length(test_bin)

  # Accuracy
  acc <- (tp + tn) / length(test_bin)

  # Likelihood ratios
  lr_pos <- sens / (1 - spec)
  lr_neg <- (1 - sens) / spec

  structure(list(
    method = "Diagnostic Test Performance",
    sensitivity = sens,
    specificity = spec,
    ppv = ppv,
    npv = npv,
    prevalence = prev,
    accuracy = acc,
    lr_positive = lr_pos,
    lr_negative = lr_neg,
    table = matrix(c(tp, fp, fn, tn), nrow = 2,
                   dimnames = list(c("Test+", "Test-"),
                                   c("Disease+", "Disease-")))
  ), class = "omni_stats")
}

# ============================================================================
# ADDITIONAL STATISTICAL TESTS
# ============================================================================

#' Fisher's Exact Test
#'
#' Perform Fisher's exact test for 2x2 tables
#' @param x Factor or table
#' @param y Factor (optional)
#' @return List with test results
#' @examples
#' data(omni_glucose)
#' stats_fisher_exact(omni_glucose$gender, omni_glucose$complications)
#' @export
stats_fisher_exact <- function(x, y = NULL) {
  if (!is.null(y)) {
    tbl <- table(x, y)
  } else {
    tbl <- x
  }

  result <- fisher.test(tbl, workspace = 2e7)

  structure(list(
    method = "Fisher's Exact Test",
    p_value = result$p.value,
    estimate = result$estimate,
    conf_int = result$conf.int,
    observed = tbl
  ), class = "omni_stats")
}

#' McNemar's Test
#'
#' Perform McNemar's test for paired proportions
#' @param x Before treatment
#' @param y After treatment
#' @return List with test results
#' @examples
#' \dontrun{
#' stats_mcnemar(before, after)
#' }
#' @export
stats_mcnemar <- function(x, y) {
  tbl <- table(x, y)

  if (nrow(tbl) != 2 || ncol(tbl) != 2) {
    stop("McNemar's test requires 2x2 table")
  }

  result <- mcnemar.test(tbl)

  structure(list(
    method = "McNemar's Test",
    statistic = result$statistic,
    p_value = result$p.value,
    observed = tbl
  ), class = "omni_stats")
}

#' Kappa Statistic
#'
#' Calculate Cohen's kappa for inter-rater agreement
#' @param rater1 Rater 1 ratings
#' @param rater2 Rater 2 ratings
#' @return List with kappa results
#' @examples
#' \dontrun{
#' stats_kappa(rater1, rater2)
#' }
#' @export
stats_kappa <- function(rater1, rater2) {
  tbl <- table(rater1, rater2)

  # Observed agreement
  p_o <- sum(diag(tbl)) / sum(tbl)

  # Expected agreement
  row_margins <- rowSums(tbl) / sum(tbl)
  col_margins <- colSums(tbl) / sum(tbl)
  p_e <- sum(row_margins * col_margins)

  # Kappa
  kappa <- (p_o - p_e) / (1 - p_e)

  # Standard error (simplified)
  se_kappa <- sqrt(p_o * (1 - p_o) / (sum(tbl) * (1 - p_e)^2))

  structure(list(
    method = "Cohen's Kappa",
    kappa = kappa,
    se = se_kappa,
    ci_lower = kappa - 1.96 * se_kappa,
    ci_upper = kappa + 1.96 * se_kappa,
    observed_agreement = p_o,
    expected_agreement = p_e,
    table = tbl
  ), class = "omni_stats")
}

#' Intraclass Correlation Coefficient
#'
#' Calculate ICC for reliability
#' @param data Data frame with measurements
#' @param method "oneway" or "twoway"
#' @return List with ICC results
#' @examples
#' \dontrun{
#' stats_icc(measurements_df)
#' }
#' @export
stats_icc <- function(data, method = "oneway") {
  # Simplified ICC calculation
  n_subjects <- nrow(data)
  n_ratings <- ncol(data)

  # Subject means
  subject_means <- rowMeans(data, na.rm = TRUE)
  grand_mean <- mean(subject_means, na.rm = TRUE)

  # Mean squares
  ms_between <- var(subject_means, na.rm = TRUE) * n_ratings
  ms_within <- mean(apply(data, 1, var, na.rm = TRUE), na.rm = TRUE)

  # ICC
  icc <- (ms_between - ms_within) / (ms_between + (n_ratings - 1) * ms_within)

  structure(list(
    method = paste("ICC (", method, ")", sep = ""),
    icc = icc,
    n_subjects = n_subjects,
    n_ratings = n_ratings,
    ms_between = ms_between,
    ms_within = ms_within
  ), class = "omni_stats")
}
