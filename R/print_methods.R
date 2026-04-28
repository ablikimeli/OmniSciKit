#' @export
print.omni_sample_size <- function(x, ...) {
  cat("\n=== Sample Size Calculation ===\n")
  cat("Method:", x$method, "\n")

  if (!is.null(x$n_per_group)) {
    cat("Sample size per group:", x$n_per_group, "\n")
    if (!is.null(x$n_total)) cat("Total sample size:", x$n_total, "\n")
  } else if (!is.null(x$n_total)) {
    cat("Total sample size:", x$n_total, "\n")
  } else if (!is.null(x$n)) {
    cat("Sample size:", x$n, "\n")
  }
  if (!is.null(x$n_disease)) {
    cat("Disease cases:", x$n_disease, "\n")
  }
  if (!is.null(x$n_healthy)) {
    cat("Healthy controls:", x$n_healthy, "\n")
  }
  if (!is.null(x$n1) && !is.null(x$n2)) {
    cat("Group 1:", x$n1, "Group 2:", x$n2, "\n")
  }
  if (!is.null(x$events_needed)) {
    cat("Events needed:", x$events_needed, "\n")
  }
  if (!is.null(x$n_events)) {
    cat("Events needed:", x$n_events, "\n")
  }
  if (!is.null(x$n_interim)) {
    cat("Interim analyses:", x$n_interim, "\n")
  }
  if (!is.null(x$n_clusters)) {
    cat("Clusters:", x$n_clusters, "\n")
  }
  if (!is.null(x$cluster_size)) {
    cat("Cluster size:", x$cluster_size, "\n")
  }
  if (!is.null(x$n_per_group)) {
    cat("Total sample size:", x$n_total, "\n")
  }
  if (!is.null(x$power)) cat("Power:", x$power, "\n")
  if (!is.null(x$alpha)) cat("Alpha:", x$alpha, "\n")

  # Formula and reference
  if (!is.null(x$formula)) cat("\nFormula:", x$formula, "\n")
  if (!is.null(x$reference)) cat("Reference:", x$reference, "\n")

  invisible(x)
}

#' @export
print.omni_stats <- function(x, ...) {
  cat("\n=== Statistical Test Results ===\n")

  if (!is.null(x$method)) cat("Method:", x$method, "\n")

  # General test statistics
  if (!is.null(x$statistic)) cat("Statistic:", round(as.numeric(x$statistic), 4), "\n")
  if (!is.null(x$p_value)) {
    if (x$p_value < 0.001) cat("P-value: < 0.001\n")
    else cat("P-value:", round(x$p_value, 4), "\n")
  }
  if (!is.null(x$estimate)) cat("Estimate:", round(as.numeric(x$estimate), 4), "\n")

  # Regression
  if (!is.null(x$r_squared)) cat("R-squared:", round(x$r_squared, 4), "\n")
  if (!is.null(x$adj_r_squared)) cat("Adjusted R-squared:", round(x$adj_r_squared, 4), "\n")
  if (!is.null(x$f_value)) cat("F-value:", round(x$f_value, 4), "\n")
  if (!is.null(x$converged)) cat("Converged:", x$converged, "\n")
  if (!is.null(x$concordance)) cat("Concordance:", round(x$concordance, 4), "\n")
  if (!is.null(x$n) && is.null(x$sensitivity)) cat("N:", x$n, "\n")
  if (!is.null(x$nevent)) cat("Events:", x$nevent, "\n")

  # Diagnostic test
  if (!is.null(x$sensitivity)) cat("Sensitivity:", round(x$sensitivity, 4), "\n")
  if (!is.null(x$specificity)) cat("Specificity:", round(x$specificity, 4), "\n")
  if (!is.null(x$ppv)) cat("PPV:", round(x$ppv, 4), "\n")
  if (!is.null(x$npv)) cat("NPV:", round(x$npv, 4), "\n")
  if (!is.null(x$prevalence)) cat("Prevalence:", round(x$prevalence, 4), "\n")
  if (!is.null(x$accuracy)) cat("Accuracy:", round(x$accuracy, 4), "\n")
  if (!is.null(x$lr_positive)) cat("LR+:", round(x$lr_positive, 4), "\n")
  if (!is.null(x$lr_negative)) cat("LR-:", round(x$lr_negative, 4), "\n")

  # Kappa / agreement
  if (!is.null(x$kappa)) cat("Kappa:", round(x$kappa, 4), "\n")
  if (!is.null(x[["se"]])) cat("SE:", round(x[["se"]], 4), "\n")
  if (!is.null(x$ci_lower) && !is.null(x$ci_upper)) {
    cat("95% CI:", round(x$ci_lower, 4), "-", round(x$ci_upper, 4), "\n")
  }
  if (!is.null(x$observed_agreement)) cat("Observed agreement:", round(x$observed_agreement, 4), "\n")
  if (!is.null(x$expected_agreement)) cat("Expected agreement:", round(x$expected_agreement, 4), "\n")

  # Bootstrap
  if (!is.null(x$bootstrap_mean)) cat("Bootstrap mean:", round(x$bootstrap_mean, 4), "\n")
  if (!is.null(x$bootstrap_sd)) cat("Bootstrap SE:", round(x$bootstrap_sd, 4), "\n")
  if (!is.null(x$n_bootstrap)) cat("Bootstrap replicates:", x$n_bootstrap, "\n")

  # Odds ratio / risk ratio
  if (!is.null(x[["or"]])) cat("Odds ratio:", round(x[["or"]], 4), "\n")
  if (!is.null(x[["odds_ratios"]]) && is.data.frame(x[["odds_ratios"]])) {
    cat("\nOdds Ratios:\n")
    print(x[["odds_ratios"]])
  }
  if (!is.null(x$relative_risk)) cat("Relative risk:", round(x$relative_risk, 4), "\n")
  if (!is.null(x$hazard_ratio)) cat("Hazard ratio:", round(x$hazard_ratio, 4), "\n")
  if (!is.null(x$nnt)) cat("NNT:", round(x$nnt, 4), "\n")

  # Survival
  if (!is.null(x$groups)) cat("Groups:", length(x$groups), "\n")

  # Alpha adjustment
  if (!is.null(x$alpha_original)) cat("Original alpha:", x$alpha_original, "\n")
  if (!is.null(x$alpha_adjusted)) cat("Adjusted alpha:", round(x$alpha_adjusted, 6), "\n")
  if (!is.null(x$k)) cat("Number of comparisons:", x$k, "\n")

  # Coefficients table for regression
  if (!is.null(x$coefficients) && is.data.frame(x$coefficients)) {
    cat("\nCoefficients:\n")
    print(x$coefficients)
  }

  # Confusion / agreement table
  if (!is.null(x$table)) {
    cat("\n")
    print(x$table)
  }

  # Formula and reference
  if (!is.null(x$formula)) cat("\nFormula:", x$formula, "\n")
  if (!is.null(x$reference)) cat("Reference:", x$reference, "\n")

  invisible(x)
}

#' @export
print.omni_power <- function(x, ...) {
  cat("\n=== Power Analysis ===\n")
  cat("Method:", x$method, "\n")
  cat("Power:", round(x$power, 4), "\n")

  if (!is.null(x$n)) cat("Sample size:", x$n, "\n")
  if (!is.null(x$n_total)) cat("Total sample size:", x$n_total, "\n")
  if (!is.null(x$n_per_group)) cat("Sample size per group:", x$n_per_group, "\n")
  if (!is.null(x$delta)) cat("Effect size (delta):", x$delta, "\n")
  if (!is.null(x$r)) cat("Correlation (r):", x$r, "\n")
  if (!is.null(x$f)) cat("Effect size (f):", x$f, "\n")
  if (!is.null(x$hr)) cat("Hazard ratio:", x$hr, "\n")
  if (!is.null(x$expected_events)) cat("Expected events:", x$expected_events, "\n")
  if (!is.null(x$p_event)) cat("Event probability:", x$p_event, "\n")
  if (!is.null(x$n_clusters)) cat("Number of clusters:", x$n_clusters, "\n")
  if (!is.null(x$alpha)) cat("Alpha:", x$alpha, "\n")
  invisible(x)
}

#' @export
print.omni_effect_size <- function(x, ...) {
  cat("\n=== Effect Size ===\n")

  if (!is.null(x$d)) {
    cat("Cohen's d:", round(x$d, 4), "\n")
    cat("|d|:", round(x$abs_d, 4), "\n")
  }
  if (!is.null(x$f)) {
    cat("Cohen's f:", round(x$f, 4), "\n")
  }
  if (!is.null(x$f2)) {
    cat("Cohen's f-squared:", round(x$f2, 4), "\n")
  }
  if (!is.null(x[["or"]])) {
    cat("Odds ratio:", round(x[["or"]], 4), "\n")
    if (!is.null(x$ci_lower) && !is.null(x$ci_upper)) {
      cat("95% CI:", round(x$ci_lower, 4), "-", round(x$ci_upper, 4), "\n")
    }
  }
  if (!is.null(x$relative_risk)) {
    cat("Relative risk:", round(x$relative_risk, 4), "\n")
  }
  if (!is.null(x$hazard_ratio)) {
    cat("Hazard ratio:", round(x$hazard_ratio, 4), "\n")
  }

  if (!is.null(x$interpretation)) cat("Interpretation:", x$interpretation, "\n")
  invisible(x)
}

#' @export
print.omni_clinical <- function(x, ...) {
  cat("\n=== Clinical Calculation ===\n")

  if (!is.null(x$bmi)) {
    cat("BMI:", x$bmi, "kg/m²\n")
    cat("Category:", x$category, "\n")
  }
  if (!is.null(x$egfr)) {
    cat("eGFR:", x$egfr, "mL/min/1.73m²\n")
    cat("CKD Stage:", x$stage, "\n")
  }

  invisible(x)
}

#' @export
print.omni_roc <- function(x, ...) {
  cat("\n=== ROC Analysis ===\n")
  cat("AUC:", round(x$auc, 4), "\n")

  if (x$auc >= 0.9)      interp <- "Excellent"
  else if (x$auc >= 0.8) interp <- "Good"
  else if (x$auc >= 0.7) interp <- "Fair"
  else if (x$auc >= 0.6) interp <- "Poor"
  else                   interp <- "Fail"

  cat("Interpretation:", interp, "\n")

  if (!is.null(x$plot)) {
    plot(x$plot)
  }

  invisible(x)
}

#' @export
print.omni_table_one <- function(x, ...) {
  cat("\n=== Table 1: Baseline Characteristics ===\n\n")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}
