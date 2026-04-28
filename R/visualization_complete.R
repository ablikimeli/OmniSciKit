#' Complete Visualization Functions for Clinical Research
#'
#' Advanced plotting functions for clinical data analysis

# ============================================================================
# SURVIVAL PLOTS
# ============================================================================

#' Enhanced Kaplan-Meier Plot
#'
#' Plot Kaplan-Meier survival curves with confidence intervals
#' @param time Survival time
#' @param status Event status (1 = event, 0 = censored)
#' @param group Optional grouping variable
#' @param title Plot title
#' @param conf_int Show confidence intervals (default TRUE)
#' @param risk_table Show risk table (default FALSE)
#' @return ggplot object
#' @examples
#' data(omni_cancer)
#' plot_km_curve(omni_cancer$time, omni_cancer$status, omni_cancer$treatment)
#' @export
plot_km_curve <- function(time, status, group = NULL,
                          title = "Kaplan-Meier Survival Curve",
                          conf_int = TRUE, risk_table = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  valid <- complete.cases(time, status)
  time <- time[valid]
  status <- status[valid]

  if (!is.null(group)) {
    group <- group[valid]
    groups <- unique(group)
    plot_data <- data.frame()

    for (g in groups) {
      g_idx <- group == g
      km <- calculate_km_for_plot(time[g_idx], status[g_idx])

      df <- data.frame(
        time = km$time,
        survival = km$survival,
        lower = pmax(0, km$survival - 1.96 * km$std_err),
        upper = pmin(1, km$survival + 1.96 * km$std_err),
        group = g
      )
      plot_data <- rbind(plot_data, df)
    }

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = survival,
                                                   color = group, fill = group)) +
      ggplot2::geom_step(linewidth = 1)

    if (conf_int) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                                     alpha = 0.2, color = NA, na.rm = TRUE)
    }
  } else {
    km <- calculate_km_for_plot(time, status)
    plot_data <- data.frame(
      time = km$time,
      survival = km$survival,
      lower = pmax(0, km$survival - 1.96 * km$std_err),
      upper = pmin(1, km$survival + 1.96 * km$std_err)
    )

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = survival)) +
      ggplot2::geom_step(color = "steelblue", linewidth = 1)

    if (conf_int) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                                     alpha = 0.2, fill = "steelblue", na.rm = TRUE)
    }
  }

  p + ggplot2::labs(title = title,
                    x = "Time",
                    y = "Survival Probability") +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_minimal()
}

#' Internal KM calculation for plotting
calculate_km_for_plot <- function(time, status) {
  ord <- order(time)
  time <- time[ord]
  status <- status[ord]

  unique_times <- unique(time[status == 1])
  unique_times <- sort(unique_times)

  survival <- numeric(length(unique_times))
  std_err <- numeric(length(unique_times))

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
  }

  list(time = unique_times, survival = survival, std_err = std_err)
}

#' Cumulative Incidence Plot
#'
#' Plot cumulative incidence function for competing risks
#' @param time Time to event
#' @param status Event status (0 = censored, 1 = event of interest, 2 = competing risk)
#' @param group Optional grouping variable
#' @param title Plot title
#' @return ggplot object
#' @examples
#' \dontrun{
#' plot_cumulative_incidence(time, status, group)
#' }
#' @export
plot_cumulative_incidence <- function(time, status, group = NULL,
                                       title = "Cumulative Incidence") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  # Simple cumulative incidence (1 - KM for competing risk)
  event_1 <- as.numeric(status == 1)

  if (!is.null(group)) {
    groups <- unique(group)
    plot_data <- data.frame()

    for (g in groups) {
      g_idx <- group == g
      km <- calculate_km_for_plot(time[g_idx], event_1[g_idx])

      df <- data.frame(
        time = km$time,
        incidence = 1 - km$survival,
        group = g
      )
      plot_data <- rbind(plot_data, df)
    }

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = incidence,
                                                   color = group)) +
      ggplot2::geom_step(linewidth = 1)
  } else {
    km <- calculate_km_for_plot(time, event_1)
    plot_data <- data.frame(
      time = km$time,
      incidence = 1 - km$survival
    )

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = incidence)) +
      ggplot2::geom_step(color = "steelblue", linewidth = 1)
  }

  p + ggplot2::labs(title = title,
                    x = "Time",
                    y = "Cumulative Incidence") +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_minimal()
}

# ============================================================================
# FOREST PLOTS
# ============================================================================

#' Forest Plot for Meta-Analysis
#'
#' Create forest plot with multiple studies
#' @param studies Study names
#' @param effects Effect sizes (e.g., log HR, log OR)
#' @param se Standard errors
#' @param subgroups Optional subgroup labels
#' @param title Plot title
#' @param summary_label Label for summary diamond
#' @return ggplot object
#' @examples
#' plot_forest_meta(
#'   studies = c("Study A", "Study B", "Study C"),
#'   effects = c(-0.5, -0.3, -0.7),
#'   se = c(0.2, 0.15, 0.25)
#' )
#' @export
plot_forest_meta <- function(studies, effects, se, subgroups = NULL,
                             title = "Forest Plot",
                             summary_label = "Summary") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  ci_lower <- effects - 1.96 * se
  ci_upper <- effects + 1.96 * se

  # Calculate summary effect (fixed effect)
  weights <- 1 / se^2
  summary_effect <- sum(weights * effects) / sum(weights)
  summary_se <- sqrt(1 / sum(weights))
  summary_ci_lower <- summary_effect - 1.96 * summary_se
  summary_ci_upper <- summary_effect + 1.96 * summary_se

  plot_data <- data.frame(
    study = factor(studies, levels = rev(studies)),
    effect = effects,
    lower = ci_lower,
    upper = ci_upper,
    se = se,
    weight = weights
  )

  if (!is.null(subgroups)) {
    plot_data$subgroup <- subgroups
  }

  # Add summary row
  summary_row <- data.frame(
    study = factor(summary_label, levels = rev(c(studies, summary_label))),
    effect = summary_effect,
    lower = summary_ci_lower,
    upper = summary_ci_upper,
    se = summary_se,
    weight = sum(weights),
    subgroup = if (!is.null(subgroups)) "Summary" else NULL
  )

  plot_data <- rbind(plot_data, summary_row)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = study)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50")

  if (!is.null(subgroups)) {
    p <- p + ggplot2::facet_grid(subgroup ~ ., scales = "free_y", space = "free")
  }

  p + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper),
                               height = 0.2, color = "steelblue") +
    ggplot2::geom_point(ggplot2::aes(x = effect, size = weight),
                        color = "steelblue") +
    ggplot2::scale_size_continuous(name = "Weight", range = c(2, 8)) +
    ggplot2::labs(title = title,
                  x = "Effect Size",
                  y = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

#' Subgroup Analysis Forest Plot
#'
#' Create forest plot for subgroup analysis
#' @param subgroups Subgroup names
#' @param effects Effect sizes within subgroups
#' @param ci_lower Lower CI bounds
#' @param ci_upper Upper CI bounds
#' @param p_interaction P-value for interaction
#' @param title Plot title
#' @return ggplot object
#' @examples
#' \dontrun{
#' plot_subgroup_forest(subgroups, effects, ci_lower, ci_upper)
#' }
#' @export
plot_subgroup_forest <- function(subgroups, effects, ci_lower, ci_upper,
                                  p_interaction = NULL,
                                  title = "Subgroup Analysis") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  plot_data <- data.frame(
    subgroup = factor(subgroups, levels = rev(subgroups)),
    effect = effects,
    lower = ci_lower,
    upper = ci_upper
  )

  subtitle <- if (!is.null(p_interaction)) {
    paste("P-interaction =", format(p_interaction, digits = 3))
  } else {
    NULL
  }

  ggplot2::ggplot(plot_data, ggplot2::aes(y = subgroup)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper),
                            height = 0.2, color = "steelblue") +
    ggplot2::geom_point(ggplot2::aes(x = effect), size = 3, color = "steelblue") +
    ggplot2::labs(title = title,
                  subtitle = subtitle,
                  x = "Effect Size",
                  y = "") +
    ggplot2::theme_minimal()
}

# ============================================================================
# DIAGNOSTIC PLOTS
# ============================================================================

#' ROC Curve - Enhanced
#'
#' Plot ROC curve with confidence interval and optimal cutoff
#' @param predictor Predictor variable
#' @param response Response variable (binary)
#' @param title Plot title
#' @param ci Show confidence interval (default TRUE)
#' @return List with AUC and plot
#' @examples
#' data(omni_diagnostic)
#' disease_binary <- ifelse(omni_diagnostic$true_disease == "阳性", 1, 0)
#' plot_roc_enhanced(omni_diagnostic$test1_result, disease_binary)
#' @export
plot_roc_enhanced <- function(predictor, response, title = "ROC Curve",
                               ci = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  valid <- complete.cases(predictor, response)
  predictor <- predictor[valid]
  response <- response[valid]

  # Convert to 0/1
  if (is.logical(response)) {
    response <- as.numeric(response)
  } else if (is.factor(response) || is.character(response)) {
    response <- as.numeric(factor(response)) - 1
  }

  # Calculate ROC
  pos <- predictor[response == 1]
  neg <- predictor[response == 0]

  auc <- wilcox.test(pos, neg, exact = FALSE)$statistic / (length(pos) * length(neg))
  if (auc < 0.5) auc <- 1 - auc

  # ROC points
  thresholds <- sort(unique(predictor))
  tpr <- numeric(length(thresholds))
  fpr <- numeric(length(thresholds))

  for (i in seq_along(thresholds)) {
    pred_pos <- predictor >= thresholds[i]
    tp <- sum(pred_pos & response == 1)
    fp <- sum(pred_pos & response == 0)
    fn <- sum(!pred_pos & response == 1)
    tn <- sum(!pred_pos & response == 0)

    tpr[i] <- tp / (tp + fn)
    fpr[i] <- fp / (fp + tn)
  }

  roc_data <- data.frame(
    FPR = c(1, fpr, 0),
    TPR = c(1, tpr, 0)
  )

  # Find optimal cutoff (Youden's index)
  youden <- tpr - fpr
  opt_idx <- which.max(youden)
  opt_threshold <- thresholds[opt_idx]
  opt_sens <- tpr[opt_idx]
  opt_spec <- 1 - fpr[opt_idx]

  p <- ggplot2::ggplot(roc_data, ggplot2::aes(x = FPR, y = TPR)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                         color = "gray50") +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("AUC = %.3f | Optimal: Sens=%.2f, Spec=%.2f",
                         auc, opt_sens, opt_spec),
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)"
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_minimal()

  structure(list(
    auc = as.numeric(auc),
    optimal_threshold = opt_threshold,
    optimal_sensitivity = opt_sens,
    optimal_specificity = opt_spec,
    plot = p,
    roc_data = roc_data
  ), class = "omni_roc")
}

#' Calibration Plot
#'
#' Plot calibration curve for prediction model
#' @param observed Observed binary outcomes
#' @param predicted Predicted probabilities
#' @param n_bins Number of bins (default 10)
#' @param title Plot title
#' @return ggplot object
#' @examples
#' \dontrun{
#' plot_calibration(observed, predicted)
#' }
#' @export
plot_calibration <- function(observed, predicted, n_bins = 10,
                             title = "Calibration Plot") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  breaks <- quantile(predicted, probs = seq(0, 1, length.out = n_bins + 1))
  breaks[1] <- 0
  breaks[length(breaks)] <- 1

  bins <- cut(predicted, breaks = breaks, include.lowest = TRUE)

  cal_data <- data.frame(
    observed = observed,
    predicted = predicted,
    bin = bins
  )

  cal_summary <- aggregate(cbind(observed, predicted) ~ bin, data = cal_data,
                           FUN = mean)
  cal_summary$n <- as.vector(table(bins))

  ggplot2::ggplot(cal_summary, ggplot2::aes(x = predicted, y = observed)) +
    ggplot2::geom_point(ggplot2::aes(size = n), color = "steelblue") +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                         color = "gray50") +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "red") +
    ggplot2::labs(title = title,
                  x = "Predicted Probability",
                  y = "Observed Proportion") +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_minimal()
}

#' Decision Curve Analysis
#'
#' Plot decision curve for clinical utility
#' @param observed Observed outcomes
#' @param predicted Predicted probabilities
#' @param title Plot title
#' @return ggplot object
#' @examples
#' \dontrun{
#' plot_decision_curve(observed, predicted)
#' }
#' @export
plot_decision_curve <- function(observed, predicted, title = "Decision Curve") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  # Threshold probabilities
  thresholds <- seq(0.01, 0.99, by = 0.01)

  net_benefit <- numeric(length(thresholds))
  treat_all <- numeric(length(thresholds))
  treat_none <- numeric(length(thresholds))

  prev <- mean(observed)

  for (i in seq_along(thresholds)) {
    pt <- thresholds[i]

    # Treat all
    treat_all[i] <- prev - (1 - prev) * pt / (1 - pt)

    # Treat none
    treat_none[i] <- 0

    # Model
    pred_pos <- predicted >= pt
    tp <- sum(pred_pos & observed == 1)
    fp <- sum(pred_pos & observed == 0)
    n <- length(observed)

    net_benefit[i] <- (tp / n) - (fp / n) * (pt / (1 - pt))
  }

  plot_data <- data.frame(
    threshold = thresholds,
    net_benefit = net_benefit,
    treat_all = treat_all,
    treat_none = treat_none
  )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = threshold)) +
    ggplot2::geom_line(ggplot2::aes(y = net_benefit, color = "Model"),
                       linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(y = treat_all, color = "Treat All"),
                       linetype = "dashed") +
    ggplot2::geom_line(ggplot2::aes(y = treat_none, color = "Treat None"),
                       linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("Model" = "steelblue",
                                           "Treat All" = "red",
                                           "Treat None" = "gray50")) +
    ggplot2::labs(title = title,
                  x = "Threshold Probability",
                  y = "Net Benefit",
                  color = "Strategy") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

# ============================================================================
# SPECIALIZED PLOTS
# ============================================================================

#' Volcano Plot
#'
#' Create volcano plot for differential analysis
#' @param effect_size Vector of effect sizes
#' @param p_values Vector of p-values
#' @param labels Optional labels
#' @param fc_threshold Fold change threshold
#' @param p_threshold P-value threshold
#' @param title Plot title
#' @return ggplot object
#' @examples
#' fc <- rnorm(100, 0, 2)
#' p <- runif(100, 0, 1)
#' plot_volcano(fc, p)
#' @export
plot_volcano <- function(effect_size, p_values, labels = NULL,
                         fc_threshold = 1, p_threshold = 0.05,
                         title = "Volcano Plot") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  neg_log_p <- -log10(p_values)

  sig <- ifelse(abs(effect_size) >= fc_threshold & p_values < p_threshold,
                "Significant",
                ifelse(abs(effect_size) >= fc_threshold,
                       "Fold Change",
                       ifelse(p_values < p_threshold,
                              "P-value",
                              "Not significant")))

  plot_data <- data.frame(
    effect_size = effect_size,
    neg_log_p = neg_log_p,
    significance = sig
  )

  if (!is.null(labels)) {
    plot_data$labels <- labels
  }

  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = effect_size, y = neg_log_p,
                                    color = significance)) +
    ggplot2::geom_point(alpha = 0.6, size = 2) +
    ggplot2::geom_vline(xintercept = c(-fc_threshold, fc_threshold),
                        linetype = "dashed", color = "gray50") +
    ggplot2::geom_hline(yintercept = -log10(p_threshold),
                        linetype = "dashed", color = "gray50") +
    ggplot2::scale_color_manual(values = c(
      "Significant" = "red",
      "Fold Change" = "orange",
      "P-value" = "blue",
      "Not significant" = "gray50"
    )) +
    ggplot2::labs(title = title,
                  x = "Effect Size",
                  y = "-log10(p-value)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  p
}

#' Bland-Altman Plot
#'
#' Create Bland-Altman plot for method comparison
#' @param method1 Measurements from method 1
#' @param method2 Measurements from method 2
#' @param title Plot title
#' @return ggplot object
#' @examples
#' \dontrun{
#' plot_bland_altman(method1, method2)
#' }
#' @export
plot_bland_altman <- function(method1, method2, title = "Bland-Altman Plot") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  mean_vals <- (method1 + method2) / 2
  diff_vals <- method1 - method2

  mean_diff <- mean(diff_vals, na.rm = TRUE)
  sd_diff <- sd(diff_vals, na.rm = TRUE)

  limits <- c(mean_diff - 1.96 * sd_diff, mean_diff + 1.96 * sd_diff)

  plot_data <- data.frame(
    mean = mean_vals,
    difference = diff_vals
  )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = mean, y = difference)) +
    ggplot2::geom_point(alpha = 0.6, color = "steelblue") +
    ggplot2::geom_hline(yintercept = mean_diff, color = "red", linewidth = 1) +
    ggplot2::geom_hline(yintercept = limits, color = "red",
                        linetype = "dashed") +
    ggplot2::labs(title = title,
                  x = "Mean of Two Methods",
                  y = "Difference (Method 1 - Method 2)") +
    ggplot2::theme_minimal()
}

#' Waterfall Plot
#'
#' Create waterfall plot for treatment response
#' @param response Change from baseline (percentage)
#' @param group Optional grouping variable
#' @param title Plot title
#' @return ggplot object
#' @examples
#' \dontrun{
#' plot_waterfall(response, group)
#' }
#' @export
plot_waterfall <- function(response, group = NULL,
                           title = "Waterfall Plot") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  plot_data <- data.frame(
    id = seq_along(response),
    response = response
  )

  if (!is.null(group)) {
    plot_data$group <- group
  }

  # Sort by response
  plot_data <- plot_data[order(plot_data$response), ]
  plot_data$id <- factor(plot_data$id, levels = plot_data$id)

  if (!is.null(group)) {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = id, y = response, fill = group))
  } else {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = id, y = response))
  }

  p + ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_hline(yintercept = 0, color = "black") +
    ggplot2::geom_hline(yintercept = c(-20, -30, 20),
                        linetype = "dashed", color = "gray50") +
    ggplot2::labs(title = title,
                  x = "Patient",
                  y = "Change from Baseline (%)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
}

#' Swimmer Plot
#'
#' Create swimmer plot for treatment duration
#' @param data Data frame with start, end, and event information
#' @param title Plot title
#' @return ggplot object
#' @examples
#' \dontrun{
#' plot_swimmer(data)
#' }
#' @export
plot_swimmer <- function(data, title = "Swimmer Plot") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  ggplot2::ggplot(data, ggplot2::aes(y = reorder(id, duration))) +
    ggplot2::geom_segment(ggplot2::aes(x = start, xend = end,
                                        yend = reorder(id, duration)),
                          linewidth = 3, color = "steelblue") +
    ggplot2::geom_point(ggplot2::aes(x = end, shape = event),
                        size = 3) +
    ggplot2::labs(title = title,
                  x = "Time",
                  y = "Patient") +
    ggplot2::theme_minimal()
}

#' Kaplan-Meier by Subgroup
#'
#' Plot KM curves with subgroup analysis
#' @param time Survival time
#' @param status Event status
#' @param treatment Treatment variable
#' @param subgroup Subgroup variable
#' @param title Plot title
#' @return ggplot object
#' @examples
#' \dontrun{
#' plot_km_subgroup(time, status, treatment, subgroup)
#' }
#' @export
plot_km_subgroup <- function(time, status, treatment, subgroup,
                              title = "Survival by Subgroup") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  subgroups <- unique(subgroup)
  plot_data <- data.frame()

  for (sg in subgroups) {
    sg_idx <- subgroup == sg
    treatments <- unique(treatment[sg_idx])

    for (trt in treatments) {
      trt_idx <- treatment == trt & sg_idx
      km <- calculate_km_for_plot(time[trt_idx], status[trt_idx])

      df <- data.frame(
        time = km$time,
        survival = km$survival,
        subgroup = sg,
        treatment = trt
      )
      plot_data <- rbind(plot_data, df)
    }
  }

  ggplot2::ggplot(plot_data,
                  ggplot2::aes(x = time, y = survival,
                               color = treatment, linetype = subgroup)) +
    ggplot2::geom_step(linewidth = 1) +
    ggplot2::labs(title = title,
                  x = "Time",
                  y = "Survival Probability") +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}
