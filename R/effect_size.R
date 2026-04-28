#' Cohen's d Effect Size
#'
#' Calculate Cohen's d for two independent groups
#' @param x Numeric vector for group 1
#' @param y Numeric vector for group 2
#' @return List with Cohen's d and interpretation
#' @examples
#' data(omni_glucose)
#' cohen_d(omni_glucose$hba1c_baseline[1:100], omni_glucose$hba1c_baseline[101:200])
#' @export
cohen_d <- function(x, y) {
  x <- na.omit(x)
  y <- na.omit(y)

  mean_x <- mean(x)
  mean_y <- mean(y)
  sd_x <- sd(x)
  sd_y <- sd(y)

  n_x <- length(x)
  n_y <- length(y)

  if (n_x < 2 || n_y < 2 || sd_x == 0 || sd_y == 0) {
    stop("Insufficient data or zero variance in one or both groups")
  }

  # Pooled standard deviation
  sd_pooled <- sqrt(((n_x - 1) * sd_x^2 + (n_y - 1) * sd_y^2) / (n_x + n_y - 2))

  d <- (mean_x - mean_y) / sd_pooled

  # Interpretation
  abs_d <- abs(d)
  if (is.na(abs_d)) {
    interpretation <- "Unknown"
  } else if (abs_d < 0.2) {
    interpretation <- "Negligible"
  } else if (abs_d < 0.5) {
    interpretation <- "Small"
  } else if (abs_d < 0.8) {
    interpretation <- "Medium"
  } else {
    interpretation <- "Large"
  }

  structure(list(
    d = d,
    abs_d = abs_d,
    interpretation = interpretation,
    sd_pooled = sd_pooled,
    n1 = n_x,
    n2 = n_y
  ), class = "omni_effect_size")
}

#' Cohen's f for ANOVA
#'
#' Calculate Cohen's f effect size for ANOVA
#' @param ss_between Between-group sum of squares
#' @param ss_within Within-group sum of squares
#' @param k Number of groups
#' @return List with Cohen's f and interpretation
#' @examples
#' cohen_f(ss_between = 100, ss_within = 400, k = 3)
#' @export
cohen_f <- function(ss_between, ss_within, k = NULL) {
  # Use eta-squared to avoid needing k/n
  eta_sq <- ss_between / (ss_between + ss_within)
  f <- sqrt(eta_sq / (1 - eta_sq))

  # Interpretation
  if (f < 0.1) {
    interpretation <- "Small"
  } else if (f < 0.25) {
    interpretation <- "Medium"
  } else {
    interpretation <- "Large"
  }

  structure(list(
    f = f,
    interpretation = interpretation,
    ss_between = ss_between,
    ss_within = ss_within,
    k = k
  ), class = "omni_effect_size")
}

#' Cohen's f-squared for Regression
#'
#' Calculate Cohen's f-squared for multiple regression
#' @param r_squared R-squared value
#' @param r_squared_base Baseline R-squared (default 0)
#' @return List with Cohen's f-squared and interpretation
#' @examples
#' cohen_f2(r_squared = 0.25)
#' @export
cohen_f2 <- function(r_squared, r_squared_base = 0) {
  f2 <- (r_squared - r_squared_base) / (1 - r_squared)

  # Interpretation
  if (f2 < 0.02) {
    interpretation <- "Small"
  } else if (f2 < 0.15) {
    interpretation <- "Medium"
  } else {
    interpretation <- "Large"
  }

  structure(list(
    f2 = f2,
    interpretation = interpretation,
    r_squared = r_squared,
    r_squared_base = r_squared_base
  ), class = "omni_effect_size")
}
