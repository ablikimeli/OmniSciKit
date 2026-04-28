#' Calculate BMI
#'
#' Calculate Body Mass Index and classify
#' @param weight Weight in kg
#' @param height Height in meters
#' @return List with BMI value and classification
#' @examples
#' calculate_bmi(weight = 70, height = 1.75)
#' @export
calculate_bmi <- function(weight, height) {
  bmi <- weight / (height^2)

  if (bmi < 18.5) {
    category <- "Underweight"
  } else if (bmi < 25) {
    category <- "Normal"
  } else if (bmi < 30) {
    category <- "Overweight"
  } else {
    category <- "Obese"
  }

  structure(list(
    bmi = round(bmi, 2),
    category = category,
    weight = weight,
    height = height
  ), class = "omni_clinical")
}

#' Calculate eGFR
#'
#' Calculate estimated Glomerular Filtration Rate using CKD-EPI 2021 equation
#' @param creatinine Serum creatinine in mg/dL
#' @param age Age in years
#' @param sex "male" or "female"
#' @param race Optional, not used in 2021 equation
#' @return List with eGFR value and CKD stage
#' @examples
#' calculate_egfr(creatinine = 1.2, age = 50, sex = "male")
#' @export
calculate_egfr <- function(creatinine, age, sex, race = NULL) {
  sex <- tolower(sex)

  if (sex == "female") {
    kappa <- 0.7
    alpha <- -0.241
    female_factor <- 1.012
  } else {
    kappa <- 0.9
    alpha <- -0.302
    female_factor <- 1
  }

  min_term <- min(creatinine / kappa, 1)
  max_term <- max(creatinine / kappa, 1)

  egfr <- 142 *
    (min_term^alpha) *
    (max_term^-1.200) *
    (0.9938^age) *
    female_factor

  # CKD Stage
  if (egfr >= 90) {
    stage <- "G1"
    description <- "Normal or high"
  } else if (egfr >= 60) {
    stage <- "G2"
    description <- "Mildly decreased"
  } else if (egfr >= 45) {
    stage <- "G3a"
    description <- "Mildly to moderately decreased"
  } else if (egfr >= 30) {
    stage <- "G3b"
    description <- "Moderately to severely decreased"
  } else if (egfr >= 15) {
    stage <- "G4"
    description <- "Severely decreased"
  } else {
    stage <- "G5"
    description <- "Kidney failure"
  }

  structure(list(
    egfr = round(egfr, 2),
    stage = stage,
    description = description,
    creatinine = creatinine,
    age = age,
    sex = sex
  ), class = "omni_clinical")
}
