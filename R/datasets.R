#' Glucose Control Study Dataset
#'
#' @description
#' A dataset containing HbA1c and clinical parameters for 500 diabetic patients.
#' Covers t-test, ANOVA, and basic statistical analyses for diabetes research.
#'
#' @format A data frame with 500 rows and 11 variables:
#' \describe{
#'   \item{patient_id}{Patient identifier (P0001-P0500)}
#'   \item{age}{Age in years}
#'   \item{gender}{Gender (男/女)}
#'   \item{bmi}{Body mass index (kg/m^2)}
#'   \item{diabetes_type}{Diabetes type (1型/2型/妊娠糖尿病)}
#'   \item{duration}{Disease duration in years}
#'   \item{hba1c_baseline}{Baseline HbA1c level (\%)}
#'   \item{treatment}{Treatment regimen (胰岛素/口服药/联合治疗)}
#'   \item{complications}{Complications (无/微血管/大血管/两者)}
#'   \item{fasting_glucose}{Fasting glucose (mmol/L)}
#'   \item{hba1c_3m}{HbA1c after 3 months (\%)}
#' }
#'
#' @source Modified from \url{https://gitee.com/ablikim-ali/data_for_learn/blob/master/clinical_stat_datasets.xlsx}
#' @keywords datasets
"omni_glucose"

#' Cardiovascular Risk Factors Dataset
#'
#' @description
#' Cardiovascular risk factors for 500 patients including blood pressure,
#' lipid profile, and lifestyle factors. Supports correlation analysis
#' and linear regression modeling.
#'
#' @format A data frame with 500 rows and 8 variables:
#' \describe{
#'   \item{patient_id}{Patient identifier (CV0001-CV0500)}
#'   \item{sbp}{Systolic blood pressure (mmHg)}
#'   \item{dbp}{Diastolic blood pressure (mmHg)}
#'   \item{ldl}{LDL cholesterol (mg/dL)}
#'   \item{age}{Age in years}
#'   \item{smoking_packyears}{Smoking pack-year history}
#'   \item{exercise_hours}{Weekly exercise hours}
#'   \item{cv_risk}{10-year cardiovascular risk score}
#' }
#'
#' @source Modified from \url{https://gitee.com/ablikim-ali/data_for_learn/blob/master/clinical_stat_datasets.xlsx}
#' @keywords datasets
"omni_cardio"

#' Cancer Treatment Outcomes Dataset
#'
#' @description
#' Cancer treatment outcomes for 300 patients including survival data.
#' Supports logistic regression and survival analysis (Kaplan-Meier, Cox PH).
#'
#' @format A data frame with 300 rows and 9 variables:
#' \describe{
#'   \item{patient_id}{Patient identifier (CA0001-CA0300)}
#'   \item{age}{Age in years}
#'   \item{cancer_stage}{Cancer stage (I期/II期/III期/IV期)}
#'   \item{treatment}{Treatment modality (手术/化疗/放疗/免疫)}
#'   \item{ecog_score}{ECOG performance status (0-3)}
#'   \item{biomarker_level}{Biomarker level}
#'   \item{response}{Treatment response (无响应/有响应)}
#'   \item{time_to_event}{Survival time in months}
#'   \item{status}{Event status (1=event, 0=censored)}
#' }
#'
#' @source Modified from \url{https://gitee.com/ablikim-ali/data_for_learn/blob/master/clinical_stat_datasets.xlsx}
#' @keywords datasets
"omni_cancer"

#' Depression Repeated Measures Dataset
#'
#' @description
#' Longitudinal depression scores for 200 patients measured at 5 time points
#' (1000 observations). Supports mixed-effects models and repeated measures ANOVA.
#'
#' @format A data frame with 1000 rows and 7 variables:
#' \describe{
#'   \item{patient_id}{Patient identifier (RM0001-RM0200)}
#'   \item{center}{Study center (中心1-中心4)}
#'   \item{time}{Time point label (第0周-第4周)}
#'   \item{time_num}{Numeric time (0-4)}
#'   \item{age}{Age in years}
#'   \item{depression_score}{Depression score (0-30)}
#'   \item{medication_adherence}{Medication adherence (0-100)}
#' }
#'
#' @source Modified from \url{https://gitee.com/ablikim-ali/data_for_learn/blob/master/clinical_stat_datasets.xlsx}
#' @keywords datasets
"omni_depression"

#' Symptom Cluster Analysis Dataset
#'
#' @description
#' Symptom scores for 600 patients across 5 symptom domains with true cluster labels.
#' Supports K-means clustering, PCA, and dimensionality reduction analyses.
#'
#' @format A data frame with 600 rows and 9 variables:
#' \describe{
#'   \item{patient_id}{Patient identifier (CL0001-CL0600)}
#'   \item{age}{Age in years}
#'   \item{gender}{Gender (男/女)}
#'   \item{fatigue}{Fatigue score (0-100)}
#'   \item{joint_pain}{Joint pain score (0-100)}
#'   \item{headache}{Headache severity (0-100)}
#'   \item{nausea}{Nausea severity (0-100)}
#'   \item{sleep_disturbance}{Sleep disturbance score (0-100)}
#'   \item{true_cluster}{True cluster label for validation (症状簇1-症状簇4)}
#' }
#'
#' @source Modified from \url{https://gitee.com/ablikim-ali/data_for_learn/blob/master/clinical_stat_datasets.xlsx}
#' @keywords datasets
"omni_symptoms"

#' Diagnostic Test Evaluation Dataset
#'
#' @description
#' Biomarker and disease status data for 450 patients with 3 diagnostic tests
#' at varying accuracy levels. Supports ROC analysis, sensitivity/specificity
#' calculation, and diagnostic test evaluation.
#'
#' @format A data frame with 450 rows and 8 variables:
#' \describe{
#'   \item{patient_id}{Subject identifier (DX0001-DX0450)}
#'   \item{true_disease}{True disease status (金标准: 有病/无病)}
#'   \item{test1_result}{Test 1 continuous result (AUC ~0.9, excellent)}
#'   \item{test2_result}{Test 2 continuous result (AUC ~0.75, moderate)}
#'   \item{test3_result}{Test 3 continuous result (AUC ~0.65, poor)}
#'   \item{test1_positive}{Test 1 binary result (阳性/阴性)}
#'   \item{test2_positive}{Test 2 binary result (阳性/阴性)}
#'   \item{test3_positive}{Test 3 binary result (阳性/阴性)}
#' }
#'
#' @source Modified from \url{https://gitee.com/ablikim-ali/data_for_learn/blob/master/clinical_stat_datasets.xlsx}
#' @keywords datasets
"omni_diagnostic"

#' Machine Learning Features Dataset
#'
#' @description
#' Clinical features for 350 patients with 20 feature variables and binary outcome.
#' Supports random forest, SVM, and other ML classification algorithms.
#'
#' @format A data frame with 350 rows and 24 variables:
#' \describe{
#'   \item{sample_id}{Sample identifier (ML0001-ML0350)}
#'   \item{age}{Age in years}
#'   \item{gender}{Gender (男/女)}
#'   \item{class}{Outcome class (稳定/进展)}
#'   \item{feature_1}{Feature variable 1}
#'   \item{feature_2}{Feature variable 2}
#'   \item{feature_3}{Feature variable 3}
#'   \item{feature_4}{Feature variable 4}
#'   \item{feature_5}{Feature variable 5}
#'   \item{feature_6}{Feature variable 6}
#'   \item{feature_7}{Feature variable 7}
#'   \item{feature_8}{Feature variable 8}
#'   \item{feature_9}{Feature variable 9}
#'   \item{feature_10}{Feature variable 10}
#'   \item{feature_11}{Feature variable 11}
#'   \item{feature_12}{Feature variable 12}
#'   \item{feature_13}{Feature variable 13}
#'   \item{feature_14}{Feature variable 14}
#'   \item{feature_15}{Feature variable 15}
#'   \item{feature_16}{Feature variable 16}
#'   \item{feature_17}{Feature variable 17}
#'   \item{feature_18}{Feature variable 18}
#'   \item{feature_19}{Feature variable 19}
#'   \item{feature_20}{Feature variable 20}
#' }
#'
#' @source Modified from \url{https://gitee.com/ablikim-ali/data_for_learn/blob/master/clinical_stat_datasets.xlsx}
#' @keywords datasets
"omni_ml"

#' Dose-Response Dataset
#'
#' @description
#' Dose-response relationship data for 250 subjects across 8 dose levels.
#' Supports non-linear regression, Emax model fitting, and dose-response analysis.
#'
#' @format A data frame with 250 rows and 3 variables:
#' \describe{
#'   \item{subject_id}{Subject identifier (DR0001-DR0250)}
#'   \item{dose_mg}{Drug dose in mg (0, 5, 10, 20, 40, 80, 160, 320)}
#'   \item{response}{Physiological response (continuous variable)}
#' }
#'
#' @source Modified from \url{https://gitee.com/ablikim-ali/data_for_learn/blob/master/clinical_stat_datasets.xlsx}
#' @keywords datasets
"omni_dose"
