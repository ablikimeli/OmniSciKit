#' OmniSciKit: Comprehensive Clinical Statistical Toolkit
#'
#' A comprehensive R package for clinical research and biostatistics.
#' Integrates sample size calculation, statistical analysis, data visualization,
#' power analysis, effect size calculation, and clinical tools.
#'
#' @section Sample Size Functions:
#' \itemize{
#'   \item \code{ss_ttest_one_sample} - One-sample t-test sample size
#'   \item \code{ss_ttest_paired} - Paired t-test sample size
#'   \item \code{ss_ttest_welch} - Welch's t-test sample size (unequal variances)
#'   \item \code{ss_mann_whitney} - Mann-Whitney U test sample size
#'   \item \code{ss_mcnemar} - McNemar's test sample size
#'   \item \code{ss_noninferiority_mean} - Non-inferiority trial (means)
#'   \item \code{ss_noninferiority_prop} - Non-inferiority trial (proportions)
#'   \item \code{ss_noninferiority_survival} - Non-inferiority trial (survival)
#'   \item \code{ss_equivalence_mean} - Equivalence trial (means)
#'   \item \code{ss_equivalence_prop} - Equivalence trial (proportions)
#'   \item \code{ss_equivalence_ratio} - Equivalence trial (ratio of means)
#'   \item \code{ss_superiority_mean} - Superiority trial (means)
#'   \item \code{ss_superiority_prop} - Superiority trial (proportions)
#'   \item \code{ss_survival_logrank} - Log-rank test sample size
#'   \item \code{ss_cox_regression} - Cox regression sample size
#'   \item \code{ss_linear_regression} - Linear regression sample size
#'   \item \code{ss_logistic_regression} - Logistic regression sample size
#'   \item \code{ss_poisson_regression} - Poisson regression sample size
#'   \item \code{ss_crossover} - Crossover design sample size
#'   \item \code{ss_group_sequential} - Group sequential design
#'   \item \code{ss_repeated_measures} - Repeated measures ANOVA
#'   \item \code{ss_cluster_randomized} - Cluster randomized trial
#'   \item \code{ss_multiple_comparisons} - Multiple comparisons (Bonferroni)
#'   \item \code{ss_prediction_binary} - Binary prediction model
#'   \item \code{ss_prediction_continuous} - Continuous prediction model
#'   \item \code{ss_prediction_survival} - Survival prediction model
#'   \item \code{ss_diagnostic_sens} - Diagnostic test (sensitivity)
#'   \item \code{ss_diagnostic_spec} - Diagnostic test (specificity)
#'   \item \code{ss_diagnostic_auc} - Diagnostic test (AUC)
#'   \item \code{ss_precision} - Precision-based (CI width)
#'   \item \code{ss_precision_prop} - Precision-based (proportion)
#' }
#'
#' @section Power Analysis:
#' \itemize{
#'   \item \code{power_ttest} - Power for t-test
#'   \item \code{power_anova} - Power for ANOVA
#'   \item \code{power_correlation} - Power for correlation
#'   \item \code{power_proportion} - Power for two proportions
#'   \item \code{power_survival} - Power for survival analysis
#'   \item \code{power_cox_regression} - Power for Cox regression
#'   \item \code{power_logistic_regression} - Power for logistic regression
#'   \item \code{power_equivalence_mean} - Power for equivalence trial
#'   \item \code{power_noninferiority_mean} - Power for non-inferiority trial
#'   \item \code{power_cluster_randomized} - Power for cluster randomized trial
#'   \item \code{plot_power_curve} - Power curve visualization
#'   \item \code{plot_effect_size_curve} - Effect size vs power curve
#' }
#'
#' @section Statistical Functions:
#' \itemize{
#'   \item \code{stats_ttest} - Two-sample t-test
#'   \item \code{stats_wilcoxon} - Wilcoxon rank-sum test
#'   \item \code{stats_anova} - One-way ANOVA
#'   \item \code{stats_correlation} - Pearson/Spearman correlation
#'   \item \code{stats_regression} - Linear regression
#'   \item \code{stats_chisq} - Chi-square test
#'   \item \code{stats_normality} - Normality test
#'   \item \code{stats_cox_regression} - Cox regression
#'   \item \code{stats_logistic_regression} - Logistic regression
#'   \item \code{stats_poisson_regression} - Poisson regression
#'   \item \code{stats_km_survival} - Kaplan-Meier survival
#'   \item \code{stats_logrank_test} - Log-rank test
#'   \item \code{stats_fisher_exact} - Fisher's exact test
#'   \item \code{stats_mcnemar} - McNemar's test
#'   \item \code{stats_kappa} - Cohen's kappa
#'   \item \code{stats_icc} - Intraclass correlation
#'   \item \code{summary_stats} - Descriptive statistics
#'   \item \code{crosstab} - Cross-tabulation with statistics
#'   \item \code{table_one} - Table 1 baseline characteristics
#' }
#'
#' @section Effect Size:
#' \itemize{
#'   \item \code{cohen_d} - Cohen's d
#'   \item \code{cohen_f} - Cohen's f for ANOVA
#'   \item \code{cohen_f2} - Cohen's f-squared for regression
#'   \item \code{odds_ratio} - Odds ratio
#'   \item \code{relative_risk} - Relative risk
#'   \item \code{hazard_ratio} - Hazard ratio
#'   \item \code{nnt} - Number Needed to Treat
#' }
#'
#' @section Diagnostic Tests:
#' \itemize{
#'   \item \code{diag_test_performance} - Sensitivity, specificity, PPV, NPV
#'   \item \code{plot_roc_enhanced} - ROC curve with AUC
#'   \item \code{plot_calibration} - Calibration plot
#'   \item \code{plot_decision_curve} - Decision curve analysis
#' }
#'
#' @section Visualization:
#' \itemize{
#'   \item \code{plot_km_curve} - Kaplan-Meier survival curve
#'   \item \code{plot_km_subgroup} - KM curve by subgroup
#'   \item \code{plot_cumulative_incidence} - Cumulative incidence
#'   \item \code{plot_forest_meta} - Forest plot for meta-analysis
#'   \item \code{plot_subgroup_forest} - Subgroup analysis forest plot
#'   \item \code{plot_volcano} - Volcano plot
#'   \item \code{plot_bland_altman} - Bland-Altman plot
#'   \item \code{plot_waterfall} - Waterfall plot
#'   \item \code{plot_swimmer} - Swimmer plot
#'   \item \code{plot_boxplot} - Box plot by group
#'   \item \code{plot_scatter} - Scatter plot with regression line
#'   \item \code{plot_correlation_matrix} - Correlation heatmap
#' }
#'
#' @section Clinical Tools:
#' \itemize{
#'   \item \code{calculate_bmi} - BMI calculation and classification
#'   \item \code{calculate_egfr} - eGFR (CKD-EPI 2021 equation)
#' }
#'
#' @section Data Conversion:
#' \itemize{
#'   \item \code{convert_wide_to_long} - Wide to long format
#'   \item \code{convert_long_to_wide} - Long to wide format
#'   \item \code{standardize_vars} - Variable standardization
#'   \item \code{merge_datasets} - Merge datasets
#'   \item \code{import_csv} - Import CSV
#'   \item \code{export_csv} - Export CSV
#'   \item \code{create_dummy_vars} - Create dummy/indicator variables
#' }
#'
#' @section Utilities:
#' \itemize{
#'   \item \code{impute_missing} - Impute missing values
#'   \item \code{multiple_imputation} - Simple multiple imputation
#'   \item \code{winsorize} - Winsorize extreme values
#'   \item \code{log_transform} - Log transformation
#'   \item \code{summary_stats} - Summary statistics
#'   \item \code{check_assumptions} - Assumption checking
#'   \item \code{multiple_testing_correction} - Multiple testing correction
#'   \item \code{fwer_control} - Family-wise error rate control
#'   \item \code{bootstrap_ci} - Bootstrap confidence interval
#'   \item \code{bootstrap_regression} - Bootstrap regression
#'   \item \code{adjust_dropout} - Adjust for dropouts
#'   \item \code{adjust_alpha} - Adjust alpha for multiple testing
#' }
#'
#' @section Built-in Datasets:
#' \itemize{
#'   \item \code{omni_glucose} - Glucose control study (500 patients)
#'   \item \code{omni_cardio} - Cardiovascular risk (400 patients)
#'   \item \code{omni_cancer} - Cancer treatment outcomes (300 patients)
#'   \item \code{omni_depression} - Depression repeated measures (200 patients)
#'   \item \code{omni_symptoms} - Symptoms cluster (600 patients)
#'   \item \code{omni_diagnostic} - Diagnostic test (450 patients)
#'   \item \code{omni_ml} - Machine learning features (350 patients)
#'   \item \code{omni_dose} - Dose-response data (250 subjects)
#' }
#'
#' @references
#' Chow, S.C., Shao, J., Wang, H., & Lokhnygina, Y. (2017). Sample Size Calculations
#' in Clinical Research (3rd ed.). Chapman and Hall/CRC.
#'
#' Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences (2nd ed.).
#' Lawrence Erlbaum Associates.
#'
#' Riley, R.D., et al. (2020). Calculating the sample size required for developing a
#' clinical prediction model. BMJ, 368, m441.
#'
#' Hsieh, F.Y., Bloch, D.A., & Larsen, M.D. (1998). A simple method of sample size
#' calculation for linear and logistic regression. Statistics in Medicine, 17(14),
#' 1623-1634.
#'
#' @keywords internal
#' @name OmniSciKit
#' @aliases OmniSciKit-package
NULL
