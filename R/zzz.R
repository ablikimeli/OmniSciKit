.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "========================================\n",
    "OmniSciKit v0.2.0 - Clinical Statistical Toolkit\n",
    "========================================\n",
    "Sample Size:  ss_ttest_one_sample(), ss_survival_logrank(), ss_noninferiority_mean()\n",
    "Effect Size:  cohen_d(), cohen_f(), cohen_f2(), odds_ratio(), relative_risk()\n",
    "Power:        power_ttest(), power_anova(), plot_power_curve()\n",
    "Statistics:   stats_ttest(), stats_anova(), stats_logistic_regression()\n",
    "Survival:     stats_km_survival(), stats_logrank_test(), plot_km_curve()\n",
    "Visualization: plot_boxplot(), plot_scatter(), plot_volcano(), plot_roc_enhanced()\n",
    "Table 1:      table_one()\n",
    "Clinical:     calculate_bmi(), calculate_egfr()\n",
    "Diagnostic:   diag_test_performance(), plot_roc_enhanced()\n",
    "Datasets:     omni_glucose, omni_cardio, omni_cancer, omni_depression,\n",
    "              omni_symptoms, omni_diagnostic, omni_ml, omni_dose\n",
    "========================================"
  )
}
