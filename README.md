# OmniSciKit: Comprehensive Clinical Statistical Toolkit

> 一个全面的临床研究和生物统计R包
# 该包由“Abel医研统计”公众号团队开发，如遇到bug或者想提供更新建议，随时欢迎联系我们。
[![License: GPL-3](https://img.shields.io/badge/License-GPL3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R >= 4.0.0](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)](https://www.r-project.org/)

## 📦 安装

```r
# 推荐：使用pak安装（自动处理依赖）
if (!requireNamespace("pak", quietly = TRUE))
    install.packages("pak")
pak::pak("ablikimeli/OmniSciKit")

# 或使用devtools
# if (!requireNamespace("devtools", quietly = TRUE))
#     install.packages("devtools")
# devtools::install_github("ablikimeli/OmniSciKit")
```

## 📋 主要功能模块

### 1. 样本量计算 (Sample Size)

涵盖临床试验主要设计类型的样本量计算：

| 函数 | 说明 |
|------|------|
| `ss_ttest_one_sample()` | 单样本t检验 |
| `ss_ttest_paired()` | 配对t检验 |
| `ss_ttest_welch()` | Welch t检验（方差不齐） |
| `ss_mann_whitney()` | Mann-Whitney U检验 |
| `ss_mcnemar()` | McNemar检验（配对比例） |
| `ss_noninferiority_mean()` | 非劣效性试验（均数） |
| `ss_noninferiority_prop()` | 非劣效性试验（率） |
| `ss_noninferiority_survival()` | 非劣效性试验（生存） |
| `ss_equivalence_mean()` | 等效性试验（均数） |
| `ss_equivalence_prop()` | 等效性试验（率） |
| `ss_equivalence_ratio()` | 等效性试验（均数比，生物等效性） |
| `ss_superiority_mean()` | 优效性试验（均数） |
| `ss_superiority_prop()` | 优效性试验（率） |
| `ss_survival_logrank()` | 生存分析（log-rank检验） |
| `ss_cox_regression()` | Cox回归 |
| `ss_linear_regression()` | 线性回归 |
| `ss_logistic_regression()` | Logistic回归 |
| `ss_poisson_regression()` | Poisson回归 |
| `ss_crossover()` | 交叉设计 |
| `ss_group_sequential()` | 成组序贯设计 |
| `ss_repeated_measures()` | 重复测量设计 |
| `ss_cluster_randomized()` | 整群随机试验 |
| `ss_multiple_comparisons()` | 多重比较校正 |
| `ss_prediction_binary()` | 临床预测模型（二分类） |
| `ss_prediction_continuous()` | 临床预测模型（连续） |
| `ss_prediction_survival()` | 临床预测模型（生存） |
| `ss_diagnostic_sens()` | 诊断试验（灵敏度） |
| `ss_diagnostic_spec()` | 诊断试验（特异度） |
| `ss_diagnostic_auc()` | 诊断试验（AUC） |
| `ss_precision()` | 精度估计（均数） |
| `ss_precision_prop()` | 精度估计（率） |

### 2. 功效分析 (Power Analysis)

| 函数 | 说明 |
|------|------|
| `power_ttest()` | t检验功效 |
| `power_anova()` | ANOVA功效 |
| `power_correlation()` | 相关分析功效 |
| `power_proportion()` | 率比较功效 |
| `power_survival()` | 生存分析功效 |
| `power_cox_regression()` | Cox回归功效 |
| `power_logistic_regression()` | Logistic回归功效 |
| `power_equivalence_mean()` | 等效性试验功效 |
| `power_noninferiority_mean()` | 非劣效性试验功效 |
| `power_cluster_randomized()` | 整群随机试验功效 |
| `plot_power_curve()` | 功效曲线可视化 |
| `plot_effect_size_curve()` | 效应量与功效曲线 |

### 3. 统计分析 (Statistics)

| 函数 | 说明 |
|------|------|
| `stats_ttest()` | 两样本t检验 |
| `stats_wilcoxon()` | Wilcoxon秩和检验 |
| `stats_anova()` | 单因素方差分析 |
| `stats_correlation()` | Pearson/Spearman相关 |
| `stats_regression()` | 线性回归 |
| `stats_logistic_regression()` | Logistic回归 |
| `stats_poisson_regression()` | Poisson回归 |
| `stats_cox_regression()` | Cox比例风险回归 |
| `stats_chisq()` | 卡方检验 |
| `stats_fisher_exact()` | Fisher精确检验 |
| `stats_mcnemar()` | McNemar检验 |
| `stats_normality()` | 正态性检验 |
| `stats_km_survival()` | Kaplan-Meier生存分析 |
| `stats_logrank_test()` | Log-rank检验 |
| `stats_kappa()` | Cohen's Kappa一致性 |
| `stats_icc()` | 组内相关系数 |

### 4. 效应量 (Effect Size)

| 函数 | 说明 |
|------|------|
| `cohen_d()` | Cohen's d |
| `cohen_f()` | Cohen's f (ANOVA) |
| `cohen_f2()` | Cohen's f² (回归) |
| `odds_ratio()` | 比值比OR |
| `relative_risk()` | 相对危险度RR |
| `hazard_ratio()` | 风险比HR |
| `nnt()` | 需治疗人数NNT |

### 5. 可视化 (Visualization)

| 函数 | 说明 |
|------|------|
| `plot_km_curve()` | Kaplan-Meier生存曲线 |
| `plot_km_subgroup()` | 亚组KM曲线 |
| `plot_cumulative_incidence()` | 累积发生率 |
| `plot_forest_meta()` | 森林图（Meta分析） |
| `plot_subgroup_forest()` | 亚组分析森林图 |
| `plot_roc_enhanced()` | ROC曲线（含最佳截断值） |
| `plot_calibration()` | 校准曲线 |
| `plot_decision_curve()` | 决策曲线DCA |
| `plot_volcano()` | 火山图 |
| `plot_bland_altman()` | Bland-Altman图 |
| `plot_waterfall()` | 瀑布图 |
| `plot_swimmer()` | 游泳图 |
| `plot_boxplot()` | 箱线图 |
| `plot_scatter()` | 散点图（含回归线） |
| `plot_correlation_matrix()` | 相关矩阵热图 |
| `plot_power_curve()` | 功效曲线 |
| `plot_effect_size_curve()` | 效应量-功效曲线 |

### 6. 诊断试验 (Diagnostic Tests)

| 函数 | 说明 |
|------|------|
| `diag_test_performance()` | 灵敏度、特异度、PPV、NPV |
| `plot_roc_enhanced()` | ROC曲线 |
| `plot_calibration()` | 校准曲线 |
| `plot_decision_curve()` | 决策曲线 |

### 7. 临床工具 (Clinical Tools)

| 函数 | 说明 |
|------|------|
| `calculate_bmi()` | BMI计算与分类 |
| `calculate_egfr()` | eGFR估算（CKD-EPI 2021） |

### 8. Table One 基线特征表

| 函数 | 说明 |
|------|------|
| `table_one()` | 基线特征表，支持分组和p值 |

### 9. 工具函数 (Utilities)

| 函数 | 说明 |
|------|------|
| `impute_missing()` | 缺失值填补 |
| `multiple_imputation()` | 多重插补 |
| `winsorize()` | 极端值处理 |
| `log_transform()` | 对数变换 |
| `standardize_vars()` | 变量标准化 |
| `check_assumptions()` | 统计假设检验 |
| `summary_stats()` | 综合描述统计 |
| `crosstab()` | 列联表分析 |
| `multiple_testing_correction()` | 多重比较校正 |
| `fwer_control()` | FWER控制 |
| `bootstrap_ci()` | Bootstrap置信区间 |
| `bootstrap_regression()` | Bootstrap回归 |
| `adjust_dropout()` | 脱落率校正 |
| `adjust_alpha()` | alpha水平校正 |
| `create_dummy_vars()` | 创建哑变量 |
| `convert_wide_to_long()` | 宽转长数据 |
| `convert_long_to_wide()` | 长转宽数据 |
| `merge_datasets()` | 数据合并 |

## 📊 内置数据集

| 数据集 | 说明 | 样本量 |
|--------|------|--------|
| `omni_glucose` | 血糖控制研究 | 500 |
| `omni_cardio` | 心血管风险因素 | 400 |
| `omni_cancer` | 癌症治疗结果 | 300 |
| `omni_depression` | 抑郁重复测量 | 200 |
| `omni_symptoms` | 症状群数据 | 600 |
| `omni_diagnostic` | 诊断试验数据 | 450 |
| `omni_ml` | 机器学习特征 | 350 |
| `omni_dose` | 剂量-反应数据 | 250 |

## 🚀 快速示例

```r
library(OmniSciKit)

# 加载数据
data(omni_glucose)

# 样本量计算：两样本t检验
ss_ttest_one_sample(mu0 = 100, mu1 = 105, sd = 15, power = 0.8)

# 非劣效性试验样本量
ss_noninferiority_mean(margin = 3, sd = 8, power = 0.8)

# 生存分析样本量
ss_survival_logrank(hr = 0.7, median_control = 12, accrual = 12, follow_up = 18)

# 功效曲线
plot_power_curve(n_range = seq(10, 200, by = 10), effect_size = 0.5)

# 统计分析
stats_ttest(omni_glucose$hba1c_3m[omni_glucose$treatment == "胰岛素"],
            omni_glucose$hba1c_3m[omni_glucose$treatment == "口服药"])

# Logistic回归
omni_glucose$poor_control <- ifelse(omni_glucose$hba1c_3m > 7, 1, 0)
stats_logistic_regression(poor_control ~ age + bmi, data = omni_glucose)

# Kaplan-Meier生存曲线
data(omni_cancer)
plot_km_curve(omni_cancer$time_to_event, omni_cancer$status, omni_cancer$treatment)

# ROC曲线
data(omni_diagnostic)
disease_binary <- ifelse(omni_diagnostic$true_disease == "阳性", 1, 0)
plot_roc_enhanced(omni_diagnostic$test1_result, disease_binary)

# Table One
table_one(omni_glucose, vars = c("age", "bmi", "hba1c_3m", "gender"),
          group = "treatment")

# 效应量
cohen_d(omni_glucose$hba1c_3m[omni_glucose$treatment == "胰岛素"],
        omni_glucose$hba1c_3m[omni_glucose$treatment == "口服药"])
```

## 📖 参考文献

- Chow, S.C., Shao, J., Wang, H., & Lokhnygina, Y. (2017). *Sample Size Calculations in Clinical Research* (3rd ed.). Chapman and Hall/CRC.
- Cohen, J. (1988). *Statistical Power Analysis for the Behavioral Sciences* (2nd ed.). Lawrence Erlbaum Associates.
- Riley, R.D., et al. (2020). Calculating the sample size required for developing a clinical prediction model. *BMJ*, 368, m441.
- Hsieh, F.Y., Bloch, D.A., & Larsen, M.D. (1998). A simple method of sample size calculation for linear and logistic regression. *Statistics in Medicine*, 17(14), 1623-1634.
- Farrington, C.P., & Manning, G. (1990). Test statistics and sample size formulae for comparative binomial trials with null hypothesis of non-zero risk difference or non-unity relative risk. *Statistics in Medicine*, 9(12), 1447-1454.
- Schoenfeld, D.A. (1983). Sample-size formula for the proportional-hazards regression model. *Biometrics*, 39(2), 499-503.

## 📄 License

GPL-3 © Ablikim Ali

## 📧 联系方式

Maintainer: Ablikim Ali <medical_stat@163.com>
