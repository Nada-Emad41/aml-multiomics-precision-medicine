# Results

## 1. Cohort and Integrated Dataset

After harmonizing patient identifiers across the clinical and molecular datasets, the predictive modeling phase included **162 matched TCGA-LAML patients**.

The machine-learning feature matrix contained **600 molecular features**:

- 500 highly variable mRNA features
- 100 highly variable miRNA features

These molecular features were integrated with clinical outcome information for survival-related classification and downstream interpretation.

---

## 2. Exploratory Multi-Omics Analysis

Exploratory analysis was used to investigate molecular variation across AML patients before predictive modeling.

Principal Component Analysis (PCA) of the integrated molecular data demonstrated substantial heterogeneity across patients rather than complete separation between prognosis groups.

The exploratory phase also included individual analysis of:

- mRNA expression
- miRNA expression
- DNA methylation

Earlier multi-omics exploration used MOFA to investigate shared sources of variation across molecular layers. This exploratory integration was separate from the later machine-learning workflow.

![Integrated Multi-Omics Feature Heatmap](/results/figures/integrated_multiomics_feature_heatmap.png)

![Integrated Multi-Omics PCA](/results/figures/pca_integrated_multiomics.png)

---

## 3. Machine-Learning Performance

Three supervised machine-learning algorithms were evaluated for survival-related classification:

- Random Forest
- Support Vector Machine (SVM)
- XGBoost

On the initial held-out test split, **XGBoost produced the highest ROC-AUC among the three evaluated models**.

| Model | Test ROC-AUC |
|---|---:|
| Random Forest | 0.579 |
| Support Vector Machine | 0.550 |
| XGBoost | **0.726** |

The initial XGBoost evaluation also achieved approximately:

- **Accuracy:** 70.7%
- **ROC-AUC:** 0.726
- **Sensitivity:** 65.0%
- **Specificity:** 76.2%

![ROC Curves — Model Comparison](/results/figures/roc_curves_model_comparison.png)

---

## 4. Cross-Validation Assessment

Because performance from a single train/test split can depend on the particular patients assigned to each subset, model performance was also evaluated using **5-fold cross-validation**.

Cross-validation produced more conservative performance estimates than the initial held-out split.

This suggests that the initial single-split XGBoost performance should not be interpreted as an externally validated estimate of clinical predictive performance.

The cross-validation analysis therefore serves as an important robustness check and highlights the uncertainty associated with predictive modeling in a relatively small patient cohort.

![Cross-Validation Model Evaluation](/results/figures/model_evaluation_cross_validation.png)

---

## 5. Candidate Biomarker Discovery

Feature-importance analysis from the XGBoost workflow was used to identify **20 candidate molecular biomarkers** associated with the prediction task.

Highly ranked mRNA features included:

- FAM127A
- NDST3
- CA1
- FBLN1
- HOXB9
- MMP9
- PTPRD

The analysis also identified miRNA features, including **hsa-miR-1266**, among the influential predictors.

These features should be interpreted as **candidate biomarkers generated from this dataset and modeling workflow**, rather than clinically validated biomarkers.

![Top 20 Candidate Biomarkers](/results/figures/top20_candidate_biomarkers.png)

The complete ranked candidate biomarker table is available here:

[View candidate biomarker results](/results/top20_candidate_biomarkers.csv)

---

## 6. Survival Analysis

Predicted prognosis groups from the XGBoost workflow were evaluated using Kaplan-Meier survival analysis.

The predicted groups showed a statistically significant difference in observed survival distributions, with a reported **log-rank p-value of 0.0044**.

This result indicates that the model-derived prognosis groups captured survival-related structure within the analyzed TCGA-LAML cohort.

However, this remains an exploratory finding and requires validation in an independent patient cohort before any clinical interpretation or application.

![Kaplan-Meier Survival Analysis](/results/figures/kaplan_meier_predicted_risk_groups.png)

---

## 7. Key Findings

The combined analysis demonstrates an end-to-end precision-medicine workflow spanning molecular exploration, multi-omics analysis, predictive modeling, candidate biomarker discovery, and survival analysis.

Key findings include:

1. **Molecular heterogeneity:** AML patients showed substantial variation across their integrated molecular profiles.

2. **Predictive modeling:** XGBoost achieved a ROC-AUC of **0.726** on the initial held-out test split.

3. **Robustness assessment:** Five-fold cross-validation produced more conservative performance estimates, emphasizing the importance of evaluating model stability beyond a single split.

4. **Candidate biomarkers:** Feature-importance analysis identified a set of potentially informative mRNA and miRNA features.

5. **Survival stratification:** Model-derived prognosis groups showed significantly different survival distributions in the analyzed cohort (**log-rank p = 0.0044**).

---

## 8. Interpretation

The project demonstrates how multiple molecular data types can be explored and combined with machine-learning methods to investigate clinically relevant patterns in AML.

Importantly, the workflow distinguishes between two related analytical stages:

- **Exploratory multi-omics analysis**, which included mRNA, miRNA, DNA methylation, and earlier MOFA-based integration.
- **Predictive precision-medicine analysis**, which used matched clinical, mRNA, and miRNA data for machine learning, biomarker prioritization, and survival analysis.

The results support the potential value of multi-omics-informed computational approaches for AML research while also demonstrating the importance of cautious interpretation, cross-validation, and independent validation.

The identified biomarkers and predictive results should therefore be considered **exploratory research findings rather than clinically validated conclusions**.
