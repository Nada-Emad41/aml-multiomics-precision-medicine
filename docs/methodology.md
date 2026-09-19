# Methodology

## 1. Study Objective

This project investigates whether molecular profiles from acute myeloid leukemia (AML) patients can support clinically relevant patient stratification and survival-related prediction. The analysis combines exploratory multi-omics profiling with a later machine-learning workflow for precision-medicine analysis.

The project was developed in two related analytical phases using TCGA-LAML data:

1. Exploratory analysis of individual omics layers and multi-omics integration.
2. Predictive modeling using matched clinical, mRNA, and miRNA data, followed by candidate biomarker and survival analyses.

## 2. Data Sources and Molecular Layers

The analysis used preprocessed TCGA-LAML molecular and clinical data.

The broader exploratory phase examined:

- mRNA gene-expression data
- miRNA-expression data
- DNA-methylation data
- Clinical and survival information

The predictive machine-learning phase focused on matched clinical, mRNA, and miRNA data. DNA methylation was part of the earlier molecular profiling and integration work but was not used as an input to the later predictive models.

## 3. Exploratory Molecular Analysis

Each molecular layer was initially examined independently to characterize its distribution, variability, and relationship with available clinical metadata.

The exploratory workflow included:

- Data inspection and preprocessing
- Expression/distribution assessment
- Normalization where required
- Heatmap-based molecular profiling
- Dimensionality reduction using PCA/MDS
- Investigation of associations between molecular patterns and clinical metadata
- Differential molecular-feature analysis

The earlier multi-omics phase also explored integration across molecular layers using MOFA to investigate shared sources of variation and molecular heterogeneity.

## 4. Predictive Multi-Omics Dataset

For the later precision-medicine analysis, patient identifiers were harmonized across the clinical, mRNA, and miRNA datasets.

After matching and preprocessing, the modeling cohort contained 162 patients.

To reduce dimensionality while retaining informative molecular variation:

- 500 highly variable mRNA features were selected.
- 100 highly variable miRNA features were selected.

This produced a combined molecular feature set of 600 features for predictive modeling.

## 5. Machine-Learning Workflow

Three supervised classification algorithms were evaluated:

- Random Forest
- Support Vector Machine (SVM)
- XGBoost

The dataset was divided using a stratified 75/25 train-test split to preserve the outcome distribution.

Model evaluation considered multiple performance measures rather than accuracy alone, including:

- Accuracy
- ROC-AUC
- Sensitivity / Recall
- Specificity
- F1 score

Five-fold cross-validation was additionally used to assess the stability and generalizability of model performance and to compare cross-validated results with the initial holdout evaluation.

## 6. Candidate Biomarker Analysis

Feature importance from the selected XGBoost model was used to rank molecular features contributing to prediction.

The top-ranked features were treated as **candidate biomarkers**, rather than validated clinical biomarkers, because the analysis did not include independent external validation or experimental confirmation.

The resulting ranked candidate features were exported for reproducibility and downstream interpretation.

## 7. Survival Analysis

Model-derived patient groups were subsequently evaluated using Kaplan-Meier survival analysis.

Survival curves were used to examine whether the predicted prognostic groups showed different observed survival patterns. This step provided an additional clinical interpretation of the predictive modeling results rather than serving as an independent validation cohort.

## 8. Reproducibility

The repository separates the analytical workflow into:

- `notebooks/` — end-to-end precision-medicine analysis
- `scripts/` — individual omics and multi-omics exploratory analyses
- `results/` — structured analytical outputs
- `results/figures/` — selected visualizations

Large source datasets are not stored directly in the repository. The repository instead focuses on the analysis workflow, selected outputs, methodology, and reproducible code.
