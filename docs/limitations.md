# Limitations

This project demonstrates an end-to-end multi-omics and precision-medicine workflow using TCGA-LAML data. The results should be interpreted as exploratory and research-oriented rather than as clinically validated predictions.

## 1. Limited Sample Size

After harmonizing the clinical, mRNA, and miRNA datasets, the predictive modeling analysis was based on a relatively small matched patient cohort.

This limits statistical power and increases the possibility that model performance may vary across different train/test splits.

## 2. No External Validation Cohort

The machine-learning models were developed and evaluated using patients from the TCGA-LAML cohort.

No independent AML cohort was used for external validation. Therefore, the reported predictive performance should not be interpreted as evidence of generalization to other patient populations.

## 3. Survival Outcome Simplification

For the machine-learning classification task, survival information was converted into a binary outcome based on patient vital status.

Although this enables comparison of classification algorithms, it does not fully represent the time-to-event nature of survival data.

Kaplan-Meier analysis was subsequently used to evaluate survival differences between model-derived patient groups.

## 4. Initial Split vs. Cross-Validation Performance

The initial train/test split produced stronger predictive performance than the more conservative cross-validation evaluation.

This difference suggests that results from a single split may overestimate model performance in a relatively small biomedical dataset.

For this reason, cross-validation results are included to provide a more cautious assessment of model stability.

## 5. Candidate Biomarkers Require Biological Validation

Feature-importance analysis was used to identify molecular features associated with model predictions.

These features should be considered **candidate biomarkers**, not validated clinical biomarkers.

Additional biological investigation, independent cohort validation, and experimental confirmation would be required before drawing clinical conclusions.

## 6. Multi-Omics Scope Differs Between Analytical Phases

The earlier exploratory phase examined **mRNA expression, miRNA expression, and DNA methylation**, including exploratory MOFA-based multi-omics integration.

The later predictive modeling phase used matched **mRNA expression and miRNA expression together with clinical information**.

Therefore, DNA methylation and MOFA-derived latent factors should not be interpreted as inputs to the later machine-learning classifiers.

## 7. Clinical Applicability

The analysis was developed as an academic and portfolio research project.

It demonstrates analytical methods for molecular profiling, multi-omics integration, machine learning, biomarker exploration, and survival analysis, but it is **not intended for clinical diagnosis, prognosis, or treatment decision-making**.
