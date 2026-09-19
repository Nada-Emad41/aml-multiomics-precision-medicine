# Data

This project uses molecular and clinical data from **The Cancer Genome Atlas – Acute Myeloid Leukemia (TCGA-LAML)** cohort.

## Data Modalities

The broader multi-omics analysis includes:

- **mRNA expression**
- **miRNA expression**
- **DNA methylation**
- **Clinical and survival information**

The later predictive modeling workflow uses matched **mRNA expression, miRNA expression, and clinical data** at the patient level.

DNA methylation was analyzed during the earlier molecular profiling and exploratory multi-omics integration phase and was not used as an input to the later machine-learning classifiers.

## Data Source

Data were obtained from the **TCGA-LAML** study and used in preprocessed form for the analytical workflows contained in this repository.

Patient identifiers were harmonized across molecular and clinical datasets before integration.

## Data Availability

Raw and intermediate datasets are **not included in this repository** because of their size and to keep the repository focused on reproducible analysis code, methodology, and results.

The repository instead provides:

- Analysis scripts for the individual molecular modalities
- Multi-omics integration code
- The predictive modeling notebook
- Processed analytical outputs and figures
- Methodology and results documentation

## Reproducibility Note

The project consists of two related analytical phases:

1. **Exploratory multi-omics profiling and integration** — mRNA, miRNA, and DNA methylation, including MOFA-based exploratory integration.
2. **Precision-medicine predictive modeling** — matched mRNA, miRNA, and clinical data used for machine-learning comparison, candidate biomarker identification, and survival analysis.

These phases are documented separately to avoid implying that DNA methylation or MOFA-derived factors were inputs to the later predictive models.
