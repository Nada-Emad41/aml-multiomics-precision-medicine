#==============================#
#        Load Libraries        #
#==============================#
library(BiocManager)
library(DESeq2)
library(NMF)
library(ggplot2)
library(purrr)
library(factoextra)
library(ggfortify)
library("MIRit")
library(survival)
library(survminer)
library(tidyverse)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(limma)
library(rlang)  # for .data pronoun
library(SummarizedExperiment)



#==============================================================#
#   Handling miRNA Expression Matrix and Clinical Metadata     #
#==============================================================#

# Load miRNA expression data as a matrix from file (adjust path and separator if needed)
mirna_data <- as.matrix(read.table("G:/integrative_course/project_AML/miRNA_DATA/mirna", quote="\"", comment.char=""))

# Display the first few rows and columns
head(mirna_data)
head(mirna_data[, 1:5])

# Check dimensions of the miRNA expression matrix
dim(mirna_data)  # Should return 705 188

# Load clinical metadata (assuming sampleID is in the second column and no row names in file)
clinical_data <- read.delim("G:/integrative_course/project_AML/miRNA_DATA/aml", row.names=1)

# Display the first few rows
head(clinical_data)

# Check dimensions of the clinical data
dim(clinical_data)  # Should return 200 97

# Extract sample IDs from clinical data
clinical_samples <- rownames(clinical_data)
clinical_samples[1]  # Example: "TCGA-AB-2802-03"

# Extract sample names from miRNA data
mirna_samples <- colnames(mirna_data)
mirna_samples[1]  # Example: "TCGA.AB.2802.03"

# Fix miRNA sample names to match clinical format (replace '.' with '-')
mirna_samples_fixed <- gsub("\\.", "-", mirna_samples)
colnames(mirna_data) <- mirna_samples_fixed

# Display sample names after replacement
head(mirna_samples_fixed)
head(clinical_samples)

# Identify common samples between clinical and miRNA datasets
common_samples <- intersect(mirna_samples_fixed, clinical_samples)
cat("Number of common samples:", length(common_samples), "\n")  # Expected: 188

# Subset miRNA expression data to only common samples
mirna_expr_data_filtered <- mirna_data[, common_samples]

# Subset clinical metadata to only common samples
clinical_data_filtered <- clinical_data[rownames(clinical_data) %in% common_samples, ]

# Check dimensions of the filtered datasets
dim(mirna_expr_data_filtered)  # Should return 705 188
dim(clinical_data_filtered)    # Should return 188 97

# Verify if sample IDs perfectly match between datasets
all(colnames(mirna_expr_data_filtered) == clinical_data_filtered$sampleID)


#==============================================================#
#      Handling Survival Data with miRNA & Clinical Metadata   #
#==============================================================#

# Load survival data
survival_data <- read.delim("G:/integrative_course/project_AML/miRNA_DATA/survival")

# Display first few rows
head(survival_data)

# Check dimensions of survival data
dim(survival_data)  # Should return 651 3

# Identify duplicated sample IDs (assumed to be in the first column)
duplicates <- survival_data[duplicated(survival_data[, 1]), 1]
length(duplicates)  # Number of duplicated samples (275)

# Display first few duplicates
head(duplicates)

# Remove duplicated rows by keeping the first occurrence
survival_data <- survival_data[!duplicated(survival_data[, 1]), ]

# Set sample ID as row names
rownames(survival_data) <- survival_data[, 1]

# Drop the first column (already set as row names)
survival_data <- survival_data[, -1]

# Check new dimensions after removing duplicates
dim(survival_data) # 376 2 

# Check format of first five sample ID to check remove duplicated
rownames(survival_data)[1:5]

#==============================================================#
#      Match Survival Data with Filtered Clinical & miRNA      #
#==============================================================#

# Identify common sample IDs across all datasets
clicn_suv_datamatch <- intersect(rownames(clinical_data_filtered), rownames(survival_data))
length(clicn_suv_datamatch)

# Subset datasets to include only matched samples
matched_surv <- survival_data[clicn_suv_datamatch, ]
matched_surv_clin <- clinical_data_filtered[clicn_suv_datamatch, ] # match survival with clinical data 
matched_surv_mir <- mirna_expr_data_filtered[, clicn_suv_datamatch]# match survival with mirna data 

# Check dimensions of the matched datasets
dim(matched_surv) # after match with clinical survival data dimensions is 188 2 from 376 2 
dim(matched_surv_clin) # 188 97 
dim(matched_surv_mir)  # 705 188

# Confirm sample ID consistency across all matched datasets
rownames(matched_surv)[1]
rownames(matched_surv_clin[1, ])
colnames(matched_surv_mir)[1]

# Save miRNA feature names for downstream analysis
mir_names <- rownames(matched_surv_mir)


sum(is.na(matched_surv_mir))
sum(is.na(matched_surv_clin))
sum(is.na(matched_surv))
colSums(is.na(matched_surv_clin))
colSums(is.na(matched_surv))
dim(matched_surv_mir)
##### Visualization of miRNA expression data in AML #####

# Histogram of raw miRNA expression
hist(matched_surv_mir, col = "pink", main = "Acute Myeloid Leukemia miRNA Expression Histogram",
     xlab = "miRNA Expression")

# Histogram of log2-transformed miRNA expression (+1 to avoid log(0))
hist(log2(matched_surv_mir + 1), col = "pink", main = "Acute Myeloid Leukemia miRNA Expression Histogram",
     xlab = "Log2 miRNA Expression")



# Boxplot for first 5 miRNAs (log2 transformed)
boxplot(log2(matched_surv_mir[1:5, ] + 1), main = "Boxplot of First 5 miRNAs Expression",
        col = rainbow(5), cex.axis = 0.8)

# Boxplot for all miRNAs
boxplot(matched_surv_mir, main = "Boxplot of miRNA Expression",
        col = rainbow(ncol(matched_surv_mir)), las = 2, cex.axis = 0.8, outline = FALSE)

# QQ plot for normality check on the 3rd miRNA
qqnorm(matched_surv_mir[3, ])
qqline(matched_surv_mir[3, ])

# Shapiro-Wilk test for normality on miRNA at column 100
shapiro.test(matched_surv_mir[, 100])

# Boxplot of first miRNA expression by Gender
miRNA_name <- rownames(matched_surv_mir)[1]
expression_df <- data.frame(Expression = as.numeric(matched_surv_mir[miRNA_name, ]),
                            Gender = matched_surv_clin$gender)

ggplot(expression_df, aes(x = Gender, y = Expression, fill = Gender)) +
    geom_boxplot() +
    labs(title = paste(miRNA_name, "Expression by Gender"), y = "Expression Level", x = "Gender") +
    scale_fill_manual(values = c("pink", "lightblue")) +
    theme_minimal()

# Top 10 miRNAs boxplots by Gender (faceted)
selected_miRNAs <- rownames(matched_surv_mir)[1:10]
long_df <- purrr::map_dfr(selected_miRNAs, function(mir) {
    data.frame(miRNA = mir, Expression = as.numeric(matched_surv_mir[mir, ]),
               Gender = matched_surv_clin$gender)
})
ggplot(long_df, aes(x = Gender, y = Expression, fill = Gender)) +
    geom_boxplot() +
    facet_wrap(~ miRNA, scales = "free_y", ncol = 5) +
    labs(title = "Top 10 miRNAs Expression by Gender", y = "Expression Level", x = "Gender") +
    scale_fill_manual(values = c("pink", "lightblue")) +
    theme_minimal() +
    theme(strip.text = element_text(size = 10))

# Boxplots of first 4 miRNAs expression by Vital Status
selected_miRNAs <- rownames(matched_surv_mir)[1:4]
long_df <- purrr::map_dfr(selected_miRNAs, function(mir) {
    data.frame(miRNA = mir, Expression = as.numeric(matched_surv_mir[mir, ]),
               ClinicalGroup = matched_surv_clin$vital_status)
})
ggplot(long_df, aes(x = ClinicalGroup, y = Expression, fill = ClinicalGroup)) +
    geom_boxplot() +
    facet_wrap(~ miRNA, scales = "free_y", ncol = 2) +
    labs(title = "First 4 miRNAs Expression by Vital Status", y = "Expression Level", x = "Vital Status") +
    theme_minimal() +
    theme(strip.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1))

# Boxplots of first 4 miRNAs expression by X_PANCAN_mirna_LAML
long_df <- purrr::map_dfr(selected_miRNAs, function(mir) {
    data.frame(miRNA = mir, Expression = as.numeric(matched_surv_mir[mir, ]),
               ClinicalGroup = matched_surv_clin$X_PANCAN_mirna_LAML)
})
ggplot(long_df, aes(x = ClinicalGroup, y = Expression, fill = ClinicalGroup)) +
    geom_boxplot() +
    facet_wrap(~ miRNA, scales = "free_y", ncol = 2) +
    labs(title = "First 4 miRNAs Expression by X_PANCAN_mirna_LAML", y = "Expression Level", x = "X_PANCAN_mirna_LAML") +
    theme_minimal() +
    theme(strip.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1))


# Boxplots of first 4 miRNAs expression by X_PANCAN_DNAMethyl_LAML
long_df <- purrr::map_dfr(selected_miRNAs, function(mir) {
    data.frame(miRNA = mir, Expression = as.numeric(matched_surv_mir[mir, ]),
               ClinicalGroup = matched_surv_clin$X_PANCAN_DNAMethyl_LAML)
})
ggplot(long_df, aes(x = ClinicalGroup, y = Expression, fill = ClinicalGroup)) +
    geom_boxplot() +
    facet_wrap(~ miRNA, scales = "free_y", ncol = 2) +
    labs(title = "First 4 miRNAs Expression by X_PANCAN_DNAMethyl_LAML ", y = "Expression Level", x = "X_PANCAN_DNAMethyl_LAML") +
    theme_minimal() +
    theme(strip.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1))


############# for X_EVENT############

# Step 1: Prepare the survival object
surv_object <- Surv(time = matched_surv$Survival, event = matched_surv$Death)

# Step 2: Clean and relabel the X_EVENT variable
matched_surv_clin_clean <- matched_surv_clin %>%
    mutate(X_EVENT = case_when(
        X_EVENT == "0" ~ "Healthy",
        X_EVENT == "1" ~ "Leukemia",
        TRUE ~ as.character(X_EVENT)
    ))

# Step 3: Fit the Kaplan-Meier survival model
km_fit <- survfit(surv_object ~ X_EVENT, data = matched_surv_clin_clean)

# Step 4: Plot the Kaplan-Meier survival curve
ggsurvplot(
    km_fit,
    data = matched_surv_clin_clean,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = TRUE,
    xlab = "Survival Days",
    ylab = "Survival Probability",
    title = "Kaplan-Meier Survival Curve by X_EVENT (Healthy vs Leukemia)",
    legend.title = "Group",
    palette = c("#8E44AD", "#27AE60"),
    risk.table.y.text.col = TRUE
)
############# for vital_status ############

# Step 1: Prepare the survival object
surv_object <- Surv(time = matched_surv$Survival, event = matched_surv$Death)

# Step 2: Clean and relabel the vital_status variable
matched_surv_clin_clean <- matched_surv_clin %>%
    mutate(vital_status = case_when(
        vital_status %in% c("Alive", "Living") ~ "Alive",
        vital_status %in% c("Dead", "Deceased") ~ "Deceased",
        TRUE ~ as.character(vital_status)
    ))

# Step 3: Fit Kaplan-Meier survival model
km_fit_vital <- survfit(surv_object ~ vital_status, data = matched_surv_clin_clean)

# Step 4: Plot Kaplan-Meier curve
ggsurvplot(
    km_fit_vital,
    data = matched_surv_clin_clean,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = TRUE,
    xlab = "Survival Days",
    ylab = "Survival Probability",
    title = "Kaplan-Meier Curve by Vital Status",
    legend.title = "Vital Status",
    palette = c("#1ABC9C", "#E74C3C"),
    risk.table.y.text.col = TRUE
)

###############For Gender###############

# Clean and relabel gender variable
matched_surv_clin_clean <- matched_surv_clin_clean %>%
    mutate(gender = case_when(
        gender %in% c("male", "Male") ~ "Male",
        gender %in% c("female", "Female") ~ "Female",
        TRUE ~ as.character(gender)
    ))

# Fit Kaplan-Meier survival model
km_fit_gender <- survfit(surv_object ~ gender, data = matched_surv_clin_clean)

# Plot Kaplan-Meier curve
ggsurvplot(
    km_fit_gender,
    data = matched_surv_clin_clean,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = TRUE,
    xlab = "Survival Days",
    ylab = "Survival Probability",
    title = "Kaplan-Meier Curve by Gender",
    legend.title = "Gender",
    palette = c("#3498DB", "#E67E22"),
    risk.table.y.text.col = TRUE
)

###############for cytogenetic risk category###############

# Step 1: Clean and relabel cytogenetic risk category
# ----------------------------------------------------
# We extract and clean the cytogenetic risk category from the clinical data.
# This variable is important in AML for predicting survival.
# Only keep entries labeled Favorable, Intermediate/Normal, or Poor.
matched_surv_clin_clean <- matched_surv_clin %>%
    mutate(cytogenetic_risk = case_when(
        acute_myeloid_leukemia_calgb_cytogenetics_risk_category == "Favorable" ~ "Favorable",
        acute_myeloid_leukemia_calgb_cytogenetics_risk_category == "Intermediate/Normal" ~ "Intermediate/Normal",
        acute_myeloid_leukemia_calgb_cytogenetics_risk_category == "Poor" ~ "Poor",
        TRUE ~ NA_character_ # Set all other values to NA
    )) %>%
    filter(!is.na(cytogenetic_risk))

# Step 2: Match the survival data to the filtered clinical data
matched_surv_clean <- matched_surv[rownames(matched_surv_clin_clean), ]

# Step 3: Create survival object AFTER filtering
surv_object_risk <- Surv(time = matched_surv_clean$Survival, event = matched_surv_clean$Death)

# Step 4: Fit Kaplan-Meier survival model using the cleaned risk categories
km_fit_risk <- survfit(surv_object_risk ~ cytogenetic_risk, data = matched_surv_clin_clean)

# Step 5: Plot Kaplan-Meier survival curve
ggsurvplot(
    km_fit_risk,
    data = matched_surv_clin_clean,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = TRUE,
    xlab = "Survival Days",
    ylab = "Survival Probability",
    title = "Kaplan-Meier Survival Curve by Cytogenetic Risk Category",
    legend.title = "Risk Category",
    palette = c("#2ECC71", "#F39C12", "#E74C3C")
)
# ==========================================================
# Kaplan-Meier Survival Analysis by Gender (Leukemia Only)
# ==========================================================

# Step 0: Clean and relabel X_EVENT in the clinical dataframe
# ------------------------------------------------------------
# This converts the event column from 0/1 to readable labels:
# - "0" → "Healthy"
# - "1" → "Leukemia"

matched_surv_clin_clean <- matched_surv_clin %>%
    mutate(X_EVENT = case_when(
        X_EVENT == "0" ~ "Healthy",
        X_EVENT == "1" ~ "Leukemia",
        TRUE ~ as.character(X_EVENT)
    ))

# Step 1: Subset clinical data to only include leukemia patients
# --------------------------------------------------------------
# We only want to study gender differences among leukemic patients,
# so we filter out the healthy ones.
matched_surv_clin_leukemic <- matched_surv_clin_clean %>%
    filter(X_EVENT == "Leukemia")

# Step 2: Match survival data to filtered clinical subset
# --------------------------------------------------------
# Ensure the survival data includes **only the same rows** as the filtered clinical data.
# This guarantees we are analyzing only leukemia patients with matching survival info.
matched_surv_leukemic <- matched_surv[rownames(matched_surv_clin_leukemic), ]

# Step 3: Create survival object
# -------------------------------
# Create the survival object using:
# - time = number of survival days
# - event = 1 if the patient died, 0 if still alive
surv_obj_leukemia_gender <- Surv(time = matched_surv_leukemic$Survival, event = matched_surv_leukemic$Death)

# Step 4: Clean the gender column
# --------------------------------
# Standardize gender labels ("male", "Male") to unified values:
# - "Male"
# - "Female"
matched_surv_clin_leukemic <- matched_surv_clin_leukemic %>%
    mutate(gender = case_when(
        gender %in% c("male", "Male") ~ "Male",
        gender %in% c("female", "Female") ~ "Female",
        TRUE ~ as.character(gender)
    ))

# Step 5: Fit the Kaplan-Meier survival model by gender
# ------------------------------------------------------
# We use the `survfit()` function to estimate survival curves
# for each gender among leukemia patients.
km_fit_gender_leukemic <- survfit(surv_obj_leukemia_gender ~ gender, data = matched_surv_clin_leukemic)



# Optional: Print the Kaplan-Meier survival model results
km_fit_gender_leukemic
# Output will show median survival time for males and females,
# number of events (deaths), and confidence intervals.
#Call: survfit(formula = surv_obj_leukemia_gender ~ gender, data = matched_surv_clin_leukemic)

#               n events median 0.95LCL 0.95UCL
#gender=FEMALE 52     52    274     212     366
#gender=MALE   60     60    244     184     365

# Step 6: Plot Kaplan-Meier survival curves using ggplot
# -------------------------------------------------------
# This creates a professional survival plot using `ggsurvplot()`:
# - risk.table = shows how many patients remain at each time point
# - pval = shows p-value for difference between curves (log-rank test)
# - conf.int = adds confidence intervals to the curves

ggsurvplot(
    km_fit_gender_leukemic,
    data = matched_surv_clin_leukemic,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = TRUE,
    xlab = "Survival Days",
    ylab = "Survival Probability",
    title = "Kaplan-Meier Curve by Gender in Leukemic Patients",
    legend.title = "Gender",
    palette = c("#2980B9", "#E67E22"),
    risk.table.y.text.col = TRUE
)
# ==========================================================
# Kaplan-Meier Survival Analysis by cytogenetics_risk_category (Leukemia Only)
# ==========================================================

# Step 0: Clean and relabel X_EVENT in clinical dataframe only
matched_surv_clin_clean <- matched_surv_clin %>%
    mutate(X_EVENT = case_when(
        X_EVENT == "0" ~ "Healthy",
        X_EVENT == "1" ~ "Leukemia",
        TRUE ~ as.character(X_EVENT)
    ))

# Step 1: Subset leukemic patients
matched_surv_clin_leukemic <- matched_surv_clin_clean %>%
    filter(X_EVENT == "Leukemia")

# Step 2: Match survival rows
matched_surv_leukemic <- matched_surv[rownames(matched_surv_clin_leukemic), ]

# Step 3: Clean cytogenetic risk
matched_surv_clin_leukemic <- matched_surv_clin_leukemic %>%
    mutate(cytogenetic_risk = case_when(
        acute_myeloid_leukemia_calgb_cytogenetics_risk_category == "Favorable" ~ "Favorable",
        acute_myeloid_leukemia_calgb_cytogenetics_risk_category == "Poor" ~ "Poor",
        acute_myeloid_leukemia_calgb_cytogenetics_risk_category == "Intermediate/Normal" ~ "Intermediate",
        TRUE ~ NA_character_
    )) %>%
    filter(!is.na(cytogenetic_risk))

# Step 4: Survival object
surv_obj_leukemia_risk <- Surv(
    time = matched_surv_leukemic[rownames(matched_surv_clin_leukemic), "Survival"],
    event = matched_surv_leukemic[rownames(matched_surv_clin_leukemic), "Death"]
)

# Step 5: Kaplan-Meier model
km_fit_risk_leukemic <- survfit(surv_obj_leukemia_risk ~ cytogenetic_risk, data = matched_surv_clin_leukemic)
#Call: survfit(formula = surv_obj_leukemia_risk ~ cytogenetic_risk, 
#data = matched_surv_clin_leukemic)
#                               n events median 0.95LCL 0.95UCL
#cytogenetic_risk=Favorable    12     12    396     151      NA
#cytogenetic_risk=Intermediate 70     70    304     242     366
#cytogenetic_risk=Poor         28     28    213      61     365
# Step 6: Plot
ggsurvplot(
    km_fit_risk_leukemic,
    data = matched_surv_clin_leukemic,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = TRUE,
    xlab = "Survival Days",
    ylab = "Survival Probability",
    title = "Kaplan-Meier Curve by Cytogenetic Risk in Leukemic Patients",
    legend.title = "Risk Group",
    palette = c("#27AE60", "#F39C12", "#C0392B"),
    risk.table.y.text.col = TRUE
)
# ==========================================================
# Kaplan-Meier Survival Analysis by age_group (Leukemia Only)
# ===========================================================

# Step 0: Clean and relabel X_EVENT
matched_surv_clin_clean <- matched_surv_clin %>%
    mutate(X_EVENT = case_when(
        X_EVENT == "0" ~ "Healthy",
        X_EVENT == "1" ~ "Leukemia",
        TRUE ~ as.character(X_EVENT)
    ))

# Step 1: Subset leukemic patients
matched_surv_clin_leukemic <- matched_surv_clin_clean %>%
    filter(X_EVENT == "Leukemia")

# Step 2: Match survival rows
matched_surv_leukemic <- matched_surv[rownames(matched_surv_clin_leukemic), ]

# Step 3: Create age group variable
matched_surv_clin_leukemic <- matched_surv_clin_leukemic %>%
    mutate(age_group = case_when(
        age_at_initial_pathologic_diagnosis < 40 ~ "<40",
        age_at_initial_pathologic_diagnosis >= 40 & age_at_initial_pathologic_diagnosis < 60 ~ "40-59",
        age_at_initial_pathologic_diagnosis >= 60 ~ "60+",
        TRUE ~ NA_character_
    )) %>%
    filter(!is.na(age_group))

# Step 4: Survival object
surv_obj_leukemia_age <- Surv(
    time = matched_surv_leukemic[rownames(matched_surv_clin_leukemic), "Survival"],
    event = matched_surv_leukemic[rownames(matched_surv_clin_leukemic), "Death"]
)

# Step 5: Kaplan-Meier model
km_fit_age_leukemic <- survfit(surv_obj_leukemia_age ~ age_group, data = matched_surv_clin_leukemic)
#Call: survfit(formula = surv_obj_leukemia_age ~ age_group, data = matched_surv_clin_leukemic)

#                 n events median 0.95LCL 0.95UCL
#age_group=<40   15     15    305     215     792
#age_group=40-59 28     28    366     306     761
#age_group=60+   69     69    184     123     273

# Step 6: Plot
ggsurvplot(
    km_fit_age_leukemic,
    data = matched_surv_clin_leukemic,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = TRUE,
    xlab = "Survival Days",
    ylab = "Survival Probability",
    title = "Kaplan-Meier Curve by Age Group in Leukemic Patients",
    legend.title = "Age Group",
    palette = c("#1ABC9C", "#F39C12", "#D35400"),
    risk.table.y.text.col = TRUE
)
################################################################################
# PCA Analysis on miRNA Expression with Clinical Variables
################################################################################

# Step 1: Transpose and Normalize the miRNA Expression Data
trans_mir <- t(matched_surv_mir)                     # Transpose matrix (samples rows, miRNAs columns)
trans_mir_scaled <- scale(trans_mir)                 # Scale (normalize) data

# Step 2: Remove miRNAs (columns) with any missing values (NA)
trans_mir_clean <- trans_mir_scaled[, colSums(is.na(trans_mir_scaled)) == 0]

# Step 3: Perform PCA on cleaned data
pca_model <- prcomp(trans_mir_clean, center = TRUE, scale. = TRUE)

# Step 4: Convert PCA results to a data frame of scores (sample projections)
pca_scores <- as.data.frame(pca_model$x)

# Step 5: Prepare data frame for plotting - first two PCs plus clinical metadata
pca_df <- cbind(pca_scores[, 1:5], matched_surv_clin)  # Include PC1 to PC5

# Step 6: Identify numerical clinical variables with valid variance and sufficient data
valid_vars <- names(matched_surv_clin)[sapply(matched_surv_clin, function(col) {
    is.numeric(col) && sum(!is.na(col)) > 1 && var(col, na.rm = TRUE) > 0
})]

# Step 7: Compute correlations between PC1 and each valid numerical clinical variable
cor_vals <- sapply(valid_vars, function(var) {
    x <- matched_surv_clin[[var]]
    if (sum(!is.na(pca_scores$PC1) & !is.na(x)) > 1) {
        cor(pca_scores$PC1, x, use = "complete.obs")
    } else {
        NA
    }
})

# Step 8: Sort correlations by absolute value descending and assign variable names
sorted_cors <- sort(abs(cor_vals), decreasing = TRUE)
names(sorted_cors) <- valid_vars


print(sorted_cors)

################################################################################
# Step 9: Define function to plot multiple PCA pairs colored by a numerical variable
################################################################################
plot_multiple_pca_numeric <- function(numeric_var_name, pca_df, pca_model) {
    pc_pairs <- list(
        c("PC1", "PC2"),
        c("PC2", "PC3"),
        c("PC3", "PC4"),
        c("PC4", "PC5")
    )
    
    pc_summary <- summary(pca_model)$importance[2, ]
    
    plot_data_list <- lapply(pc_pairs, function(pair) {
        pc_x <- pair[1]
        pc_y <- pair[2]
        
        temp_df <- pca_df[!is.na(pca_df[[numeric_var_name]]), ]
        temp_df$PCX <- temp_df[[pc_x]]
        temp_df$PCY <- temp_df[[pc_y]]
        temp_df$pair <- paste0(
            pc_x, " (", round(100 * pc_summary[as.numeric(sub("PC", "", pc_x))], 2), "%) vs ",
            pc_y, " (", round(100 * pc_summary[as.numeric(sub("PC", "", pc_y))], 2), "%)"
        )
        return(temp_df)
    })
    
    all_data <- do.call(rbind, plot_data_list)
    
    ggplot(all_data, aes(x = PCX, y = PCY, color = .data[[numeric_var_name]], shape = gender)) +
        geom_point(size = 3.5, alpha = 0.85) +
        facet_wrap(~pair, scales = "free", ncol = 2) +
        labs(
            title = paste("PCA Plots by", numeric_var_name),
            x = NULL,
            y = NULL,
            color = numeric_var_name
        ) +
        scale_color_gradientn(colors = c("blue", "green", "yellow", "red")) +
        theme_minimal(base_size = 13) +
        theme(
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            strip.text = element_text(size = 12)
        )
}

################################################################################
# Step 10: Plot top 3 correlated numerical variables across 4 PCA combinations
################################################################################
num_var1 <- names(sorted_cors)[1]
print(plot_multiple_pca_numeric(num_var1, pca_df, pca_model))

num_var2 <- names(sorted_cors)[2]
print(plot_multiple_pca_numeric(num_var2, pca_df, pca_model))

num_var3 <- names(sorted_cors)[3]
print(plot_multiple_pca_numeric(num_var3, pca_df, pca_model))

num_var4 <- names(sorted_cors)[4]
print(plot_multiple_pca_numeric(num_var4, pca_df, pca_model))

num_var5 <- names(sorted_cors)[5]
print(plot_multiple_pca_numeric(num_var5, pca_df, pca_model))

# Step 11: Identify categorical variables (factor or character)
categorical_vars <- names(matched_surv_clin)[
    sapply(matched_surv_clin, function(col) is.factor(col) || is.character(col))
]

# Step 12: Filter categorical variables with less than 50% NA and at least 2 unique levels
categorical_vars <- categorical_vars[sapply(categorical_vars, function(var) {
    na_ratio <- sum(is.na(matched_surv_clin[[var]])) / nrow(matched_surv_clin)
    levels_count <- length(unique(na.omit(matched_surv_clin[[var]])))
    (na_ratio < 0.5) && (levels_count > 1)
})]

print("Filtered categorical variables for ANOVA:")
print(categorical_vars)

# Step 13: Perform ANOVA of PC1 ~ each categorical variable
pca_with_metadata <- cbind(pca_scores, matched_surv_clin)

anova_results <- lapply(categorical_vars, function(var) {
    tryCatch({
        aov_model <- aov(PC1 ~ get(var), data = pca_with_metadata)
        summary(aov_model)
    }, error = function(e) NULL)
})

names(anova_results) <- categorical_vars

# Step 14: Extract p-values from ANOVA results
anova_pvals <- sapply(anova_results, function(result) {
    if (!is.null(result) && !is.null(result[[1]]) && "Pr(>F)" %in% colnames(result[[1]])) {
        result[[1]][["Pr(>F)"]][1]
    } else {
        NA
    }
})

print("ANOVA p-values for categorical variables:")
print(anova_pvals)

# Step 15: Filter significant categorical vars (p < 0.05) and sort
significant_vars <- anova_pvals[anova_pvals < 0.05]
sorted_sig_vars <- sort(significant_vars)
print("Significant categorical variables sorted by p-value:")
print(sorted_sig_vars)

################################################################################
# Step 16: Define function to plot PCA vs categorical variables
################################################################################
plot_multiple_pca_categorical <- function(cat_var_name, pca_df, pca_model) {
    pc_pairs <- list(
        c("PC1", "PC2"),
        c("PC2", "PC3"),
        c("PC3", "PC4"),
        c("PC4", "PC5")
    )
    
    pc_summary <- summary(pca_model)$importance[2, ]
    
    plot_data_list <- lapply(pc_pairs, function(pair) {
        pc_x <- pair[1]
        pc_y <- pair[2]
        
        temp_df <- pca_df[!is.na(pca_df[[cat_var_name]]), ]
        temp_df$PCX <- temp_df[[pc_x]]
        temp_df$PCY <- temp_df[[pc_y]]
        temp_df$pair <- paste0(
            pc_x, " (", round(100 * pc_summary[as.numeric(sub("PC", "", pc_x))], 2), "%) vs ",
            pc_y, " (", round(100 * pc_summary[as.numeric(sub("PC", "", pc_y))], 2), "%)"
        )
        return(temp_df)
    })
    
    all_data <- do.call(rbind, plot_data_list)
    
    categories <- unique(na.omit(all_data[[cat_var_name]]))
    n_colors <- max(length(categories), 3)
    color_palette <- RColorBrewer::brewer.pal(min(n_colors, 9), "Set3")
    if (length(categories) > length(color_palette)) {
        color_palette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(categories))
    }
    
    ggplot(all_data, aes(x = PCX, y = PCY, color = .data[[cat_var_name]], shape = gender)) +
        geom_point(size = 3.5, alpha = 0.85) +
        facet_wrap(~pair, scales = "free", ncol = 2) +
        labs(
            title = paste("PCA Plots by", cat_var_name),
            x = NULL,
            y = NULL,
            color = cat_var_name
        ) +
        scale_color_manual(values = color_palette) +
        theme_minimal(base_size = 13) +
        theme(
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            strip.text = element_text(size = 12)
        )
}


# Step 17: Plot PCA pair combinations for significant categorical variables
print(plot_multiple_pca_categorical("leukemia_french_american_british_morphology_code", pca_with_metadata, pca_model))
print(plot_multiple_pca_categorical("X_PANCAN_mirna_LAML", pca_with_metadata, pca_model))
print(plot_multiple_pca_categorical("acute_myeloid_leukemia_calgb_cytogenetics_risk_category", pca_with_metadata, pca_model))
print(plot_multiple_pca_categorical("X_PANCAN_DNAMethyl_LAML", pca_with_metadata, pca_model))


################################################################################
# Differential miRNA Expression Analysis (DESeq2)
################################################################################

# Check if column names of expression matrix match clinical data
all(colnames(matched_surv_mir) == rownames(matched_surv_clin))
all(colnames(matched_surv_mir) == rownames(matched_surv))

# Store miRNA gene names
miRNA_names <- rownames(matched_surv_mir)

# Round and convert expression values to integers
matched_surv_mir <- apply(round(matched_surv_mir), 2, as.integer)

# Set row names as miRNA gene names
rownames(matched_surv_mir) <- miRNA_names

# Count NA values
cat("Number of NA values in miRNA expression matrix: ", sum(is.na(matched_surv_mir)), "\n")

# Select binary clinical variables
sel_matched_surv_clin <- names(matched_surv_clin)[sapply(matched_surv_clin, function(column) {
    length(unique(column[!is.na(column)])) == 2
})]

# Print selected variables and frequency tables
for (var in sel_matched_surv_clin) {
    print(var)
    print(table(matched_surv_clin[var]))
}

################################################################################

# DESeq2 Analysis by X_EVENT

# Remove samples with missing X_EVENT
cleaned_X_Event_clin <- matched_surv_clin[!is.na(matched_surv_clin$X_EVENT), ]
cleaned_X_Event_mir <- matched_surv_mir[, rownames(cleaned_X_Event_clin)]

# Recode X_EVENT factor levels
table(cleaned_X_Event_clin$X_EVENT)
cleaned_X_Event_clin$X_EVENT <- factor(ifelse(cleaned_X_Event_clin$X_EVENT == 1, "Leukemia", "Healthy"))

# Create DESeq2 dataset
dds_X_Event <- DESeqDataSetFromMatrix(
    countData = cleaned_X_Event_mir,
    colData = cleaned_X_Event_clin,
    design = ~ X_EVENT
)

# Run DESeq2
dds_1 <- DESeq(dds_X_Event)

# Extract DESeq2 results
res_1 <- results(dds_1, contrast = c("X_EVENT", "Leukemia", "Healthy"))

# Remove NA results
res_1 <- as.data.frame(res_1[complete.cases(res_1), ])

# Filter significant DE miRNAs
sig_mir_1 <- res_1[res_1$padj < 0.05 & abs(res_1$log2FoldChange) > 1.2, ]

# Export significant results
write.csv(sig_mir_1, file = "deseq_significant_DEGs_X_EVENT.csv", quote = FALSE, row.names = TRUE)

# Plot Volcano Plot
par(mfrow = c(1, 1))
plot(
    res_1$log2FoldChange,
    -log10(res_1$padj),
    pch = 20,
    col = ifelse(res_1$padj < 0.05 & res_1$log2FoldChange > 1.2, "blue",
                 ifelse(res_1$padj < 0.05 & res_1$log2FoldChange < -1.2, "red", "grey")),
    main = "Volcano Plot: Leukemia vs Healthy DEmiRNAs",
    xlab = "Log2 Fold Change",
    ylab = "-Log10 Adjusted P-value",
    cex = 1.2
)
abline(h = -log10(0.05), col = "black", lty = 2)
abline(v = c(-1.2, 1.2), col = "black", lty = 2)
legend("topright", legend = c("Upregulated", "Downregulated", "Not Significant"),
       col = c("blue", "red", "grey"), pch = 20, bty = "n", cex = 0.8)

# Plot heatmap
dds_1_heatmap <- estimateSizeFactors(dds_1)
normalized_counts1 <- as.data.frame(counts(dds_1_heatmap, normalized = TRUE))
DEmiRNA1 <- as.matrix(normalized_counts1[rownames(normalized_counts1) %in% rownames(sig_mir_1), ])
aheatmap(log2(DEmiRNA1 + 1), annCol = cleaned_X_Event_clin$X_EVENT, main = "Heatmap of Significant DEmiRNAs")

# Plot PCA
vsd_1 <- varianceStabilizingTransformation(dds_1, blind = TRUE)
pca_model <- prcomp(t(assay(vsd_1)))
pca_scores <- as.data.frame(pca_model$x)
percentVar <- round(100 * summary(pca_model)$importance[2, 1:5], 2)
names(percentVar) <- paste0("PC", 1:5)
sample_info <- colData(vsd_1)
pca_scores$X_EVENT <- sample_info$X_EVENT
make_pca_plot <- function(df, x_pc, y_pc, percentVar, group_var = "X_EVENT") {
    ggplot(df, aes_string(x = x_pc, y = y_pc, color = group_var)) +
        geom_point(size = 3) +
        labs(
            title = paste(x_pc, "vs", y_pc),
            x = paste0(x_pc, ": ", percentVar[x_pc], "% variance"),
            y = paste0(y_pc, ": ", percentVar[y_pc], "% variance"),
            color = group_var
        ) +
        theme_minimal() +
        theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))
}
p1 <- make_pca_plot(pca_scores, "PC1", "PC2", percentVar)
p2 <- make_pca_plot(pca_scores, "PC2", "PC3", percentVar)
p3 <- make_pca_plot(pca_scores, "PC3", "PC4", percentVar)
p4 <- make_pca_plot(pca_scores, "PC4", "PC5", percentVar)
combined_plot <- (p1 | p2) / (p3 | p4)
print(combined_plot)

################################################################################

# DESeq2 Analysis by Gender

cleaned_gender_clin <- matched_surv_clin[!is.na(matched_surv_clin$gender), ]
cleaned_gender_mir <- matched_surv_mir[, rownames(cleaned_gender_clin)]
cleaned_gender_clin$gender <- as.factor(cleaned_gender_clin$gender)

dds_gender <- DESeqDataSetFromMatrix(
    countData = cleaned_gender_mir,
    colData = cleaned_gender_clin,
    design = ~ gender
)

dds_2 <- DESeq(dds_gender)
res_2 <- results(dds_2, contrast = c("gender", "FEMALE", "MALE"))
res_2 <- as.data.frame(res_2[complete.cases(res_2), ])
sig_mir_2 <- res_2[res_2$padj < 0.05 & abs(res_2$log2FoldChange) > 0.05, ]
write.csv(sig_mir_2, file = "deseq_significant_DEGs_Gender.csv", quote = FALSE, row.names = TRUE)

# Plot Volcano Plot
par(mfrow = c(1, 1))
plot(
    res_2$log2FoldChange,
    -log10(res_2$padj),
    pch = 20,
    col = ifelse(res_2$padj < 0.05 & res_2$log2FoldChange > 1.2, "blue",
                 ifelse(res_2$padj < 0.05 & res_2$log2FoldChange < -1.2, "red", "grey")),
    main = "Volcano Plot: Male vs Female DEmiRNAs",
    xlab = "Log2 Fold Change",
    ylab = "-Log10 Adjusted P-value",
    cex = 1.2
)
abline(h = -log10(0.05), col = "black", lty = 2)
abline(v = c(-1.2, 1.2), col = "black", lty = 2)
legend("topright", legend = c("Upregulated", "Downregulated", "Not Significant"),
       col = c("blue", "red", "grey"), pch = 20, bty = "n", cex = 0.8)

# Plot heatmap
dds_2_heatmap <- estimateSizeFactors(dds_2)
normalized_counts2 <- as.data.frame(counts(dds_2_heatmap, normalized = TRUE))
DEmiRNA2 <- as.matrix(normalized_counts2[rownames(normalized_counts2) %in% rownames(sig_mir_2), ])
aheatmap(log2(DEmiRNA2 + 1), annCol = cleaned_gender_clin$gender, main = "Heatmap of Significant DEmiRNAs")


# Plot PCA
vsd_2 <- varianceStabilizingTransformation(dds_2, blind = TRUE)
pca_model <- prcomp(t(assay(vsd_2)))
pca_scores <- as.data.frame(pca_model$x)
percentVar <- round(100 * summary(pca_model)$importance[2, 1:5], 2)
names(percentVar) <- paste0("PC", 1:5)
sample_info <- colData(vsd_2)
pca_scores$gender <- sample_info$gender
make_pca_plot <- function(df, x_pc, y_pc, percentVar, group_var = "gender") {
    ggplot(df, aes_string(x = x_pc, y = y_pc, color = group_var)) +
        geom_point(size = 3) +
        labs(
            title = paste(x_pc, "vs", y_pc),
            x = paste0(x_pc, ": ", percentVar[x_pc], "% variance"),
            y = paste0(y_pc, ": ", percentVar[y_pc], "% variance"),
            color = group_var
        ) +
        theme_minimal() +
        theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))
}
p1 <- make_pca_plot(pca_scores, "PC1", "PC2", percentVar)
p2 <- make_pca_plot(pca_scores, "PC2", "PC3", percentVar)
p3 <- make_pca_plot(pca_scores, "PC3", "PC4", percentVar)
p4 <- make_pca_plot(pca_scores, "PC4", "PC5", percentVar)
combined_plot <- (p1 | p2) / (p3 | p4)
print(combined_plot)

#  DESeq2 Analysis by vital_status
cleaned_vital_status_clin <- matched_surv_clin[!is.na(matched_surv_clin$vital_status), ]
cleaned_vital_status_mir <- matched_surv_mir[, rownames(cleaned_vital_status_clin)]
cleaned_vital_status_clin$vital_status <- as.factor(cleaned_vital_status_clin$vital_status)

dds_vital_status <- DESeqDataSetFromMatrix(
    countData = cleaned_vital_status_mir,
    colData = cleaned_vital_status_clin,
    design = ~ vital_status
)

dds_3 <- DESeq(dds_vital_status)
res_3 <- results(dds_3, contrast = c("vital_status", "LIVING", "DECEASED"))
res_3 <- as.data.frame(res_3[complete.cases(res_3), ])
sig_mir_3 <- res_3[res_3$padj < 0.05 & abs(res_3$log2FoldChange) > 0.05, ]
# Export significant DEGs to CSV
write.csv(sig_mir_3, file = "deseq_significant_DEGs_vital_status.csv", quote = FALSE, row.names = TRUE)

# Draw DEGs volcano plot
par(mfrow = c(1, 1))
# Create the volcano plot
plot(
    res_3$log2FoldChange, 
    -log10(res_3$padj), 
    pch = 20, 
    col = ifelse(res_3$padj < 0.05 & res_3$log2FoldChange > 1.2, "blue", 
                 ifelse(res_3$padj < 0.05 & res_3$log2FoldChange < -1.2, "red", "grey")),
    main = "Volcano Plot: Living vs Deceased DEmiRNAs",
    xlab = "Log2 Fold Change",
    ylab = "-Log10 Adjusted P-value",
    cex = 1.2  # Adjust point size
)
# Add significance threshold lines
abline(h = -log10(0.05), col = "black", lty = 2)  # Horizontal line for p-value threshold
abline(v = c(-1.2, 1.2), col = "black", lty = 2)  # Vertical lines for fold change thresholds

# Add a legend
legend("topleft", 
       legend = c("Upregulated", "Downregulated", "Not Significant"), 
       col = c("blue", "red", "grey"), 
       pch = 20, 
       bty = "n",  # No box around the legend
       cex = 0.8)  # Adjust legend text size

#Draw heatmap
dds_3_heatmap <- estimateSizeFactors(dds_3)
normalized_counts3 <- as.data.frame(counts(dds_3_heatmap, normalized = TRUE))
DEmiRNA3 <- as.matrix(normalized_counts3[rownames(normalized_counts3) %in% rownames(sig_mir_3), ])
aheatmap(log2(DEmiRNA3 + 1), annCol = cleaned_vital_status_clin$vital_status, main = "Heatmap of Significant DEmiRNAs")



# Plot PCA
vsd_3 <- varianceStabilizingTransformation(dds_3, blind = TRUE)
pca_model <- prcomp(t(assay(vsd_3)))
pca_scores <- as.data.frame(pca_model$x)
percentVar <- round(100 * summary(pca_model)$importance[2, 1:5], 2)
names(percentVar) <- paste0("PC", 1:5)
sample_info <- colData(vsd_3)
pca_scores$vital_status <- sample_info$vital_status
make_pca_plot <- function(df, x_pc, y_pc, percentVar, group_var = "vital_status") {
    ggplot(df, aes_string(x = x_pc, y = y_pc, color = group_var)) +
        geom_point(size = 3) +
        labs(
            title = paste(x_pc, "vs", y_pc),
            x = paste0(x_pc, ": ", percentVar[x_pc], "% variance"),
            y = paste0(y_pc, ": ", percentVar[y_pc], "% variance"),
            color = group_var
        ) +
        theme_minimal() +
        theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))
}
p1 <- make_pca_plot(pca_scores, "PC1", "PC2", percentVar)
p2 <- make_pca_plot(pca_scores, "PC2", "PC3", percentVar)
p3 <- make_pca_plot(pca_scores, "PC3", "PC4", percentVar)
p4 <- make_pca_plot(pca_scores, "PC4", "PC5", percentVar)
combined_plot <- (p1 | p2) / (p3 | p4)
print(combined_plot)

# =======================
# Phase II: Run DESeq2 
# =======================

# Phase II: Re-run DESeq2 on the most significant categorical variables and run PCA

# Example: for Morphology

sorted_sig_vars

# Filter and prepare morphology data
cleaned_morphology_clin <- matched_surv_clin[!is.na(matched_surv_clin$leukemia_french_american_british_morphology_code), ]
cleaned_morphology_mir <- matched_surv_mir[, rownames(cleaned_morphology_clin)]

# Clean factor levels
cleaned_morphology_clin$leukemia_french_american_british_morphology_code[cleaned_morphology_clin$leukemia_french_american_british_morphology_code == "M0 Undifferentiated"] <- "M0_Undifferentiated"
cleaned_morphology_clin$leukemia_french_american_british_morphology_code[cleaned_morphology_clin$leukemia_french_american_british_morphology_code == "Not Classified"] <- "Not_Classified"
cleaned_morphology_clin$leukemia_french_american_british_morphology_code <- as.factor(cleaned_morphology_clin$leukemia_french_american_british_morphology_code)

# Step 1: Create DESeq2 dataset and run DESeq
dds_morphology <- DESeqDataSetFromMatrix(
    countData = cleaned_morphology_mir,
    colData = cleaned_morphology_clin,
    design = ~ leukemia_french_american_british_morphology_code
)

dds_4 <- DESeq(dds_morphology)

# Step 2: Extract results for M0 vs M1
res_4 <- results(dds_4, contrast = c("leukemia_french_american_british_morphology_code", "M0_Undifferentiated", "M1"))
res_4 <- as.data.frame(res_4[complete.cases(res_4), ])
sig_mir_4 <- res_4[res_4$padj < 0.05 & abs(res_4$log2FoldChange) > 0.05, ]

# Step 3: Export significant miRNAs
write.csv(sig_mir_4, file = "deseq_significant_DEGs_morphology.csv", quote = FALSE, row.names = TRUE)

# Step 4: Volcano plot
par(mfrow = c(1, 1))
plot(
    res_4$log2FoldChange,
    -log10(res_4$padj),
    pch = 20,
    col = ifelse(res_4$padj < 0.05 & res_4$log2FoldChange > 1.2, "blue",
                 ifelse(res_4$padj < 0.05 & res_4$log2FoldChange < -1.2, "red", "grey")),
    main = "Volcano Plot: Morphology M0 vs M1",
    xlab = "Log2 Fold Change",
    ylab = "-Log10 Adjusted P-value",
    cex = 1.2
)
abline(h = -log10(0.05), col = "black", lty = 2)
abline(v = c(-1.2, 1.2), col = "black", lty = 2)
legend("topright", legend = c("Upregulated", "Downregulated", "Not Significant"), col = c("blue", "red", "grey"), pch = 20, bty = "n", cex = 0.8)

# Step 5: Heatmap for significant miRNAs
dds_4_heatmap <- estimateSizeFactors(dds_4)
normalized_counts4 <- as.data.frame(counts(dds_4_heatmap, normalized = TRUE))
DEmiRNA4 <- as.matrix(normalized_counts4[rownames(normalized_counts4) %in% rownames(sig_mir_4), ])

pdf("heatmap_DEmiRNAs.pdf", width = 20, height = 10)
aheatmap(
    log2(DEmiRNA4 + 1),
    color = "-RdYlBu2:100",
    scale = "none",
    Rowv = TRUE,
    Colv = TRUE,
    annCol = data.frame(Morphology = cleaned_morphology_clin$leukemia_french_american_british_morphology_code),
    main = "Heatmap of miRNA M0 vs M1",
    fontsize = 14,
    cexRow = 1.2,
    cexCol = 1.0,
    labRow = rownames(DEmiRNA4)
)
dev.off()

#pca
vsd_4 <- varianceStabilizingTransformation(dds_4, blind = TRUE)
pca_model <- prcomp(t(assay(vsd_4)))
pca_scores <- as.data.frame(pca_model$x)
percentVar <- round(100 * summary(pca_model)$importance[2, 1:5], 2)
names(percentVar) <- paste0("PC", 1:5)
sample_info <- colData(vsd_4)
pca_scores$leukemia_french_american_british_morphology_code <- sample_info$leukemia_french_american_british_morphology_code
make_pca_plot <- function(df, x_pc, y_pc, percentVar, group_var = "leukemia_french_american_british_morphology_code") {
    ggplot(df, aes_string(x = x_pc, y = y_pc, color = group_var)) +
        geom_point(size = 3) +
        labs(
            title = paste(x_pc, "vs", y_pc),
            x = paste0(x_pc, ": ", percentVar[x_pc], "% variance"),
            y = paste0(y_pc, ": ", percentVar[y_pc], "% variance"),
            color = group_var
        ) +
        theme_minimal() +
        theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))
}
p1 <- make_pca_plot(pca_scores, "PC1", "PC2", percentVar)
p2 <- make_pca_plot(pca_scores, "PC2", "PC3", percentVar)
p3 <- make_pca_plot(pca_scores, "PC3", "PC4", percentVar)
p4 <- make_pca_plot(pca_scores, "PC4", "PC5", percentVar)
combined_plot <- (p1 | p2) / (p3 | p4)
print(combined_plot)


# ================================
# Phase III: Run DESeq2 for X_PANCAN_mirna_LAML Clusters
# ================================

# Step 1: Clean clinical data
cleaned_PANCAN_mirna_clin <- matched_surv_clin[!is.na(matched_surv_clin$X_PANCAN_mirna_LAML), ]
cleaned_PANCAN_mirna_mir <- matched_surv_mir[, rownames(cleaned_PANCAN_mirna_clin)]
cleaned_PANCAN_mirna_clin$X_PANCAN_mirna_LAML <- gsub(" ", "_", cleaned_PANCAN_mirna_clin$X_PANCAN_mirna_LAML)
cleaned_PANCAN_mirna_clin$X_PANCAN_mirna_LAML <- as.factor(cleaned_PANCAN_mirna_clin$X_PANCAN_mirna_LAML)

# Step 2: Create DESeq2 dataset and run DESeq
dds_PANCAN_mirna <- DESeqDataSetFromMatrix(
    countData = cleaned_PANCAN_mirna_mir,
    colData = cleaned_PANCAN_mirna_clin,
    design = ~ X_PANCAN_mirna_LAML
)

dds_5 <- DESeq(dds_PANCAN_mirna)

# Step 3: Extract results for cluster_1 vs cluster_2
res_5 <- results(dds_5, contrast = c("X_PANCAN_mirna_LAML", "cluster_1", "cluster_2"))
res_5 <- as.data.frame(res_5[complete.cases(res_5), ])
sig_mir_5 <- res_5[res_5$padj < 0.05 & abs(res_5$log2FoldChange) > 0.05, ]
write.csv(sig_mir_5, file = "deseq_significant_DEGs_X_PANCAN_mirna_LAML.csv", quote = FALSE, row.names = TRUE)

# Step 4: Volcano plot
par(mfrow = c(1, 1))
plot(
    res_5$log2FoldChange,
    -log10(res_5$padj),
    pch = 20,
    col = ifelse(res_5$padj < 0.05 & res_5$log2FoldChange > 1.2, "blue",
                 ifelse(res_5$padj < 0.05 & res_5$log2FoldChange < -1.2, "red", "grey")),
    main = "Volcano Plot: Cluster 1 vs Cluster 2 DEmiRNAs",
    xlab = "Log2 Fold Change",
    ylab = "-Log10 Adjusted P-value",
    cex = 1.2
)
abline(h = -log10(0.05), col = "black", lty = 2)
abline(v = c(-1.2, 1.2), col = "black", lty = 2)
legend("topright", legend = c("Upregulated", "Downregulated", "Not Significant"), col = c("blue", "red", "grey"), pch = 20, bty = "n", cex = 0.8)

# Step 5: Heatmap for clusters
dds_5_heatmap <- estimateSizeFactors(dds_5)
normalized_counts5 <- as.data.frame(counts(dds_5_heatmap, normalized = TRUE))
DEmiRNA5 <- as.matrix(normalized_counts5[rownames(normalized_counts5) %in% rownames(sig_mir_5), ])

pdf("heatmap_DEmiRNAsX_PANCAN_mirna_LAML.pdf", width = 20, height = 10)
aheatmap(
    log2(DEmiRNA5 + 1),
    color = "-RdYlBu2:100",
    scale = "none",
    Rowv = TRUE,
    Colv = TRUE,
    annCol = data.frame(Cluster = cleaned_PANCAN_mirna_clin$X_PANCAN_mirna_LAML),
    main = "Heatmap of Cluster 1 vs 2",
    fontsize = 14,
    cexRow = 1.2,
    cexCol = 1.0,
    labRow = rownames(DEmiRNA5)
)
dev.off()


#pca
vsd_5 <- varianceStabilizingTransformation(dds_5, blind = TRUE)
pca_model <- prcomp(t(assay(vsd_5)))
pca_scores <- as.data.frame(pca_model$x)
percentVar <- round(100 * summary(pca_model)$importance[2, 1:5], 2)
names(percentVar) <- paste0("PC", 1:5)
sample_info <- colData(vsd_5)
pca_scores$X_PANCAN_mirna_LAML <- sample_info$X_PANCAN_mirna_LAML
make_pca_plot <- function(df, x_pc, y_pc, percentVar, group_var = "X_PANCAN_mirna_LAML") {
    ggplot(df, aes_string(x = x_pc, y = y_pc, color = group_var)) +
        geom_point(size = 3) +
        labs(
            title = paste(x_pc, "vs", y_pc),
            x = paste0(x_pc, ": ", percentVar[x_pc], "% variance"),
            y = paste0(y_pc, ": ", percentVar[y_pc], "% variance"),
            color = group_var
        ) +
        theme_minimal() +
        theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))
}
p1 <- make_pca_plot(pca_scores, "PC1", "PC2", percentVar)
p2 <- make_pca_plot(pca_scores, "PC2", "PC3", percentVar)
p3 <- make_pca_plot(pca_scores, "PC3", "PC4", percentVar)
p4 <- make_pca_plot(pca_scores, "PC4", "PC5", percentVar)
combined_plot <- (p1 | p2) / (p3 | p4)
print(combined_plot)


###########DNA Methylation subtype analysis#############

# Filter out samples with NA values in DNA methylation subtype
cleaned_PANCAN_DNAMethyl_clin <- matched_surv_clin[!is.na(matched_surv_clin$X_PANCAN_DNAMethyl_LAML), ]
cleaned_PANCAN_DNAMethyl_mir <- matched_surv_mir[, rownames(cleaned_PANCAN_DNAMethyl_clin)]

# Step 2: Check levels of DNA methylation subtype and convert to factor (replace spaces with underscores)
table(cleaned_PANCAN_DNAMethyl_clin$X_PANCAN_DNAMethyl_LAML)
cleaned_PANCAN_DNAMethyl_clin$X_PANCAN_DNAMethyl_LAML <- gsub(" ", "_", cleaned_PANCAN_DNAMethyl_clin$X_PANCAN_DNAMethyl_LAML)
cleaned_PANCAN_DNAMethyl_clin$X_PANCAN_DNAMethyl_LAML  <- as.factor(cleaned_PANCAN_DNAMethyl_clin$X_PANCAN_DNAMethyl_LAML)

# Create DESeqDataSet object for DNA methylation subtype groups
dds_PANCAN_DNAMethyl <- DESeqDataSetFromMatrix(
    countData = cleaned_PANCAN_DNAMethyl_mir,
    colData = cleaned_PANCAN_DNAMethyl_clin,
    design = ~ X_PANCAN_DNAMethyl_LAML
)

# Run DESeq2 differential expression analysis
dds_6 <- DESeq(dds_PANCAN_DNAMethyl)
res_6 <- results(dds_6, contrast = c("X_PANCAN_DNAMethyl_LAML","cluster_1", "cluster_2"))

# Remove NA values from results and convert to data frame
res_6 <- as.data.frame(res_6[complete.cases(res_6), ])

# Filter significant miRNAs with adjusted p-value < 0.05 and absolute log2 fold change > 0.05
sig_mir_6 <- res_6[res_6$padj < 0.05 & abs(res_6$log2FoldChange) > 0.05, ]

# Export significant miRNAs to CSV
write.csv(sig_mir_6, file = "deseq_significant_DEGs_X_PANCAN_DNAMethyl_LAML.csv", quote = FALSE, row.names = TRUE)

# Plot volcano plot of differential expression results
par(mfrow = c(1, 1))
plot(
    res_6$log2FoldChange, 
    -log10(res_6$padj), 
    pch = 20, 
    col = ifelse(res_6$padj < 0.05 & res_6$log2FoldChange > 1.2, "blue", 
                 ifelse(res_6$padj < 0.05 & res_6$log2FoldChange < -1.2, "red", "grey")),
    main = "Volcano Plot: Cluster 1 vs Cluster 2 Significant miRNAs",
    xlab = "Log2 Fold Change",
    ylab = "-Log10 Adjusted P-value",
    cex = 1.2
)
abline(h = -log10(0.05), col = "black", lty = 2)
abline(v = c(-1.2, 1.2), col = "black", lty = 2)

legend("topleft", 
       legend = c("Upregulated", "Downregulated", "Not Significant"), 
       col = c("blue", "red", "grey"), 
       pch = 20, 
       bty = "n",
       cex = 0.8)

# Heatmap for significant miRNAs
# Normalize counts using size factor estimation
dds_6_heatmap <- estimateSizeFactors(dds_6)
normalized_counts6 <- as.data.frame(counts(dds_6_heatmap, normalized = TRUE))

# Extract counts for significant miRNAs
DEmiRNA6 <- as.matrix(normalized_counts6[rownames(normalized_counts6) %in% rownames(sig_mir_6), ])

# Plot heatmap of log2 normalized counts for significant miRNAs
pdf("heatmap_DEmiRNAs_X_PANCAN_DNAMethyl_LAML.pdf", width = 20, height = 10)
aheatmap(
    log2(DEmiRNA6 + 1),
    color = "-RdYlBu2:100",
    scale = "none",
    Rowv = TRUE,
    Colv = TRUE,
    annCol = data.frame(DNAMethylationSubtype = cleaned_PANCAN_DNAMethyl_clin$X_PANCAN_DNAMethyl_LAML),
    main = "Heatmap of Significant miRNAs: Cluster 1 vs Cluster 2",
    fontsize = 16,
    cexRow = 1.2,
    cexCol = 1.0,
    labRow = rownames(DEmiRNA6)
)
dev.off()

#pca
vsd_6 <- varianceStabilizingTransformation(dds_6, blind = TRUE)
pca_model <- prcomp(t(assay(vsd_6)))
pca_scores <- as.data.frame(pca_model$x)
percentVar <- round(100 * summary(pca_model)$importance[2, 1:5], 2)
names(percentVar) <- paste0("PC", 1:5)
sample_info <- colData(vsd_6)
pca_scores$X_PANCAN_DNAMethyl_LAML <- sample_info$X_PANCAN_DNAMethyl_LAML
make_pca_plot <- function(df, x_pc, y_pc, percentVar, group_var = "X_PANCAN_DNAMethyl_LAML") {
    ggplot(df, aes_string(x = x_pc, y = y_pc, color = group_var)) +
        geom_point(size = 3) +
        labs(
            title = paste(x_pc, "vs", y_pc),
            x = paste0(x_pc, ": ", percentVar[x_pc], "% variance"),
            y = paste0(y_pc, ": ", percentVar[y_pc], "% variance"),
            color = group_var
        ) +
        theme_minimal() +
        theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))
}
p1 <- make_pca_plot(pca_scores, "PC1", "PC2", percentVar)
p2 <- make_pca_plot(pca_scores, "PC2", "PC3", percentVar)
p3 <- make_pca_plot(pca_scores, "PC3", "PC4", percentVar)
p4 <- make_pca_plot(pca_scores, "PC4", "PC5", percentVar)
combined_plot <- (p1 | p2) / (p3 | p4)
print(combined_plot)


##############for Cytogenetic Risk Category################

# Filter samples without NA in risk category
cleaned_risk_category_clin <- matched_surv_clin[!is.na(matched_surv_clin$acute_myeloid_leukemia_calgb_cytogenetics_risk_category), ]
cleaned_risk_category_mir <- matched_surv_mir[, rownames(cleaned_risk_category_clin)]

# Replace '/' in risk category labels and convert to factor
table(cleaned_risk_category_clin$acute_myeloid_leukemia_calgb_cytogenetics_risk_category)
cleaned_risk_category_clin$acute_myeloid_leukemia_calgb_cytogenetics_risk_category <- 
    gsub("/", "_", cleaned_risk_category_clin$acute_myeloid_leukemia_calgb_cytogenetics_risk_category)
cleaned_risk_category_clin$acute_myeloid_leukemia_calgb_cytogenetics_risk_category <- 
    as.factor(cleaned_risk_category_clin$acute_myeloid_leukemia_calgb_cytogenetics_risk_category)

# Create DESeqDataSet for cytogenetic risk groups
dds_risk_category <- DESeqDataSetFromMatrix(
    countData = cleaned_risk_category_mir,
    colData = cleaned_risk_category_clin,
    design = ~ acute_myeloid_leukemia_calgb_cytogenetics_risk_category
)

# Run DESeq2 differential expression analysis
dds_7 <- DESeq(dds_risk_category)
res_7 <- results(dds_7, contrast = c("acute_myeloid_leukemia_calgb_cytogenetics_risk_category","Favorable", "Intermediate_Normal"))

# Remove NA and convert to data frame
res_7 <- as.data.frame(res_7[complete.cases(res_7), ])

# Filter significant miRNAs with padj < 0.05 and abs log2 fold change > 0.05
sig_mir_7 <- res_7[res_7$padj < 0.05 & abs(res_7$log2FoldChange) > 0.05, ]

# Export significant miRNAs to CSV
write.csv(sig_mir_7, file = "deseq_significant_DEGs_risk.csv", quote = FALSE, row.names = TRUE)

# Volcano plot for cytogenetic risk comparison
par(mfrow = c(1, 1))
plot(
    res_7$log2FoldChange, 
    -log10(res_7$padj), 
    pch = 20, 
    col = ifelse(res_7$padj < 0.05 & res_7$log2FoldChange > 1.2, "blue", 
                 ifelse(res_7$padj < 0.05 & res_7$log2FoldChange < -1.2, "red", "grey")),
    main = "Volcano Plot: Favorable vs Intermediate/Normal Risk miRNAs",
    xlab = "Log2 Fold Change",
    ylab = "-Log10 Adjusted P-value",
    cex = 1.2
)
abline(h = -log10(0.05), col = "black", lty = 2)
abline(v = c(-1.2, 1.2), col = "black", lty = 2)

legend("topright", 
       legend = c("Upregulated", "Downregulated", "Not Significant"), 
       col = c("blue", "red", "grey"), 
       pch = 20, 
       bty = "n",
       cex = 0.8)

# Heatmap for significant miRNAs in cytogenetic risk
dds_7_heatmap <- estimateSizeFactors(dds_7)
normalized_counts7 <- as.data.frame(counts(dds_7_heatmap, normalized = TRUE))
DEmiRNA7 <- as.matrix(normalized_counts7[rownames(normalized_counts7) %in% rownames(sig_mir_7), ])

pdf("heatmap_DEmiRNAs_risk_category.pdf", width = 20, height = 10)
aheatmap(
    log2(DEmiRNA7 + 1),
    color = "-RdYlBu2:100",
    scale = "none",
    Rowv = TRUE,
    Colv = TRUE,
    annCol = data.frame(CytogeneticRiskCategory = cleaned_risk_category_clin$acute_myeloid_leukemia_calgb_cytogenetics_risk_category),
    main = "Heatmap of Significant miRNAs: Favorable vs Intermediate/Normal Risk",
    fontsize = 16,
    cexRow = 1.2,
    cexCol = 1.0,
    labRow = rownames(DEmiRNA7)
)
dev.off()

#pca
vsd_7 <- varianceStabilizingTransformation(dds_7, blind = TRUE)
pca_model <- prcomp(t(assay(vsd_7)))
pca_scores <- as.data.frame(pca_model$x)
percentVar <- round(100 * summary(pca_model)$importance[2, 1:5], 2)
names(percentVar) <- paste0("PC", 1:5)
sample_info <- colData(vsd_7)
pca_scores$acute_myeloid_leukemia_calgb_cytogenetics_risk_category <- sample_info$acute_myeloid_leukemia_calgb_cytogenetics_risk_category
make_pca_plot <- function(df, x_pc, y_pc, percentVar, group_var = "acute_myeloid_leukemia_calgb_cytogenetics_risk_category") {
    ggplot(df, aes_string(x = x_pc, y = y_pc, color = group_var)) +
        geom_point(size = 3) +
        labs(
            title = paste(x_pc, "vs", y_pc),
            x = paste0(x_pc, ": ", percentVar[x_pc], "% variance"),
            y = paste0(y_pc, ": ", percentVar[y_pc], "% variance"),
            color = group_var
        ) +
        theme_minimal() +
        theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))
}
p1 <- make_pca_plot(pca_scores, "PC1", "PC2", percentVar)
p2 <- make_pca_plot(pca_scores, "PC2", "PC3", percentVar)
p3 <- make_pca_plot(pca_scores, "PC3", "PC4", percentVar)
p4 <- make_pca_plot(pca_scores, "PC4", "PC5", percentVar)
combined_plot <- (p1 | p2) / (p3 | p4)
print(combined_plot)



# VST normalization and export for DESeq2 results from different clinical groups
normalize_and_export <- function(dds, filename) {
    vst_result <- varianceStabilizingTransformation(dds, blind = FALSE)
    norm_data <- assay(vst_result)
    boxplot(norm_data, las = 2, outline = FALSE)
    write.csv(as.matrix(norm_data), file = filename, quote = FALSE, row.names = TRUE)
    return(norm_data)
}
# Normalize datasets and export normalized miRNA expression matrices
norm_gender <- normalize_and_export(dds_2, "Normalized_gender_DEmiRNAs.csv")
norm_vital <- normalize_and_export(dds_3, "Normalized_Vital_Status_DEmiRNAs.csv")
norm_morphology <- normalize_and_export(dds_4, "Normalized_Leukemia_Morphology_DEmiRNAs.csv")
norm_miRNA <- normalize_and_export(dds_5, "Normalized_PanCan_miRNA_DEmiRNAs.csv")
norm_dnAmethyl <- normalize_and_export(dds_6, "Normalized_DNAmethyl_LAML_DEmiRNAs.csv")


# Organize DEG (Differentially Expressed Genes) results and sort them by adjusted p-value
ordered_results <- lapply(list(sig_mir_1, sig_mir_2, sig_mir_3, sig_mir_4, sig_mir_5, sig_mir_6, sig_mir_7), 
                          function(res) res[order(res$padj), ])

# Find common and union miRNAs among results 4 to 6
common_456 <- Reduce(intersect, lapply(ordered_results[4:6], rownames))
union_456 <- Reduce(union, lapply(ordered_results[4:6], rownames))

# Extract normalized miRNA data for common miRNAs
miR_norm_morphology <- norm_morphology[common_456, ]
miR_norm_miRNA <- norm_miRNA[common_456, ]
miR_norm_dnAmethyl <- norm_dnAmethyl[common_456, ]

# Find common samples across all three datasets
shared_samples <- Reduce(intersect, list(
    colnames(miR_norm_morphology),
    colnames(miR_norm_miRNA),
    colnames(miR_norm_dnAmethyl)
))

# Transpose normalized data matrices to have samples as rows
miR_norm_morphology <- t(miR_norm_morphology[, shared_samples])
miR_norm_miRNA <- t(miR_norm_miRNA[, shared_samples])
miR_norm_dnAmethyl <- t(miR_norm_dnAmethyl[, shared_samples])
norm_vital <- t(norm_vital[, shared_samples])
norm_gender <- t(norm_gender[, shared_samples])

# Combine clinical metadata with survival data into one dataframe
combined_meta <- cbind(matched_surv_clin, matched_surv)

# Convert some clinical variables to numeric or factor types for analysis
combined_meta$Survival <- as.numeric(combined_meta$Survival)
combined_meta$Death <- as.numeric(combined_meta$Death)
combined_meta$leukemia_french_american_british_morphology_code <- factor(combined_meta$leukemia_french_american_british_morphology_code)
combined_meta$X_PANCAN_mirna_LAML <- factor(combined_meta$X_PANCAN_mirna_LAML)
combined_meta$X_PANCAN_DNAMethyl_LAML <- factor(combined_meta$X_PANCAN_DNAMethyl_LAML)

# Keep only samples present in both miRNA data and clinical metadata
final_samples <- intersect(rownames(miR_norm_morphology), rownames(combined_meta))
combined_meta <- combined_meta[final_samples, ]

# Merge clinical metadata with normalized miRNA expression data for each dataset
merged_morphology <- cbind(combined_meta, miR_norm_morphology)
merged_miRNA <- cbind(combined_meta, miR_norm_miRNA)
merged_dnAmethyl <- cbind(combined_meta, miR_norm_dnAmethyl)

# Define a common ggplot2 theme for consistent styling of survival plots
common_ggtheme <- theme_minimal() + theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent")
)

# Function to plot Kaplan-Meier survival curves with customized colors and styling
plot_km_curve <- function(fit_obj, data, title) {
    strata_names <- names(fit_obj$strata)
    n_strata <- length(strata_names)
    
    # Choose color palette for the survival groups
    palette_colors <- brewer.pal(min(n_strata, 8), "Dark2")
    if (n_strata > 8) {
        palette_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(n_strata)
    }
    
    # Generate the survival plot with risk table and p-value
    ggsurvplot(
        fit_obj, data = data,
        palette = palette_colors,
        pval = TRUE,
        conf.int = TRUE,
        risk.table = TRUE,
        risk.table.col = "strata",
        risk.table.height = 0.25,
        surv.median.line = "none",
        xlab = "Time (days)",
        ylab = "Survival Probability",
        title = title,
        legend.title = "Group",
        legend.labs = strata_names,
        fontsize = 4,
        font.title = c(16, "bold"),
        font.x = c(14, "plain"),
        font.y = c(14, "plain"),
        font.tickslab = c(13, "plain"),
        risk.table.fontsize = 3.5,
        ggtheme = common_ggtheme,
        na.rm = TRUE
    )
}

# Remove rows with missing Survival or Death values before survival analysis
merged_morphology_clean <- merged_morphology[!is.na(merged_morphology$Survival) & !is.na(merged_morphology$Death), ]
merged_miRNA_clean <- merged_miRNA[!is.na(merged_miRNA$Survival) & !is.na(merged_miRNA$Death), ]
merged_dnAmethyl_clean <- merged_dnAmethyl[!is.na(merged_dnAmethyl$Survival) & !is.na(merged_dnAmethyl$Death), ]

# Fit Kaplan-Meier survival curve by morphology groups
fit_morphology <- survfit(Surv(Survival, Death) ~ leukemia_french_american_british_morphology_code, 
                          data = merged_morphology_clean)
print(plot_km_curve(fit_morphology, merged_morphology_clean, "Survival Based on Morphology Code"))

# Fit Kaplan-Meier survival curve by miRNA group
fit_miRNA <- survfit(Surv(Survival, Death) ~ X_PANCAN_mirna_LAML, 
                     data = merged_miRNA_clean)
print(plot_km_curve(fit_miRNA, merged_miRNA_clean, "Survival Based on Leukemia miRNA"))

# Fit Kaplan-Meier survival curve by DNAmethylation group
fit_dnAmethyl <- survfit(Surv(Survival, Death) ~ X_PANCAN_DNAMethyl_LAML, 
                         data = merged_dnAmethyl_clean)
print(plot_km_curve(fit_dnAmethyl, merged_dnAmethyl_clean, "Survival Based on Leukemia DNAmethylation"))


# Define miRNAs to analyze (from morphology normalized data)
top_mir_4 <- colnames(miR_norm_morphology)

# Step 1: Perform univariate Cox regression for each miRNA
cox_results_4 <- data.frame(Gene = character(), HR = numeric(), CI_lower = numeric(),
                            CI_upper = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

for (gene in top_mir_4) {
    # Create Cox model formula dynamically
    formula_cox <- as.formula(paste("Surv(Survival, Death) ~", paste0("`", gene, "`")))
    
    # Fit Cox proportional hazards model
    cox_model <- coxph(formula_cox, data = merged_morphology)
    summary_cox <- summary(cox_model)
    
    # Store results: Hazard Ratio, Confidence Interval, p-value
    cox_results_4 <- rbind(cox_results_4, data.frame(
        Gene = gene,
        HR = summary_cox$coefficients[1, "exp(coef)"],
        CI_lower = summary_cox$conf.int[1, "lower .95"],
        CI_upper = summary_cox$conf.int[1, "upper .95"],
        p_value = summary_cox$coefficients[1, "Pr(>|z|)"]
    ))
}

# Step 2: Adjust p-values for multiple testing (False Discovery Rate)
cox_results_4$Adjusted_P_Value <- p.adjust(cox_results_4$p_value, method = "fdr")

# Step 3: Select significant miRNAs with adjusted p-value < 0.05
sig_cox_mir_4 <- cox_results_4[cox_results_4$Adjusted_P_Value < 0.05, ]

# Stop if no significant miRNAs found
if (nrow(sig_cox_mir_4) == 0) {
    stop("No significant miRNAs found with adjusted p-value < 0.05")
}

# Choose top significant miRNA
selected_mir_4 <- sig_cox_mir_4$Gene[1]

# Step 4: Create High/Low expression groups based on median split
group_name <- paste0(selected_mir_4, "_group")
merged_morphology[[group_name]] <- ifelse(
    merged_morphology[[selected_mir_4]] > median(merged_morphology[[selected_mir_4]], na.rm = TRUE),
    "High", "Low"
)

# Step 5: Fit Kaplan-Meier curve by High/Low expression groups
km_fit_4 <- survfit(Surv(Survival, Death) ~ get(group_name), data = merged_morphology)

# Step 6: Plot Kaplan-Meier survival curve
ggsurvplot(km_fit_4, data = merged_morphology,
           risk.table = TRUE,
           pval = TRUE,
           title = paste("Kaplan-Meier Survival Curve for", selected_mir_4, "\n~ leukemia_french_american_british_morphology_code"),
           legend.labs = c("Low Expression", "High Expression"),
           palette = c("blue", "red")  # Customize line colors here
)

# Repeat the same Cox and KM analysis steps for miRNA normalized data (norm_miRNA)
top_mir_5 <- colnames(miR_norm_morphology)

cox_results_5 <- data.frame(Gene = character(), HR = numeric(), CI_lower = numeric(),
                            CI_upper = numeric(), p_value = numeric(), stringsAsFactors=FALSE)

for (gene in top_mir_5) {
    cox_model <- coxph(as.formula(paste("Surv(Survival, Death) ~", gene)), data = merged_miRNA)
    summary_cox <- summary(cox_model)
    
    cox_results_5 <- rbind(cox_results_5, data.frame(
        Gene = gene,
        HR = summary_cox$coefficients[1, "exp(coef)"],
        CI_lower = summary_cox$conf.int[1, "lower .95"],
        CI_upper = summary_cox$conf.int[1, "upper .95"],
        p_value = summary_cox$coefficients[1, "Pr(>|z|)"]
    ))
}

cox_results_5$Adjusted_P_Value <- p.adjust(cox_results_5$p_value, method = "fdr")
sig_cox_mir_5 <- cox_results_5[cox_results_5$Adjusted_P_Value < 0.05, ]

if (nrow(sig_cox_mir_5) == 0) {
    stop("No significant miRNAs found with adjusted p-value < 0.05 in miRNA dataset")
}

selected_mir_5 <- sig_cox_mir_5$Gene[1]
group_name_5 <- paste0(selected_mir_5, "_group")
merged_miRNA[[group_name_5]] <- ifelse(
    merged_miRNA[[selected_mir_5]] > median(merged_miRNA[[selected_mir_5]], na.rm = TRUE),
    "High", "Low"
)

km_fit_5 <- survfit(Surv(Survival, Death) ~ get(group_name_5), data = merged_miRNA)
ggsurvplot(km_fit_5, data = merged_miRNA,
           risk.table = TRUE,
           pval = TRUE,
           title = paste("Kaplan-Meier Survival Curve for", selected_mir_5, "\n~ X_PANCAN_mirna_LAML"),
           legend.labs = c("Low Expression", "High Expression"),
           palette = c("blue", "red")
)

# Repeat the same Cox and KM analysis for DNAmethylation normalized data
top_mir_6 <- colnames(miR_norm_dnAmethyl)

cox_results_6 <- data.frame(Gene = character(), HR = numeric(), CI_lower = numeric(),
                            CI_upper = numeric(), p_value = numeric(), stringsAsFactors=FALSE)

for (gene in top_mir_6) {
    cox_model <- coxph(as.formula(paste("Surv(Survival, Death) ~", gene)), data = merged_dnAmethyl)
    summary_cox <- summary(cox_model)
    
    cox_results_6 <- rbind(cox_results_6, data.frame(
        Gene = gene,
        HR = summary_cox$coefficients[1, "exp(coef)"],
        CI_lower = summary_cox$conf.int[1, "lower .95"],
        CI_upper = summary_cox$conf.int[1, "upper .95"],
        p_value = summary_cox$coefficients[1, "Pr(>|z|)"]
    ))
}

cox_results_6$Adjusted_P_Value <- p.adjust(cox_results_6$p_value, method = "fdr")
sig_cox_mir_6 <- cox_results_6[cox_results_6$Adjusted_P_Value < 0.05, ]

if (nrow(sig_cox_mir_6) == 0) {
    stop("No significant miRNAs found with adjusted p-value < 0.05 in DNAmethylation dataset")
}

selected_mir_6 <- sig_cox_mir_6$Gene[1]
group_name_6 <- paste0(selected_mir_6, "_group")
merged_dnAmethyl[[group_name_6]] <- ifelse(
    merged_dnAmethyl[[selected_mir_6]] > median(merged_dnAmethyl[[selected_mir_6]], na.rm = TRUE),
    "High", "Low"
)

km_fit_6 <- survfit(Surv(Survival, Death) ~ get(group_name_6), data = merged_dnAmethyl)
ggsurvplot(km_fit_6, data = merged_dnAmethyl,
           risk.table = TRUE,
           pval = TRUE,
           title = paste("Kaplan-Meier Survival Curve for", selected_mir_6, "\n~ X_PANCAN_DNAMethyl_LAML"),
           legend.labs = c("Low Expression", "High Expression"),
           palette = c("blue", "red")
)

