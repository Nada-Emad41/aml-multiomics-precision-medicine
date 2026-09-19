# -------------------------------------------------
# Load Required Libraries
# -------------------------------------------------
library(MOFA2)
library(basilisk)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(miRBaseConverter)
library(readr)
library(psych)

# -------------------------------------------------
# Define file paths (modify these to match your system)
# -------------------------------------------------
base_path <- "G:/integrative_course/project_AML"
aml_file <- file.path(base_path, "miRNA_DATA/aml")
methylation_file <- file.path(base_path, "methyl.csv")
expression_file <- file.path(base_path, "exp")
mirna_file <- file.path(base_path, "mirna")

# -------------------------------------------------
# Read Data Files with Error Handling
# -------------------------------------------------

# Check if files exist before reading
if (!file.exists(aml_file)) {
    stop("AML metadata file not found at: ", aml_file)
}
if (!file.exists(methylation_file)) {
    stop("Methylation data file not found at: ", methylation_file)
}
if (!file.exists(expression_file)) {
    stop("Expression data file not found at: ", expression_file)
}
if (!file.exists(mirna_file)) {
    stop("miRNA data file not found at: ", mirna_file)
}

# Read the data files
cat("Reading data files...\n")
aml <- read.delim(aml_file, row.names = 1)
methylation_data <- read.csv(methylation_file, header = TRUE, row.names = 1)
expression_data <- as.matrix(read.table(expression_file, quote="\"", comment.char=""))
mirna_data <- as.matrix(read.table(mirna_file, quote="\"", comment.char=""))

# Display data dimensions for quality check
cat("Data dimensions after loading:\n")
cat("AML metadata:", dim(aml), "\n")
cat("Expression data:", dim(expression_data), "\n")
cat("Methylation data:", dim(methylation_data), "\n")
cat("miRNA data:", dim(mirna_data), "\n")

# -------------------------------------------------
# Normalize Expression and miRNA Data (log2 + 1)
# -------------------------------------------------
cat("Applying log2(x+1) normalization to expression and miRNA data...\n")

# Check for negative values before log transformation
if (any(expression_data < 0, na.rm = TRUE)) {
    warning("Negative values found in expression data. Consider reviewing data preprocessing.")
}
if (any(mirna_data < 0, na.rm = TRUE)) {
    warning("Negative values found in miRNA data. Consider reviewing data preprocessing.")
}

expression_data <- log2(expression_data + 1)
mirna_data <- log2(mirna_data + 1)

cat("Normalization completed.\n")

# -------------------------------------------------
# Harmonize Omics Data Matrices
# -------------------------------------------------
cat("Harmonizing omics data matrices...\n")

omics_data_list <- list(methylation_data, expression_data, mirna_data)
sample_ids <- Reduce(union, lapply(omics_data_list, colnames))

cat("Total unique samples across all datasets:", length(sample_ids), "\n")

omics_data_list <- lapply(omics_data_list, function(mat) {
    mat <- as.matrix(mat)
    missing_cols <- setdiff(sample_ids, colnames(mat))
    if (length(missing_cols) > 0) {
        cat("Adding", length(missing_cols), "missing columns with NA values\n")
        na_mat <- matrix(NA, nrow = nrow(mat), ncol = length(missing_cols))
        colnames(na_mat) <- missing_cols
        rownames(na_mat) <- rownames(mat)
        mat <- cbind(mat, na_mat)
    }
    mat <- mat[, sample_ids, drop = FALSE]
    return(mat)
})

methylation_data <- omics_data_list[[1]]
expression_data <- omics_data_list[[2]]
mirna_data <- omics_data_list[[3]]

# Display final dimensions
cat("Final harmonized data dimensions:\n")
cat("Expression data:", dim(expression_data), "\n")
cat("Methylation data:", dim(methylation_data), "\n")
cat("miRNA data:", dim(mirna_data), "\n")

# -------------------------------------------------
# Process Metadata
# -------------------------------------------------
cat("Processing metadata...\n")

samples <- gsub("-", ".", rownames(aml))
diff_meta <- setdiff(samples, sample_ids)

if (length(diff_meta) > 0) {
    cat("Removing", length(diff_meta), "metadata entries not present in omics data\n")
}

samples_data <- aml[!(samples %in% diff_meta), ]
rownames(samples_data) <- gsub("-", ".", rownames(samples_data))
samples_data$sample <- rownames(samples_data)

cat("Final metadata dimensions:", dim(samples_data), "\n")

# -------------------------------------------------
# Prepare the MOFA Data Object
# -------------------------------------------------
cat("Preparing MOFA data object...\n")

omics_data <- list(
    "Gene Expression" = data.matrix(expression_data),
    "Methylation"     = data.matrix(methylation_data),
    "miRNA"           = data.matrix(mirna_data)
)

# Data quality checks
cat("Data quality summary:\n")
for (view_name in names(omics_data)) {
    view_data <- omics_data[[view_name]]
    cat(view_name, "- Features:", nrow(view_data), "Samples:", ncol(view_data), "\n")
    cat("  Missing values:", sum(is.na(view_data)), "\n")
    cat("  Data range:", round(range(view_data, na.rm = TRUE), 3), "\n")
}

# Optional: Check histograms
cat("Generating data distribution plots...\n")
par(mfrow = c(1, 3))
hist(expression_data, main = "Gene Expression", xlab = "Log2(Expression + 1)")
hist(methylation_data, main = "Methylation", xlab = "Beta Value")
hist(mirna_data, main = "miRNA", xlab = "Log2(miRNA + 1)")

# Create MOFA object
Omics_MOFA_Obj <- create_mofa(omics_data)

# Assign metadata BEFORE training (critical step)
samples_metadata(Omics_MOFA_Obj) <- samples_data

# Overview of training data
cat("Plotting data overview...\n")
plot_data_overview(Omics_MOFA_Obj)

# -------------------------------------------------
# Define MOFA Training Options
# -------------------------------------------------
cat("Setting up MOFA training options...\n")

# Data Options
data_opts <- get_default_data_options(Omics_MOFA_Obj)
data_opts$scale_views <- TRUE
data_opts$use_float32 <- TRUE

# Model Options
model_opts <- get_default_model_options(Omics_MOFA_Obj)
model_opts$num_factors <- 15  # Conservative choice for sample size

# Training Options
train_opts <- get_default_training_options(Omics_MOFA_Obj)
train_opts$verbose <- TRUE
train_opts$seed <- 42
train_opts$gpu_mode <- TRUE
train_opts$convergence_mode <- "slow"
train_opts$maxiter <- 5000

# Prepare the model
cat("Preparing MOFA model...\n")
Omics_MOFA_Obj <- prepare_mofa(
    object = Omics_MOFA_Obj,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
)

# -------------------------------------------------
# Run the MOFA Model
# -------------------------------------------------
cat("Starting MOFA training (this may take 10-20 minutes)...\n")
start_time <- Sys.time()

Omics_MOFA_Obj <- run_mofa(Omics_MOFA_Obj, use_basilisk = TRUE)

end_time <- Sys.time()
cat("MOFA training completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

# Reset plotting layout
par(mfrow = c(1, 1))

# -------------------------------------------------
# Post-Training Analysis
# -------------------------------------------------
cat("Performing post-training analysis...\n")

# Plot correlation between factors
plot_factor_cor(Omics_MOFA_Obj)

# Calculate and display variance explained by each factor
r2 <- get_variance_explained(Omics_MOFA_Obj)
cat("Total variance explained by all factors:\n")
print(r2$r2_total)

cat("Variance explained per factor (first 10):\n")
print(head(r2$r2_per_factor, 10))

# Plot the variance explained
plot_variance_explained(Omics_MOFA_Obj)

# -------------------------------------------------
# Association Testing Between Factors and Covariates
# -------------------------------------------------
cat("Testing associations with clinical covariates...\n")

# Check if covariates exist in metadata
available_covariates <- c("vital_status", "gender", "X_PANCAN_mirna_LAML", "age_at_initial_pathologic_diagnosis")
existing_covariates <- available_covariates[available_covariates %in% colnames(samples_data)]

if (length(existing_covariates) > 0) {
    correlate_factors_with_covariates(
        Omics_MOFA_Obj, 
        covariates = existing_covariates, 
        plot = "log_pval", 
        return_data = FALSE
    )
} else {
    cat("Warning: No specified covariates found in metadata\n")
    cat("Available columns:", colnames(samples_data), "\n")
}

# -------------------------------------------------
# Factor Visualization
# -------------------------------------------------
cat("Creating factor visualizations...\n")

# Get actual number of factors
n_factors <- Omics_MOFA_Obj@dimensions$K
cat("Number of factors learned:", n_factors, "\n")

# Use actual number of factors for plotting
plot_factors <- min(n_factors, 15)  # Plot up to 15 factors

if ("vital_status" %in% colnames(samples_data)) {
    plot_factor(
        Omics_MOFA_Obj,
        factors = 1:plot_factors,
        color_by = "vital_status"
    )
}

if ("gender" %in% colnames(samples_data)) {
    plot_factor(
        Omics_MOFA_Obj,
        factors = 1:plot_factors,
        color_by = "gender"
    )
}

if ("X_PANCAN_mirna_LAML" %in% colnames(samples_data)) {
    plot_factor(
        Omics_MOFA_Obj,
        factors = 1:plot_factors,
        color_by = "X_PANCAN_mirna_LAML"
    )
}
# Detailed visualization for top factors
if (n_factors >= 8 && "gender" %in% colnames(samples_data)) {
    plot_factor(
        Omics_MOFA_Obj, 
        factors = 8, 
        color_by = "gender",
        add_violin = TRUE,
        dodge = TRUE
    )
}

if (n_factors >= 1 && "gender" %in% colnames(samples_data)) {
    plot_factor(
        Omics_MOFA_Obj, 
        factors = 1, 
        color_by = "gender",
        dodge = TRUE,
        add_violin = TRUE
    )
}

# -------------------------------------------------
# Feature Weight Visualization
# -------------------------------------------------
cat("Creating feature weight visualizations...\n")

# Heatmap of weights across factors for "Gene Expression" view
plot_weights_heatmap(
    Omics_MOFA_Obj, 
    view = "Gene Expression", 
    factors = 1:plot_factors,
    show_colnames = FALSE
)

# Heatmap of weights across factors for "Methylation" view
plot_weights_heatmap(
    Omics_MOFA_Obj, 
    view = "Methylation", 
    factors = 1:plot_factors,
    show_colnames = FALSE
)

# Heatmap of weights across factors for "Methylation" view
plot_weights_heatmap(
    Omics_MOFA_Obj, 
    view = "miRNA", 
    factors = 1:plot_factors,
    show_colnames = FALSE
)
# Plot top features for specific factors
if (n_factors >= 8) {
    plot_top_weights(
        Omics_MOFA_Obj, 
        view = "miRNA", 
        factor = 8, 
        nfeatures = 35
    )
}
# Plot top features for specific factors
if (n_factors >= 8) {
        plot_top_weights(
            Omics_MOFA_Obj, 
            view = "Gene Expression", 
            factor = 8, 
            nfeatures = 35
        )
    }
# Plot top features for specific factors
if (n_factors >= 8) {
    plot_top_weights(
        Omics_MOFA_Obj, 
        view = "Methylation", 
        factor = 8, 
        nfeatures = 35
    )
}

# -------------------------------------------------
# Sample Ordination and Data Heatmaps
# -------------------------------------------------
cat("Creating sample ordination plots...\n")

# Scatter plot ordination
if (n_factors >= 8 && "gender" %in% colnames(samples_data)) {
    plot_data_scatter(
        Omics_MOFA_Obj,
        factor = 8,
        color_by = "gender"
    )
}

# Gene expression scatter plot
if ("gender" %in% colnames(samples_data)) {
    plot_data_scatter(
        Omics_MOFA_Obj, 
        view = "Gene Expression",
        factor = 1,  
        features = 15,
        color_by = "gender"
    ) + labs(y = "RNA expression")
}

# Methylation heatmap
if ("gender" %in% colnames(samples_data)) {
    plot_data_heatmap(
        Omics_MOFA_Obj, 
        view = "Methylation", 
        factor = 1, 
        features = 19, 
        show_rownames = TRUE,
        denoise = TRUE, 
        annotation_samples = "gender"
    )
}
