#Load library
library(knitr)
library(limma)
library(pheatmap)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(DMRcatedata)
library(lumi)
library(labdsv)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(magick)
BiocManager::install("magick")
##########################################

### Loading the Data ###
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
meta_data <- aml
beta_data <- methy
color_palette <- brewer.pal(8,"Set1")

##########################################

### Pre-processing the Data & Quick Check###
beta_data <- as.matrix(beta_data)
m_data <- beta2m(beta_data)
GR_Set <- makeGenomicRatioSetFromMatrix(beta_data)


# Identify Categorical Columns & Replace NA with "Unknown"
categorical_cols <- sapply(meta_data, function(col) is.character(col) | is.factor(col))

# Convert to Character to Handle Factors & Convert Back to Factor
meta_data[categorical_cols] <- lapply(aml[categorical_cols], function(col) {
  ifelse(is.na(col), "Unknown", as.character(col))
}) 
meta_data[categorical_cols] <- lapply(meta_data[categorical_cols], as.factor)


# Rename the SampleID Column to Match the Beta Values
meta_data$sampleID <- gsub("\\-", ".", (meta_data$sampleID))


# Give the Rows Names & Remove the SampleID Column
row.names(meta_data) <- meta_data$sampleID
meta_data <- meta_data[,-1]


# Check Samples Arrangement
all(rownames(meta_data) == colnames(beta_data)) # FALSE

# Get Intersection in Data
common_data <- intersect(colnames(m_data), rownames(meta_data))
# Subset Data
m_data <- m_data[, common_data, drop = FALSE]
meta_data <- meta_data[common_data, , drop = FALSE]

all(rownames(meta_data) == colnames(beta_data)) # TRUE

# Check for Duplicate Rows
sum(duplicated(meta_data))

# Point of Interest
table(meta_data$vital_status)

##########################################

# Get the age 
meta_data$age <- round(-meta_data$days_to_birth / 365.25)
hist(meta_data$age, main = "Age Distribution", xlab = "Age", col = "lightblue", breaks = 50)

meta_data$age_group <- cut(
  meta_data$age,
  breaks = c(-Inf, 18, 30, 45, 65, Inf),
  labels = c("Young", "Young Adult", "Adult", "Middle-Age", "Elder"),
  right = FALSE
)

table(meta_data$age_group)

# Plot the Age Distribution
ggplot(meta_data, aes(x = age_group, y = age, fill = age_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  labs(title = "Age Distribution by Manual Bins", x = "Cluster", y = "Age", fill = "Age Group") +
  theme_minimal() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

##########################################


### Visualizing ###
par(mfrow=c(1,1))
# 1. Beta Values Distribution
hist(beta_data, 
     main="Histogram of Beta Values", 
     xlab="Beta Value", 
     ylab="Frequency", 
     col="lightblue", 
     breaks=100)

# 2. M Values Distribution
hist(m_data, 
     main="Histogram of M Values", 
     xlab="M Value", 
     ylab="Frequency", 
     col="lightgreen", 
     breaks=100)

# 3. Density Plot of Beta Values
densityPlot(beta_data, main="Beta values",
            sampGroups=meta_data$vital_status,
            legend=FALSE, xlab="Beta values")

# 4. Density Plot of M Values
densityPlot(m_data, main="M values",
            sampGroups=meta_data$vital_status,
            legend=FALSE, xlab="M values")

# 5. Density Plot of Beta Values
pheatmap(beta_matrix, scale = "row", show_rownames = FALSE, show_colnames = FALSE)

# Use ggplot2 for Density Plot
# Convert to data frame
beta_data_long <- as.data.frame(beta_data)
beta_data_long$CpG <- rownames(beta_data_long) # Add row names as a column
# Melt into long format
beta_data_long <- melt(beta_data_long, id.vars = "CpG", variable.name = "Sample", value.name = "BetaValue")
# Plot
ggplot(beta_data_long, aes(x = BetaValue, fill = Sample)) +
  geom_density(alpha = 0.5) +
  labs(title = "Beta Value Distribution", x = "Beta Value", y = "Density") +
  theme_minimal() +
  theme(legend.position = "none")



##########################################

### Data exploration ###

# MDS Plots
par(mfrow=c(1,2))

plotMDS(m_data, top=10000, pch=16,
        gene.selection="common",
        col=color_palette[meta_data$age_group], main="MDS plot of M-values")

legend("topleft", legend=levels(meta_data$age_group), 
       bg="white", pch=16, cex=0.6, ncol=2, inset=0.02, 
       text.col=color_palette, col = color_palette)


# Examine higher dimensions to look at other sources of variation
par(mfrow=c(2,3))

#make list() features 
# Updated list of top features for clustering
top_features_list <- colnames(meta_data)

# Define a color palette for visualization
color_palette <- rainbow(length(unique(meta_data$vital_status)))

# Loop through each feature and generate MDS plots
for (feature in top_features_list) {
  # Ensure the feature exists in meta_data
  if (feature %in% colnames(meta_data)) {
    # Plot MDS for different dimensions
    plotMDS(m_data, top=1000, pch = 16,
            gene.selection="common",
            col=color_palette[as.numeric(meta_data[[feature]])], dim=c(1,2))
    title(main = paste("MDS Plot - Dim 1 vs 2 (Colored by", feature, ")"))
    
    plotMDS(m_data, top=1000, pch = 16,
            gene.selection="common",
            col=color_palette[as.numeric(meta_data[[feature]])], dim=c(2,3))
    title(main = paste("MDS Plot - Dim 2 vs 3 (Colored by", feature, ")"))
    
    plotMDS(m_data, top=1000, pch = 16,
            gene.selection="common",
            col=color_palette[as.numeric(meta_data[[feature]])], dim=c(3,4))
    title(main = paste("MDS Plot - Dim 3 vs 4 (Colored by", feature, ")"))
    
    plotMDS(m_data, top=1000, pch = 16,
            gene.selection="common",
            col=color_palette[as.numeric(meta_data[[feature]])], dim=c(4,1))
    title(main = paste("MDS Plot - Dim 4 vs 1 (Colored by", feature, ")"))
  } else {
    warning(paste("Feature", feature, "not found in meta_data. Skipping."))
  }
}



plotMDS(m_data, top=1000, pch = 16,
        gene.selection="common",
        col=color_palette[meta_data$gender], dim=c(1,2))
plotMDS(m_data, top=1000, pch = 16,
        gene.selection="common",
        col=color_palette[meta_data$gender], dim=c(2,3))
plotMDS(m_data, top=1000, pch = 16,
        gene.selection="common",
        col=color_palette[meta_data$gender], dim=c(3,4))
plotMDS(m_data, top=1000, pch = 16,
        gene.selection="common",
        col=color_palette[meta_data$gender], dim=c(4,1))
plotMDS(m_data, top=1000, pch = 16,
        gene.selection="common",
        col=color_palette[meta_data$gender], dim=c(4,2))
plotMDS(m_data, top=1000, pch = 16,
        gene.selection="common",
        col=color_palette[meta_data$gender], dim=c(3,1))

legend("topleft", legend=levels(meta_data$gender), 
       bg="white", pch=16, cex=0.6, ncol=2, 
       text.col=color_palette, col = color_palette)

meta_data$gender



##########################################
# Filtration

xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)


# Exclude Cross Reactive Probes
keep <- !(featureNames(GR_Set) %in% xReactiveProbes$TargetID)
GR_Set_Filtered <- GR_Set[keep,]

# Exclude Probes on X and Y Chromosomes
keep <- !(featureNames(GR_Set_Filtered) %in% ann450k$Name[ann450k$chr %in% 
                                                  c("chrX","chrY")])
GR_Set_Filtered <- GR_Set_Filtered[keep,]

# Remove Probes with SNPs at CpG site
GR_Set_Filtered <- dropLociWithSnps(GR_Set)

# Visualize the Filtration
par(mfrow=c(1,2))
densityPlot(getBeta(GR_Set_Filtered), main="Beta values",
            sampGroups=meta_data$vital_status,
            legend=FALSE, xlab="Beta values")

# Density Plot of M Values
densityPlot(getM(GR_Set_Filtered), main="M values",
            sampGroups=meta_data$vital_status,
            legend=FALSE, xlab="M values")

##########################################
# Linear Models Differential Methylation Analysis


# Factor of Interest
vital_status <- factor(meta_data$vital_status)
age_group <- factor(meta_data$age_group)
gender <- factor(meta_data$gender)

# Individual Effect
individual <- factor(row.names(meta_data))

# Design Matrix
design <- model.matrix(~0+vital_status, data=meta_data)
colnames(design) <- c(levels(vital_status))

# Fit The Linear Model
fit <- lmFit(getM(GR_Set_Filtered), design)

# Contrast Matrix For Specific Comparisons
contMatrix <- makeContrasts(DECEASED-LIVING,
                            levels=design)
contMatrix

# Fit The Contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)


# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))


# Get Table of Results
ann450kSub <- ann450k[match(rownames(getM(GR_Set_Filtered)),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
head(DMPs)

write.table(DMPs, file="DMPs.csv", sep=",", row.names=FALSE)

# Add a column for -log10(p-value) for better visualization
DMPs$negLogPValue <- -log10(DMPs$P.Value)

# Create the volcano plot
volcano_plot <- ggplot(DMPs, aes(x = logFC, y = negLogPValue)) +
  geom_point(aes(color = ifelse(P.Value < 0.05 & abs(logFC) > 0.5, "Significant", "Not Significant")), alpha = 0.6) +
  scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 P-Value",
    title = "Volcano Plot",
    color = "Significance"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") + # FDR cutoff line
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") # LogFC cutoff lines

# Display the plot
print(volcano_plot)


# Plot Top 4 Most Significantly Differentially Methylated CpGs
par(mfrow=c(2,3))
sapply(rownames(DMPs)[1:6], function(cpg){
  plotCpg(getBeta(GR_Set_Filtered), cpg=cpg, pheno=meta_data$vital_status, ylab = "Beta Values")
})


################################################
# Differential Methylation Analysis of Regions
myAnnotation <- cpg.annotate(object = getM(GR_Set_Filtered), datatype = "array", what = "M",
                             analysis.type = "differential", design = design,
                             contrasts = TRUE, cont.matrix = contMatrix,
                             coef = "DECEASED - LIVING", arraytype = "450K")

DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges

# Setup Grouping Variables & Colors

groups <- color_palette[1:length(unique(meta_data$vital_status))]
names(groups) <- levels(factor(meta_data$vital_status))
cols <- groups[as.character(factor(meta_data$vital_status))]

# draw the plot for the top DMR

par(mfrow=c(1,1))


# Plot the top DMR

DMR.plot(ranges = results.ranges, 
         dmr = 1,  # Plot the first DMR
         CpGs = getBeta(GR_Set_Filtered[, 1:20]), 
         phen.col = cols[1:20],
         what = "Beta", 
         arraytype = "450K", 
         genome = "hg19"
)

# make it better 
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges, 
         dmr = 1, 
         CpGs = getBeta(GR_Set_Filtered[, 1:5]), 
         phen.col = cols[1:5],
         what = "Beta",
         arraytype = "450K", 
         genome = "hg19", 
)



# Save to a larger PNG file
png("dmr_plot.png", width = 3000, height = 2000, res = 300)  # Adjust width/height as needed
DMR.plot(ranges = results.ranges, 
         dmr = 1, 
         CpGs = getBeta(GR_Set_Filtered), 
         phen.col = cols,
         what = "Beta", 
         arraytype = "450K", 
         genome = "hg19"
)
dev.off()  # Close the device


