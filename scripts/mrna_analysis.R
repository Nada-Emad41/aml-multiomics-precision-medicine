#==============================#
#        Load Libraries        #
#==============================#
library(BiocManager)
library(DESeq2)
library(NMF)
library(ggplot2)
library(limma)
library(RColorBrewer)
#==============================================================#
#   Handling RNA Expression Matrix and Clinical Metadata     #
#==============================================================#

# Load RNA expression data as a matrix from file (adjust path and separator if needed)
RNA_data <- as.matrix(read.table("D:/integrative_course/project AML/RNA_seq/exp", quote="\"", comment.char=""))

# Display the first few rows and columns
head(RNA_data)
head(RNA_data[, 1:5])

# Check dimensions of the RNA expression matrix
dim(RNA_data)  # Should return 705 188

# Load clinical metadata (assuming sampleID is in the second column and no row names in file)
clinical_data <- read.delim("D:/integrative_course/project AML/RNA_seq/aml", row.names=1)

# Display the first few rows
head(clinical_data)

# Check dimensions of the clinical data
dim(clinical_data)  # Should return 200 97

# Extract sample IDs from clinical data
clinical_samples <- rownames(clinical_data)
clinical_samples[1]  # Example: "TCGA-AB-2802-03"

# Extract sample names from RNA data
RNA_samples <- colnames(RNA_data)
RNA_samples[1]  # Example: "TCGA.AB.2802.03"

# Fix RNA sample names to match clinical format (replace '.' with '-')
RNA_samples_fixed <- gsub("\\.", "-", RNA_samples)
colnames(RNA_data) <- RNA_samples_fixed

# Display sample names after replacement
head(RNA_samples_fixed)
head(clinical_samples)

# Identify common samples between clinical and RNA datasets
common_samples <- intersect(RNA_samples_fixed, clinical_samples)
cat("Number of common samples:", length(common_samples), "\n")  # Expected: 188

# Subset RNA expression data to only common samples
RNA_expr_data_filtered <- RNA_data[, common_samples]

# Subset clinical metadata to only common samples
clinical_data_filtered <- clinical_data[rownames(clinical_data) %in% common_samples, ]

# Check dimensions of the filtered datasets
dim(RNA_expr_data_filtered)  # Should return 20531   173
dim(clinical_data_filtered)    # Should return 173  97

# Verify if sample IDs perfectly match between datasets
all(colnames(RNA_expr_data_filtered) == clinical_data_filtered$sampleID)


#==============================================================#
#      Handling Survival Data with RNA & Clinical Metadata   #
#==============================================================#

# Load survival data
survival_data <- read.delim("D:/integrative_course/project AML/RNA_seq/survival")

# Display first few rows
head(survival_data)

# Check dimensions of survival data
dim(survival_data)  # Should return 363   2

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


# Filter Out NA in Death (or other variable of interest)

survival_data <- survival_data[!is.na(survival_data$Death), ]


###Factorizing Death variable in survival data

survival_data$Death <- factor(survival_data$Death, levels = c(0, 1),labels = c("Alive", "Dead"))
str(survival_data)
table(survival_data$Death)#Alive 130 & Dead   233


#==============================================================#
#      Match Survival Data with Filtered Clinical & RNA      #
#==============================================================#

# Identify common sample IDs across all datasets
clicn_suv_datamatch <- intersect(rownames(clinical_data_filtered), rownames(survival_data))
length(clicn_suv_datamatch)

# Subset datasets to include only matched samples
matched_surv <- survival_data[clicn_suv_datamatch, ]
matched_surv_clin <- clinical_data_filtered[clicn_suv_datamatch, ] # match survival with clinical data 
matched_surv_rna <- RNA_expr_data_filtered[, clicn_suv_datamatch]# match survival with RNA data 

# Check dimensions of the matched datasets
dim(matched_surv) # 162   2 
dim(matched_surv_clin) # 162  97
dim(matched_surv_rna)  # 20531   162



# Confirm sample ID consistency across all matched datasets
rownames(matched_surv)[1]
rownames(matched_surv_clin[1, ])
colnames(matched_surv_rna)[1]



# Save RNA feature names for downstream analysis
rna_names <- rownames(matched_surv_rna)


sum(is.na(matched_surv_rna))
sum(is.na(matched_surv_clin))
sum(is.na(matched_surv))
colSums(is.na(matched_surv_clin))
colSums(is.na(matched_surv))
dim(matched_surv_rna)
##### Visualization of RNA expression data in AML #####


#explore the data distribution using the histogram plot
hist(matched_surv_rna, col = "orange", main="Histogram")

#scaling the data using log2 transformation to better visulization
# we use (+1) to avoid the infinity character when we log zero valus 
hist(log2(matched_surv_rna+1), col = "orange", main="Histogram")
####boxplot for the first four genes

boxplot(log2(matched_surv_rna[1:4,]+1))


# QQ plot for the normality
qqnorm(matched_surv_rna[3,])
qqline(matched_surv_rna[3,])
####Q-Q plot for the first four genes

par(mfrow = c(2, 2))  # Layout for 4 plots

for (i in 1:4) {
  qqnorm(matched_surv_rna[i, ], main = paste("Q-Q plot Gene", i))
  qqline(matched_surv_rna[i, ])
}


#####Shapiro test for the first four genes
shapiro.test(matched_surv_rna[1,])
shapiro.test(matched_surv_rna[2,])
shapiro.test(matched_surv_rna[3,])
shapiro.test(matched_surv_rna[4,])

#shapiro.test(matched_surv_rna[100,])

#explore if is there any missing expression value (empty cell)
sum(is.na(matched_surv_rna))
is.na(matched_surv_rna)
is.null(matched_surv_rna)
sum(is.nan(matched_surv_rna))





# Check if column names of expression matrix match pheno (survival) data
pheno=matched_surv[colnames(matched_surv_rna),]


#The deseq2 package require the count data values to be integers 
#save the gene names in a variable
genes=row.names(matched_surv_rna)
#convert the data values to integers
matched_surv_rna=apply(round(matched_surv_rna),2,as.integer)

head(matched_surv_rna)
#rename the rows of the data
row.names(matched_surv_rna)=genes




################PCA###########
#specify how many conditions do you want to compare according to 
#the phenotypic table
cond1="Dead" 
cond2="Alive"




#creat a deseq dataset object
dds= DESeqDataSetFromMatrix( countData = matched_surv_rna , colData = pheno, design = ~ Death)

###Variance stabilizing transformation
vsd <- vst(dds, blind=TRUE)
vsd_mat <- assay(vsd)  # Extract the matrix for plotting



# MDS plots for major sources of variation

# First plot: by gender
group_gender <- factor(matched_surv_clin$gender)
pal_gender <- brewer.pal(n = max(3, length(levels(group_gender))), "Set1")[1:length(levels(group_gender))]

par(mfrow=c(1,2))
plotMDS(vsd_mat, top=1000, gene.selection="common",
        col=pal_gender[group_gender])
legend("top", legend=levels(group_gender), text.col=pal_gender,
       bg="white", cex=0.7)

# Second plot: by survival (Death)
group_death <- factor(matched_surv$Death)
pal_death <- brewer.pal(n = max(3, length(levels(group_death))), "Dark2")[1:length(levels(group_death))]

plotMDS(vsd_mat, top=1000, gene.selection="common",
        col=pal_death[group_death])
legend("top", legend=levels(group_death), text.col=pal_death,
       bg="white", cex=0.7)



# Higher dimensions to look at hidden variation

par(mfrow=c(1,3))
plotMDS(vsd_mat, top=1000, gene.selection="common",
        col=pal_death[group_death], dim=c(1,3))
legend("top", legend=levels(group_death), text.col=pal_death,
       cex=0.7, bg="white")


# Get all unique dim combinations from 1 to 5 (without repeats)
dim_combos <- combn(1:5, 2, simplify = FALSE)

# MDS plots for gender
par(mfrow=c(2,2))  # adjust as needed based on output size
for (dims in dim_combos) {
  plotMDS(vsd_mat, top=1000, gene.selection="common",
          col=pal_gender[group_gender], dim=dims,
          main=paste0("Gender: Dim ", dims[1], " vs ", dims[2]))
  legend("topright", legend=levels(group_gender), text.col=pal_gender,
         cex=0.7, bg="white")
}

# MDS plots for death
par(mfrow=c(2,2))  # reset plotting window
for (dims in dim_combos) {
  plotMDS(vsd_mat, top=1000, gene.selection="common",
          col=pal_death[group_death], dim=dims,
          main=paste0("Death: Dim ", dims[1], " vs ", dims[2]))
  legend("topright", legend=levels(group_death), text.col=pal_death,
         cex=0.7, bg="white")
}





###### DO the differential EXP analysis using DeSeq2

#specify how many conditions do you want to compare according to 
#the phenotypic table
cond1="Dead" 
cond2="Alive"




#creat a deseq dataset object
dds= DESeqDataSetFromMatrix( countData = matched_surv_rna , colData = pheno, design = ~ Death)


#run the deseq2 worflow
dds.run = DESeq(dds)
#specifying teh contrast (to make a res object based on two specific conditions)
res=results(dds.run, contrast = c("Death",cond1 ,cond2))

# remove nulls
res=as.data.frame(res[complete.cases(res), ])

#chose the statistical significant differentialy expressed genes (DEGs) based
#on the p adjusted value less than 0.05 and biological significance  based
#on the fold change more than 2
deseq.deg0=res[res$padj < 0.05 & abs(res$log2FoldChange)>0.5,]#####project reporting requirement
deseq.deg1=res[res$padj < 0.05 & abs(res$log2FoldChange)>1.2,]
deseq.deg2=res[res$padj < 0.01 & abs(res$log2FoldChange)>1.2,]
deseq.deg=res[res$padj < 0.01 & abs(res$log2FoldChange)>2,]

####N.B:
####deseg.degs0 were 431 with(p<0.05&log2fold change of >0.5)
####deseg.degs1 were 247 with(p<0.05&log2fold change of >1.2)
####deseg.degs2 were 153 with(p<0.05&log2fold change of >2)
####deseg.degs were 59 with(p<0.01&log2fold change of >2)


#export the Degs into your current folder for further analysthis
write.csv(as.matrix(deseq.deg0),file="deseq.deg.csv", quote=F,row.names=T)


write.csv(as.matrix(deseq.deg),file="deseq.deg_filtered.csv", quote=F,row.names=T)

#drow DEGs volcano plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Alive vs Dead DEGs"))
with(subset(res, padj<.01 & (log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<.01 & (log2FoldChange)< -2), points(log2FoldChange, -log10(padj), pch=20, col="red"))
legend(x=-6.2,y=11.2,c("upregulated","downgulated"), cex=.8, bty="n", col=c("blue","red"),pch=19)


####draw heatmap####
#normalize the data
dds2 <- estimateSizeFactors(dds)
normalized_counts <- as.data.frame(counts(dds2, normalized=TRUE))

#extract counts values of DEGs only for each stage
exp.degs=as.matrix(normalized_counts[rownames(normalized_counts) %in% rownames(deseq.deg), ])
heatmap(log2(exp.degs+1))
aheatmap(log2(exp.degs+1), annCol =pheno$Death, main="mRNA Alive vs Dead")



