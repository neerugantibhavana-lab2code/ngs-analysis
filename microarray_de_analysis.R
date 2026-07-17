#!/usr/bin/env Rscript
# ==============================================================================
# Microarray Differential Expression Analysis
# Dataset : E-MTAB-6095 (RNase I-deficient E. coli vs wild-type BW25113)
# Array   : Affymetrix GeneChip E.coli Genome 2.0 Array
#
# Aim: Detect differentially expressed genes between WT and mutant samples
#      using PM/MM probe intensities, RMA normalization, and limma.
#
# Data source: ArrayExpress -> search "E-MTAB-6095" -> Data files -> Download all -> unzip
# ==============================================================================

suppressPackageStartupMessages({
  library(oligo)
  library(limma)
  library(dplyr)
})

# ---- 1. Setup ----------------------------------------------------------------
data_dir <- "data/E-MTAB-6095"     # directory containing .CEL files
anno_file <- "data/anno.txt"        # gene ID <-> annotation mapping file

celFiles <- list.celfiles(data_dir, full.names = TRUE)
rawData  <- read.celfiles(celFiles)

# ---- 2. Extract PM / MM probe intensities ------------------------------------
pmdata <- data.frame(pm(rawData))
mmdata <- data.frame(mm(rawData))

# ---- 3. Exploratory QC: density plots of PM/MM intensities -------------------
plot(density(mmdata[, 1]), xlim = c(1, 1000), main = "Density plot for MM data")
for (i in 2:ncol(mmdata)) lines(density(mmdata[, i]), col = i)

plot(density(pmdata[, 1]), xlim = c(1, 1000), main = "Density plot for PM data")
for (i in 2:ncol(pmdata)) lines(density(pmdata[, i]), col = i)

# ---- 4. Positive (true signal) vs negative (background noise) ---------------
positive <- pmdata[, 1] - mmdata[, 1]   # PM > MM  -> true signal
negative <- mmdata[, 1] - pmdata[, 1]   # MM > PM  -> background noise

plot(positive, main = "PM-MM signal separation")
points(negative, col = "red")

# ---- 5. Log2 transform PM/MM and compare distributions -----------------------
pmlog2 <- log2(pmdata)
mmlog2 <- log2(mmdata)

par(mfrow = c(1, 2))
hist(pmdata[, 1], main = "PM data (raw)", xlab = "PM intensity", col = "pink")
hist(pmlog2[, 1], main = "PM data (log2)", xlab = "log2(PM)", col = "light blue")
par(mfrow = c(1, 1))

# ---- 6. MA plot: raw data -----------------------------------------------------
A_raw <- 0.5 * log2(pmdata[, 1] + pmdata[, 4])
M_raw <- log2(pmdata[, 1]) - log2(pmdata[, 4])
plot(A_raw, M_raw, xlab = "Average log2 intensity", ylab = "log2 fold change",
     main = "MA plot (raw data)")

# ---- 7. RMA normalization ------------------------------------------------------
normalizedData <- rma(rawData)
normalizedData <- exprs(normalizedData)
rawLog2 <- log2(exprs(rawData))

# MA plot after normalization
A_norm <- 0.5 * log2(normalizedData[, 1] + normalizedData[, 4])
M_norm <- log2(normalizedData[, 1]) - log2(normalizedData[, 4])
plot(A_norm, M_norm, xlab = "Average log2 intensity", ylab = "log2 fold change",
     main = "MA plot (normalized data)")

# ---- 8. Density & boxplots: raw vs normalized ---------------------------------
sample_names <- c("WT1", "WT2", "WT3", "M1", "M2", "M3")

par(mfrow = c(1, 2))
plot(density(rawLog2[, 1]), xlim = c(1, 15), main = "Raw Data")
for (i in 2:ncol(rawLog2)) lines(density(rawLog2[, i]), col = i)

plot(density(normalizedData[, 1]), xlim = c(1, 15), main = "Normalized Data")
for (i in 2:ncol(normalizedData)) lines(density(normalizedData[, i]), col = i)
par(mfrow = c(1, 1))

par(mfrow = c(1, 2))
boxplot(rawLog2, names = sample_names, main = "Boxplot: log2(rawData)", col = "light green")
boxplot(normalizedData, names = sample_names, main = "Boxplot: Normalized data", col = "pink")
par(mfrow = c(1, 1))

# ---- 9. Differential expression with limma ------------------------------------
Groups <- factor(c("WT", "WT", "WT", "MUT", "MUT", "MUT"))
Groups <- relevel(Groups, ref = "WT")

design <- model.matrix(~Groups)
fit <- lmFit(normalizedData, design)
fit <- eBayes(fit)

DiffExp <- topTable(fit, number = Inf, adjust.method = "BH")
DiffExp$ID <- rownames(DiffExp)

# ---- 10. Filter up/down-regulated genes ---------------------------------------
upregulated_0.6   <- subset(DiffExp, adj.P.Val < 0.05 & logFC > 0.6)
upregulated_1     <- subset(DiffExp, adj.P.Val < 0.05 & logFC > 1)
downregulated_1   <- subset(DiffExp, adj.P.Val < 0.05 & logFC < -1)

# ---- 11. Volcano plot -----------------------------------------------------------
no_DiffExp   <- -log10(DiffExp$adj.P.Val)
up_DiffExp   <- -log10(upregulated_1$adj.P.Val)
down_DiffExp <- -log10(downregulated_1$adj.P.Val)

plot(DiffExp$logFC, no_DiffExp, xlab = "log2FC", ylab = "-log10(adj.P.Val)",
     main = "Volcano Plot")
points(upregulated_1$logFC, up_DiffExp, col = "red")
points(downregulated_1$logFC, down_DiffExp, col = "blue")

# ---- 12. Annotate DE results ----------------------------------------------------
anno <- read.delim(anno_file)
write.table(DiffExp, "DiffExp_Table.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df <- inner_join(anno, DiffExp, by = "ID")
write.table(df, "DiffExp_Annotated.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("Analysis complete. Annotated results written to DiffExp_Annotated.tsv\n")
