#!/usr/bin/env Rscript
# ==============================================================================
# RNA-Seq Quantification and Differential Expression Analysis
# Aim: Generate a gene x sample count matrix, compute RPKM, and identify
#      differentially expressed genes between two conditions (A vs B) using edgeR.
#
# Inputs : Sorted/aligned SAM files for 4 samples (conA_rep1/2, conB_rep1/2)
#          + a reference GTF annotation file
# ==============================================================================

suppressPackageStartupMessages({
  library(Rsubread)  # read alignment & feature counting
  library(dplyr)      # data manipulation
  library(tidyr)      # reshaping
  library(edgeR)       # normalization & DE analysis
})

gtf_file <- "GCF_000146045.2_R64_genomic.gtf"

# ---- 1. Feature counting (reads per gene) --------------------------------------
samples <- c("conA_rep1", "conA_rep2", "conB_rep1", "conB_rep2")

feat_counts <- lapply(samples, function(s) {
  featureCounts(files = paste0(s, ".sam"),
                annot.ext = gtf_file,
                isGTFAnnotationFile = TRUE,
                GTF.featureType = "exon",
                GTF.attrType = "gene_id")
})
names(feat_counts) <- samples

Feat_Count <- data.frame(lapply(feat_counts, function(x) x$counts))
colnames(Feat_Count) <- paste0(samples, ".sam")
Feat_Count <- cbind(ID = rownames(Feat_Count), Feat_Count)

# ---- 2. RPKM computation --------------------------------------------------------
lib_sizes <- sapply(feat_counts, function(x) sum(x$counts))
pmsf <- lib_sizes / 1e6                       # per-million scaling factor

rpm <- mapply(function(fc, p) fc$counts / p, feat_counts, pmsf)
colnames(rpm) <- paste0(samples, ".sam")

# gene coordinates extracted from GTF via awk (run once, outside R):
#   awk -F'\t' '{match($9,/gene_id "([^"]+)"/,a); if(a[1]!="") print a[1], $4, $5}' \
#       GCF_000146045.2_R64_genomic.gtf > gene_start_end.txt
gene_data <- read.table("gene_start_end.txt", header = TRUE, sep = "\t") %>%
  separate(ID.Start.End, into = c("ID", "Start", "End"), sep = "\\s+") %>%
  distinct(ID, .keep_all = TRUE)

df <- inner_join(Feat_Count, gene_data, by = "ID")
write.table(df, "Feat_Count_Gene_length.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

gene_length_kb <- abs(as.numeric(df$End) - as.numeric(df$Start)) / 1000

RPKM <- rpm / gene_length_kb
RPKM <- data.frame(ID = df$ID, RPKM)
RPKM_log <- data.frame(ID = df$ID, log2(RPKM[, -1]))

# ---- 3. edgeR differential expression (condition A vs B) ------------------------
group <- c(rep("A", 2), rep("B", 2))

y <- DGEList(counts = Feat_Count[, -1], group = group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~group)
y <- estimateDisp(y, design)

# BCV plot
plotBCV(y)

DEet <- exactTest(y, pair = c("A", "B"))
et_table <- DEet$table
et_table <- cbind(ID = rownames(et_table), et_table)

up_et   <- subset(et_table, logFC >= 1  & PValue <= 0.05)
down_et <- subset(et_table, logFC <= -1 & PValue <= 0.05)

# ---- 4. Volcano plot --------------------------------------------------------------
plot(et_table$logFC, -log10(et_table$PValue),
     xlab = "logFC", ylab = "-log10(PValue)", main = "Volcano Plot")
points(up_et$logFC,   -log10(up_et$PValue),   col = "red")
points(down_et$logFC, -log10(down_et$PValue), col = "blue")

# ---- 5. Heatmaps of up/down-regulated genes --------------------------------------
up_rpkm   <- semi_join(RPKM, data.frame(up_et),   by = "ID")
down_rpkm <- semi_join(RPKM, data.frame(down_et), by = "ID")

up_mat   <- as.matrix(up_rpkm[, -1])
down_mat <- as.matrix(down_rpkm[, -1])

heatmap(up_mat, Rowv = NA, Colv = NA, main = "Upregulated Genes", scale = "row",
        col = colorRampPalette(c("blue", "white", "red"))(100))
heatmap(down_mat, Rowv = NA, Colv = NA, main = "Downregulated genes", scale = "row",
        col = colorRampPalette(c("blue", "white", "red"))(100))

write.table(et_table, "DE_results_edgeR.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("Differential expression analysis complete: DE_results_edgeR.tsv\n")
