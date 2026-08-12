#!/usr/bin/env Rscript
# ==============================================================================
# Analysis of Alternative Splicing Events using SGSeq
# Aim: Identify and visualize alternative splicing events (exons, introns,
#      junctions) from sorted BAM files and a GTF annotation.
#
# Inputs: Sorted/indexed BAM files produced by 01_alignment_pipeline.sh
# ==============================================================================

suppressPackageStartupMessages(library(SGSeq))

gtf_file <- "GCF_000146045.2_R64_genomic.gtf"

# ---- 1. Prepare sample sheet -----------------------------------------------------
file_bam <- c("A_1_sorted.bam", "A_2_sorted.bam", "B_1_sorted.bam", "B_2_sorted.bam")
sample_name <- c("WT1", "WT2", "Mut1", "Mut2")
v <- data.frame(sample_name, file_bam)

# ---- 2. Extract BAM file info (library size, read length, strand) --------------
bam_info <- getBamInfo(v)
print(bam_info)

# ---- 3. Import and convert annotation -------------------------------------------
transcripts  <- importTranscripts(file = gtf_file)
tx_features  <- convertToTxFeatures(transcripts)

# ---- 4. Predict and quantify splicing features -----------------------------------
AnaFeat <- analyzeFeatures(bam_info, features = tx_features)

# FPKM values for each feature
fpkm_values <- FPKM(AnaFeat)
write.table(fpkm_values, "SGSeq_FPKM.tsv", sep = "\t", quote = FALSE)

# ---- 5. View feature metadata -----------------------------------------------------
feature_ranges <- rowRanges(AnaFeat)
print(feature_ranges)

# ---- 6. Visualize splice features for a gene of interest ------------------------
gene_of_interest <- "YAL030W"
plotFeatures(AnaFeat, geneName = gene_of_interest, include = "both")

# ---- 7. Highlight novel junctions -------------------------------------------------
Ana_Feat_novel <- analyzeFeatures(bam_info)
annotated <- annotate(Ana_Feat_novel, tx_features)
plotFeatures(annotated, geneName = gene_of_interest, include = "both", color_novel = "red")

# ---- 8. Splice graph + per-sample coverage plots ----------------------------------
gene_id_of_interest <- 500   # example internal SGSeq gene ID

par(mfrow = c(5, 1), mar = c(1, 3, 1, 1))
plotSpliceGraph(rowRanges(Ana_Feat_novel), geneID = gene_id_of_interest,
                toscale = "none", color_novel = "red")
for (j in 1:4) {
    plotCoverage(Ana_Feat_novel[, j], geneID = gene_id_of_interest, toscale = "none")
}
par(mfrow = c(1, 1))

cat("SGSeq splicing analysis complete. FPKM table written to SGSeq_FPKM.tsv\n")
