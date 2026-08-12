#!/usr/bin/env bash
# ==============================================================================
# Metagenomic Analysis using QIIME2
# Aim: Assess microbial community structure, diversity, and taxonomic
#      composition from 16S amplicon sequencing data.
#
# Requires: QIIME 2 (conda environment activated), a pre-trained taxonomy
#           classifier (e.g. gg-13-8-99-515-806-nb-classifier.qza)
# ==============================================================================
set -euo pipefail

METADATA="sample-metadata.tsv"
CLASSIFIER="gg-13-8-99-515-806-nb-classifier.qza"
SAMPLING_DEPTH=1103   # chosen from the feature-table summary (table.qzv)

# ---- 1. Import sequence data ---------------------------------------------------------
qiime tools import \
    --type 'EMPSingleEndSequences' \
    --input-path input_qiime \
    --output-path single-end-sequences.qza

# ---- 2. Demultiplex -------------------------------------------------------------------
qiime demux emp-single \
    --i-seqs single-end-sequences.qza \
    --m-barcodes-file "${METADATA}" \
    --m-barcodes-column barcode-sequence \
    --o-per-sample-sequences demux.qza \
    --o-error-correction-details demux-details.qza

qiime demux summarize \
    --i-data demux.qza \
    --o-visualization demux.qzv

# ---- 3. Denoise with DADA2 --------------------------------------------------------------
qiime dada2 denoise-single \
    --i-demultiplexed-seqs demux.qza \
    --p-trim-left 0 \
    --p-trunc-len 120 \
    --o-representative-sequences rep-seqs.qza \
    --o-table table.qza \
    --o-denoising-stats stats.qza

qiime metadata tabulate \
    --m-input-file stats.qza \
    --o-visualization stats.qzv

# ---- 4. Feature table / representative sequence summaries --------------------------------
qiime feature-table summarize \
    --i-table table.qza \
    --m-sample-metadata-file "${METADATA}" \
    --o-visualization table.qzv

qiime feature-table tabulate-seqs \
    --i-data rep-seqs.qza \
    --o-visualization rep-seqs.qzv

# ---- 5. Phylogenetic tree construction (MAFFT + FastTree) ---------------------------------
qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences rep-seqs.qza \
    --o-alignment aligned-rep-seqs.qza \
    --o-masked-alignment masked-aligned-rep-seqs.qza \
    --o-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza

qiime tools export --input-path rooted-tree.qza --output-path rooted-tree

# ---- 6. Alpha and beta diversity -----------------------------------------------------------
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted-tree.qza \
    --i-table table.qza \
    --p-sampling-depth "${SAMPLING_DEPTH}" \
    --m-metadata-file "${METADATA}" \
    --output-dir diversity-core-metrics-phylogenetic

# Alpha diversity significance (Faith's PD, evenness)
qiime diversity alpha-group-significance \
    --i-alpha-diversity diversity-core-metrics-phylogenetic/faith_pd_vector.qza \
    --m-metadata-file "${METADATA}" \
    --o-visualization faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
    --i-alpha-diversity diversity-core-metrics-phylogenetic/evenness_vector.qza \
    --m-metadata-file "${METADATA}" \
    --o-visualization evenness-group-significance.qzv

# Beta diversity significance (unweighted UniFrac, by body-site / subject)
qiime diversity beta-group-significance \
    --i-distance-matrix diversity-core-metrics-phylogenetic/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file "${METADATA}" \
    --m-metadata-column body-site \
    --p-pairwise \
    --o-visualization unweighted-unifrac-body-site-group-significance.qzv

qiime diversity beta-group-significance \
    --i-distance-matrix diversity-core-metrics-phylogenetic/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file "${METADATA}" \
    --m-metadata-column subject \
    --p-pairwise \
    --o-visualization unweighted-unifrac-subject-group-significance.qzv

# ---- 7. PCoA / Emperor plots ------------------------------------------------------------------
qiime emperor plot \
    --i-pcoa diversity-core-metrics-phylogenetic/unweighted_unifrac_pcoa_results.qza \
    --m-metadata-file "${METADATA}" \
    --p-custom-axes days-since-experiment-start \
    --o-visualization unweighted-unifrac-emperor-days-since-experiment-start.qzv

# ---- 8. Alpha rarefaction (sequencing depth vs diversity) --------------------------------------
qiime diversity alpha-rarefaction \
    --i-table table.qza \
    --i-phylogeny rooted-tree.qza \
    --p-max-depth 4000 \
    --m-metadata-file "${METADATA}" \
    --o-visualization alpha-rarefaction.qzv

# ---- 9. Taxonomic classification --------------------------------------------------------------
qiime feature-classifier classify-sklearn \
    --i-classifier "${CLASSIFIER}" \
    --i-reads rep-seqs.qza \
    --o-classification taxonomy.qza

qiime metadata tabulate \
    --m-input-file taxonomy.qza \
    --o-visualization taxonomy.qzv

qiime taxa barplot \
    --i-table table.qza \
    --i-taxonomy taxonomy.qza \
    --m-metadata-file "${METADATA}" \
    --o-visualization taxa-bar-plots.qzv

echo "QIIME2 pipeline complete. Key outputs: table.qzv, taxa-bar-plots.qzv,"
echo "  alpha-rarefaction.qzv, unweighted-unifrac-*-group-significance.qzv"
