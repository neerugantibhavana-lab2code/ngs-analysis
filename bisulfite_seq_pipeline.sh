#!/usr/bin/env bash
# ==============================================================================
# Bisulfite Sequencing (DNA Methylation) Analysis Pipeline
# Aim: Map DNA methylation at single-base resolution.
#
# Principle: Sodium bisulfite converts unmethylated cytosines (C) -> uracil (T),
#            while methylated cytosines (5mC) remain unconverted. Comparing
#            converted vs unconverted reads reveals per-base methylation status.
#
# Requires: Bismark (bismark_genome_preparation, bismark, bismark_methylation_extractor),
#           bowtie2, samtools
# ==============================================================================
set -euo pipefail

REF_DIR="BS_seq_WGS"          # directory containing the reference genome FASTA
READ1="${REF_DIR}/R1.fastq.gz"
READ2="${REF_DIR}/R2.fastq.gz"

# ---- 1. Genome preparation for bisulfite mapping --------------------------------------
bismark_genome_preparation "${REF_DIR}"

# ---- 2. Align bisulfite-treated reads (internally performs C->T / G->A conversion) -----
bismark -bowtie2 "${REF_DIR}" "${READ1}"
bismark -bowtie2 "${REF_DIR}" "${READ2}"

# Output SAM files, e.g. R1_bismark_bt2.sam / R2_bismark_bt2.sam
# Methylation call tags within the SAM file:
#   x = methylated, h = hemimethylated, . = unmethylated

# ---- 3. Convert SAM -> BAM and sort ----------------------------------------------------
for label in R1 R2; do
    samtools view -bS "${label}_bismark_bt2.sam" > "${label}_bismark_bt2.bam"
    samtools sort -o "${label}_bismark_bt2_sorted.bam" "${label}_bismark_bt2.bam"
done

# ---- 4. Extract methylation calls (CpG / CHG / CHH context) ----------------------------
bismark_methylation_extractor --bedGraph --counts R1_bismark_bt2.bam
bismark_methylation_extractor --bedGraph --counts R2_bismark_bt2.bam

# Generates *.bedGraph.gz files and per-sample methylation summary reports
# (check the splitting report for % methylation at CpG, CHG, and CHH contexts).

# ---- 5. Visualization in IGV ------------------------------------------------------------
echo "Load in IGV for visualization:"
echo "  Reference genome : ref.fa"
echo "  Tracks           : R1_bismark_bt2.bedGraph.gz, R2_bismark_bt2.bedGraph.gz"
echo "  (Bar height/color intensity reflects % methylation at each cytosine, 0-100%)"

echo "Bisulfite sequencing pipeline complete."
