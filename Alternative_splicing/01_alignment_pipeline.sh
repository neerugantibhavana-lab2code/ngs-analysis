#!/usr/bin/env bash
# ==============================================================================
# Alternative Splicing - Alignment & BAM Preparation Pipeline
# Aim: Align RNA-Seq reads with a splice-aware aligner (HISAT2), then convert,
#      sort, and index the resulting BAM files for downstream splice analysis
#      with SGSeq (see 02_sgseq_analysis.R).
#
# Requires: hisat2, samtools
# ==============================================================================
set -euo pipefail

REF_GENOME="Copy_of_GCF_000146045.2_R64_genomic.fna"
INDEX_PREFIX="Hisat_index"

declare -A SAMPLES=(
  [A_1]="Copy_of_conA_rep1.fq"
  [A_2]="Copy_of_conA_rep2.fq"
  [B_1]="Copy_of_conB_rep1.fq"
  [B_2]="Copy_of_conB_rep2.fq"
)

# ---- 1. Build the genome index (Burrows-Wheeler Transform / FM-index) --------
hisat2-build "${REF_GENOME}" "${INDEX_PREFIX}"

# ---- 2. Align each single-end sample -----------------------------------------
for label in "${!SAMPLES[@]}"; do
    fastq="${SAMPLES[$label]}"
    echo ">>> Aligning ${label} (${fastq})"
    hisat2 -x "${INDEX_PREFIX}" -U "${fastq}" -S "${label}.sam"
done

# ---- 3. Convert SAM -> BAM ------------------------------------------------------
for label in "${!SAMPLES[@]}"; do
    samtools view -bS "${label}.sam" > "${label}.bam"
done

# ---- 4. Sort BAM files ------------------------------------------------------------
for label in "${!SAMPLES[@]}"; do
    samtools sort -o "${label}_sorted.bam" "${label}.bam"
done

# ---- 5. Index sorted BAM files -----------------------------------------------------
for label in "${!SAMPLES[@]}"; do
    samtools index "${label}_sorted.bam"
done

echo "Alignment, sorting, and indexing complete for: ${!SAMPLES[*]}"
