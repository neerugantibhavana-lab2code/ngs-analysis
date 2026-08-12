#!/usr/bin/env bash
# ==============================================================================
# RNA-Seq Alignment Pipeline
# Aim: Align single-end RNA-Seq reads to a reference genome using Bowtie2
#      and produce sorted, indexed alignments ready for feature counting.
#
# Usage: ./01_alignment_pipeline.sh
# Requires: bowtie2, samtools
# ==============================================================================
set -euo pipefail

REF_GENOME="Copy_of_GCF_000146045.2_R64_genomic.fna"
INDEX_PREFIX="Sach_index"

SAMPLES=("conA_rep1" "conA_rep2" "conB_rep1" "conB_rep2")

# ---- 1. Build genome index ---------------------------------------------------
bowtie2-build "${REF_GENOME}" "${INDEX_PREFIX}"

# ---- 2. Align each sample (single-end reads) ---------------------------------
for sample in "${SAMPLES[@]}"; do
    echo ">>> Aligning ${sample}"
    bowtie2 -x "${INDEX_PREFIX}" -U "Copy_of_${sample}.fq.gz" -S "${sample}.sam"
done

echo "Alignment complete. SAM files generated for: ${SAMPLES[*]}"
