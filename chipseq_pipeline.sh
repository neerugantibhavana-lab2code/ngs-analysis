#!/usr/bin/env bash
# ==============================================================================
# ChIP-Seq Analysis Pipeline
# Aim: Identify protein-DNA interaction sites (peaks) genome-wide, annotate
#      them relative to genes, and discover enriched DNA motifs.
#
# Requires: sra-tools (prefetch, fastq-dump), fastqc, bowtie2, samtools,
#           macs2, HOMER (mergePeaks, annotatePeaks.pl), bedtools, meme
# ==============================================================================
set -euo pipefail

REF_GENOME="E.coli_BW25113.fasta"
REF_GTF="E.coli_BW25113_annotation.gff"
REF_BED="E.coli_BW25113_annotation.bed"
INDEX_PREFIX="Ecoli_index"

CHIP_SAMPLE="SO_4933_C_CHR1_R1"   # ChIP (treatment) sample
INPUT_SAMPLE="SO_4933_C_INR1_R1" # input control sample
OUTDIR="SampleC"
PEAK_PREFIX="Sorted_C"

# ---- 1. (Optional) Download raw reads from SRA --------------------------------------
# prefetch SRR_ID
# fastq-dump SRR_ID
# fastq-dump -x 10000 SRR_ID   # subset of reads only

# ---- 2. Quality control --------------------------------------------------------------
fastqc "${CHIP_SAMPLE}.fastq"
fastqc "${INPUT_SAMPLE}.fastq"
# Duplicate reads can be removed later with samtools rmdup or Picard MarkDuplicates.

# ---- 3. Build genome index and align reads --------------------------------------------
bowtie2-build "${REF_GENOME}" "${INDEX_PREFIX}"

bowtie2 -x "${INDEX_PREFIX}" -U "${CHIP_SAMPLE}.fastq"  -S align_C_CHR1.sam
bowtie2 -x "${INDEX_PREFIX}" -U "${INPUT_SAMPLE}.fastq" -S align_C_INR1.sam

# ---- 4. SAM -> BAM conversion and sorting ------------------------------------------------
samtools view -bS align_C_CHR1.sam > align_C_CHR1.bam
samtools view -bS align_C_INR1.sam > align_C_INR1.bam

samtools sort align_C_CHR1.bam -o align_C_CHR1_sorted.bam
samtools sort align_C_INR1.bam -o align_C_INR1_sorted.bam

# ---- 5. Peak calling with MACS2 -----------------------------------------------------------
macs2 callpeak \
    -t align_C_CHR1_sorted.bam \
    -c align_C_INR1_sorted.bam \
    -f BAM -g 4.6e6 \
    --outdir "${OUTDIR}" -n "${PEAK_PREFIX}" \
    -B -p 0.001

# Visualize in IGV:
#   1. Load reference genome (REF_GENOME)
#   2. Load ${OUTDIR}/${PEAK_PREFIX}_peaks.narrowPeak
#   3. Right-click -> Autoscale

# ---- 6. Merge nearby peaks -------------------------------------------------------------------
mergePeaks "${OUTDIR}/${PEAK_PREFIX}_peaks.narrowPeak" -d 500 > merge500.txt

# ---- 7. Annotate peaks (HOMER) -----------------------------------------------------------------
annotatePeaks.pl "${OUTDIR}/${PEAK_PREFIX}_peaks.narrowPeak" "${REF_GENOME}" \
    -gff "${REF_GTF}" > annotated_peaks.tsv

# ---- 8. Identify closest genes to each peak (BEDTools) ------------------------------------------
closestBed -a "${OUTDIR}/${PEAK_PREFIX}_peaks.narrowPeak" -b "${REF_BED}" \
    | awk -F "\t" '{if($18 == "gene"){print}}' > C_closest_genes.bed

# ---- 9. Extract peak sequences for motif discovery -----------------------------------------------
bedtools getfasta -fi "${REF_GENOME}" -bed "${OUTDIR}/${PEAK_PREFIX}_peaks.narrowPeak" \
    -fo seq.fasta

# ---- 10. DNA motif discovery (MEME Suite) -----------------------------------------------------------
meme seq.fasta -dna -nmotifs 3

echo "ChIP-Seq pipeline complete."
echo "  Peaks             : ${OUTDIR}/${PEAK_PREFIX}_peaks.narrowPeak"
echo "  Annotated peaks    : annotated_peaks.tsv"
echo "  Closest genes      : C_closest_genes.bed"
echo "  Motif results (MEME): meme_out/"
