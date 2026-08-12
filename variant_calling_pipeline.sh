#!/usr/bin/env bash
# ==============================================================================
# Variant Calling Pipeline
# Aim: Identify and annotate genetic variants (SNPs/Indels) from raw
#      sequencing reads using BWA-MEM alignment + GATK best practices.
#
# Requires: fastqc, bwa, samtools, picard, gatk, awk
# ==============================================================================
set -euo pipefail

REF_GENOME="GCF_000146045.2_R64_genomic.fa"
SAMPLE="conA"
FASTQ="${SAMPLE}.fastq"

# ---- 1. Quality control -------------------------------------------------------------
fastqc "${FASTQ}"
# Per-base quality > 30 and no adapter contamination -> trimming not required.

# ---- 2. Alignment with BWA-MEM -------------------------------------------------------
bwa index "${REF_GENOME}"
bwa mem "${REF_GENOME}" "${FASTQ}" > "${SAMPLE}.bam"

# ---- 3. Sort and index alignment -----------------------------------------------------
samtools sort "${SAMPLE}.bam" -o "${SAMPLE}_sorted.bam"
samtools index "${SAMPLE}_sorted.bam"

# View in IGV manually:
#   Genome -> Load genome -> REF_GENOME
#   File   -> Load from file -> ${SAMPLE}_sorted.bam

# ---- 4. Mark duplicates (Picard) -------------------------------------------------------
picard MarkDuplicates \
    I="${SAMPLE}_sorted.bam" \
    O="${SAMPLE}_markdup.bam" \
    M=marked_dup_metrics.txt

# ---- 5. Add read group information -----------------------------------------------------
picard AddOrReplaceReadGroups \
    I="${SAMPLE}_markdup.bam" \
    O="${SAMPLE}_grpadded.bam" \
    RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM="${SAMPLE}" RGCN=bi

# ---- 6. Prepare reference genome for GATK ------------------------------------------------
gatk CreateSequenceDictionary -R "${REF_GENOME}"
samtools faidx "${REF_GENOME}"

# ---- 7. Initial variant calling (used as "known sites" for BQSR) ------------------------
gatk HaplotypeCaller \
    -I "${SAMPLE}_grpadded.bam" \
    -R "${REF_GENOME}" \
    -O "${SAMPLE}_variants.vcf"

# ---- 8. Base Quality Score Recalibration (BQSR) ------------------------------------------
gatk BaseRecalibrator \
    -I "${SAMPLE}_grpadded.bam" \
    -R "${REF_GENOME}" \
    --known-sites "${SAMPLE}_variants.vcf" \
    -O "BQSR_out_${SAMPLE}.table"

gatk ApplyBQSR \
    -I "${SAMPLE}_grpadded.bam" \
    -R "${REF_GENOME}" \
    -bqsr "BQSR_out_${SAMPLE}.table" \
    -O "${SAMPLE}_BQSR_applied_reads.bam"

# ---- 9. Final, high-confidence variant calling -------------------------------------------
gatk HaplotypeCaller \
    -R "${REF_GENOME}" \
    -I "${SAMPLE}_BQSR_applied_reads.bam" \
    -O "${SAMPLE}_variants_final.vcf"

# ---- 10. Variant filtration -----------------------------------------------------------------
gatk VariantFiltration \
    -R "${REF_GENOME}" \
    -V "${SAMPLE}_variants_final.vcf" \
    --filter-expression "QD < 10.0"          --filter-name "LowQD" \
    --filter-expression "FS > 60.0"          --filter-name "HighFS" \
    --filter-expression "MQ < 40.0"          --filter-name "LowMQ" \
    --filter-expression "MQRankSum < -12.5"  --filter-name "LowMQRankSum" \
    -O "${SAMPLE}_variants_final.filtered.vcf"

# ---- 11. Prepare for annotation ---------------------------------------------------------------
awk '($1 ~ "^NC")' "${SAMPLE}_variants_final.filtered.vcf" > "${SAMPLE}_variants_final_no_header.txt"
cut -f 1 "${SAMPLE}_variants_final_no_header.txt" | sort | uniq > gene_id.txt

# Reference chromosome info (downloaded separately from NCBI: S. cerevisiae R64 genome)
# cut -f 3,9 S_cerevisiae_refseq.tsv > chr_name_refseq_id_of_variants.txt

echo "Variant calling complete: ${SAMPLE}_variants_final.filtered.vcf"
echo "Next step: annotate ${SAMPLE}_variants_final.filtered.vcf with Ensembl VEP"
echo "  (species: your organism; expect consequence categories such as"
echo "   missense_variant, synonymous_variant, frameshift_variant, etc.)"
