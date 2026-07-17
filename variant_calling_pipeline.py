#!/usr/bin/env python3
"""
Variant Calling Pipeline
Aim: Identify and annotate genetic variants (SNPs/Indels) from raw
     sequencing reads using BWA-MEM alignment + GATK best practices.

Requires: fastqc, bwa, samtools, picard, gatk - all installed and on $PATH
"""

import subprocess

REF_GENOME = "GCF_000146045.2_R64_genomic.fa"
SAMPLE = "conA"
FASTQ = f"{SAMPLE}.fastq"


def run(cmd):
    print(f">>> {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def quality_control():
    run(["fastqc", FASTQ])
    # Per-base quality > 30 and no adapter contamination -> trimming not required.


def align_with_bwa():
    run(["bwa", "index", REF_GENOME])
    with open(f"{SAMPLE}.bam", "wb") as out:
        print(">>> bwa mem", REF_GENOME, FASTQ)
        subprocess.run(["bwa", "mem", REF_GENOME, FASTQ], stdout=out, check=True)


def sort_and_index():
    run(["samtools", "sort", f"{SAMPLE}.bam", "-o", f"{SAMPLE}_sorted.bam"])
    run(["samtools", "index", f"{SAMPLE}_sorted.bam"])
    # View in IGV manually:
    #   Genome -> Load genome -> REF_GENOME
    #   File   -> Load from file -> {SAMPLE}_sorted.bam


def mark_duplicates():
    run([
        "picard", "MarkDuplicates",
        f"I={SAMPLE}_sorted.bam",
        f"O={SAMPLE}_markdup.bam",
        "M=marked_dup_metrics.txt",
    ])


def add_read_groups():
    run([
        "picard", "AddOrReplaceReadGroups",
        f"I={SAMPLE}_markdup.bam",
        f"O={SAMPLE}_grpadded.bam",
        "RGID=4", "RGLB=lib1", "RGPL=illumina", "RGPU=unit1",
        f"RGSM={SAMPLE}", "RGCN=bi",
    ])


def prepare_reference_for_gatk():
    run(["gatk", "CreateSequenceDictionary", "-R", REF_GENOME])
    run(["samtools", "faidx", REF_GENOME])


def initial_variant_calling():
    run([
        "gatk", "HaplotypeCaller",
        "-I", f"{SAMPLE}_grpadded.bam",
        "-R", REF_GENOME,
        "-O", f"{SAMPLE}_variants.vcf",
    ])


def base_quality_recalibration():
    run([
        "gatk", "BaseRecalibrator",
        "-I", f"{SAMPLE}_grpadded.bam",
        "-R", REF_GENOME,
        "--known-sites", f"{SAMPLE}_variants.vcf",
        "-O", f"BQSR_out_{SAMPLE}.table",
    ])
    run([
        "gatk", "ApplyBQSR",
        "-I", f"{SAMPLE}_grpadded.bam",
        "-R", REF_GENOME,
        "-bqsr", f"BQSR_out_{SAMPLE}.table",
        "-O", f"{SAMPLE}_BQSR_applied_reads.bam",
    ])


def final_variant_calling():
    run([
        "gatk", "HaplotypeCaller",
        "-R", REF_GENOME,
        "-I", f"{SAMPLE}_BQSR_applied_reads.bam",
        "-O", f"{SAMPLE}_variants_final.vcf",
    ])


def filter_variants():
    run([
        "gatk", "VariantFiltration",
        "-R", REF_GENOME,
        "-V", f"{SAMPLE}_variants_final.vcf",
        "--filter-expression", "QD < 10.0",         "--filter-name", "LowQD",
        "--filter-expression", "FS > 60.0",         "--filter-name", "HighFS",
        "--filter-expression", "MQ < 40.0",         "--filter-name", "LowMQ",
        "--filter-expression", "MQRankSum < -12.5", "--filter-name", "LowMQRankSum",
        "-O", f"{SAMPLE}_variants_final.filtered.vcf",
    ])


def prepare_for_annotation():
    filtered_vcf = f"{SAMPLE}_variants_final.filtered.vcf"
    no_header = f"{SAMPLE}_variants_final_no_header.txt"

    with open(filtered_vcf) as f_in, open(no_header, "w") as f_out:
        for line in f_in:
            if line.startswith("NC"):
                f_out.write(line)

    chrom_ids = set()
    with open(no_header) as f_in:
        for line in f_in:
            chrom_ids.add(line.split("\t")[0])

    with open("gene_id.txt", "w") as f_out:
        for chrom in sorted(chrom_ids):
            f_out.write(chrom + "\n")

    # Reference chromosome info (downloaded separately from NCBI: e.g. S. cerevisiae R64
    # genome). Once you have S_cerevisiae_refseq.tsv, extract columns 3 and 9:
    #   with open("S_cerevisiae_refseq.tsv") as f, open("chr_name_refseq_id_of_variants.txt", "w") as out:
    #       for line in f:
    #           cols = line.rstrip("\n").split("\t")
    #           out.write(f"{cols[2]}\t{cols[8]}\n")


def main():
    quality_control()
    align_with_bwa()
    sort_and_index()
    mark_duplicates()
    add_read_groups()
    prepare_reference_for_gatk()
    initial_variant_calling()
    base_quality_recalibration()
    final_variant_calling()
    filter_variants()
    prepare_for_annotation()

    print(f"Variant calling complete: {SAMPLE}_variants_final.filtered.vcf")
    print(f"Next step: annotate {SAMPLE}_variants_final.filtered.vcf with Ensembl VEP")
    print("  (expect consequence categories such as missense_variant,")
    print("   synonymous_variant, frameshift_variant, etc.)")


if __name__ == "__main__":
    main()
