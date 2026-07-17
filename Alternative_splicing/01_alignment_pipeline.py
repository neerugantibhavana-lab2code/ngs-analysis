#!/usr/bin/env python3
"""
Alternative Splicing - Alignment & BAM Preparation Pipeline
Aim: Align RNA-Seq reads with a splice-aware aligner (HISAT2), then convert,
     sort, and index the resulting BAM files for downstream splice analysis
     with SGSeq (see 02_sgseq_analysis.R).

Requires: hisat2, samtools installed and on $PATH
"""

import subprocess

REF_GENOME = "Copy_of_GCF_000146045.2_R64_genomic.fna"
INDEX_PREFIX = "Hisat_index"

SAMPLES = {
    "A_1": "Copy_of_conA_rep1.fq",
    "A_2": "Copy_of_conA_rep2.fq",
    "B_1": "Copy_of_conB_rep1.fq",
    "B_2": "Copy_of_conB_rep2.fq",
}


def run(cmd):
    print(f">>> {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def build_index():
    run(["hisat2-build", REF_GENOME, INDEX_PREFIX])


def align(label, fastq):
    run(["hisat2", "-x", INDEX_PREFIX, "-U", fastq, "-S", f"{label}.sam"])


def sam_to_sorted_bam(label):
    run(["samtools", "view", "-bS", f"{label}.sam", "-o", f"{label}.bam"])
    run(["samtools", "sort", "-o", f"{label}_sorted.bam", f"{label}.bam"])
    run(["samtools", "index", f"{label}_sorted.bam"])


def main():
    build_index()

    for label, fastq in SAMPLES.items():
        align(label, fastq)

    for label in SAMPLES:
        sam_to_sorted_bam(label)

    print(f"Alignment, sorting, and indexing complete for: {', '.join(SAMPLES)}")


if __name__ == "__main__":
    main()
