#!/usr/bin/env python3
"""
Bisulfite Sequencing (DNA Methylation) Analysis Pipeline
Aim: Map DNA methylation at single-base resolution.

Principle: Sodium bisulfite converts unmethylated cytosines (C) -> uracil (T),
           while methylated cytosines (5mC) remain unconverted. Comparing
           converted vs unconverted reads reveals per-base methylation status.

Requires: Bismark (bismark_genome_preparation, bismark,
          bismark_methylation_extractor), bowtie2, samtools - all on $PATH
"""

import subprocess

REF_DIR = "BS_seq_WGS"   # directory containing the reference genome FASTA
READ1 = f"{REF_DIR}/R1.fastq.gz"
READ2 = f"{REF_DIR}/R2.fastq.gz"


def run(cmd):
    print(f">>> {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def prepare_genome():
    run(["bismark_genome_preparation", REF_DIR])


def align_reads():
    # Internally performs C->T / G->A conversion before mapping
    run(["bismark", "-bowtie2", REF_DIR, READ1])
    run(["bismark", "-bowtie2", REF_DIR, READ2])
    # Output SAM files, e.g. R1_bismark_bt2.sam / R2_bismark_bt2.sam
    # Methylation call tags within the SAM file:
    #   x = methylated, h = hemimethylated, . = unmethylated


def sam_to_sorted_bam(label):
    run(["samtools", "view", "-bS", f"{label}_bismark_bt2.sam", "-o", f"{label}_bismark_bt2.bam"])
    run(["samtools", "sort", "-o", f"{label}_bismark_bt2_sorted.bam", f"{label}_bismark_bt2.bam"])


def extract_methylation(label):
    run(["bismark_methylation_extractor", "--bedGraph", "--counts", f"{label}_bismark_bt2.bam"])
    # Generates *.bedGraph.gz files and per-sample methylation summary reports
    # (check the splitting report for % methylation at CpG, CHG, and CHH contexts).


def main():
    prepare_genome()
    align_reads()

    for label in ["R1", "R2"]:
        sam_to_sorted_bam(label)
        extract_methylation(label)

    print("Load in IGV for visualization:")
    print("  Reference genome : ref.fa")
    print("  Tracks           : R1_bismark_bt2.bedGraph.gz, R2_bismark_bt2.bedGraph.gz")
    print("  (Bar height/color intensity reflects % methylation at each cytosine, 0-100%)")
    print("Bisulfite sequencing pipeline complete.")


if __name__ == "__main__":
    main()
