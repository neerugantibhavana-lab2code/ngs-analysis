#!/usr/bin/env python3
"""
RNA-Seq Alignment Pipeline
Aim: Align single-end RNA-Seq reads to a reference genome using Bowtie2
     and produce SAM files ready for feature counting.

Requires: bowtie2 installed and on $PATH
"""

import subprocess

REF_GENOME = "Copy_of_GCF_000146045.2_R64_genomic.fna"
INDEX_PREFIX = "Sach_index"
SAMPLES = ["conA_rep1", "conA_rep2", "conB_rep1", "conB_rep2"]


def run(cmd):
    """Print and execute a shell command, raising on failure."""
    print(f">>> {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def build_index():
    run(["bowtie2-build", REF_GENOME, INDEX_PREFIX])


def align_sample(sample):
    fastq = f"Copy_of_{sample}.fq.gz"
    sam_out = f"{sample}.sam"
    run(["bowtie2", "-x", INDEX_PREFIX, "-U", fastq, "-S", sam_out])


def main():
    build_index()
    for sample in SAMPLES:
        align_sample(sample)
    print(f"Alignment complete. SAM files generated for: {', '.join(SAMPLES)}")


if __name__ == "__main__":
    main()
