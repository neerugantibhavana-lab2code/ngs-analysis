#!/usr/bin/env python3
"""
De Novo Genome Assembly Pipeline
Aim: Assemble short paired-end reads into contigs (Velvet), assess quality
     (QUAST), predict genes (Augustus), scaffold with a reference (RagTag),
     and assess completeness (BUSCO).

Requires: fastqc, velveth, velvetg, quast.py, augustus, gffread,
          ragtag.py, busco - all installed and on $PATH
"""

import subprocess
from pathlib import Path

READ1 = "SR1.fastq"
READ2 = "SR2.fastq"
REF_GENOME = "GCA_014131755.1_ASM1413175v1_genomic.fna"
KMER = "31"

ASSEMBLY_DIR = "Assembly_step1"
CONTIGS = f"{ASSEMBLY_DIR}/contigs.fa"
SCAFFOLD = "ragtag_output/ragtag.scaffold.fasta"


def run(cmd):
    print(f">>> {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def quality_control():
    run(["fastqc", READ1])
    run(["fastqc", READ2])
    # Phred quality > 30 observed in this run -> no trimming required.
    # If trimming is needed, call Trimmomatic here before proceeding.


def assemble_with_velvet():
    run(["velveth", ASSEMBLY_DIR, KMER, "-shortPaired", "-fastq", "-separate", READ1, READ2])
    run(["velvetg", ASSEMBLY_DIR, "-cov_cutoff", "5", "-unused_reads", "yes"])


def assess_contigs():
    run(["quast.py", "-o", "quast_contigs", CONTIGS])


def predict_genes():
    with open("contig_1.gff", "w") as out:
        print(">>> augustus --species=E_coli_K12", CONTIGS)
        subprocess.run(["augustus", "--species=E_coli_K12", CONTIGS], stdout=out, check=True)
    run(["gffread", "-w", "multi.fasta", "-g", "contig_1.gff"])


def scaffold_with_ragtag():
    run(["ragtag.py", "correct", REF_GENOME, CONTIGS])
    run(["ragtag.py", "scaffold", REF_GENOME, "ragtag_output/ragtag.correct.fasta"])


def assess_scaffold():
    run(["quast.py", "-o", "quast_scaffold", SCAFFOLD])


def assess_completeness():
    run(["busco", "-i", SCAFFOLD, "-m", "geno", "-o", "busco_output", "--auto-lineage-prok"])


def main():
    quality_control()
    assemble_with_velvet()
    assess_contigs()
    predict_genes()
    scaffold_with_ragtag()
    assess_scaffold()
    assess_completeness()

    print("Genome assembly pipeline complete.")
    print("  Contig QUAST report  : quast_contigs/")
    print("  Scaffold QUAST report: quast_scaffold/")
    print("  BUSCO report         : busco_output/")


if __name__ == "__main__":
    main()
