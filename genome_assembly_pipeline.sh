#!/usr/bin/env bash
# ==============================================================================
# De Novo Genome Assembly Pipeline
# Aim: Assemble short paired-end reads into contigs (Velvet), assess quality
#      (QUAST), predict genes (Augustus), scaffold with a reference (RagTag),
#      and assess completeness (BUSCO).
#
# Requires: fastqc, velvet (velveth/velvetg), quast.py, augustus, gffread,
#           ragtag.py, busco
# ==============================================================================
set -euo pipefail

READ1="SR1.fastq"
READ2="SR2.fastq"
REF_GENOME="GCA_014131755.1_ASM1413175v1_genomic.fna"
KMER=31

# ---- 1. Quality control ------------------------------------------------------------
fastqc "${READ1}"
fastqc "${READ2}"
# Phred quality > 30 observed -> no trimming required in this run.
# If trimming is needed: run Trimmomatic here before proceeding.

# ---- 2. Assembly with Velvet (De Bruijn graph) --------------------------------------
velveth Assembly_step1 "${KMER}" -shortPaired -fastq -separate "${READ1}" "${READ2}"
velvetg Assembly_step1 -cov_cutoff 5 -unused_reads yes

CONTIGS="Assembly_step1/contigs.fa"

# ---- 3. Assess contig quality with QUAST --------------------------------------------
quast.py -o quast_contigs "${CONTIGS}"

# ---- 4. Ab initio gene prediction with Augustus -------------------------------------
augustus --species=E_coli_K12 "${CONTIGS}" > contig_1.gff

# ---- 5. Extract coding sequences ------------------------------------------------------
gffread -w multi.fasta -g contig_1.gff

# ---- 6. Scaffolding with RagTag ---------------------------------------------------------
ragtag.py correct "${REF_GENOME}" "${CONTIGS}"
ragtag.py scaffold "${REF_GENOME}" ragtag_output/ragtag.correct.fasta

SCAFFOLD="ragtag_output/ragtag.scaffold.fasta"

# ---- 7. Assess scaffold quality with QUAST -----------------------------------------------
quast.py -o quast_scaffold "${SCAFFOLD}"

# ---- 8. Assess assembly completeness with BUSCO ------------------------------------------
busco -i "${SCAFFOLD}" -m geno -o busco_output --auto-lineage-prok

echo "Genome assembly pipeline complete."
echo "  Contig QUAST report : quast_contigs/"
echo "  Scaffold QUAST report: quast_scaffold/"
echo "  BUSCO report         : busco_output/"
