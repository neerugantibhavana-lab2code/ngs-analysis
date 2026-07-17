#!/usr/bin/env python3
"""
ChIP-Seq Analysis Pipeline
Aim: Identify protein-DNA interaction sites (peaks) genome-wide, annotate
     them relative to genes, and discover enriched DNA motifs.

Requires: sra-tools (prefetch, fastq-dump), fastqc, bowtie2, samtools,
          macs2, HOMER (mergePeaks, annotatePeaks.pl), bedtools, meme
          - all installed and on $PATH
"""

import subprocess

REF_GENOME = "E.coli_BW25113.fasta"
REF_GTF = "E.coli_BW25113_annotation.gff"
REF_BED = "E.coli_BW25113_annotation.bed"
INDEX_PREFIX = "Ecoli_index"

CHIP_SAMPLE = "SO_4933_C_CHR1_R1"     # ChIP (treatment) sample
INPUT_SAMPLE = "SO_4933_C_INR1_R1"    # input control sample
OUTDIR = "SampleC"
PEAK_PREFIX = "Sorted_C"
PEAKS_FILE = f"{OUTDIR}/{PEAK_PREFIX}_peaks.narrowPeak"


def run(cmd):
    print(f">>> {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def download_from_sra(srr_id, n_reads=None):
    """Optional: fetch raw reads from SRA (skip if FASTQs already exist)."""
    run(["prefetch", srr_id])
    cmd = ["fastq-dump"]
    if n_reads:
        cmd += ["-x", str(n_reads)]
    cmd.append(srr_id)
    run(cmd)


def quality_control():
    run(["fastqc", f"{CHIP_SAMPLE}.fastq"])
    run(["fastqc", f"{INPUT_SAMPLE}.fastq"])
    # Duplicate reads can be removed later with samtools rmdup or Picard MarkDuplicates.


def align_reads():
    run(["bowtie2-build", REF_GENOME, INDEX_PREFIX])
    run(["bowtie2", "-x", INDEX_PREFIX, "-U", f"{CHIP_SAMPLE}.fastq", "-S", "align_C_CHR1.sam"])
    run(["bowtie2", "-x", INDEX_PREFIX, "-U", f"{INPUT_SAMPLE}.fastq", "-S", "align_C_INR1.sam"])


def sam_to_sorted_bam(label):
    run(["samtools", "view", "-bS", f"align_C_{label}.sam", "-o", f"align_C_{label}.bam"])
    run(["samtools", "sort", f"align_C_{label}.bam", "-o", f"align_C_{label}_sorted.bam"])


def call_peaks():
    run([
        "macs2", "callpeak",
        "-t", "align_C_CHR1_sorted.bam",
        "-c", "align_C_INR1_sorted.bam",
        "-f", "BAM", "-g", "4.6e6",
        "--outdir", OUTDIR, "-n", PEAK_PREFIX,
        "-B", "-p", "0.001",
    ])
    # Visualize in IGV:
    #   1. Load reference genome (REF_GENOME)
    #   2. Load PEAKS_FILE
    #   3. Right-click -> Autoscale


def merge_peaks():
    with open("merge500.txt", "w") as out:
        print(">>> mergePeaks", PEAKS_FILE, "-d 500")
        subprocess.run(["mergePeaks", PEAKS_FILE, "-d", "500"], stdout=out, check=True)


def annotate_peaks():
    with open("annotated_peaks.tsv", "w") as out:
        print(">>> annotatePeaks.pl", PEAKS_FILE, REF_GENOME, "-gff", REF_GTF)
        subprocess.run(
            ["annotatePeaks.pl", PEAKS_FILE, REF_GENOME, "-gff", REF_GTF],
            stdout=out, check=True,
        )


def closest_genes():
    # closestBed output is piped through an awk filter keeping only "gene" features
    cmd = f"closestBed -a {PEAKS_FILE} -b {REF_BED} | awk -F '\\t' '{{if($18 == \"gene\"){{print}}}}' > C_closest_genes.bed"
    print(f">>> {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def extract_peak_sequences():
    run(["bedtools", "getfasta", "-fi", REF_GENOME, "-bed", PEAKS_FILE, "-fo", "seq.fasta"])


def discover_motifs():
    run(["meme", "seq.fasta", "-dna", "-nmotifs", "3"])


def main():
    quality_control()
    align_reads()
    sam_to_sorted_bam("CHR1")
    sam_to_sorted_bam("INR1")
    call_peaks()
    merge_peaks()
    annotate_peaks()
    closest_genes()
    extract_peak_sequences()
    discover_motifs()

    print("ChIP-Seq pipeline complete.")
    print(f"  Peaks              : {PEAKS_FILE}")
    print("  Annotated peaks    : annotated_peaks.tsv")
    print("  Closest genes      : C_closest_genes.bed")
    print("  Motif results (MEME): meme_out/")


if __name__ == "__main__":
    main()
