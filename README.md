# NGS Analysis & Bioinformatics Coursework

This repository contains bioinformatics analysis workflows and practical exercises completed as part of my **M.Sc. Biotechnology & Bioinformatics coursework at the Institute of Bioinformatics and Applied Biotechnology (IBAB), Bengaluru**.

The work covers different areas of genomics, transcriptomics, epigenomics, metagenomics, and sequence data analysis. These projects were carried out to gain practical experience with commonly used bioinformatics tools and to understand the major steps involved in NGS data analysis.

## Analyses Included

| Analysis | What I Practiced | Tools |
|---|---|---|
| **Microarray Analysis** | Differential gene expression analysis | R, `oligo`, `limma` |
| **RNA-seq Analysis** | Read alignment, read quantification and differential expression analysis | Bowtie2, Rsubread, edgeR |
| **Alternative Splicing** | Splice-aware alignment and analysis of splice junctions | HISAT2, SGSeq |
| **Genome Assembly** | De novo assembly, assembly quality assessment and gene prediction | Velvet, QUAST, Augustus, RagTag, BUSCO |
| **Variant Calling** | Identification of SNPs and small insertions/deletions from sequencing data | BWA-MEM, Picard, GATK |
| **ChIP-seq Analysis** | Read alignment, peak calling and downstream analysis | Bowtie2, MACS2, HOMER, BEDTools, MEME |
| **QIIME2 / Metagenomics** | 16S amplicon sequencing and microbial community analysis | QIIME 2 |
| **Bisulfite Sequencing** | Analysis of DNA methylation sequencing data | Bismark, Bowtie2, samtools |

## General Workflow

Depending on the analysis, the practical workflows involved steps such as:

**Raw sequencing data → Quality assessment → Pre-processing → Alignment / Assembly → Quantification or Variant/Peak Analysis → Statistical or Downstream Analysis → Interpretation**

The exact steps and tools vary between projects.

## Repository Structure

- `RNAseq_Analysis/` – RNA-seq analysis
- `Alternative_splicing/` – Alternative splicing analysis
- `microarray_de_analysis.R` – Microarray differential expression analysis
- `genome_assembly_pipeline.sh` – Genome assembly workflow
- `variant_calling_pipeline.sh` – Variant calling workflow
- `chipseq_pipeline.sh` – ChIP-seq workflow
- `qiime2_pipeline.sh` – QIIME2 / metagenomics workflow
- `bisulfite_seq_pipeline.sh` – Bisulfite sequencing workflow

## About the Work

These analyses were completed as **academic coursework and practical exercises** during my M.Sc. program.

The purpose of this repository is to document my **hands-on learning and exposure to NGS and bioinformatics analysis workflows**, rather than to present them as production or industry-developed pipelines.

The scripts contain comments and workflow steps that helped me understand how different tools are used together in a typical bioinformatics analysis.

## Skills Practiced

- NGS data analysis
- Genomics and transcriptomics
- RNA-seq analysis
- Differential expression analysis
- Genome assembly
- Variant calling
- ChIP-seq analysis
- Alternative splicing analysis
- DNA methylation analysis
- Metagenomics / 16S analysis
- R programming
- Linux / shell scripting
- Bioinformatics tools and databases

## Note

The scripts were developed and practiced using **course/lab datasets**. File paths, sample names and input data may need to be modified when using the workflows with other datasets.

This repository represents my **coursework-based practical experience and learning in bioinformatics**, and I am continuing to build on these skills through further projects and applications.