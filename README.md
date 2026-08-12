# NGS Analysis & Bioinformatics Coursework

A collection of NGS and bioinformatics analysis workflows completed as part of my M.Sc. Biotechnology and Bioinformatics coursework at the Institute of Bioinformatics and Applied Biotechnology (IBAB).
The repository contains practical analyses of sequencing and genomic data analysis workflows, including transcriptomics, genome assembly, variant analysis, ChIP-seq, metagenomics, bisulfite sequencing, and microarray analysis.

## Analyses Included

| Analysis | Tools / Methods |
|---|---|
| RNA-seq Analysis | FastQC, Trimmomatic, Bowtie2, SAMtools, featureCounts/Rsubread, edgeR |
| Alternative Splicing | HISAT2, SAMtools, SGSeq |
| Microarray Analysis | R, oligo, RMA, limma |
| Genome Assembly | Velvet, Canu, QUAST, Augustus, gffread, RagTag, BUSCO |
| Variant Calling & Annotation | BWA-MEM, SAMtools, Picard, GATK, IGV, Ensembl VEP |
| ChIP-seq Analysis | SRA Toolkit, FastQC, Bowtie2, SAMtools, MACS2, IGV, HOMER, BEDTools, MEME |
| Bisulfite Sequencing | Bismark, Bowtie2, SAMtools |
| Metagenomics | QIIME 2 |

## RNA-seq Analysis

The RNA-seq workflow covers processing of raw sequencing reads through differential expression analysis.

### Workflow

- Quality assessment of raw FASTQ reads using **FastQC**
- Read trimming using **Trimmomatic**
- Reference genome indexing and read alignment using **Bowtie2**
- Alignment file processing using **SAMtools**
- Generation of gene-level counts using **featureCounts/Rsubread**
- Read-count normalization and filtering
- Differential gene expression analysis using **edgeR**
- Generation and interpretation of differential expression plots

## Alternative Splicing Analysis

- Quality assessment of RNA-seq reads
- Splice-aware alignment using **HISAT2**
- BAM processing using **SAMtools**
- Transcript feature and splice-event analysis using **SGSeq**
- Analysis of exon, intron and splice-junction features
- Visualization of splice features and coverage

## Microarray Analysis

Microarray differential expression analysis was performed using R.

### Workflow

- Microarray data retrieval
- Processing of Affymetrix CEL files using **oligo**
- Log transformation and quality assessment
- RMA normalization
- Density plots, boxplots and MA plots
- Experimental design and statistical modeling
- Differential expression analysis using **limma**
- Identification of differentially expressed genes
- Volcano plot generation and interpretation

## De novo Genome Assembly

The genome assembly workflow includes:

- Quality assessment of sequencing reads
- Read preprocessing
- De novo assembly using **Velvet** and **Canu**
- Assembly quality assessment using **QUAST**
- Gene prediction using **Augustus**
- Coding sequence extraction using **gffread**
- Assembly scaffolding/correction using **RagTag**
- Assembly completeness assessment using **BUSCO**

## Variant Calling & Annotation

The variant analysis workflow includes:

- Read quality assessment
- Reference genome alignment using **BWA-MEM**
- BAM processing using **SAMtools**
- Alignment visualization using **IGV**
- Duplicate/read-group processing using **Picard**
- Variant calling using **GATK HaplotypeCaller**
- Base Quality Score Recalibration
- Variant filtering
- Variant annotation using **Ensembl VEP**
- Interpretation of predicted variant consequences

## ChIP-seq Analysis

The ChIP-seq workflow includes:

- Retrieval of sequencing data using **SRA Toolkit**
- FASTQ generation and quality assessment
- Read alignment using **Bowtie2**
- BAM processing using **SAMtools**
- Peak calling using **MACS2**
- Peak visualization using **IGV**
- Peak annotation using **HOMER**
- Genomic interval analysis using **BEDTools**
- DNA motif analysis using **MEME Suite**

## Bisulfite Sequencing

The repository also contains a bisulfite sequencing workflow covering:

- Read preprocessing and quality assessment
- Bisulfite-aware alignment using **Bismark**
- Alignment processing using **Bowtie2/SAMtools**
- Downstream methylation analysis

## Metagenomics

The metagenomics workflow uses **QIIME 2** for processing and analysis of sequencing data.

The workflow covers sequence preprocessing, quality assessment, feature analysis and downstream interpretation of microbial community data.

## Tools & Technologies

### Programming & Scripting
- R
- Bash / Linux command line
- AWK

### NGS & Transcriptomics
- FastQC
- Trimmomatic
- Bowtie2
- HISAT2
- SAMtools
- featureCounts / Rsubread
- edgeR
- limma
- SGSeq

### Genome Assembly
- Velvet
- Canu
- QUAST
- Augustus
- gffread
- RagTag
- BUSCO

### Variant Analysis
- BWA-MEM
- Picard
- GATK
- IGV
- Ensembl VEP

### ChIP-seq
- SRA Toolkit
- MACS2
- HOMER
- BEDTools
- MEME Suite

### Microarray
- R
- oligo
- RMA
- limma

### Other NGS Workflows
- Bismark
- QIIME 2

## Repository Structure

```text
ngs-analysis/
│
├── Alternative_splicing/
├── RNAseq_Analysis/
│
├── bisulfite_seq_pipeline.sh
├── chipseq_pipeline.sh
├── genome_assembly_pipeline.sh
├── microarray_de_analysis.R
├── qiime2_pipeline.sh
├── variant_calling_pipeline.sh
│
└── README.md
