# NGS Analysis Pipelines

This repository contains bioinformatics workflows and analysis pipelines developed as part of my M.Sc. Biotechnology & Bioinformatics coursework and academic projects at the Institute of Bioinformatics and Applied Biotechnology (IBAB), Bengaluru. The repository includes reproducible pipelines for genomics, transcriptomics, epigenomics, metagenomics, and sequence analysis using R, and widely used bioinformatics tools.

## Contents

| Folder | Analysis | Key Tools |
|---|---|---|
| `01_Microarray_Analysis` | Differential expression from Affymetrix microarray data (E-MTAB-6095: RNase I-deficient *E. coli* vs WT) | `oligo`, `limma` (R) |
| `02_RNAseq_Analysis` | Alignment, RPKM quantification, and DE analysis of RNA-Seq reads | Bowtie2, Rsubread, edgeR (R) |
| `03_Alternative_Splicing` | Splice-aware alignment and splice graph/junction analysis | HISAT2, SGSeq (R) |
| `04_Genome_Assembly` | De novo genome assembly, QC, gene prediction, scaffolding | Velvet, QUAST, Augustus, RagTag, BUSCO |
| `05_Variant_Calling` | SNP/Indel calling pipeline following GATK best practices | BWA-MEM, Picard, GATK |
| `06_ChIP_Seq` | Peak calling, annotation, and motif discovery | Bowtie2, MACS2, HOMER, BEDTools, MEME |
| `07_QIIME2_Metagenomics` | 16S amplicon microbial community & diversity analysis | QIIME 2 |
| `08_Bisulfite_Sequencing` | DNA methylation mapping at single-base resolution | Bismark, Bowtie2, samtools |

## Usage

Each script/folder assumes:
- Raw input files (FASTQ/CEL/BAM, reference genome, annotation) are placed in
  the working directory or a `data/` subfolder as indicated in the script.
- Required tools are installed and on `$PATH` (or, for R scripts, the listed
  Bioconductor/CRAN packages are installed).

Edit the variables at the top of each script (sample names, file paths,
reference genome) before running.

```bash
# Example: run the ChIP-Seq pipeline
cd 06_ChIP_Seq
chmod +x chipseq_pipeline.sh
./chipseq_pipeline.sh
```

```r
# Example: run the microarray DE analysis in R
setwd("01_Microarray_Analysis")
source("microarray_de_analysis.R")
```

## Notes

- Scripts are adapted from analyses originally run on lab/course datasets;
  paths have been generalized for portability. Substitute your own reference
  genome, annotation, and sample FASTQ/CEL files as needed.
- Each pipeline was documented step-by-step to support learning and
  reproducibility, not just execution — comments explain the purpose of each
  step (QC, alignment, normalization, statistical testing, visualization).
