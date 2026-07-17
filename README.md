# NGS Analysis Pipelines

This repository contains bioinformatics workflows and analysis pipelines developed 
as part of my M.Sc. Biotechnology & Bioinformatics coursework and academic 
projects at the Institute of Bioinformatics and Applied Biotechnology (IBAB), Bengaluru.
The repository includes reproducible pipelines for genomics, transcriptomics, 
epigenomics, metagenomics, and sequence analysis using Python, R, and widely used bioinformatics tools.


## Contents

| Folder | Analysis | Language | Key Tools |
|---|---|---|---|
| `01_Microarray_Analysis` | Differential expression from Affymetrix microarray data (E-MTAB-6095: RNase I-deficient *E. coli* vs WT) | R | `oligo`, `limma` |
| `02_RNAseq_Analysis` | (1) Alignment, (2) RPKM quantification and DE analysis | Python + R | Bowtie2 · Rsubread, edgeR |
| `03_Alternative_Splicing` | (1) Splice-aware alignment, (2) splice graph/junction analysis | Python + R | HISAT2 · SGSeq |
| `04_Genome_Assembly` | De novo genome assembly, QC, gene prediction, scaffolding | Python | Velvet, QUAST, Augustus, RagTag, BUSCO |
| `05_Variant_Calling` | SNP/Indel calling pipeline following GATK best practices | Python | BWA-MEM, Picard, GATK |
| `06_ChIP_Seq` | Peak calling, annotation, and motif discovery | Python | Bowtie2, MACS2, HOMER, BEDTools, MEME |
| `07_QIIME2_Metagenomics` | 16S amplicon microbial community & diversity analysis | Python | QIIME 2 |
| `08_Bisulfite_Sequencing` | DNA methylation mapping at single-base resolution | Python | Bismark, Bowtie2, samtools |

### Why Python for some pipelines and R for others?

Tools like Bowtie2, samtools, GATK, MACS2, HOMER, bedtools, MEME, QIIME2, and
Bismark are standalone command-line programs, there's no R or Python package
that reimplements them, so any script that runs them is just calling out to
the same binaries. Rather than raw shell scripts, those pipelines are written
in **Python**, using the `subprocess` module to call each tool, this keeps
the orchestration logic readable, testable, and easy to extend (e.g. adding
loops, config files, logging) compared to bash.

The statistics, heavy steps (differential expression, splice graph analysis)
use Bioconductor packages: `limma`, `edgeR`, `SGSeq`, that only exist in
**R**. There's no equivalent Python library that reproduces the same
statistical methods, so those stay in R rather than being awkwardly
reimplemented.

## Usage

Each script/folder assumes:
- Raw input files (FASTQ/CEL/BAM, reference genome, annotation) are placed in
  the working directory or a `data/` subfolder as indicated in the script.
- Required command-line tools are installed and on `$PATH`.
- For R scripts, the listed Bioconductor/CRAN packages are installed.

Edit the variables at the top of each script (sample names, file paths,
reference genome) before running.

```bash
# Example: run the ChIP-Seq pipeline
cd 06_ChIP_Seq
python3 chipseq_pipeline.py
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
