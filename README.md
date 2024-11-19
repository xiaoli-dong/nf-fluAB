# nf-fluAB Pipeline

## Introduction

**nf-fluAB** is a bioinformatics pipeline for assembling and classifying influenza NGS sequence data generated using Illumina or Nanopore platforms. The pipeline is implemented using **Nextflow**, a workflow tool that enables running tasks across multiple compute infrastructures in a highly portable manner. It leverages **Docker** and **Singularity** containers, ensuring a trivial installation process and highly reproducible results.

## Pipeline Summary

### nf-fluAB Pipeline Overview

The **nf-fluAB** pipeline processes NGS influenza A and B viral genome sequence data through the following steps:

### 1. Quality Control (QC) and De-hosting

Perform quality control checks on the raw sequencing data to assess its quality and remove low-quality reads, ensuring that only high-quality data is used for downstream analysis. The following tools are used for sequence quality control:

#### Short Reads:
- **Short Read Statistics**: `seqkit stats`  
  Collects basic statistics of the short read sequences to assess quality.
  
- **Short Read Quality Control**: `fastp` or `bbduk`  
  Performs adapter trimming, quality filtering, and base correction of the short reads.

#### Long Reads:
- **Nanopore Long Read Adapter Trimming**: `Porechop`  
  Removes adapter sequences from Nanopore long reads.
  
- **Nanopore Long Read Quality and Length Filter**: `chopper`  
  Filters long reads based on quality and length thresholds.

- **Nanopore Long Read Statistics**: `seqkit stats`  
  Provides summary statistics on long read data for quality assessment.

#### De-hosting:
- **Remove Host Sequences**: `Hostile`  
  Removes contaminating host DNA sequences from both short and long read datasets.

### 2. Reference Genome Search (Mash)

Quality-controlled reads are "screened" against the influenza database using **Mash** to identify the most closely related reference genomes from the database.

### 3. NGS Reads Mapping and Mapping File Preprocessing

Align quality controlled short or long reads to the selected reference genome using a suitable read aligner and Process the alignment files to generate a clean mapping file.

- **Short Reads Mapping**: `bwa`, `minimap2`, `samtools`, `Picard`
- **Long Reads Mapping**: `minimap2`, `samtools`, `Picard`

### 4. Variant Calling and VCF File Preprocessing

 Variants (mutations) are then called from the aligned reads, identifying SNPs (single nucleotide polymorphisms) and indels (insertions and deletions) in the viral genome. normalize the vcf files and filter out the variants which can cause frameshift

- **Short Read Variant Calling**: `freebayes`, `bcftools`, `snpEff`
- **Long Read Variant Calling**: `Clair3`, `bcftools`, `snpEff`

### 5. Consensus Sequence Generation

Generate consensus sequences for each viral segment based on the variant-calling step, allowing the reconstruction of the full viral genome or specific segments.

- **Consensus Generation**: `bcftools`

### 6. Viral Segment Classification

Classify and annotate the individual viral segments (e.g., HA, NA, PB2, PB1, etc.)

- **Flu Typing**: `blastn` search against a flu typing database
- **Clade Assignment**: `nextclade`

### 7. Report Generation

Summarize the analysis and generate the following reports:
- **Analysis report**
- **Software version control report**

### Pipeline Reference Databases
- **Flu database**
- **Typing database**
- **Nextclade database**

---

## Quick Start

### 1. Install Nextflow (>=21.10.3)

Install **Nextflow** and any of the following for full pipeline reproducibility:
- Docker
- Singularity
- Podman
- Shifter
- Charliecloud

You can follow the [Nextflow installation guide](https://www.nextflow.io/docs/latest/getstarted.html). Conda can also be used to install Nextflow and manage software within pipelines. However, using Conda is recommended only as a last resort; see the [docs](https://www.nextflow.io/docs/latest/usage.html#containerization) for more details.

### 2. Download the Pipeline and Test It

To download and test the pipeline on a minimal dataset, run the following command:

```bash
nextflow run xiaoli-dong/nf-fluAB -profile test,YOURPROFILE --outdir <OUTDIR>
