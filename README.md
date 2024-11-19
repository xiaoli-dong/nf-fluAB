# nf-fluAB Pipeline

## Introduction

**nf-fluAB** is a bioinformatics pipeline for assembling and classifying influenza NGS sequence data generated using Illumina or Nanopore platforms. The pipeline is implemented with **Nextflow**, a workflow tool designed to run tasks across multiple compute infrastructures in a portable manner. It leverages **Docker** and **Singularity** containers, ensuring an easy installation process and highly reproducible results.

## Pipeline Summary

### nf-fluAB Pipeline Overview

The **nf-fluAB** pipeline processes NGS influenza A and B viral genome sequence data through the following steps:

### 1. Quality Control (QC) and De-hosting

The raw sequencing data undergoes quality control checks to assess its quality and remove low-quality reads, ensuring that only high-quality data is used for downstream analysis. The following tools are used for sequence quality control:

#### Short Reads:
- **Short Read Statistics**: `seqkit stats`  
  Collects basic statistics of short read sequences to assess quality.
  
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

Align quality-controlled short or long reads to the selected reference genome using a suitable read aligner. The resulting alignment files are processed to generate a clean mapping file.

- **Short Reads Mapping**: `bwa`, `minimap2`, `samtools`, `Picard`
- **Long Reads Mapping**: `minimap2`, `samtools`, `Picard`

### 4. Variant Calling and VCF File Preprocessing

Variants (mutations) are called from the aligned reads, identifying SNPs (single nucleotide polymorphisms) and indels (insertions and deletions) in the viral genome. The VCF files are normalized and filtered to remove variants that may cause frameshifts.

- **Short Read Variant Calling**: `freebayes`, `bcftools`, `snpEff`
- **Long Read Variant Calling**: `Clair3`, `bcftools`, `snpEff`

### 5. Consensus Sequence Generation

Generate consensus sequences using `bcftools` for each viral segment based on the variant-calling step, allowing the reconstruction of the full viral genome or specific segments.

### 6. Viral Segment Classification

Classify and annotate the individual viral segments (e.g., HA, NA, PB2, PB1, etc.) to help characterize the influenza strain.

- **Flu Typing**: `blastn` search against a flu typing database
- **Clade Assignment**: `nextclade`

### 7. Report Generation

Summarize the analysis and generate the following reports:
- **Analysis report**
- **Software version control report**

## Pipeline Reference Databases

### Build Influenza A and B fasta sequence database, mash sketch database, and snp database for snpEff 
1. **Download fasta format sequecne data** and save it as sequences.fasta: Public available fasta sequence data downloaded from: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=taxid:197911&VirusLineage_ss=taxid:197912&VirusLineage_ss=taxid:197913&VirusLineage_ss=taxid:1511083&LabHost_s=include
1. **Downlaod metadata** and save it as BVBRC_genome.csv from https://www.bv-brc.org/view/Taxonomy/11308#view_tab=genomes&filter=false 

### Typing database
### Nextclade database

---

## Quick Start

### 1. Install Nextflow (>=21.10.3)

Install **Nextflow** along with one of the following for full pipeline reproducibility:
- Docker
- Singularity
- Podman
- Shifter
- Charliecloud

You can follow the [Nextflow installation guide](https://www.nextflow.io/docs/latest/getstarted.html). Conda can also be used to install Nextflow and manage software within pipelines. However, it is recommended to use Conda as a last resort. For more details, see the [Nextflow documentation](https://www.nextflow.io/docs/latest/usage.html#containerization).

### 2. Download the Pipeline and Test It

To download and test the pipeline on a minimal dataset, run the following command:

```bash
nextflow run xiaoli-dong/nf-fluAB -profile test,YOURPROFILE --outdir <OUTDIR>
```
### 3. Customize Configuration

Configuration is required for Nextflow to know how to fetch the necessary software. This is typically done through a config profile (`YOURPROFILE`). You can chain multiple config profiles in a comma-separated string.

The pipeline includes several pre-configured profiles:

- `docker`
- `singularity`
- `podman`
- `shifter`
- `charliecloud`
- `conda`

For example, to run the pipeline with Docker, use the following command:

```bash
nextflow run xiaoli-dong/nf-fluAB -profile test,docker
```
### 4. Start Running Your Analysis

To run the pipeline with your own input data, use the following command:

```bash
nextflow run xiaoli-dong/nf-fluAB --input samplesheet.csv --outdir <OUTDIR> -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
```
## Credits

- **nf-fluAB** was written by **Xiaoli Dong**.
- The Illumina part of the pipeline was primarily based on Dr. **Matthew Croxen**'s **flu pipeline**.
- Extensive support was provided by the **ProvLab Research Team** in Calgary for genrating testing dataset, technical inputs.
