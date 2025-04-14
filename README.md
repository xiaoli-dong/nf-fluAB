# nf-fluAB Pipeline

## Introduction

**nf-fluAB** is a bioinformatics pipeline for the assembly, typing, and lineage assignment of influenza NGS data generated from Illumina or Nanopore platforms. Built with **Nextflow**,it enables portable and scalable execution across a range of computing environments. The use of Docker and Singularity containers ensures easy installation and highly reproducible results.


## Pipeline Summary
The pipeline takes a samplesheet and corresponding FASTQ files as input. It performs quality control (QC), identifies the closest publicly available reference(s), maps reads, calls variants, and generates a consensus contig for each segment. In cases of co-infection, multiple references are selected, and multiple consensus contigs are generated per segment. Segments are typed against a local database, and lineage is determined using Nextclade. At the end of the workflow, a comprehensive master summary report is generated for each sample, along with a consolidated analysis overview.

![Pipeline Diagram](assets/nf-fluab-drawio.svg)


* QC
  * Illumina QC ([fastp](https://github.com/OpenGene/fastp) or[BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) -> [seqkit stats](https://bioinf.shenwei.me/seqkit/usage/#stats))
  * Nanopore QC ([Porechop](https://github.com/rrwick/Porechop) -> [chopper](https://github.com/wdecoster/chopper) -> [seqkit stats](https://bioinf.shenwei.me/seqkit/usage/#stats))
* dehost ([hostile](https://github.com/bede/hostile) -> [seqkit stats](https://bioinf.shenwei.me/seqkit/usage/#stats))
* seek references ([mash screen](https://github.com/marbl/Mash) -> filter mash screen output)
* mapping and post-processing bam files
  * Illumina mapping ([bwa](https://github.com/lh3/bwa) or [minimap2](https://github.com/lh3/minimap2) -> [samtools sort, index](https://www.htslib.org/doc/samtools.html) -> [picard MARKDUPLICATES](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) -> [samtools coverage](https://www.htslib.org/doc/samtools-coverage.html) -> bedtools genomecov)
  * Nanopore mapping (minimap2 -> samtools sort, index -> samtools coverage -> samtools coverage -> [bedtools genomecov](https://bedtools.readthedocs.io/en/latest/content/overview.html) )
* variant calling and post-processing vcf files
  * Illumina data variant calling and post-processing ([freebayes](https://github.com/freebayes/freebayes) or [bcftools](https://samtools.github.io/bcftools/bcftools.html) -> bcftools sort, index -> bcftools norm -> bcftools filter (low quality, low depth) -> [snpeff](https://pcingola.github.io/SnpEff/) -> bcftools filter (frameshift)
  * Nanopore data variant calling and post-processing (clair3 -> bcftools sort, index -> bcftools filter (low quality, low depth) -> snpeff -> bcftools filter (frameshift)
* consensus calling (bcftools consensus -> SEQKIT fx2tab (stats) 
* consensus typing (blastn against typing database)
* Lineage determination (nextclade)
* Summary report

## Pipeline required reference sequences and databases
1. [Influenza A and B  primer sequences](assets/flu-primers.fa): used for sequence data qulaity control process
2. [flu typing database](assets/typing.fa): used for assembled segment typing 
3. [nextcade dataset](https://github.com/nextstrain/nextclade_data/tree/master/data/nextstrain/flu): used for flu A&B lineage determination
4. Influenza A and B reference databases: used for seeking the closest related public available sequences and used as the references in the assembly. [Go to reference database building guide](#build-influenza-a-and-b-reference-databases)

### Build Influenza A and B reference databases
The following procedure will build the influenza A and B reference database. Then by using hte reference fasta sequences generated, the script also generated the mash sketch database (msh file) and also the snpEff flu database

 1. Public available fasta sequence data is downloaded from [Influenza Virus Data Hub](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=taxid:197911&VirusLineage_ss=taxid:197912&VirusLineage_ss=taxid:197913&VirusLineage_ss=taxid:1511083&LabHost_s=include) and saved as sequences.fasta
 1. Metadata is downloaded from [BV-BRC](https://www.bv-brc.org/view/Taxonomy/11308#view_tab=genomes&filter=false) and saved as BVBRC_genome.csv
 1. run the following command to build the database
   ```
   #create conda environment for creating database
   mamba create -n '$env_name' mash=2.3 snpeff=5.2 vadr=1.6.4 biopython=1.84 entrez-direct=22.4 diamond=2.1.11 cd-hit=4.8.1 -y
   #running the script
   path_to_bin_directory/make_db.sh -i path_to/sequences.fasta -o outdir -c number_of_cpus -g path_to/BVBRC_genome.csv -d output_database_prefix
   
   ```

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
