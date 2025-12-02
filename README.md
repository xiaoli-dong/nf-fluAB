# nf-fluAB Pipeline

> 🕒 **Last updated:** April 14, 2025

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

>### Build Influenza A and B reference databases
>The following procedure will build the influenza A and B reference database. Then by using the reference fasta sequences generated, the script also generated the mash sketch database (msh file) and also the snpEff flu database
> ```
 >#create conda environment for creating database
 >mamba create -n '$env_name' mash=2.3 snpeff=5.2 vadr=1.6.4 biopython=1.84 entrez-direct=22.4 diamond=2.1.11 cd-hit=4.8.1 -y
 >#running the script
 >path_to_bin_directory/make_db.sh -i path_to/sequences.fasta -o outdir -c number_of_cpus -g path_to/BVBRC_genome.csv -d output_database_prefix
 >
 >```

---

## Quick Start
>If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with -profile test before running the workflow on actual data.

### Check pipeline command line options
You can download the nf-fluab from github to local computer or you can directly run the pipeline from github remotely. The following is an example how to check the command line options without downloading the pipeline locally:

```
# running directly from github without downloading or cloning
nextflow run xiaoli-dong/nf-fluab -r 7f72d6c --help
```
To run the analysis with your data, prepare a csv format samplesheet, which contains the sequenence information for each of your samples, as input. The samplesheet can contain only Illumina data or only nanopore data, it cannot not accept data from both Illumina and nanopore data in one analysis. See below for what the samplesheet looks like:

Illumina data analysis sample sheet example
```
sample,fastq_1,fastq_2,long_fastq
# comment lines will be ignoreed
sample1,sample1_R1.fastq.gz,sample1_R2.fastq.gz,NA
sample2,sample2_R1.fastq.gz,sample2_R2.fastq.gz,NA
sample3,sample3_R1.fastq.gz,sample3_R2.fastq.gz,NA
```

Nanopore data analysis sample sheet example
```
sample,fastq_1,fastq_2,long_fastq
# comment lines will be ignoreed
sample1,NA,NA,sample1.fastq.gz
sample2,NA,NA,sample2.fastq.gz
sample3,NA,NA,sample3.fastq.gz
```

Now, you can run the pipeline using:

```bash
nextflow run xiaoli-dong/nf-fluAB --input samplesheet.csv --outdir <OUTDIR> -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --platform <illumina or nanopore>
```
## Credits

- **nf-fluAB** was written by **Xiaoli Dong**.
- The Illumina part of the pipeline was primarily based on **Dr. Matthew Croxen**'s **flu pipeline**.
- Extensive support was provided by the **ProvLab Research Team** in Calgary for generating sequencing data and conducting the pipeline testing:
  >- **Linda Lee**
  >- **Johanna M Thayer**
  >- **Fitsum Getachew**
  >- **Petya Kolva**
  >- **Anita Wong**
  >- **Kanti Pabbaraju**
- **Dr. Tarah Lynch**, and **Dr. Matthew Croxen** for extensive technical inputs.
