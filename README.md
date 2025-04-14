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

>### Build Influenza A and B reference databases
>The following procedure will build the influenza A and B reference database. Then by using hte reference fasta sequences generated, the script also generated the mash sketch database (msh file) and also the snpEff flu database
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
You can clone or download the nf-fluab from github to local computer or you can directly run the pipeline from github. To check the pipeline command line options:
```
# running directly from github without downloading or cloning
nextflow run xiaoli-dong/nf-fluab -r 7f72d6c --help
N E X T F L O W  ~  version 23.04.1
Launching `https://github.com/xiaoli-dong/nf-fluab` [pensive_hugle] DSL2 - revision: 7f72d6cefbace3f02b7c52a2421cbd08f15692f5


------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
  nf-core/influenza v1.0.2-g7f72d6c
------------------------------------------------------
Typical pipeline command:

  nextflow run nf-core/influenza --input samplesheet.csv --genome GRCh37 -profile docker

Input/output options
  --input                      [string]  Path to comma-separated file containing information about the samples in the experiment.
  --outdir                     [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud 
                                         infrastructure. 
  --platform                   [string]  Specifies the sequencing platform of the input reads - available options are 'illumina|nanopore'. [default: 
                                         illumina] 
  --email                      [string]  Email address for completion summary.

qc_options
  --skip_illumina_reads_qc     [boolean] skip illumina read qc step [default: false]
  --illumina_reads_qc_tool     [string]  illumina read quality processing tool, the available options are 'fastp|bbduk'. [default: fastp]
  --flu_primers                [string]  flu sequencing primer [default: /nfs/APL_Genomics/db/prod/fluAB/flu-primers.fa]
  --hostile_human_ref_bowtie2  [string]  hostile human genome index file [default: 
                                         /nfs/APL_Genomics/db/prod/hostile/bowtie2_indexes/human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401] 
  --hostile_human_ref_minimap2 [string]  hostile human reference genome for minimap2 [default: 
                                         /nfs/APL_Genomics/db/prod/hostile/minimap2_ref/human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401.fa.gz] 

retrieve_reference_options
  --flu_db_msh                 [string]  flu database mash sketch file [default: /nfs/APL_Genomics/db/prod/fluAB/influenzaDB/sequences.msh]
  --flu_db_fasta               [string]  Fasta format flu database file [default: /nfs/APL_Genomics/db/prod/fluAB/influenzaDB/sequences.msh]
  --mashthreshold              [number]  mash screen minimum identity to report [default: 0.9]
  --max_p_value                [number]  mash screen Maximum p-value to report. [default: 0.1]

reference_based_assembly_options
  --mapping_tool               [string]  illumina read mapping tool. The available options are: 'minimap2|bwa'. [default: minimap2]
  --variant_caller             [string]  variant caller. The available options are 'bcftools|freebayes|clair3' [default: bcftools]
  --mindepth                   [integer] require at least this depth to process a site for variants. [default: 10]
  --lower_ambiguity_freq       [number]  lowest alt frequency for a site to be considered in generating consensus. [default: 0.25]
  --upper_ambiguity_freq       [number]  upper alt frequency for a site to be considered in generating consensus site as iupac. [default: 0.75]

annotation_options
  --typing_db                  [string]  typing database [default: /nfs/APL_Genomics/db/prod/fluAB/typing.fa]
  --nextclade_dataset_base     [string]  nextclade database base directory [default: /nfs/APL_Genomics/db/prod/fluAB/nextclade]
  --minblastident              [integer] blastn search typing_db minimum percent identity. [default: 70]
  --minblastcov                [integer] blastn search agaist typing_db percent minimum  query coverage per hsp. [default: 75]

!! Hiding 19 params, use --show_hidden_params to show them !!
------------------------------------------------------
If you use nf-core/influenza for your analysis please cite:

* The nf-core framework
  https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
  https://github.com/nf-core/influenza/blob/master/CITATIONS.md
------------------------------------------------------
```
To download and test the pipeline on a minimal dataset, run the following command:

```bash
nextflow run xiaoli-dong/nf-fluAB -profile test,YOURPROFILE --outdir <OUTDIR>
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
