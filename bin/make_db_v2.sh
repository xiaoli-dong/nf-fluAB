#!/bin/bash

# ---------------------------------------------
# fluA, B sequence data, metadata and pipeline
# Sequence data, XML, and CSV files were downloaded from:
# https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=taxid:197911&VirusLineage_ss=taxid:197912&VirusLineage_ss=taxid:197913&VirusLineage_ss=taxid:1511083&LabHost_s=include on 2024-08-12
# Meta data downloaded from: https://www.bv-brc.org/view/Taxonomy/11308#view_tab=genomes&filter=false
# BVBRC_genome.csv (2024-08-23)

# ---------------------------------------------
# Usage function to display help message
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -i <input_fasta>    Input FASTA file (default: sequences.fasta)"
    echo "  -o <outdir>         Output directory (default: ./output)"
    echo "  -c <cpus>           Number of CPUs to use (default: 8)"
    echo "  -g <genome_csv>     Path to BVBRC genome CSV file (default: BVBRC_genome.csv)"
    echo "  -d <outdb_prefix>   Output database prefix (default: influenzaDB-20240823)"
    echo "  -h                  Show this help message"
    exit 1
}

# Default values
input_fasta="sequences.fasta"
outdir="./output"
cpus=8
genome_csv="BVBRC_genome.csv"
prefix="sequence"
outdb_prefix="influenzaDB-20240823"
bindir=$(dirname "$0")

# Parse command-line arguments
while getopts "i:o:c:g:d:h" opt; do
    case $opt in
        i) input_fasta="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        c) cpus="$OPTARG" ;;
        g) genome_csv="$OPTARG" ;;
        d) outdb_prefix="$OPTARG" ;;
        h|*) usage ;;
    esac
done

# ---------------------------------------------
# Conda Environment Handling

env_name="nf-fluab"

# Get environment paths using `conda info`
env_exists=$(conda info --envs | awk -v env="$env_name" '$1 == env {print "true"}')

# Check if the environment exists
if [ "$env_exists" == "true" ]; then
    echo "Environment '$env_name' exists."
else
    echo "Error: Environment '$env_name' does not exist."
    echo "Please using the following command to creating it before running the script..."
    echo "mamba create -n '$env_name' mmseqs2=15.6f452 mash=2.3 snpeff=5.2 vadr=1.6.4 biopython=1.84 -y"
    exit 0
fi
#pip install requests Bio
# Check if the environment is activated
if [[ "$CONDA_DEFAULT_ENV" != "$env_name" ]]; then
    echo "Activating the '$env_name' Conda environment..."
    eval "$(conda shell.bash hook)"
    conda activate "$env_name"
else
    echo "The '$env_name' Conda environment is already active."
fi

# ---------------------------------------------
# Directory and Output Setup

mkdir -p "$outdir"

# ---------------------------------------------
# Step 1: Filter and reformat FASTA Sequences
#filter sequeces: too short, too long, has ambiguour bases, has not metadata
echo "Filtering and reformatting FASTA sequences..."
if ! [ -f "${outdir}/${prefix}.reformat.fasta" ]; then
    
    python ${bindir}/reformatSeqs.py \
        --fasta sequences.fasta \
        --minlen 700 \
        --maxlen 3000 \
        --maxambigs 0 \
        --bvbrc BVBRC_genome.csv \
        > "${outdir}/${prefix}.reformat.fasta" 2> reformatSeq.log.txt
fi

# ---------------------------------------------
# Step 2: MMseqs2 Clustering
echo "Running MMseqs2 clustering..."
if ! [ -f "${outdir}/${prefix}.cluster0.99_rep_seq.fasta" ]; then
    mmseqs easy-cluster \
        --threads "$cpus" \
        -c 0.7 \
        --cov-mode 0 \
        --alignment-mode 3 \
        --min-seq-id 0.99 \
        --similarity-type 2 \
        --cluster-mode 2 \
        --cluster-reassign 1 \
        -v 2 \
        "${outdir}/${prefix}.reformat.fasta" \
        "${outdir}/${prefix}.cluster0.99" \
        tmp \
        >& "${outdir}/${prefix}.cluster0.99.log.txt"
fi

# python ${bindir}/filter_singleton.py \
#     --fasta output/sequence.cluster0.99_rep_seq.fasta \
#     --tsv output/sequence.cluster0.99_cluster.tsv \
#     > output/sequence.cluster0.99_rep_seq.singularity_out.fasta \
#     2> filter_singleton.log.txt

# ---------------------------------------------
# Step 3: Trim Terminal Ambiguities
echo "Trimming terminal ambiguities..."
if ! [ -f "${outdir}/${prefix}.vadr_trim.fasta" ]; then
    fasta-trim-terminal-ambigs.pl \
        "${outdir}/${prefix}.cluster0.99_rep_seq.fasta" \
        > "${outdir}/${prefix}.vadr_trim.fasta"
fi

# ---------------------------------------------
# Step 4: VADR Annotation
echo "Annotating sequences with VADR..."
if ! [ -f "${outdir}/${prefix}.annotate/${prefix}.annotate.vadr.pass.fa" ]; then
    v-annotate.pl \
        --cpu "$cpus" \
        --split -r --atgonly --xnocomp --nomisc \
        --alt_fail extrant5,extrant3 \
        --mkey flu \
        --mdir /nfs/APL_Genomics/db/prod/vadr/vadr-models-flu-1.6.3-2 \
        "${outdir}/${prefix}.cluster0.99_rep_seq.fasta" \
        "${outdir}/${prefix}.annotate"
fi

# ---------------------------------------------
# Step 5: Extract CDS Sequences
echo "Extracting CDS sequences..."
if ! [ -f "${outdir}/${prefix}.annotate.cds.fasta" ]; then
    python "${bindir}/extract_cds.py" \
        --fasta "${outdir}/${prefix}.annotate/${prefix}.annotate.vadr.pass.fa" \
        --bvbrc "$genome_csv" \
        --sgm "${outdir}/${prefix}.annotate/${prefix}.annotate.vadr.sgm" \
        > "${outdir}/${prefix}.annotate.cds.fasta"
fi

# ---------------------------------------------
# Step 6: Increase Clustering Coverage for CDS
echo "Increasing clustering coverage for CDS sequences..."
if ! [ -f "${outdir}/${outdb_prefix}_rep_seq.fasta" ]; then
    mmseqs easy-cluster \
        --threads "$cpus" \
        -c 0.9 \
        --cov-mode 0 \
        --alignment-mode 3 \
        --min-seq-id 0.97 \
        --similarity-type 2 \
        --cluster-mode 2 \
        --cluster-reassign 1 \
        -v 2 \
        "${outdir}/${prefix}.annotate.cds.fasta" \
        "${outdir}/${outdb_prefix}" \
        tmp
fi

mv ${outdir}/${outdb_prefix}_rep_seq.fasta ${outdir}/${outdb_prefix}.fasta
# ---------------------------------------------
# Step 7: Mash Sketch
echo "Performing Mash sketch..."
if ! [ -f "${outdir}/${outdb_prefix}.msh" ]; then
    mash sketch \
        -i "${outdir}/${outdb_prefix}.fasta" \
        -o "${outdir}/${outdb_prefix}.msh"
fi

if ! [ -d "${outdir}/${outdb_prefix}" ]; then
    mkdir ${outdir}/${outdb_prefix}
fi

cp ${outdir}/${outdb_prefix}.fasta ${outdir}/${outdb_prefix}/sequences.fasta
cp ${outdir}/${outdb_prefix}.msh ${outdir}/${outdb_prefix}/sequences.msh

# start to build database for snpEff
# ---------------------------------------------
# Step 8: Extract Genome IDs
echo "Extracting genome IDs..."
#if ! [ -f "${outdir}/genome_ids.txt" ]; then
grep ">" "${outdir}/${outdb_prefix}.fasta" | cut -f1 -d$' ' | cut -c 2- > "${outdir}/genome_ids.txt"
#fi

# ---------------------------------------------
# Step 9: Download GenBank Files
echo "Downloading GenBank files..."
if ! [ -d "${outdir}/gb_dir" ]; then
    mkdir "${outdir}/gb_dir"
fi
python "${bindir}/download_genbank.py" \
    --input "${outdir}/genome_ids.txt" \
    --outdir "${outdir}/gb_dir"


# ---------------------------------------------
# Step 10: Concatenate GenBank Files
echo "Concatenating GenBank files..."
if [ -f "${outdir}/data/fluab/genes.gbk" ]; then
    rm -rf "${outdir}/data/fluab/"
fi

#if ! [ -f "${outdir}/data/fluab/genes.gbk" ]; then
mkdir -p "${outdir}/data/fluab"
cat "${outdir}/gb_dir"/*.gb > "${outdir}/data/fluab/genes.gbk"
#fi

# ---------------------------------------------
# Step 11: Prepare snpEff Config
echo "Preparing snpEff config file..."
perl -ne '
    BEGIN { @ids = (); }
    chomp; push(@ids, $_);
    END {
        print "fluab.genome : influenza AB\n";
        print "\tfluab.chromosomes: ";
        print join(", ", @ids);
        print "\n";
    }
' "${outdir}/genome_ids.txt" > "${outdir}/snpEff.config"

# ---------------------------------------------
# Step 12: Convert Coding Sequences to GTF Format
echo "Converting coding sequences to GTF format..."
python "${bindir}/coding2gtf22.py" \
    --input "${outdir}/${outdb_prefix}.fasta" \
    --output "${outdir}/data/fluab/genes.gtf"

cp ${outdir}/${outdb_prefix}/sequences.fasta  ${outdir}/data/fluab/sequences.fa 
# ---------------------------------------------
# Step 13: Build SNP Database with snpEff
echo "Building SNP database with snpEff..."
snpEff build \
    -gtf22 -v \
    -c "${outdir}/snpEff.config" \
    -noCheckCds -noCheckProtein \
    -maxErrorRate 0.2 \
    fluab &> "${outdir}/log.gtf.txt"

# ---------------------------------------------
# Deactivate Conda Environment
echo "Deactivating nf-fluab environment..."
conda deactivate

# Final message
echo "Pipeline completed successfully."
