#!/bin/bash

# ---------------------------------------------
# fluA, B sequence data, metadata and pipeline
# Sequence data, XML, and CSV files were downloaded from:
# https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=taxid:197911&VirusLineage_ss=taxid:197912&VirusLineage_ss=taxid:197913&VirusLineage_ss=taxid:1511083&LabHost_s=include on 2024-08-12
# Meta data downloaded from: https://www.bv-brc.org/view/Taxonomy/11308#view_tab=genomes&filter=false
# example cmd: sh make_db.sh -i sequences_nt.20250310.fasta -g BVBRC_genome.20250310.csv -d influenzaDB-20250310
#sh ../bin/make_db.sh -i sequences_nt.20250331.fasta -o ./output -c 8 -g BVBRC_genome.20250331.csv -d influenzaDB-20250331
#sh ../bin/make_db.sh -i sequences_nt.20250310.fasta -o ./output -c 8 -g BVBRC_genome.20250310.csv -d influenzaDB-20250310 >& log.txt &
# ---------------------------------------------
# Usage function to display help message
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -i <input_fasta>    Input FASTA file (required)"
    echo "  -o <outdir>         Output directory (default: ./output)"
    echo "  -c <cpus>           Number of CPUs to use (default: 8)"
    echo "  -g <genome_csv>     Path to BVBRC genome CSV file (required)"
    echo "  -d <outdb_prefix>   Output database prefix (default:influenzaDB)"
    echo "  -h                  Show this help message"
    echo "Example:"
    echo "  $0 -i sequences_nt.20250310.fasta -o ./output -c 8 -g BVBRC_genome.20250310.csv -d influenzaDB-20250310"
    exit 1
}

# Default values
input_fasta="sequences.fasta"
outdir="./output"
cpus=8
genome_csv="BVBRC_genome.csv"
prefix="sequences"
outdb_prefix="influenzaDB"
gb_dir="./gb_dir"

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

# Check if both arguments are provided
if [ -z "$input_fasta" ] || [ -z "$genome_csv" ]; then
    echo "Both -a and -b options are required."
    usage
fi
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
    echo "mamba create -n '$env_name' mash=2.3 snpeff=5.2 vadr=1.6.4 biopython=1.84 entrez-direct=22.4 diamond=2.1.11 cd-hit=4.8.1 -y"
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
# Step 1: Filter and reformat FASTA Sequences filter out:
# parital sequecnes
# sequences containing ambiguour bases
# do not have metadata
echo "Filtering and reformatting FASTA sequences..."
cmd="python ${bindir}/reformatSeqs.py \
        --fasta $input_fasta \
        --bvbrc $genome_csv \
        > ${prefix}.reformat.fasta 2> reformatSeq.log.txt"
echo $cmd

if ! [ -f "${prefix}.reformat.fasta" ]; then
    cmd_to_run=$(echo $cmd)
    eval $cmd_to_run
fi

# ---------------------------------------------
# Step 2.1: MMseqs2 Clustering
rep_fasta="${outdir}/${prefix}.cluster0.99_rep_seq.fasta"
rep_clstr="${outdir}/${prefix}.cluster0.99_rep_seq.fasta.clstr"
echo "Running cd-hit-est clustering..."
cmd="cd-hit-est \
        -i ${prefix}.reformat.fasta \
        -o ${rep_fasta} \
        -c  0.99 -d 0  -g 1 -M 0 -T 8 -p 1 -sc 1 \
        2> ${outdir}/${prefix}.cluster0.99.log.txt"
echo $cmd
if ! [ -f "${rep_fasta}" ]; then

    cmd_to_run=$(echo $cmd)
    eval $cmd_to_run

fi

# Step 2.2:
echo "Running pick repseq..."
rep_picked_fasta="${outdir}/${prefix}.cluster0.99_rep_picked.fasta"
rep_picked_clstr="${rep_picked_fasta}.clstr"

cmd="python ${bindir}/pick_rep_from_cdhit.py \
        -f ${prefix}.reformat.fasta \
        -c ${rep_clstr} \
        --out_clstr ${rep_picked_clstr} \
        --out_fasta ${rep_picked_fasta} \
        2> ${outdir}/pick_rep.0.99.stderr.txt"
echo $cmd

if ! [ -f "${outdir}/${prefix}.cluster0.99_rep_picked.fasta" ]; then
    cmd_to_run=$(echo $cmd)
    eval $cmd_to_run
fi

# ---------------------------------------------
# Step 3: Trim Terminal Ambiguities
echo "Trimming terminal ambiguities..."
cmd="fasta-trim-terminal-ambigs.pl --minlen 500 --maxlen 10000 ${rep_picked_fasta} > ${outdir}/${prefix}.vadr_trim.fasta"
echo $cmd

if ! [ -f "${outdir}/${prefix}.vadr_trim.fasta" ]; then
    cmd_to_run=$(echo $cmd)
    eval $cmd_to_run
fi

# ---------------------------------------------
# Step 4: VADR Annotation
echo "Annotating sequences with VADR..."
cmd="v-annotate.pl \
        --cpu $cpus \
        --split -r --atgonly --xnocomp --nomisc \
        --alt_fail extrant5,extrant3 \
        --mkey flu \
        --mdir /nfs/APL_Genomics/db/prod/vadr/vadr-models-flu-1.6.3-2 \
        ${outdir}/${prefix}.vadr_trim.fasta \
        ${outdir}/${prefix}.annotate"
echo $cmd

if ! [ -f "${outdir}/${prefix}.annotate/${prefix}.annotate.vadr.pass.fa" ]; then
    cmd_to_run=$(echo $cmd)
    eval $cmd_to_run
fi

# ---------------------------------------------
# Step 5: Extract CDS Sequences
echo "Extracting CDS sequences..."
cmd="python ${bindir}/extract_cds.py \
        --fasta ${outdir}/${prefix}.annotate/${prefix}.annotate.vadr.pass.fa \
        --bvbrc "$genome_csv" \
        --sgm ${outdir}/${prefix}.annotate/${prefix}.annotate.vadr.sgm \
        > ${outdir}/${prefix}.annotate.cds.fasta"
echo $cmd

if ! [ -f "${outdir}/${prefix}.annotate.cds.fasta" ]; then
    cmd_to_run=$(echo $cmd)
    eval $cmd_to_run
fi

# ---------------------------------------------
# Step 6: Increase Clustering Coverage for CDS
echo "Increasing clustering coverage for CDS sequences..."
rep_fasta="${outdir}/cds.cluster0.97_rep_seq.fasta"
rep_clstr="${outdir}/cds.cluster0.97_rep_seq.fasta.clstr"

cmd="cd-hit-est -c 0.97 -d 0 -g 1 -M 0 -T 8 -s 0.9 -aL 0.9 -aS 0.9 -p 1 -sc 1 \
        -i ${outdir}/${prefix}.annotate.cds.fasta \
        -o ${rep_fasta}"
echo $cmd

if ! [ -f "${rep_fasta}" ]; then
    cmd_to_run=$(echo $cmd)
    eval $cmd_to_run
fi

echo "Running pick repseq..."
rep_picked_fasta="${outdir}/cds.cluster0.97_rep_picked.fasta"
rep_picked_clstr="${rep_picked_fasta}.clstr"

cmd="python ${bindir}/pick_rep_from_cdhit.py \
        -f ${outdir}/${prefix}.annotate.cds.fasta \
        -c ${rep_clstr} \
        --out_clstr ${rep_picked_clstr} \
        --out_fasta ${rep_picked_fasta} \
        2> ${outdir}/pick_rep.0.97.stderr.txt"
echo $cmd

if ! [ -f "${rep_picked_fasta}" ]; then
    cmd_to_run=$(echo $cmd)
    eval $cmd_to_run
fi


# downloading refseq viral protein db
if [ ! -d "./viral_protein_dir" ]; then
  echo "Directory does not exist."
  mkdir -p ./viral_protein_dir
else
  echo "./viral_protein_dir Directory exists."
fi

echo "Getting influenza protein sequences from refseq to correct some segment misassignments in NCBI..."
cmd="wget -P ./viral_protein_dir https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz"
echo $cmd
if ! [ -f "./viral_protein_dir/viral.1.protein.faa.gz" ]; then
    cmd_to_run=$(echo $cmd)
    eval $cmd_to_run
fi

cmd="zcat ./viral_protein_dir/viral.1.protein.faa.gz | grep -E \"^(>.*(Influenza A|Influenza B).*)\" > ./viral_protein_dir/fluab_protein_ids.txt"
echo $cmd
if ! [ -f "./viral_protein_dir/fluab_protein_ids.txt" ]; then
    cmd_to_run=$(echo $cmd)
    eval $cmd_to_run
fi

cmd="zcat ./viral_protein_dir/viral.1.protein.faa.gz \
    | perl ${bindir}/getSeqs.pl \
    -s - -t file \
    -q ./viral_protein_dir/fluab_protein_ids.txt > ./viral_protein_dir/fluab_refseq_protein.fasta"
echo $cmd
if ! [ -f "./viral_protein_dir/fluab_refseq_protein.fasta" ]; then
    cmd_to_run=$(echo $cmd)
    eval $cmd_to_run
fi

cmd="diamond makedb --in ./viral_protein_dir/fluab_refseq_protein.fasta -d ./viral_protein_dir/fluab_refseq_protein"
echo $cmd
if ! [ -f "./viral_protein_dir/fluab_refseq_protein.dmnd" ]; then
    cmd_to_run=$(echo $cmd)
    eval $cmd_to_run
fi

echo "correcting segment misassignments in NCBI..."

cmd="diamond blastx -d ./viral_protein_dir/fluab_refseq_protein \
    -q  ${rep_picked_fasta} \
    -o  ${rep_picked_fasta}.blastx_fmt6.tsv \
    --evalue 1e-5 -k 1 --outfmt 6"
echo $cmd

if ! [ -f " ${rep_picked_fasta}.blastx_fmt6.tsv" ]; then
    cmd_to_run=$(echo $cmd)
    eval $cmd_to_run
fi
echo "validating desc ..."
cmd="python ${bindir}/validate_fludb_desc.py \
    ./viral_protein_dir/fluab_refseq_protein.fasta \
    ${rep_picked_fasta}.blastx_fmt6.tsv \
    ${bindir}/flu_protein_to_segment.csv \
    ${rep_picked_fasta} \
    >  ${outdir}/${outdb_prefix}.fasta"


echo $cmd

if ! [ -f "${outdir}/${outdb_prefix}.fasta" ]; then
    cmd_to_run=$(echo $cmd)
    eval $cmd_to_run
fi


#mv ${outdir}/${outdb_prefix}_rep_seq.fasta ${outdir}/${outdb_prefix}.fasta
# ---------------------------------------------
# Step 7: Mash Sketch
echo "Performing Mash sketch..."
cmd="mash sketch \
        -i ${outdir}/${outdb_prefix}.fasta \
        -o ${outdir}/${outdb_prefix}.msh"
echo $cmd

if ! [ -f "${outdir}/${outdb_prefix}.msh" ]; then
    cmd_to_run=$(echo $cmd)
    eval $cmd_to_run
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
if ! [ -d "${gb_dir}" ]; then
    mkdir "./${gb_dir}"
fi

cmd="python ${bindir}/download_genbank.py \
    --input ${outdir}/genome_ids.txt \
    --outdir ./${gb_dir}"

cmd_to_run=$(echo $cmd)
eval $cmd_to_run

# ---------------------------------------------
# Step 10: Concatenate GenBank Files
echo "Concatenating GenBank files..."
if [ -f "${outdir}/snpeff/data/fluab/genes.gbk" ]; then
    rm -rf "${outdir}/snpeff/data/fluab/"
fi

#if ! [ -f "${outdir}/data/fluab/genes.gbk" ]; then
mkdir -p "${outdir}/snpeff/data/fluab"
cat "./gb_dir"/*.gb > "${outdir}/snpeff/data/fluab/genes.gbk"
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
' "${outdir}/genome_ids.txt" > "${outdir}/snpeff/snpEff.config"

# ---------------------------------------------
# Step 12: Convert Coding Sequences to GTF Format
echo "Converting coding sequences to GTF format..."
cmd="python ${bindir}/coding2gtf22.py \
    --input ${outdir}/${outdb_prefix}.fasta \
    --output ${outdir}/snpeff/data/fluab/genes.gtf"

cmd_to_run=$(echo $cmd)
eval $cmd_to_run

cp ${outdir}/${outdb_prefix}/sequences.fasta  ${outdir}/snpeff/data/fluab/sequences.fa
# ---------------------------------------------
# Step 13: Build SNP Database with snpEff
echo "Building SNP database with snpEff..."
cmd="snpEff build \
    -gtf22 -v \
    -c ${outdir}/snpeff/snpEff.config \
    -noCheckCds -noCheckProtein \
    -maxErrorRate 0.2 \
    fluab 2> ${outdir}/log.gtf.txt"

cmd_to_run=$(echo $cmd)
eval $cmd_to_run
# ---------------------------------------------
# Deactivate Conda Environment
echo "Deactivating nf-fluab environment..."
conda deactivate

# Final message
echo "Pipeline completed successfully."
