#!/bin/bash

<<comment
# fluA, B
sequence data, xml, and csv file downloaded from: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=taxid:197911&VirusLineage_ss=taxid:197912&VirusLineage_ss=taxid:197913&VirusLineage_ss=taxid:1511083&LabHost_s=include on 2024-08-12

#meta data downloaded from: 2024-08-23
https://www.bv-brc.org/view/Taxonomy/11308#view_tab=genomes&filter=false
BVBRC_genome.csv
comment

#example:
#./makedb.sh -i sequences.fasta -o ./my_outdir -c 16 -g BVBRC_genome.csv -d custom_outdb_prefix
# Function to display help message
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
        i)
            input_fasta="$OPTARG"
            ;;
        o)
            outdir="$OPTARG"
            ;;
        c)
            cpus="$OPTARG"
            ;;
        g)
            genome_csv="$OPTARG"
            ;;
        d)
            outdb_prefix="$OPTARG"
            ;;
        h|*)
            usage
            ;;
    esac
done

# Create output directory if it doesn't exist
mkdir -p "$outdir"

# Define paths for Singularity containers and commands
mmseq2="singularity run /nfs/APL_Genomics/apps/production/singularity_home/downloaded_images/mmseqs2:15.6f452--pl5321h6a68c12_3 mmseqs"
mash="singularity run /nfs/APL_Genomics/apps/production/NXF_HOME/NXF_SINGULARITY_CACHEDIR/depot.galaxyproject.org-singularity-mash-2.3--he348c14_1.img mash"
snpEff="singularity run /nfs/APL_Genomics/apps/production/NXF_HOME/NXF_SINGULARITY_CACHEDIR/depot.galaxyproject.org-singularity-snpeff%3A5.2--hdfd78af_1.img snpEff"

# Filter FASTA sequences: Minimum length 700, Maximum length 3000, Max ambiguities 0
echo "Filtering FASTA sequences..."
if ! [ -f "${outdir}/${prefix}.filter.fasta" ]; then
    python ${bindir}/filter_fasta.py \
        --fasta "$input_fasta" \
        --minlen 700 \
        --maxlen 3000 \
        --maxambigs 0 \
        > "${outdir}/${prefix}.filter.fasta"
fi

# Reduce dataset for running MMseqs2 clustering
echo "Running MMseqs2 clustering..."
if ! [ -f "${outdir}/${prefix}.cluster0.99_rep_seq.fasta" ]; then
    ${mmseq2} easy-cluster \
        --threads "${cpus}" \
        -c 0.7 \
        --cov-mode 0 \
        --alignment-mode 3 \
        --min-seq-id 0.99 \
        --similarity-type 2 \
        --cluster-mode 2 \
        --cluster-reassign 1 \
        -v 2 \
        "${outdir}/${prefix}.filter.fasta" \
        "${outdir}/${prefix}.cluster0.99" \
        tmp \
        >& "${outdir}/${prefix}.cluster0.99.log.txt"
fi

# Activate VADR Conda environment
echo "Activating VADR environment..."
# ACTIVATE conda
eval "$(conda shell.bash hook)"
conda activate

conda activate vadr

# Trim terminal ambiguities from the clustered sequences
echo "Trimming terminal ambiguities..."
if ! [ -f "${outdir}/${prefix}.vadr_trim.fasta" ]; then
    fasta-trim-terminal-ambigs.pl \
        "${outdir}/${prefix}.cluster0.99_rep_seq.fasta" \
        > "${outdir}/${prefix}.vadr_trim.fasta"
fi
# Annotate sequences with VADR
echo "Annotating sequences with VADR..."
if ! [ -f "${outdir}/${prefix}.annotate/${prefix}.annotate.vadr.pass.fa" ]; then
    v-annotate.pl \
        --cpu "${cpus}" \
        --split -r --atgonly --xnocomp --nomisc \
        --alt_fail extrant5,extrant3 \
        --mkey flu \
        --mdir /nfs/APL_Genomics/db/prod/vadr/vadr-models-flu-1.6.3-2 \
        "${outdir}/${prefix}.cluster0.99_rep_seq.fasta" \
        "${outdir}/${prefix}.annotate"
fi

# Deactivate Conda environment
echo "Deactivating VADR environment..."
conda deactivate

# Extract CDS from the sequences annotated by VADR
echo "Extracting CDS sequences..."
if ! [ -f "${outdir}/${prefix}.annotate.cds.fasta" ]; then
    python ${bindir}/extract_cds.py \
        --fasta "${outdir}/${prefix}.annotate/${prefix}.annotate.vadr.pass.fa" \
        --bvbrc "$genome_csv" \
        --sgm "${outdir}/${prefix}.annotate/${prefix}.annotate.vadr.sgm" \
        > "${outdir}/${prefix}.annotate.cds.fasta"
fi

# Increase clustering coverage for CDS sequences using MMseqs2
echo "Increasing clustering coverage for CDS sequences..."
if ! [ -f "${outdir}/${outdb_prefix}_rep_seq.fasta" ]; then
    ${mmseq2} easy-cluster \
        --threads "${cpus}" \
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

# Perform Mash sketch on the clustered representative sequences
echo "Performing Mash sketch..."
if ! [ -f "${outdir}/${outdb_prefix}_rep_seq.msh" ]; then
    ${mash} sketch \
        -i "${outdir}/${outdb_prefix}_rep_seq.fasta" \
        -o "${outdir}/${outdb_prefix}_rep_seq.msh"
fi

# Get genome IDs of the representative sequences
echo "Extracting genome IDs..."
if ! [ -f "${outdir}/genome_ids.txt" ]; then
    grep ">" "${outdir}/${outdb_prefix}_rep_seq.fasta" | cut -f1 -d$' ' | cut -c 2- > "${outdir}/genome_ids.txt"
fi

# Download GenBank files for the genomes using the extracted IDs
echo "Downloading GenBank files..."
if ! [ -d "${outdir}/gb_dir" ]; then
    mkdir "${outdir}/gb_dir"
    python ${bindir}/download_genbank.py \
        --input "${outdir}/genome_ids.txt" \
        --outdir "${outdir}/gb_dir"
fi

# Concatenate GenBank files into one
echo "Concatenating GenBank files..."
if ! [ -f "${outdir}/data/influenza_a_b/genes.gbk" ]; then
    cat "${outdir}/gb_dir"/*.gb > "${outdir}/data/influenza_a_b/genes.gbk"
fi

# Prepare snpEff config file for annotation
echo "Preparing snpEff config file..."
perl -ne '
    BEGIN {
        @ids = ();
    }
    chomp;
    push(@ids, $_);
    END {
        print "influenza_a_b.genome : influenza AB\n";
        print "\tinfluenza_a_b.chromosomes: ";
        print join(", ", @ids);
        print "\n";
    }
' "${outdir}/genome_ids.txt" > "${outdir}/snpEff.config"

# Convert coding sequences to GTF format
echo "Converting coding sequences to GTF format..."
python ${bindir}/coding2gtf22.py \
    --input "${outdir}/${outdb_prefix}_rep_seq.fasta" \
    --output "${outdir}/data/influenza_a_b/genes.gtf"

# Build the SNP database with snpEff
echo "Building SNP database with snpEff..."
${snpEff} build \
    -gtf22 -v \
    -c "${outdir}/snpEff.config" \
    -noCheckCds -noCheckProtein \
    -maxErrorRate 0.2 \
    influenza_a_b &> "${outdir}/log.gtf.txt"

echo "Pipeline completed successfully."
