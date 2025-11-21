#!/bin/bash

# Function to display usage instructions
usage() {
  cat <<EOF
Usage: $0 -i <nanopore sequence directory> -r <run_name> [-o <output_directory>] [-b <barcode_range>]

This script processes raw sequencing data files for the specified 'run'.

Required arguments:
  -i <input_dir>          Nanopore sequence run directory or fastq_pass dir
  -r <run_name>           The name of the run to be processed

Optional arguments:
  -o <output_directory>   Directory to store processed files (default: ./nf-fluAB_analysis)
  -b <barcode_range>      Comma-separated list of barcodes or ranges (e.g., 01-05,13,23; default: 01-96)

Examples:
  Process all barcodes (default):
    $0 -i /path/to/Run123/sequence_run_dir -r Run123

  Process barcodes 01 to 12 only:
    $0 -i /nfs/APL_Genomics/raw/ncov-raw/241227_S_N_186/no_sample_id -r 241227_S_N_186 -b 01-12


  Process specific barcodes 01, 05, and 08:
    $0 -i  /path/to/Run123/sequence_run_dir -r Run123 -o /path/to/output_dir -b 01,05,08

  Process mixed ranges:
    $0  -i /nfs/APL_Genomics/raw/ncov-raw/241010_S_N_155/no_sample -r 241010_S_N_155 -b 01-05,23,13
EOF
  exit 1
}

# Default values
input_dir="."
output_dir="./nf-fluAB_analysis"
barcode_range="01-96"
script_dir=$(dirname "$(realpath "$0")")
echo ${script_dir}
# Parse command-line arguments
while getopts ":i:r:o:b:h" opt; do
  case $opt in
    i) input_dir=$OPTARG ;;
    r) run=$OPTARG ;;
    o) output_dir=$OPTARG ;;
    b) barcode_range=$OPTARG ;;
    h) usage ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Check required params
if [ -z "$run" ] || [ -z "$input_dir" ]; then
  echo "Error: Both -i <input_dir> and -r <run_name> are required." >&2
  usage
fi

# Create output directories
mkdir -p fastq
mkdir -p "$output_dir/raw_data"

echo "Processing run: $run"
echo "Input directory: $input_dir"
echo "Output directory: $output_dir"
echo "Barcode range(s): $barcode_range"

# Function to expand barcode ranges (e.g. "01-03,05,08-10")
expand_barcodes() {
  local IFS=',' ranges=($1)
  local result=()
  for range in "${ranges[@]}"; do
    if [[ "$range" =~ ^([0-9]{2})-([0-9]{2})$ ]]; then
      start=${BASH_REMATCH[1]}
      end=${BASH_REMATCH[2]}
      for ((i=10#$start; i<=10#$end; i++)); do
        result+=("$(printf "%02d" $i)")
      done
    elif [[ "$range" =~ ^[0-9]{2}$ ]]; then
      result+=("$range")
    else
      echo "Invalid barcode range format: $range" >&2
      exit 1
    fi
  done
  echo "${result[@]}"
}

# Expand barcode list
barcodes=($(expand_barcodes "$barcode_range"))

# Concatenate fastq files per barcode
for dir in $(find "$input_dir" -type d -name fastq_pass); do
  echo "Found fastq_pass directory: $dir"
  for bc in "${barcodes[@]}"; do
    bc_dir="$dir/barcode$bc"
    if ls "$bc_dir"/*.fastq.gz 1>/dev/null 2>&1; then
      echo "Concatenating barcode $bc..."
      cat "$bc_dir"/*.fastq.gz > "fastq/barcode${bc}.fastq.gz"
    else
      echo "Warning: No fastq.gz files found for barcode $bc in $bc_dir"
    fi
  done
done

# Copy concatenated files to output raw_data folder with run prefix
for x in fastq/barcode*.fastq.gz; do
  fname=${x##*/}  # strip path
  echo "Copying $x to raw_data/${run}-${fname}"
  cp "$x" "$output_dir/raw_data/${run}-${fname}"
done

# Generate samplesheet.csv
cd "$output_dir/raw_data" || { echo "Failed to cd to $output_dir/raw_data"; exit 1; }
echo "sample,fastq_1,fastq_2,long_fastq" > ../samplesheet.csv

for file in ${run}-barcode*.fastq.gz; do
  if [[ $file =~ ${run}-barcode([0-9]{2}).fastq.gz ]]; then
    bc="${BASH_REMATCH[1]}"
    echo "${run}-S$bc,NA,NA,./raw_data/${run}-barcode$bc.fastq.gz" >> ../samplesheet.csv
  fi
done

cd ..

# Copy batch and config files from script directory using realpath
#script_dir=$(realpath "$(dirname "$0")")


cp "${script_dir}/slurm_nanopore.batch" . 2>/dev/null || echo "Warning: slurm_nanopore.batch not found."
cp "${script_dir}/fluab_routine.config" . 2>/dev/null || echo "Warning: fluab_routine.config not found."

echo "✅ Processing complete."
echo "Output directory: $output_dir"
echo "Samplesheet: $output_dir/samplesheet.csv"
