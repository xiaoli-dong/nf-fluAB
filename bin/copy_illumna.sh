#!/bin/bash

# Function to display usage instructions
usage() {
  cat <<EOF
Usage: $0 -i <input_directory> -r <run_name> -s <sample_spec> [-o <output_directory>]

This script processes specified sequencing samples from an Illumina run.

Required arguments:
  -i <input_dir>        Illumina sequence run directory (e.g., path/to/Alignment_1)
  -r <run_name>         The name of the run to be processed
  -s <sample_spec>      Sample specification (e.g., "51,56,61" or "5,7,12-18,76")

Optional arguments:
  -o <output_directory> Directory to store processed files (default: ./nf-fluAB_analysis)

Examples:
  $0 -i /path/to/Alignment_1 -r RUN2023 -s "51,56,61"
  $0 -i /path/to/Alignment_1 -r RUN2023 -s "5-8,12-15" -o /path/to/output
EOF
  exit 1
}

# Function to expand sample ranges
expand_ranges() {
  local ranges=$1
  local result=()

  IFS=',' read -ra parts <<< "$ranges"
  for part in "${parts[@]}"; do
    if [[ $part =~ ^([0-9]+)-([0-9]+)$ ]]; then
      for ((i=${BASH_REMATCH[1]}; i<=${BASH_REMATCH[2]}; i++)); do
        result+=("$i")
      done
    else
      result+=("$part")
    fi
  done

  echo "${result[@]}"
}

# Initialize variables
output_dir="./nf-fluAB_analysis"
input_dir=""
run=""
sample_spec=""
script_dir=$(dirname "$(readlink -f "$0")")

# Parse command-line arguments
while getopts ":i:r:s:o:h" opt; do
  case $opt in
    i) input_dir=$(readlink -f "$OPTARG") ;;
    r) run=$OPTARG ;;
    s) sample_spec=$OPTARG ;;
    o) output_dir=$(readlink -f "$OPTARG") ;;
    h) usage ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Validate required arguments
if [ -z "$run" ]; then
  echo "Error: 'run' parameter is required." >&2
  usage
fi

if [ -z "$input_dir" ]; then
  echo "Error: 'input_dir' parameter is required." >&2
  usage
fi

if [ -z "$sample_spec" ]; then
  echo "Error: 'sample_spec' parameter is required." >&2
  usage
fi

if [ ! -d "$input_dir" ]; then
  echo "Error: Input directory '$input_dir' does not exist." >&2
  exit 1
fi

# Expand sample specification
echo "Expanding sample specification: $sample_spec"
sample_array=($(expand_ranges "$sample_spec"))
if [ ${#sample_array[@]} -eq 0 ]; then
  echo "Error: No valid samples found in specification '$sample_spec'" >&2
  exit 1
fi

# Create output directories
mkdir -p "$output_dir/raw_data" || { echo "Error creating output directory"; exit 1; }
mkdir -p "$output_dir/fastq" || { echo "Error creating fastq directory"; exit 1; }

# Initialize samplesheet
samplesheet="$output_dir/samplesheet.csv"
echo "sample,fastq_1,fastq_2,long_fastq" > "$samplesheet"

# Main processing loop
for sample_num in "${sample_array[@]}"; do
  echo "Processing sample: $sample_num"

  # Initialize arrays
  r1_files=()
  r2_files=()

  # Find files using temporary file to avoid process substitution
  tmp_file=$(mktemp)
  find "$input_dir" -type f -name "*${sample_num}_S${sample_num}_*.fastq.gz" > "$tmp_file"

  # Read found files into array
  while IFS= read -r file; do
    filename=$(basename "$file")
    if [[ "$filename" =~ _R1_ ]]; then
      r1_files+=("$file")
    elif [[ "$filename" =~ _R2_ ]]; then
      r2_files+=("$file")
    fi
  done < "$tmp_file"
  rm "$tmp_file"

  # Verify file pairs
  if [ ${#r1_files[@]} -eq 0 ]; then
    echo "Warning: No R1 files found for sample $sample_num" >&2
    continue
  fi
  if [ ${#r2_files[@]} -eq 0 ]; then
    echo "Warning: No R2 files found for sample $sample_num" >&2
    continue
  fi

  # Sort files by lane number
  IFS=$'\n' sorted_r1=($(sort -V <<<"${r1_files[*]}"))
  IFS=$'\n' sorted_r2=($(sort -V <<<"${r2_files[*]}"))
  unset IFS

  # Process first pair
  r1_file="${sorted_r1[0]}"
  r2_file="${sorted_r2[0]}"
  r1_basename=$(basename "$r1_file")
  r2_basename=$(basename "$r2_file")

  # Copy files
  cp "$r1_file" "$output_dir/fastq/"
  cp "$r2_file" "$output_dir/fastq/"
  cp "$r1_file" "$output_dir/raw_data/${run}-${r1_basename}"
  cp "$r2_file" "$output_dir/raw_data/${run}-${r2_basename}"

  # Add to samplesheet
  echo "${run}-${sample_num},./raw_data/${run}-${r1_basename},./raw_data/${run}-${r2_basename},NA" >> "$samplesheet"
done

# Copy configuration files
echo "Copying configuration files..."
cp "$script_dir/slurm_illumina.batch" "$output_dir/"
cp "$script_dir/fluab_routine.config" "$output_dir/"

echo "Processing complete. Results in: $output_dir"
echo "Samplesheet created at: $samplesheet"
