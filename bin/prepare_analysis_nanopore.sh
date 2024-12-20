#!/bin/bash

# Function to display usage instructions
usage() {
  echo "Usage: $0 -r <run_name> -o <output_directory> -b <barcode_range>"
  echo
  echo "This script processes raw sequencing data files for the specified 'run'."
  echo "-r <run_name>           The name of the run to be processed."
  echo "-o <output_directory>   The directory to store the processed files (default: ./analysis_apl)"
  echo "-b <barcode_range>      Range of barcodes to process (default: 1-96)"
  echo
  echo "Example:"
  echo "  $0 -r my_run_name -o /path/to/output -b 1-96"
  exit 1
}

# Default values for the options
output_dir="./analysis_apl"
barcode_range="1-96"

# Parse command-line arguments using getopts
while getopts ":r:o:b:" opt; do
  case $opt in
    r)
      run=$OPTARG
      ;;
    o)
      output_dir=$OPTARG
      ;;
    b)
      barcode_range=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      ;;
  esac
done

# Ensure the 'run' parameter is provided
if [ -z "$run" ]; then
  echo "Error: 'run' parameter is required."
  usage
fi

# Ensure the output directory exists
mkdir -p "$output_dir/fastq"
mkdir -p "$output_dir/raw_data"

echo "Processing files for run: $run"
echo "Output directory: $output_dir"
echo "Barcode range: $barcode_range"

# Concatenate the fastq files for each sample
# Use a loop based on the barcode range (1 to 96 by default)
IFS="-" read -r start_barcode end_barcode <<< "$barcode_range"
for var in $(seq "$start_barcode" "$end_barcode"); do
  # If the barcode has a leading zero (e.g., barcode01, barcode02, ...), we need to handle this for both single and double-digit barcodes
  echo "Processing barcode $var"
  
  # Check if the barcode has a leading zero (1-9) or not (10-96) and process accordingly
  if [ $var -lt 10 ]; then
    cat */*/fastq_pass/barcode0$var/*.fastq.gz > "$output_dir/fastq/barcode0$var.fastq.gz"
  else
    cat */*/fastq_pass/barcode$var/*.fastq.gz > "$output_dir/fastq/barcode$var.fastq.gz"
  fi
done

# Loop through fastq.gz files in the fastq directory
for x in "$output_dir"/fastq/barcode*.fastq.gz; do
  fname=${x/*\//}  # Extract file name from the full path
  
  # Print debugging information
  echo "Processing sample: $x"
  echo "Run: $run"

  # Copy the file to the new directory
  cp "$x" "$output_dir/raw_data/${run}-${fname}"
done

# Change to the directory where the files were copied
cd "$output_dir/raw_data"

# Print the header once before processing any files
echo "sample,fastq_1,fastq_2,long_fastq" > "../samplesheet.csv"

# Loop through the fastq files and process each one with Perl
ls -l *.fastq.gz | while read -r line; do
  # Extract the filename from the ls -l output
  file=$(echo "$line" | awk '{print $NF}')
  
  # Run the Perl one-liner with the updated 'run' variable
  perl -e "
    if (\$ARGV[0] =~ /.*?${run}-barcode0?(\d+).fastq.gz/) {
      print \"${run}-S\$1,NA,NA,./raw_data/${run}-barcode\$1.fastq.gz\n\";
    }
  " "$file" >> "../samplesheet.csv"
done

cd ..
# Copy the batch and config files
cp /nfs/APL_Genomics/apps/production/influenza/slurm_nanopore.batch . 
cp /nfs/APL_Genomics/apps/production/influenza/fluab_routine.config .

cd -
