#!/bin/bash

# Function to display usage instructions
usage() {
  echo "Usage: $0 -i <nanopore sequence directory> -r <run_name> -o <output_directory> -b <barcode_range>"
  echo
  echo "This script processes raw sequencing data files for the specified 'run'."
  echo "-i <input_dir>          Nanopore sequence run directory or fastq_pass dir"
  echo "-r <run_name>           The name of the run to be processed."
  echo "-o <output_directory>   The directory to store the processed files (default: ./analysis_apl)"
  echo "-b <barcode_range>      Range of barcodes to process (default: 01-96)"
  echo
  echo "Example:"
  echo "  $0 -i nanopore_seq_dir -r my_run_name -o /path/to/output -b 01-96"
  exit 1
}

# Default values for the options
input_dir="."
output_dir="./analysis_apl"
barcode_range="01-96"

# Parse command-line arguments using getopts
while getopts ":i:r:o:b:" opt; do
  case $opt in
    i)
      input_dir=$OPTARG
      ;;
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

# Ensure the output directory 
if [ ! -d fastq ]; then
    echo "Directory does not exist. Creating it now."
    mkdir -p fastq  # The -p flag ensures that parent directories are created if they don't exist.
fi

mkdir -p "$output_dir/raw_data"

echo "Processing files for run: $run"
echo "Input nanopore sequence directory: $input_dir"
echo "Output directory: $output_dir"
echo "Barcode range: $barcode_range"

# Concatenate the fastq files for each sample
# Use a loop based on the barcode range (01 to 96 by default)
IFS="-" read -r start_barcode end_barcode <<< "$barcode_range"
echo $start_barcode

for dir in $(find $input_dir -type d -name fastq_pass); do 
  echo $dir;
  for (( i=$start_barcode; i<=$end_barcode; i++ )); do
    var=$(printf "%02d\n" $i)
    echo "Processing barcode $dir/barcode$var"
    cat $dir/barcode${var}/*.fastq.gz > "fastq/barcode${var}.fastq.gz";
  done 
done


# Loop through fastq.gz files in the fastq directory
for x in fastq/*barcode*.fastq.gz; do
  echo $x
  fname=${x/*\//}  # Extract file name from the full path
  fname=${fname//${run}-/}  # get rid of run name from the fastq file
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
    if (\$ARGV[0] =~ /.*?${run}-barcode(\d+).fastq.gz/) {
      print \"${run}-S\$1,NA,NA,./raw_data/${run}-barcode\$1.fastq.gz\n\";
    }
  " "$file" >> "../samplesheet.csv"
done

cd ..
# Copy the batch and config files
script_dir=$(dirname "$0")
cp ${script_dir}/slurm_nanopore.batch . 
cp ${script_dir}/fluab_routine.config .

cd -
