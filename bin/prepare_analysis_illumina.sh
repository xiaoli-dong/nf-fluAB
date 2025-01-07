#!/bin/bash

# Function to display usage instructions
usage() {
  echo "Usage: $0 -i <illumina sequence Alignment_* directory> -r <run_name> -o <output_directory>"
  echo
  echo "This script processes raw sequencing data files for the specified 'run'."
  echo "-i <input_dir>        Illumina sequence run directory, when you have multiple copies of the basecalls, provide the basecall directory. eg dir_to/Alignment_2 dir"
  echo "-r <run_name>         The name of the run to be processed."
  echo "-o <output_directory> The directory to store the processed files (default: ./analysis_apl)"
  echo
  echo "Example:"
  echo "  $0 -i my_alignment_dir -r my_run_name -o /path/to/output"
  exit 1
}

# Default output directory
output_dir="./analysis_apl"
input_dir="."

# Parse command-line arguments using getopts
while getopts ":i:r:o:" opt; do
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
    echo "copy fastq files"
    for dir in $(find $input_dir -type d -name Fastq); do 
      find $dir -type f -name "*.fastq.gz" ! -name 'Undetermined*' -exec cp {} fastq \;
    done
fi
# Create directories for analysis (based on output_dir)
mkdir -p "$output_dir/raw_data"

echo "Processing files for run: $run"
echo "Output directory: $output_dir"

# Loop through fastq.gz files in the ./fastq directory
for x in $(find ./fastq -name "*.fastq.gz"); do
  # Extract the file name without the path
  fname=${x/*\//}

  # Print debugging information
  echo "Processing file: $x"
  echo "Run: $run"
  echo "File name: $fname"
  
  # Remove the 'run' prefix from the file name (if present)
  fname=${fname/$run\_/}
  
  # Copy the file to the new directory (output_dir/raw_data)
  cp "$x" "$output_dir/raw_data/$run-$fname"
done

# Change to the directory where the files were copied
cd "$output_dir/raw_data"

# Print the header once before processing any files
echo "sample,fastq_1,fastq_2,long_fastq" > "../samplesheet.csv"

# Loop through the *R1* files and process each one with Perl
ls -l *R1_* | while read -r line; do
  # Extract the filename from the ls -l output
  file=$(echo "$line" | awk '{print $NF}')
  
  # Run the Perl one-liner with updated 'run' variable
  perl -e "
    if (\$ARGV[0] =~ /.*?${run}-(\d+)(\S+?)R1(\S+)/) {
      print \"${run}-S\$1,./raw_data/${run}-\$1\$2R1\$3,./raw_data/${run}-\$1\$2R2\$3,NA\n\";
    }
  " "$file" >> "../samplesheet.csv"
done

cd ..

# Copy the batch and config files
script_dir=$(dirname "$0")
cp ${script_dir}/slurm_illumina.batch . 
cp ${script_dir}/fluab_routine.config .

cd -
