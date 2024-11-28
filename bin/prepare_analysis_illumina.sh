#!/bin/bash

# Create directories for analysis
mkdir -p analysis_apl/raw_data

# Initial 'run' value from current working directory
run=$(basename "$PWD")  # Set 'run' variable to the current directory name

echo "Processing files in directory: $run"

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
  
  # Copy the file to the new directory
  cp "$x" "./analysis_apl/raw_data/$run-$fname"
done


# Change to the directory where the files were copied
cd analysis_apl/raw_data

# Print the header once before processing any files
echo "sample,fastq_1,fastq_2,long_fastq" > ../samplesheet.csv

# Loop through the *R1* files and process each one with Perl
ls -l *R1_* | while read -r line; do
  # Extract the filename from the ls -l output
  file=$(echo "$line" | awk '{print $NF}')
  
  # Run the Perl one-liner with updated 'run' variable
  perl -e "
    if (\$ARGV[0] =~ /.*?${run}-(\d+)(\S+?)R1(\S+)/) {
      print \"${run}-S\$1,./raw_data/${run}-\$1\$2R1\$3,./raw_data/${run}-\$1\$2R2\$3,NA\n\";
    }
  " "$file" >> ../samplesheet.csv
done
cp /nfs/APL_Genomics/apps/production/influenza/slurm_illumina.batch . 

cd -