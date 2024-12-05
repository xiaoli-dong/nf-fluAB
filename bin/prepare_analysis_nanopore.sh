#!/bin/bash
mkdir fastq
#concatnate the fastq for each samples
for var in {1..9}; do  
  echo $var; 
  cat */*/fastq_pass/barcode0$var/*.fastq.gz > fastq/barcode$var.fastq.gz; 
done

for var in {10..96}; do  
  echo $var; 
  cat */*/fastq_pass/barcode$var/*.fastq.gz > fastq/barcode$var.fastq.gz; 
done
#pipeline launcher has prepred fastq direcotry
# Create directories for analysis
mkdir -p analysis_apl/raw_data

# Initial 'run' value from current working directory
run=$(basename "$PWD")  # Set 'run' variable to the current directory name

echo "Processing files in directory: $run"

#cd ./fastq
# Loop through fastq.gz files in the ./fastq directory
for x in fastq/barcode*.fastq.gz; do
  fname=${x/*\//}  
  # Print debugging information
  echo "Processing sample: $x"
  echo "Run: $run"

  # Copy the file to the new directory
  cp $x ./analysis_apl/raw_data/${run}-${fname}
done

# Change to the directory where the files were copied
cd analysis_apl/raw_data

# Print the header once before processing any files
echo "sample,fastq_1,fastq_2,long_fastq" > ../samplesheet.csv

# Loop through the *R1* files and process each one with Perl
ls -l *.fastq.gz | while read -r line; do
  # Extract the filename from the ls -l output
  file=$(echo "$line" | awk '{print $NF}')
  
  # Run the Perl one-liner with updated 'run' variable
  perl -e "
    if (\$ARGV[0] =~ /.*?${run}-barcode0?(\d+).fastq.gz/) {
      print \"${run}-S\$1,NA,NA,./raw_data/${run}-barcode\$1.fastq.gz\n\";
    }
  " "$file" >> ../samplesheet.csv
done
cp /nfs/APL_Genomics/apps/production/influenza/slurm_nanopore.batch . 
cp /nfs/APL_Genomics/apps/production/influenza/fluab_routine.config .

cd -

