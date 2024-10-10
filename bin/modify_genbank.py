from Bio import SeqIO

# Define the input and output file paths
input_file = "input.gb"    # Path to the input GenBank file
output_file = "output.gb"  # Path to the output GenBank file

def remove_mat_peptide_feature(record):
    # Filter out 'mat_peptide' features
    record.features = [feature for feature in record.features if feature.type != 'mat_peptide']
    return record

def main():
    # Read the GenBank file
    records = list(SeqIO.parse(input_file, "genbank"))

    # Process each record to remove 'mat_peptide' features
    for record in records:
        record = remove_mat_peptide_feature(record)

    # Write the modified records to a new GenBank file
    with open(output_file, "w") as output_handle:
        SeqIO.write(records, output_handle, "genbank")

if __name__ == "__main__":
    main()