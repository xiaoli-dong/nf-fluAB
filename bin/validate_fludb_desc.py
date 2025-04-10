import sys
import csv
from collections import defaultdict
import traceback
# Parse FASTA file and return sequences and descriptions in two dictionaries
def parse_protein_fasta(fasta_file):
    sequences = defaultdict(str)
    descriptions = {}

    # Correct description mapping (example)
    desc_correct = {
        "YP_308664.1": "polymerase PB2",
        "YP_308666.1": "polymerase PA",
        "YP_308665.1": "polymerase PB1"
    }

    try:
        with open(fasta_file, 'r') as file:
            seq_id = None
            description = None
            seq = ''

            for line in file:
                line = line.strip()

                if line.startswith('>'):  # Header line starts with '>'
                    if seq_id:  # Save the previous sequence before moving to the next
                        sequences[seq_id] = seq
                        descriptions[seq_id] = desc_correct.get(seq_id, description)

                    seq_id = line[1:].split()[0]  # Get the sequence ID (first part of the header)

                    # Extract description part after sequence ID and up to the first square bracket '['
                    description = line.split(' ', 1)[1].split('[')[0].strip()

                    seq = ''  # Reset sequence string for the next sequence
                else:
                    seq += line  # Append sequence data

            # Save the last sequence in the file
            if seq_id:
                sequences[seq_id] = seq
                descriptions[seq_id] = desc_correct.get(seq_id, description)

    except Exception as e:
        print(f"Error reading FASTA file: {e}", file=sys.stderr)
        sys.exit(1)

    return sequences, descriptions

def parse_fasta_header(header):
    # Remove the leading '>' character from the header
    header = header.lstrip('>')
    
    # Split at the first space to separate the sequence ID from the description
    parts = header.split(' ', 1)  # Split into at most 2 parts
    
    # The first part is the sequence ID, the rest is the description
    seq_id = parts[0]
    description = parts[1] if len(parts) > 1 else ""  # Handle case where there's no description
    
    return seq_id, description

# Parse FASTA file and return sequences and descriptions in two dictionaries
def parse_fasta(fasta_file):
    sequences = defaultdict(str)
    descriptions = {}
    
    try:
        with open(fasta_file, 'r') as file:
            seq_id = None
            description = None
            seq = ''

            for line in file:
                line = line.strip()
               
                if line.startswith('>'):  # Header line starts with '>'
                    if seq_id:  # Save the previous sequence before moving to the next
                        sequences[seq_id] = seq
                        descriptions[seq_id] = description
                    
                    seq_id, description = parse_fasta_header(line)
                    seq = ''  # Reset sequence string for the next sequence
                else:
                    seq += line  # Append sequence data
                
            # Save the last sequence in the file
            if seq_id:
                sequences[seq_id] = seq
                descriptions[seq_id] = description

    except Exception as e:
        print(f"Error reading FASTA file: {e}", file=sys.stderr)
        # Capture the detailed exception
        error_details = traceback.format_exc()
        
        print(f"Detailed traceback: {error_details}", file=sys.stderr)
        sys.exit(1)
    
    return sequences, descriptions

# Parse TSV file and return a dictionary of key-value pairs
def parse_tsv(tsv_file):
    tsv_dict = {}

    try:
        with open(tsv_file, 'r') as file:
            for line in file:
                columns = line.strip().split('\t')
                if len(columns) >= 2:
                    key = columns[0].strip()
                    value = columns[1].strip()
                    tsv_dict[key] = value
    except Exception as e:
        print(f"Error reading TSV file: {e}", file=sys.stderr)
        sys.exit(1)

    return tsv_dict

# Parse CSV file and return two dictionaries for descriptions and gene names
def parse_csv(csv_file):
    desc2segid_dict = {}
    desc2genename_dict = {}

    try:
        with open(csv_file, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                if len(row) >= 3:
                    key = row[0].strip().lower()  # First column as key
                    value1 = row[1].strip()  # Second column as segment ID
                    value2 = row[2].strip()  # Third column as gene name
                    desc2segid_dict[key] = value1
                    desc2genename_dict[key] = value2
    except Exception as e:
        print(f"Error reading CSV file: {e}", file=sys.stderr)
        sys.exit(1)

    return desc2segid_dict, desc2genename_dict

# Main function to process the files and produce output
def main():
    if len(sys.argv) != 5:
        print("Usage: python script.py <protein_fluab_protein_fasta_file> <diamond_blastx_format_6_tsv_file> <flu_protein_desc_to_segment_csv_file> <fluab_database_fasta_file>")
        sys.exit(1)

    fluab_protein_fasta_file = sys.argv[1]
    diamon_blastx_fmt6_tsv_file = sys.argv[2]
    flu_protein_desc2segment_file = sys.argv[3]
    fludb_fasta = sys.argv[4]

    # Parse the input files
    # >YP_009118470.1 polymerase PB2 [Influenza A virus (A/Shanghai/02/2013(H7N9))]
    # >YP_009118471.1 polymerase PB1 [Influenza A virus (A/Shanghai/02/2013(H7N9))]
    aa_sequences, aa_descriptions = parse_protein_fasta(fluab_protein_fasta_file)
    seq2proteinhits_dict = parse_tsv(diamon_blastx_fmt6_tsv_file)
    desc2segid_dict, desc2genename_dict = parse_csv(flu_protein_desc2segment_file)

    # Print the extracted sequences and descriptions from the FASTA file
    print("Protein FASTA file contents:", file=sys.stderr)
    """  for aa_seqid, aa_description in aa_descriptions.items():
        print(f"ID: {aa_seqid}, Description: {aa_description}", file=sys.stderr)
        #print(f"Sequence: {aa_sequences[seq_id]}\n") """

    # Print the extracted key-value pairs from the TSV file
    print("Diamon blastx fmt6 TSV file contents:", file=sys.stderr)
    for dna_seqid, aa_seqid in seq2proteinhits_dict.items():
        #print(f"Key: {dna_seqid}, Value: {aa_seqid}", file=sys.stderr)
        # Retrieve the description for the sequence from the descriptions dictionary
        aa_seq_description = aa_descriptions.get(aa_seqid, "Description not found")
        segid = desc2segid_dict.get(aa_seq_description.lower(), "Segment ID not found")
        segname = desc2genename_dict.get(aa_seq_description.lower(), "Segment ID not found")
        #print(f"New Segment ID: {segid}|{segname}")

    # Print unique descriptions from the sequences dict
    """ unique_descriptions = set(aa_descriptions.values())  # Set to remove duplicates
    print("\nUnique Descriptions:")
    for desc in unique_descriptions:
        print(desc) """

    # Parse the additional FASTA file
    dna_sequences, dna_descriptions = parse_fasta(fludb_fasta)
    #print(dna_descriptions)
    print(f"Parsed descriptions from {fludb_fasta}:", file=sys.stderr)
    for dna_seqid, dna_description in dna_descriptions.items():
        aa_sseqid = seq2proteinhits_dict.get(dna_seqid, "Hit ID not found")
        aa_description = aa_descriptions.get(aa_sseqid, "Description not found")

        segid = desc2segid_dict.get(aa_description.lower(), "Segment ID not found")
        segname = desc2genename_dict.get(aa_description.lower(), "Gene Name not found")
        import re
        #print(f"New Segment ID: {segid}|{segname}")
        new_value = segid + '|' + segname
        if new_value not in dna_description: 
            modified_desc = re.sub(r'(\|)[^|]+(\|)[^|]+(\|)', f'|{new_value}|', dna_description)
            print(f"ID: {dna_seqid}, Description: {dna_description} change to {modified_desc}", file=sys.stderr)
            print(f">{dna_seqid} {modified_desc}\n{dna_sequences[dna_seqid]}")
        else:
            print(f">{dna_seqid} {dna_description}\n{dna_sequences[dna_seqid]}")

if __name__ == "__main__":
    main()
