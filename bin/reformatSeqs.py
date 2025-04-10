#!/usr/bin/env python
import sys
import argparse
from collections import defaultdict
import re
from csv import DictReader
import os

__author__ = 'Xiaoli Dong'
__version__ = '0.1.0'
__maintainer__ = 'Xiaoli Dong'
__email__ = 'xiaoli.dong@albertaprecisionlabs.ca'
__status__ = 'Dev'

def count_ambiguous_bases(sequence):
    ambiguous_bases = "RYSMKWBVDHNX"  # Set of ambiguous base codes
    count = 0
    
    # Loop through each base in the sequence and check if it's ambiguous
    for base in sequence.upper():  # Convert to uppercase to handle both cases
        if base in ambiguous_bases:
            count += 1
    return count
def remove_whitespace_regex(input_string):
    return re.sub(r'\s+', '', input_string)  # '\s+' matches all whitespace characters

"""
Filter fasta:
    seq len < seq_minlen 
    seq len > seq_maxlen
    The ambiguous number is > maxambigs 
    add seq length to the end of header
"""    
def validate_sequences(file_path, output_file, seq_minlen, seq_maxlen, maxambigs):
    seqids_above_cutoff = [] 
    with open(file_path, 'r') as input_file, open(output_file, 'w') as output:
        sequence = ""
        header = None
        for line in input_file:
            if line.startswith('>'):
                if header:
                    sequence = remove_whitespace_regex(sequence)
                    seqlen = len(sequence)
                    ambiguous_count = count_ambiguous_bases(sequence)

                    if (seq_minlen <= seqlen <= seq_maxlen) and (ambiguous_count <= maxambigs):
                        
                        # Write the sequence with the previous header, including its length
                        output.write(f"{header} len={seqlen}\n")
                        output.write(f"{sequence}\n")
                        # Extract the sequence ID without the version number
                        seqid = header.split(' ')[0][1:]

                        seqids_above_cutoff.append(seqid)
                # Write the current header (without any sequence yet)
                header = line.strip()  # Remove newline characters
                # Use regex to remove the version number (e.g., ".1", ".2", etc.)
                header = re.sub(r'\.\d+', '', header)
                sequence = ""  # Reset sequence
            else:
                sequence += line.strip()  # Add sequence lines, stripping newlines

        # After the loop, write the last sequence with its header and length
        if header:
            sequence = remove_whitespace_regex(sequence)
            seqlen = len(sequence)
            ambiguous_count = count_ambiguous_bases(sequence)
            if seq_minlen <= seqlen <= seq_maxlen and ambiguous_count < maxambigs:
                # Write the sequence with the previous header, including its length
                output.write(f"{header} len={seqlen}\n")
                output.write(f"{sequence}\n")
                # Extract the sequence ID without the version number
                seqid = header.split(' ')[0][1:]
                seqids_above_cutoff.append(seqid)

    return  seqids_above_cutoff


"""
Extract and reformat fasta:
    seq must passed all the cutoff 
    seq has metadata available
    seq has subtype/lineage available
"""    
def extract_refromat_fasta(file_path, metadata_dict):
    
    #filter sequence by length 
    #https://pmc.ncbi.nlm.nih.gov/articles/PMC3074182/
    segment_length = {
        "1": 2341, 
        "2": 2341, 
        "3": 2233, 
        "4": 1778, 
        "5": 1565,  
        "6": 1413, 
        "7": 1027, 
        "8": 890
    }
    with open(file_path, 'r') as file:
        header = None
        sequence_lines = [] 
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                #>PQ264394 |Influenza A virus (A/Massachusetts/47/2024(H3N2)) segment 6 neuraminidase (NA) gene, complete cds len=1442
                if header:
                    seqid = header.split()[0]
                    seqlen = int(header.split()[-1].split("=")[1])
                    sequence = remove_whitespace_regex(''.join(sequence_lines))
                    #print(f"{header} before checking metadata", file=sys.stderr)
                    if seqid in metadata_dict:
            
                        metadata = metadata_dict[seqid]
                        segid = metadata.split("|")[1]
                        segid_avg_len = segment_length[segid]
                        if segid_avg_len * 0.90 <= seqlen <= segid_avg_len * 1.10:
                            print(f">{seqid} {metadata} len={seqlen}")
                            print(sequence)
                    else:
                        # Handle case where sequence ID is not in metadata
                        print(f"{seqid} has no metadata available", file=sys.stderr)
                       
                header = line[1:]  # Remove the '>' character
                sequence_lines = []
            else:
                sequence_lines.append(line)
        
        # Save the last sequence
        if header:
            seqid = header.split()[0]
            sequence = remove_whitespace_regex(''.join(sequence_lines))
            if seqid in metadata_dict:
                metadata = metadata_dict[seqid]
                segid = metadata.split("|")[1]
                segid_avg_len = segment_length[segid]
                if segid_avg_len * 0.90 <= seqlen <= segid_avg_len * 1.10:
                    print(f">{seqid} {metadata} len={seqlen}")
                    print(sequence)
            else:
                # Handle case where sequence ID is not in metadata
                print(f"{seqid} has no metadata available", file=sys.stderr)
        
# Function to check if both 'H' and 'N' are in the string, ignoring case
def contains_h_and_n(str):
    str_lower = str.lower()
    return 'h' in str_lower and 'n' in str_lower

# Function to check if the string contains 'Yamagata' or 'Victoria', ignoring case
def contains_yamagata_or_victoria(str):
    str_lower = str.lower()
    return 'yamagata' in str_lower or 'victoria' in str_lower 
def starts_with_h(str):
    return str.lower().startswith('h')
"""
extract metadata for all the sequences ids included in ids_above_cutoff
exclude the ids:
    no lineage/subtype availabe
    genome is Deprecated
    quality is poor
    bad segment ids 
"""
def extract_metadata_from_bvbrc_csv(bvbrc, output_file):

    metadata = defaultdict(str)
    segid2gname = {
        "1": "PB2", 
        "2": "PB1", 
        "3": "PA", 
        "4": "HA", 
        "5": "NP",  
        "6": "NA", 
        "7": "M", 
        "8": "NS"
    }  
    gname2segid = {
        "PB2": "1",
        "PB1": "2",
        "PA": "3", 
        "HA": "4", 
        "NP": "5",  
        "NA": "6", 
        "M": "7", 
        "NS": "8"
    }

    
    keys_to_include = [
            'GenBank Accessions', 
            'Host Common Name', 
            'Segment', 
            'Protein',
            'Subtype', 
            'Isolation Country', 
            'Strain',  
            'H1 Clade Global',
            'H1 Clade US',
            'H5 Clade'
    ]

    #print(*keys_to_include, sep="\t")
    with open(bvbrc, 'r') as input_file, open(output_file, 'w') as output:
    #with open(bvbrc, 'r') as f: 

        dict_reader = DictReader(input_file) 
        count = 0
        for row in dict_reader:
            extracted_columns = []
            lineage = None
            is_flua = False
            is_flub = False
            valid_sub_type = None
            

            if row['Species'] == "Betainfluenzavirus influenzae":
                is_flub = True
            if row['Species'] == "Alphainfluenzavirus influenzae":
                is_flua = True

            if (not is_flua) and (not is_flub):
                continue
           
            if row["Genome Status"] in ["Deprecated"]:
                continue
            if row["Genome Quality"] in ["Poor"]:
                continue
            if row['Segment'] not in list(segid2gname.keys()) + list(gname2segid.keys()):
                continue

            # extract lineage:  ";lineage:Victorial;"
            lineage_pattern = r";lineage:(\S+?);"
            lineage_match = re.search(lineage_pattern, str(row))
            # flu B has no subtype, use lineage
            if is_flub and lineage_match:
                # Lineage is the part after "lineage:" and before the semicolon
                lineage = lineage_match.group(1)
                row['Subtype'] = lineage

            subtype = row['Subtype'] 
            if (not (contains_h_and_n(subtype) and starts_with_h(subtype))) and (not contains_yamagata_or_victoria(subtype)):
                continue
            # #exclude items whose subtype/lineage is not available or not valid
            # valid_sub_type = re.match(r"^H|^N|^Yamagata|^Victoria", subtype)
            # if valid_sub_type is None:
            #     continue
            #print(row, file=sys.stderr)
            subset_row_dict = dict(filter(lambda item: item[0] in keys_to_include, row.items()))

            # in bvrc file, some of the entries is using gene name as segment id in 'Segment' column, fix it
            if row['Segment'] in list(gname2segid.keys()):
                subset_row_dict['Segment'] = gname2segid[row['Segment']]
            
            subset_row_dict['Protein'] = segid2gname[subset_row_dict['Segment']]
            
            for x in keys_to_include:
                extracted_columns.append(subset_row_dict[x])

            extracted_columns = [x.replace(' ', "_") for x in extracted_columns]
            extracted_columns = ['na' if x == '' else x for x in extracted_columns]
            seqid = extracted_columns[0]
            desc = '|'.join(extracted_columns[1:])
            if desc:
                metadata[seqid] =  desc
                #print(f"gggggggggggggggggg {seqid} = {desc}", file=sys.stderr)
                # Write the sequence with the previous header, including its length
                output.write(f"{seqid} {desc}\n")

            #print(*extracted_columns, sep="\t")
            count += 1
            if count %100 ==0:
                print(f"************************* count={count}", file=sys.stderr)

    return metadata
def main():
    script_name = os.path.basename(__file__)
    # Initialize parser
    parser = argparse.ArgumentParser(
        description=f"Example of {script_name} usage.",
        epilog=(f"Example commands:\npython {script_name}\n\t--fasta seq_nt.fasta\n\t--bvbrc BVBRC_genome.csv\n\t> reformat.fasta 2> log.txt\n")
        ,formatter_class=argparse.RawTextHelpFormatter
    )

    # Adding optional argument
    parser.add_argument("-f", "--fasta", help = "fasta", required=True)
    parser.add_argument("--minlen", type=int, default=800, help = "min sequence length", required=False)
    parser.add_argument("--maxlen", type=int, default=2500, help = "max sequence length", required=False)
    parser.add_argument("--maxambigs", type=int, default=0, help = "max number of ambigs bases", required=False)
    parser.add_argument("--bvbrc", help = "bvbrc csv file", required=True)
    
    #python ../bin/reformatSeqs.py  --fasta sequences_nt.20250331.fasta --bvbrc BVBRC_genome.20250331.csv > sequences_nt.20250331.reformat.fasta  2> log.txt
    # Read arguments from command line
    args = parser.parse_args()
    fasta_base_name = os.path.splitext(args.fasta)[0]  # Removes the suffix/extension
    filtered_fasta = fasta_base_name + "_filter.fasta"
    ids_above_cutoff = validate_sequences(args.fasta, filtered_fasta, args.minlen, args.maxlen, args.maxambigs)
    #read bvbrc-genome metadata
    print(f"*******************Start parsing bvbrc csv file", file=sys.stderr)
    bvbrc_base_name = os.path.splitext(args.bvbrc)[0]  # Removes the suffix/extension
    reformat_and_filtered_bvbrc_file = bvbrc_base_name + "_filter_reformat.csv"
    metadata_dict =extract_metadata_from_bvbrc_csv(args.bvbrc, reformat_and_filtered_bvbrc_file)
    print(f"*******************Finish parsing bvbrc csv file", file=sys.stderr)
    extract_refromat_fasta(filtered_fasta, metadata_dict)

if __name__ == "__main__":
    main()
