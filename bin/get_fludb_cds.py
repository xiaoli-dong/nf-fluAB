#!/usr/bin/env python
import argparse
import re
from collections import defaultdict
import sys


# Function to process sequences and calculate the smallest start and largest end for those with multiple ranges
def extract_cds_coordinates_from_fasta(cds_fasta_file):
    sequence_coordinates = {}

    with open(cds_fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            # Skip empty lines and non-header lines
            if not line or not line.startswith(">"):
                continue
            """
            >join(NC_026431.1:1..26,NC_026431.1:715..982) |matrix protein 2 [Influenza A virus (A/California/07/2009(H1N1))]
            >NC_026431.1:1..759 |matrix protein 1 [Influenza A virus (A/California/07/2009(H1N1))]
            >complement(DQ098265.1:127..777) |NSP [Influenza A virus (A/Moscow/346/2003(H3N2))]
            """
            pattern = r'>(join\(([^)]+)\))|(complement\(([^)]+)\))'
            # Check if 'join' exists in the line
            if '>join' in line or '>complement' in line:
                # If it's a 'join' sequence
                match = re.match(pattern, line)
                if match:
                    # Extract the coordinates part inside the join() using group 1
                    coordinates_str = match.group(1)
                    seq_id = coordinates_str.split(":")[0]  # Extract sequence ID (e.g., NC_007376.1)
                    coordinates = []
                    for coord in coordinates_str.split(','):
                        #Regular expression pattern to capture the ranges in the format: <start..end>, start..end, or start..>end
                        pattern = re.compile(r"[<]*(\d+)\.\.[>]*(\d+)")
                        # Search the pattern in the string
                        matches = pattern.findall(line)
                        # Output the matches
                        for match in matches:
                            start = int(match[0])
                            end = int(match[1])
                            coordinates.append((start, end))
                        
                    
                    # Add coordinates to dictionary
                    if seq_id in sequence_coordinates:
                        sequence_coordinates[seq_id].extend(coordinates)
                    else:
                        sequence_coordinates[seq_id] = coordinates

            else:
                # If it's not a 'join' sequence, directly extract coordinates
                # >NC_026431.1:1..759 |matrix protein 1 [Influenza A virus (A/California/07/2009(H1N1))]
                # >MN875140.1:<1..>1445 |nucleocapsid protein, partial [Influenza A virus]
                #print(line)
                pattern = re.compile(r"^>(\S+):[<]*(\d+)\.\.[>]*(\d+)")
                #seq_id = line.split()[0][1:].split(":")[0]  # Extract sequence ID
                # Regular expression to match sequences with valid ranges (number..number or number..>number)
                
                #print(f"{line}",file=sys.stderr)
                match = pattern.search(line)
                if match:
                    seq_id = match.group(1)
                    start = int(match.group(2))
                    end = int(match.group(3)) 
                #coord_part = line.split()[0][1:].split(":")[1]
                #start, end = map(int, coord_part.split('..'))
                
                    # Add coordinates to dictionary
                    if seq_id in sequence_coordinates:
                        sequence_coordinates[seq_id].append((start, end))
                    else:
                        sequence_coordinates[seq_id] = [(start, end)]
   
    # Now, for each sequence ID, find the smallest start and largest end if there are multiple ranges
    seq_range_dict_multi = {}
    seq_range_dict_single = {}
    for seq_id, coordinates in sequence_coordinates.items():
        print(f"{seq_id}",file=sys.stderr)
        print(f"{coordinates}",file=sys.stderr)
        if len(coordinates) > 1:  # If there are multiple ranges
            min_start = min(coordinates, key=lambda x: x[0])[0]
            max_end = max(coordinates, key=lambda x: x[1])[1]
            seq_range_dict_multi[seq_id] = (min_start, max_end)
        else:
            print(seq_id)
            print(coordinates)
            start, end = coordinates[0]
            seq_range_dict_single[seq_id] = (start, end)

    return seq_range_dict_multi, seq_range_dict_single      

def extract_id_desc(genome_file):
    seqid_to_desc = {}
    
    with open(genome_file, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            
            if line.startswith(">"):  # Header line
                # Use regex to capture sequence ID and description after the vertical bar
                match = re.match(r'>(\S+)\s+\|(.+)', line)
                if match:
                    seqid = match.group(1)  # Sequence ID
                    description = match.group(2).strip()  # Description after the vertical bar
                    seqid_to_desc[seqid] = description  # Store in dictionary
    
    return seqid_to_desc

# Function to process the CDS file and genome file based on seq_range_dict
def extract_cds_from_genome(seq_range_dict_multi, genome_file_path, genomeid_to_desc):
    with open(genome_file_path, 'r') as genome_file:
        current_id = None
        current_sequence = []

        for line in genome_file:
            line = line.strip()  # Remove leading/trailing whitespaces

            # If the line starts with '>', it's a new sequence header
            if line.startswith(">"):
                if (current_id is not None) and (current_id in seq_range_dict_multi):
                    # Store the last sequence under its ID
                    sequence = ''.join(current_sequence)
                    start, end = seq_range_dict_multi[current_id]
                    subseq = sequence[start-1:end]
                    desc = genomeid_to_desc[current_id]
                    print(">" + current_id + ":" + str(start) + ".." + str(end) + " |" + desc)
                    print(subseq)
                # Start a new sequence
                current_id = line.split()[0][1:]  # Remove '>'
                current_sequence = []
            else:
                # This is part of the sequence (not a header)
                current_sequence.append(line)
        # After finishing the file, store the last sequence
        if (current_id is not None) and (current_id in seq_range_dict_multi):
            sequence = ''.join(current_sequence)
            start, end = seq_range_dict_multi[current_id]
            subseq = sequence[start-1:end]
            desc = genomeid_to_desc[current_id]
            print(">" + current_id + ":" + str(start) + ".." + str(end) + " |" + desc)
            print(subseq)
                    
# Function to process the CDS file and genome file based on seq_range_dict
def extract_cds_without_alternative_splicing(seq_range_dict_single, cds_file_path, genomeid_to_desc):
    with open(cds_file_path, 'r') as cds_file:
        current_id = None
        current_sequence = []

        for line in cds_file_path:
            line = line.strip()  # Remove leading/trailing whitespaces
            
            """
            >join(NC_007367.1:26..51,NC_007367.1:740..1007) |matrix protein 2 [Influenza A virus (A/New York/392/2004(H3N2))]
            >join(NC_007370.1:27..56,NC_007370.1:529..864) |nonstructural protein 2 [Influenza A virus (A/New York/392/2004(H3N2))]
            """
            # those sequences have alternative splicng
            if line.startswith(">join("):
                continue
            # Skip lines containing ':<' or '..>'
            """
            >PV390950.1:<1..71 |matrix protein 1, partial [Influenza A virus]
            >PV391022.1:<1..>1090 |neuraminidase, partial [Influenza A virus]
            >PV391023.1:<1..689 |nonstructural protein 1, partial [Influenza A virus]
            """
            #those are partial sequences
            # if ':<' in line or '..>' in line:
            #     continue
            
            # If the line starts with '>', it's a new sequence header
            if line.startswith(">"):
                if (current_id is not None) and (seq_id in seq_range_dict_single):
                    # Store the last sequence under its ID
                    sequence = ''.join(current_sequence)
                    
                    #print(">" + current_id )
                    if 'N' not in sequence.upper():
                        print(">" + current_id + " |" + genomeid_to_desc[seq_id])
                        print(sequence)
                    else:
                        print(f"{current_id} is has Ns",file=sys.stderr)
                # Start a new sequence
                current_id = line.split()[0][1:]  # Remove '>'
                seq_id = current_id.split(":")[0]
                current_sequence = []
            else:
                # This is part of the sequence (not a header)
                current_sequence.append(line)

        # After finishing the file, store the last sequence
        if (current_id is not None) and (seq_id in seq_range_dict_single):
            # Store the last sequence under its ID
            sequence = ''.join(current_sequence)
            #print(">" + current_id)
            if 'N' not in sequence.upper():
                print(">" + current_id + " |" + genomeid_to_desc[seq_id])
                print(sequence)
            else: 
                print(f"{current_id} has Ns",file=sys.stderr)       
def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="""extract fluAB cds subset:
            for the alternative splicing gene, get the biggest rang contain mutiple splicing location;
            get rid of partial cds;
            get rid of sequences containing any Ns"""
    )
    # Adding optional argument
    parser.add_argument("-c", "--cds_file", help = "cds fasta file", required=True)
    parser.add_argument("-g", "--genome_file", help = "genome fasta file", required=True)
    parser.add_argument("--maxambigs", type=int, default=0, help = "max number of ambigs bases", required=False)
   
    # Parse command-line arguments
    args = parser.parse_args()
    genomeid_to_desc = extract_id_desc(args.genome_file)
    seq_range_dict_multi, seq_range_dict_single = extract_cds_coordinates_from_fasta(args.cds_file)
    extract_cds_without_alternative_splicing(seq_range_dict_single, args.cds_file, genomeid_to_desc)
    extract_cds_from_genome(seq_range_dict_multi, args.genome_file, genomeid_to_desc)

if __name__ == "__main__":
    main()
