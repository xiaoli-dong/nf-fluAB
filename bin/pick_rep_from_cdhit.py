#!/usr/bin/env python
import argparse
import re
import sys
import logging
from collections import defaultdict
import traceback


def parse_cd_hit_clstr(clstr_file):
    """
    Parses a cd-hit cluster file and returns a list of representative sequences for each cluster.
    """
    clusters = []
    current_cluster = None

    try:
        with open(clstr_file, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith(">Cluster"):
                    # New cluster found, append the previous one if it exists
                    if current_cluster:
                        clusters.append(current_cluster)
                    current_cluster = {'cluster_id': line, 'sequences': []}
                else:
                    # Sequence line in the cluster
                    current_cluster['sequences'].append(line)

            # Append the last cluster if it exists
            if current_cluster:
                clusters.append(current_cluster)
    except FileNotFoundError:
        logging.error(f"The file {clstr_file} does not exist.")
        exit(1)

    return clusters


def parse_fasta_header(header):
    # Remove the leading '>' character from the header
    header = header.lstrip('>')

    # Split at the first space to separate the sequence ID from the description
    parts = header.split(' ', 1)  # Split into at most 2 parts

    # The first part is the sequence ID, the rest is the description
    seq_id = parts[0]
    # Handle case where there's no description
    description = parts[1] if len(parts) > 1 else ""

    return seq_id, description


def extract_subtype_from_fasta(fasta_file):
    id2subtype = {}
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()  # Remove leading/trailing whitespace

            if line.startswith('>'):  # Header line (starts with '>')

                seqid, desc = parse_fasta_header(line)
                subtype = desc.split('|')[3]
                id2subtype[seqid] = subtype  # Add to dictionary
    return id2subtype


def extract_sequence_details(seq):
    """
    Extracts sequence ID, length, and subtype using a regular expression.
    """
    # 0       2233nt, >NC_007371.1... *
    # 1       2233nt, >NC_007376.1... at 1:2233:1:2233/+/93.01%
    # pattern = r"(\d+)nt, >(\S+)(?:--.*?\|.*?\|.*?\|([^|]+))"
    pattern = r'(\d+)(?=nt).*?>([^...]+)'
    match = re.search(pattern, seq)

    if match:
        length = int(match.group(1))
        seqid = match.group(2)
        # subtype = match.group(3)
        return length, seqid  # , subtype
    return None


def get_representative_sequence(cluster, id2subtype):
    """
    For a given cluster, selects the representative sequence based on specific criteria.
    """
    # Initialize variables
    rep_seq = None
    seq_with_subtype = []
    subtype_count = defaultdict(int)
    default_rep = None
    abundant_subtype = None
    abundant_count = 0

    # Process sequences in the cluster
    for seq in cluster['sequences']:
        if seq.endswith('... *'):
            default_rep = seq

        seq_details = extract_sequence_details(seq)
        if seq_details:
            length, seqid = seq_details
            subtype = id2subtype[seqid]
            subtype_count[subtype] += 1
            seq_with_subtype.append((length, seq, subtype))

    # Sort by subtype frequency (decreasing order)
    sorted_subtype_count = sorted(
        subtype_count.items(), key=lambda x: x[1], reverse=True)
    total_count = sum(subtype_count.values())

    # Find the most abundant subtype
    if sorted_subtype_count:
        abundant_subtype, abundant_count = sorted_subtype_count[0]

    # Sort sequences by length (descending)
    sorted_seq_with_subtype = sorted(
        seq_with_subtype, key=lambda x: x[0], reverse=True)

    # If sequences exist, select the longest as the initial representative
    if sorted_seq_with_subtype:
        rep_seq = sorted_seq_with_subtype[0][1]

    # Check if the cluster is mixed and needs inspection
    if abundant_count / total_count < 0.7:
        logging.info(
            "Cluster is mixed, may need further inspection and will be excluded for now.")
        rep_seq = None
        for subtype, count in sorted_subtype_count:
            proportion = count / total_count
            logging.info(
                f"Subtype: {subtype}, Proportion: {proportion:.2f}, Count: {count} out of {total_count}")

    # Ensure that the representative sequence matches the abundant subtype
    elif rep_seq and abundant_count / total_count >= 0.7:
        length, seqid = extract_sequence_details(rep_seq)
        rep_subtype = id2subtype[seqid]

        if abundant_subtype != rep_subtype:
            logging.info(
                f"Representative sequence {seqid} (subtype={rep_subtype}) does not match the abundant subtype ({abundant_subtype}).")
            logging.info(f"Subtype counts: {sorted_subtype_count}")

            # Iterate through sorted sequences to find one with the abundant subtype
            for length, seq, subtype in sorted_seq_with_subtype:
                if subtype == abundant_subtype:
                    logging.info(
                        f"Found representative with abundant subtype ({abundant_subtype}), count: {abundant_count} out of {total_count}.")
                    logging.info(seq)
                    rep_seq = seq
                    break

    # Return the final representative sequence (or default if no representative found)
    return rep_seq if rep_seq is not None else default_rep


def parse_fasta_header(header):
    # Remove the leading '>' character from the header
    header = header.lstrip('>')

    # Split at the first space to separate the sequence ID from the description
    parts = header.split(' ', 1)  # Split into at most 2 parts

    # The first part is the sequence ID, the rest is the description
    seq_id = parts[0]
    # Handle case where there's no description
    description = parts[1] if len(parts) > 1 else ""

    return seq_id, description


def extract_fasta_subset(input_file, seqids, output_file=None):
    """
    Extract a subset of sequences from a FASTA file based on the provided sequence IDs.

    :param input_file: str, path to the input FASTA file.
    :param seqids: list of str, list of sequence IDs to extract.
    :param output_file: str, optional, path to save the extracted sequences. If None, returns as a string.
    
    :return: None (if output_file is provided), otherwise the extracted sequences as a string.
    """
    extracted_sequences = []
    current_seqid = None
    current_sequence = []
    try:
        with open(input_file, 'r') as file:
            for line in file:
                line = line.strip()

                if line.startswith('>'):  # Header line (sequence ID)
                    # print(line)
                    # If the previous sequence ID was in the list, store the sequence
                    if current_seqid and current_seqid in seqids:
                        extracted_sequences.append(
                            f">{current_seqid} {current_sequence_desc}\n{''.join(current_sequence)}")

                    # Reset for the new sequence
                    current_seqid, current_sequence_desc = parse_fasta_header(
                        line)
                    current_sequence = []  # Reset sequence content
                else:
                    current_sequence.append(line)  # Append sequence lines

            # Handle the last sequence in the file
            if current_seqid and current_seqid in seqids:
                # extracted_sequences.append(f">{current_seqid}\n{''.join(current_sequence)}")
                extracted_sequences.append(
                    f">{current_seqid} {current_sequence_desc}\n{''.join(current_sequence)}")
    except Exception as e:
        print(f"Error reading input file: {e}")
        # Capture the detailed exception
        error_details = traceback.format_exc()

        print(f"Detailed traceback: {error_details}")
        sys.exit(1)

    # If an output file is specified, write the result to it
    if output_file:
        with open(output_file, 'w') as out_file:
            out_file.write("\n".join(extracted_sequences))
    else:
        # Otherwise, return the extracted sequences as a string
        return "\n".join(extracted_sequences)


def parse_clusters_and_select_representatives(clusters, id2subtype, fasta_file, output_file):
    """
    Parses the clusters and selects the representative sequences for each cluster.
    """
    representative_sequences = []
    representative_seqids = []

    for cluster in clusters:
        rep_seq = get_representative_sequence(cluster, id2subtype)
        if rep_seq is not None:
            representative_sequences.append(
                {'cluster_id': cluster['cluster_id'], 'seq': rep_seq})
            len, id = extract_sequence_details(rep_seq)
            representative_seqids.append(id)

    extract_fasta_subset(fasta_file, representative_seqids, output_file)
    return representative_sequences


def write_output_file(output_file, representative_sequences):
    """
    Writes the representative sequences to an output file.
    """
    try:
        with open(output_file, 'w') as out_file:
            for rep in representative_sequences:
                out_file.write(
                    f"Cluster ID: {rep['cluster_id']}, Sequence: {rep['seq']}\n")
        logging.info(f"Representative sequences saved to {output_file}")
    except Exception as e:
        logging.error(f"Error writing to file {output_file}: {e}")
        exit(1)


def print_representative_sequences(representative_sequences):
    """
    Prints the representative sequences to the console.
    """
    for rep in representative_sequences:
        print(f"Cluster ID: {rep['cluster_id']}, Sequence: {rep['seq']}")


def main():
    description = """
    This program filters out rows from the mash screen:
    It keeps only the row with the best identity for each named segment,
    and for each type (e.g., H1, H3), it keeps one reference that passed the identity cutoff.
    """
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-c", "--clstr", required=True, help="Input CLSTR file from cd-hit-est output."
    )
    parser.add_argument(
        "-f", "--fasta", required=True, help="Input fasta file used to generated cd-hit-est results"
    )
    parser.add_argument(
        "--out_clstr", required=True, help="Output file for representatives."
    )
    parser.add_argument(
        "--out_fasta", required=True, help="Output file for representatives."
    )

    parser.add_argument(
        "--print", action="store_true", help="Print the representative sequences to the console."
    )
    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(level=logging.INFO)
    id2subtype = extract_subtype_from_fasta(args.fasta)

    # Parse the CLSTR file and get the clusters
    clusters = parse_cd_hit_clstr(args.clstr)

    # Select representative sequences for each cluster
    representative_sequences = parse_clusters_and_select_representatives(
        clusters, id2subtype, args.fasta, args.out_fasta)

    # Write the representative sequences to the output file

    write_output_file(args.out_clstr, representative_sequences)

    # Optionally print the sequences to the console
    if args.print:
        print_representative_sequences(representative_sequences)


if __name__ == "__main__":
    main()
