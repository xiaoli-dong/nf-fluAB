import argparse

def extract_na_ids_from_fasta(file_path):
    ids_with_na = []
    id2header = {}
    
    # Open the FASTA file
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        # Process each line in the file
        for line in lines:
            # Remove leading/trailing whitespace
            line = line.strip()
            
            # Only process header lines (those starting with '>')
            if line.startswith('>'):
                # Get the sequence ID (part after '>' and before space)
                sequence_id = line[1:].split(' ')[0]  # Strip '>' and split by space to get the ID
                id2header[sequence_id] = line[1:].split(' ')[1] 
                # Split by '|' and check if the subtype is 'na' (4th element in the parts list)
                parts = line.split('|')
                if parts[3].lower() == 'na':
                    # Add sequence ID to the list if the subtype is 'na'
                    ids_with_na.append(sequence_id)
    
    return ids_with_na, id2header

def parse_mmseq2_tsv(tsv, headers):
    
    # Initialize an empty dictionary to store the clusters
    cluster_dict = {}

    # Open and read the input file
    with open(tsv, 'r') as file:
        for line in file:
            # Remove leading/trailing whitespace and split the line into cluster_id and member_id
            cluster_id, member_id = line.strip().split()
            # If the cluster ID is not in the dictionary, add it with an empty list
            if cluster_id not in cluster_dict:
                cluster_dict[cluster_id] = []
            
            # Append the member_id to the list of the corresponding cluster_id
            cluster_dict[cluster_id].append(member_id)
    
    for cluster_id, members in cluster_dict.items():
        segid2seqs = {}
        rep_header = headers[cluster_id]
        parts_rep = rep_header.split('|')
        segid_rep = parts_rep[1]
        if segid_rep not in segid2seqs:
            segid2seqs[segid_rep] = []
            
        segid2seqs[segid_rep].append(cluster_id)
        for member in members:
            if member == cluster_id:
                continue
            member_header = headers[member]
            parts_member = member_header.split('|')
            segid_member = parts_member[1]
            if segid_member not in segid2seqs:
                segid2seqs[segid_member] = []
            segid2seqs[segid_member].append(member)
        
        #print(f"cluster id {cluster_id}")
        # Loop through each cluster's values (the lists of members)
        if len(segid2seqs) != 1:
            print(f"cluster id {cluster_id} segid {segid_rep}")
            for segid, seqs in segid2seqs.items():
                print(f"segid {segid} {len(seqs)}")
                if segid != segid_rep:
                    print(f"segid {segid} {seqs}")
                # # Loop through each member in the list of members for this cluster
                # for seq in seqs:
                #     print(f"  seq: {seq}")
            print("*********************************")      
        #print()  # Print a blank line between clusters for readability
            
    return cluster_dict
    """ # Print the resulting hash table (cluster_id -> list of member_ids)
    for cluster_id, members in cluster_dict.items():
        print(f"Cluster ID: {cluster_id}, Members: {members}") """


def main():
    # Set up argparse to handle command-line input
    parser = argparse.ArgumentParser(description="Extract sequence IDs with 'na' as the subtype from a FASTA file.")
    #parser.add_argument("fasta_file", help="Path to the input FASTA file")
    parser.add_argument("--fasta_rep", help = "rep fasta after mmseq2 cluster", required=True)
    parser.add_argument("--fasta_all", help = "input fasta to mmseq2 cluster", required=True)
    parser.add_argument("--tsv", help = "tsv from mmseq2 cluster program", required=True)
    # Parse the command-line arguments
    args = parser.parse_args()
    
    # Extract IDs with 'na' as subtype
    na_ids_rep, headers_rep = extract_na_ids_from_fasta(args.fasta_rep)
    na_ids_all, headers_all = extract_na_ids_from_fasta(args.fasta_all)
    cluster_dict = parse_mmseq2_tsv(args.tsv, headers_all)
    for cluster_id, members in cluster_dict.items():
        
        if cluster_id not in na_ids_rep:
            continue
        print(f"Cluster ID: {cluster_id} missing subtype, check ...")
        # Loop through each member in the list of members for this cluster
        for member in members:
            if member not in na_ids_all:
                print(f"replace {cluster_id} with the Member: {member}")
                break

    """ # Print the result
    if na_ids_rep:
        print("Sequence IDs with 'na' as subtype:")
        for seq_id in na_ids_rep:
            print(seq_id)
    else:
        print("No sequences with 'na' as subtype found.")
 """
if __name__ == "__main__":
    main()
