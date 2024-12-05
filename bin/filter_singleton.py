#!/usr/bin/env python
import sys
import argparse
import csv
from collections import defaultdict
import re
from csv import DictReader

__author__ = 'Xiaoli Dong'
__version__ = '0.1.0'
__maintainer__ = 'Xiaoli Dong'
__email__ = 'xiaoli.dong@albertaprecisionlabs.ca'
__status__ = 'Dev'

def filterFasta(file_path, singletons):
    with open(file_path, 'r') as file:
        
        header = None
        sequence_lines = []
        
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    # Save the last sequence
                    # remove all the whitespace characters (space, tab, newlines)
                    sequence = ''.join(sequence_lines)
                    filter_seq(header, sequence, singletons)
                       
                header = line[1:]  # Remove the '>' character
                sequence_lines = []
            else:
                sequence_lines.append(line)
        
        # Save the last sequence
        if header:
            sequence = ''.join(sequence_lines)
            filter_seq(header, sequence, singletons)
           
        
    #return sequences

def filter_seq(header, sequence, singletons):
    
    seqid = header.split()[0]
    if seqid in singletons:
        print(f"{seqid} is a singleton", file=sys.stderr)
        
    else:
        # If sequence passes all filters
        print(f">{header}")
        print(sequence)

def read_cluster_csv(file_path_tsv):

    from collections import Counter

    # Read the file and process the clusters
    with open(file_path_tsv, "r") as f:
        # Assuming each line contains a cluster ID in the first column, separated by tabs
        clusters = [line.strip().split("\t")[0] for line in f]

    # Count occurrences of each cluster
    cluster_counts = Counter(clusters)

    # Find singletons (clusters with only one member)
    singletons = [cluster for cluster, count in cluster_counts.items() if count == 1]

    # Output results
    print(
        f"Singleton clusters: {singletons}",
        file=sys.stderr
    )
    
    print(
        f"Total number of singleton clusters: {len(singletons)}",
        file=sys.stderr
    )
    return singletons

# Initialize parser
parser = argparse.ArgumentParser()

# Adding optional argument
parser.add_argument("-f", "--fasta", help = "fasta", required=True)
parser.add_argument("-t", "--tsv", help = "mmseqs cluster tsv file", required=True)


# Read arguments from command line
args = parser.parse_args()

#read mmseqs cluster tsv file in the format of 
#cluster-representative cluster-member
#Q0KJ32 Q0KJ32
#Q0KJ32 C0W539
#Q0KJ32 D6KVP9
#E3HQM9 E3HQM9
#E3HQM9 F0YHT8
singletons = read_cluster_csv(args.tsv)

fasta = filterFasta(args.fasta, singletons)



