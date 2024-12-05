#!/usr/bin/env python
import argparse
import csv
from collections import defaultdict
import re
import glob

__author__ = 'Xiaoli Dong'
__version__ = '0.1.0'
__maintainer__ = 'Xiaoli Dong'
__email__ = 'xiaoli.dong@albertaprecisionlabs.ca'
__status__ = 'Dev'


def filterFasta(file_path, minlen, maxlen, maxambigs):
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
                    filter(header, sequence, int(minlen), int(maxlen), int(maxambigs))
                       
                header = line[1:]  # Remove the '>' character
                sequence_lines = []
            else:
                sequence_lines.append(line)
        
        # Save the last sequence
        if header:
            sequence = ''.join(sequence_lines)
            filter(header, sequence, int(minlen), int(maxlen), int(maxambigs))
           
        
    #return sequences

def filter(header, sequence, minlen, maxlen, maxambigs):
    
    allowed_bases=['A','T','G','C']
    total_atcg_count = len([b for b in sequence if b in allowed_bases])
    if len(sequence) - total_atcg_count <= maxambigs  and len(sequence) >= minlen and len(sequence) <= maxlen:
        print(">" + header)
        print(sequence)
    

# Initialize parser
parser = argparse.ArgumentParser()

# Adding optional argument
parser.add_argument("-f", "--fasta", help = "fasta", required=True)
parser.add_argument("--minlen", default=700, help = "min sequence length (default: 700)")
parser.add_argument("--maxlen", default=3000, help = "max sequence length (default: 3000)")
parser.add_argument("--maxambigs", default=0, help = "max number of ambigs bases (default: 0)")

# Read arguments from command line
args = parser.parse_args()

fasta = filterFasta(args.fasta, args.minlen, args.maxlen, args.maxambigs)
