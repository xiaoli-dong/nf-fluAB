#!/usr/bin/env python
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

def extract_cds_from_fasta(fasta, cds_coords, minlen, maxlen, maxambigs):

    with open(fasta, 'r') as file:    
        header = None
        sequence_lines = []
        
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    # Save the last sequence
                    # remove all the whitespace characters (space, tab, newlines)
                    sequence = ''.join(sequence_lines)
                    extract(header, sequence, cds_coords, minlen, maxlen, maxambigs)
                   
                       
                header = line[1:]  # Remove the '>' character
                sequence_lines = []
            else:
                sequence_lines.append(line)
        
        # Save the last sequence
        if header:
            sequence = ''.join(sequence_lines)
            extract(header, sequence, cds_coords, minlen, maxlen, maxambigs)
def extract(header, sequence, cds_coords, minlen, maxlen, maxambigs):
    
    allowed_bases=['A','T','G','C']
    seqid = header.split()[0]
    seqid_no_ver = seqid.split('.')[0]
    subseq=''

    if seqid in cds_coords:
        
        if cds_coords[seqid]['strand'] == '-':
            # start_index = seq_to -1
            # end_index = seq_from
            start_index = cds_coords[seqid]['seq_to'] -1 
            end_index = cds_coords[seqid]['seq_from']
            
        elif cds_coords[seqid]['strand'] == '+':
            # start_index = seq_from -1
            # end_index = seq_to
            start_index = cds_coords[seqid]['seq_from'] -1 
            end_index = cds_coords[seqid]['seq_to']

        subseq = sequence[start_index:end_index]
        total_atcg_count = len([b for b in subseq if b in allowed_bases])
        subseq_len = len(subseq)

        if (subseq_len - total_atcg_count) <= maxambigs  and subseq_len >= minlen and subseq_len <= maxlen:
            #print(">" + seqid_no_ver + " " + metadata_dict[seqid_no_ver] + "|" + str(subseq_len))
            print(">" + header + "|" + str(subseq_len))
            print(subseq)


def parseSGM(sgm):
    
    cds_coord_pass = defaultdict(dict)
    #pattern = r'^\d+\.1\.1\s+.*'
   
    with open (sgm, 'r+') as in_f:
      
        for line in in_f:

            # match = re.search(pattern, line)
            # if match:
            if not line.startswith("#"):
                #print("Pattern found:", match.group())
                line_list = line.split()
                seq_from = int(line_list[10])
                seq_to = int(line_list[11])
                sgm_len = int(line_list[14])

                strand = line_list[15]
                trc = line_list[16]
                pf = line_list[3]
                seq_name = line_list[1]

                if trc == 'no' and pf == 'PASS' :
                    #print(line)
                    if seq_name not in cds_coord_pass :
                        #cds_coord_pass[seq_name] = seq_from + '..' + seq_to + ':' + strand
                        cds_coord_pass[seq_name]['seq_from'] = seq_from
                        cds_coord_pass[seq_name]['seq_to'] = seq_to
                        cds_coord_pass[seq_name]['strand'] = strand
                        cds_coord_pass[seq_name]['sgm_len'] = sgm_len

                    elif cds_coord_pass[seq_name]['sgm_len'] < sgm_len:
                        #cds_coord_pass[seq_name] = seq_from + '..' + seq_to + ':' + strand
                        cds_coord_pass[seq_name]['seq_from'] = seq_from
                        cds_coord_pass[seq_name]['seq_to'] = seq_to
                        cds_coord_pass[seq_name]['strand'] = strand
                        cds_coord_pass[seq_name]['sgm_len'] = sgm_len   
    return cds_coord_pass
               

# Initialize parser
parser = argparse.ArgumentParser()

# Adding optional argument
parser.add_argument("-f", "--fasta", help = "fasta", required=True)
parser.add_argument("--minlen", type=int, default=700, help = "min sequence length", required=False)
parser.add_argument("--maxlen", type=int, default=3000, help = "max sequence length", required=False)
parser.add_argument("--maxambigs", type=int, default=0, help = "max number of ambigs bases", required=False)
parser.add_argument("--sgm", help = "sgm output from vadr program", required=True)
parser.add_argument("--bvbrc", help = "bvbrc csv file", required=True)


# Read arguments from command line
args = parser.parse_args()

cds_coord_pass = parseSGM(args.sgm)

extract_cds_from_fasta(args.fasta, cds_coord_pass, args.minlen, args.maxlen, args.maxambigs)

