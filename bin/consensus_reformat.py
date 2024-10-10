#!/usr/bin/env python
import argparse
from collections import defaultdict
import re

#>OQ367231 Human|3|PA|H1N1|Germany|A/Rheinland-Pfalz/USAFSAM-13946/2022|na|na|na|2151 
def readFasta(fasta):
    id2seq = defaultdict(str)
    id2desc = defaultdict(str)
    with open (fasta, 'r+') as in_f:
        seqid = ""
        for line in in_f:   
            if line[0] == ">":
                header = line[1:].split()
                seqid = header[0]
                id2desc[seqid] = header[1]
                id2seq[seqid] = ""  
            else:
                # Remove all whitespace characters (including spaces, tabs, newlines)
                id2seq[seqid] += re.sub(r'\s+', '', line)
    return id2seq, id2desc

def reheader(id2desc, seqPrefix):
    id2newid = defaultdict(str)
    # seqprefix example: runid-mapping_basecaller: 230721_S_I_314-S87-bwa_freebayes
    #Human|3|PA|H1N1|Germany|A/Rheinland-Pfalz/USAFSAM-13946/2022|na|na|na|2151 
    for seqid, desc in id2desc.items():
        # Split the string by the pipe character '|'
        parts = desc.split('|')
        segid = parts[1]
        newid = f"{seqPrefix}-seg_{segid}-ref_{seqid}"
        id2newid[seqid] = newid
    return id2newid


# remove contigs with too may Ns, default is 20%
def filterSeq(id2seq, maxambigs):
    allowed_bases=['A','T','G','C']
    for seqid in list(id2seq):
    #for seqid, seq in id2seq.items():
        seq = id2seq[seqid]
        total_atcg_count = len([b for b in seq if b in allowed_bases])
        seqlen = len(seq)
        if ((seqlen - total_atcg_count)/seqlen) > maxambigs:
            del id2seq[seqid]
    return id2seq

def main():
    # Initialize parser
    parser = argparse.ArgumentParser()

    # Adding optional argument
    parser.add_argument("-f", "--fasta", help = "fasta", required=True)
    parser.add_argument("-p", "--prefix", required=True, help=f"Prefix of the sequence id\n")
    parser.add_argument("--maxambigs", type=float, default=0.25, help = "max ratio of the ambiguity bases (default: 0.25)")

    # Read arguments from command line
    args = parser.parse_args()

    id2seq, id2desc = readFasta(args.fasta)
    id2newid = reheader(id2desc, args.prefix)
    final_id2seq = filterSeq(id2seq, args.maxambigs)

    # Sort the dictionary by its values in ascending order
    sorted_id2newid = sorted(id2newid.items(), key=lambda item: item[1])
    for key, value in sorted_id2newid:
        if key in final_id2seq:
            print(f">{value}")
            print(final_id2seq[key])


if __name__ == "__main__":
    main()
