#!/usr/bin/env python

import argparse
import csv

def main():

    description = "this program is change the chrom id in the bed file"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-b",
        "--bed",
        required=True,
        help=f"bed format file \n",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help=f"fasta header\n",
    )
    parser.add_argument("-o", "--output", required=True, help=f"Output file name\n")
    args = parser.parse_args()
    
    seqname_dict = {}

    #MH356668 Human|7|M|H1N1|Kenya|A/Kenya/035/2018|A|na|na|na
    with open(args.input) as file:
        
        reader = csv.reader(file, delimiter=' ')
        for line in reader:
            print(line[0])
            print(line[1])
            seqname_dict[line[0]] = line[1]
        
    file.close()
   
    fout = open(args.output, "w", encoding="utf-8")
    with open(args.bed) as file:
        reader = csv.reader(file, delimiter="\t")
        #JQ396193        1       3215
        #JQ396193        2       3224
        for line in reader:
            refid = line[0]
            print(".......")
            print(line[0])
            refdesc = seqname_dict[refid]
            segment_number = refdesc.split("|")[1]
            line[0] = f"seg{segment_number}-{refid}"
            fout.write("\t".join(line) + "\n")

    file.close()
    fout.close()       

if __name__ == "__main__":
    main()
