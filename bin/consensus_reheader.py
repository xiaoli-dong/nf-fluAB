#!/usr/bin/env python

import argparse
import csv

def main():

    description = "this program is reformat the flu segment sequence header"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help=f"seqkit fx2tab --length -C N -H output \n",
    )
    parser.add_argument("-o", "--output", required=True, help=f"Output file name\n")
    parser.add_argument("-p", "--prefix", required=True, help=f"Prefix of the sequence id\n")
    parser.add_argument("-c", "--cov", required=True, default=0, help=f"minimum ratio of the non_Ns\n")



    args = parser.parse_args()
    

    fout = open(args.output, "w", encoding="utf-8")
    # Simple Way to Read TSV Files in Python using csv
    ##name   seq     length  N
    #MH356668 Human|7|M|H1N1|Kenya|A/Kenya/035/2018|A|na|na|na       ATGAGTCTTCTAACCGAGG       689     169
    # open .tsv file
    with open(args.input) as file:
        
        reader = csv.reader(file, delimiter="\t")
        
        # printing data line by line
        header = next(reader)
        fout.write("\t".join(header) + "\n")
        #fout.write(header)

        for l in reader:
          
            if int(l[3])/int(l[2]) > 1 - float(args.cov):
                #print(l)
                continue
            
            else: 
                refid = l[0].split(" ")[0]
                segment_number, segment_name = l[0].split("|")[1:3]
                l[0] = f"{args.prefix}-segment_{segment_number}-ref_accession_{refid}"
                fout.write("\t".join(l) + "\n")

    file.close()
    fout.close()       

if __name__ == "__main__":
    main()
