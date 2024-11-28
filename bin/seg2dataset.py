#!/usr/bin/env python

import argparse
import csv
from pathlib import Path

from itertools import groupby


def get_flu_typedata(type_str):
    switcher = {
        "H1": "flu_h1n1pdm_ha",
        "H3": "flu_h3n2_ha",
        "Victoria": "flu_vic_ha",
        "Yamagata": "flu_yam_ha",
        #h5nx input dataset is downloaded from https://github.com/nextstrain/nextclade_data/tree/master/data/community/moncla-lab/iav-h5/ha
        "H5": "flu_h5nx_ha"
    }
    # get() method of dictionary data type returns
    # value of passed argument if it is present
    # in dictionary otherwise second argument will
    # be assigned as default value of passed argument
    return switcher.get(type_str, "NA")


def fasta_iter(fasta_name):
    """
    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    given a fasta file. yield tuples of header, sequence

    Fasta iterator from https://www.biostars.org/p/710/#120760
    """
    with open(fasta_name) as fh:
        # ditch the boolean (x[0]) and just keep the header or sequence since
        # we know they alternate.
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            headerStr = header.__next__()[1:].strip()

            # join all sequence lines to one.
            seq = "".join(s.strip() for s in faiter.__next__())

            yield (headerStr, seq)


def splitfasta_by_id(fasta, input_dir):
    fiter = fasta_iter(fasta)

    for ff in fiter:
        fasta_header, seq = ff
        # print(fasta_header)
        # print(seq)
        # Use first ID as separated by spaces as the "sequence name"
        # (equivalent to "chromosome" in other cases)
        seqname = fasta_header.split()[0]

        with open(f"{input_dir}/{seqname}.segcontig.fa", "w") as contig_file:
            contig_file.write(f">{fasta_header}\n{seq}")
        contig_file.close()


def main():
    description = "program is used filter out the rows from the mash screen out, whose access ids are in the exclusion list; then output the best record for each flu segments"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-t",
        "--tsv",
        required=True,
        help=f"blastn format 6 output '6 std qlen slen qcovs qcovhsp'\n",
    )

    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        help=f"Contig fasta file\n",
    )
    parser.add_argument(
        "-d",
        "--dbdir",
        default="./",
        help=f"Typedata directory\n",
    )
    parser.add_argument(
        "-i",
        "--sid",
        default="./",
        help=f"Sample id\n",
    )

    args = parser.parse_args()

    input_dir = Path(args.fasta).parent.absolute()
    # print(input_dir)

    # split fasta contig to its own file. The file named as sequence id
    splitfasta_by_id(args.fasta, input_dir)

    #print(f"seqid\tfasta_path\ttypedata")
    empty = True

    with open(args.tsv) as tsv_file:
        reader = csv.reader(tsv_file, delimiter="\t")
        next(reader, None)  # skip the headers
        # filter out the exclusions
        for line in reader:
            
            # print(line)
            # S10_T1_segment_6        GQ377078~~6~~N1 96.373  1158    39      2       205     1359    252     1409    0.0     1903    1360    1410    85      85

            (target, segid, typing) = line[1].split("~~")
            typedata = get_flu_typedata(typing)
            if typedata != "NA":
                print(f"{line[0]}\t{input_dir}/{line[0]}.segcontig.fa\t{args.dbdir}/{typedata}")
                empty = False

    """ if empty:
        print(f"{args.sid}_segment_1\t{input_dir}/{args.sid}_segment_1.segcontig.fa\t{args.dbdir}/NA") """

    tsv_file.close()


if __name__ == "__main__":
    main()
