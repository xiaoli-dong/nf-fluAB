#!/usr/bin/env python

import argparse
import csv, operator
import json


def main():

    description = (
        "program is used filter out the rows from the mash screen out, whose access ids are in the exclusion list; "
    )
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-i", "--input", required=True, help=f"The output file name from the mash screen command\n")
    parser.add_argument("-p", "--prefix", required=True, help=f"Prefix of the output file\n")
    parser.add_argument(
        "-e",
        "--excludes",
        help=f"a file contains a list of ncbi accession ids, which are needed to be excluded\n",
    )

    args = parser.parse_args()
    ex_list = []

    if args.excludes is not None:
        # Open the file in read mode
        with open(args.excludes, "r") as file:
            # Iterate over the lines of the file
            for line in file:
                # Remove the newline character at the end of the line
                line = line.strip()
                # Append the line to the list
                ex_list.append(line)
        # Print the list of lines
        # print(ex_list)

    # [identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment]
    # 0.901759        114/1000        907     0       AB239126        Human|6|NA|H5N1|Viet_Nam|A/Hanoi/30408/2005|N1|na|na|na
    # 0.931235        224/1000        1046    0       AB434290        Pig|6|NA|H1N1|Thailand|A/swine/Ratchaburi/NIAH1481/2000|N1|na|na|na
    # open .tsv file

    f = open(f"{args.prefix}.screen", "w", encoding="utf-8")
    segments = {"PB1": False, "PB2": False, "PA": False, "HA": False, "NP": False, "NA": False, "M": False, "NS": False}

    with open(args.input) as file:

        csvData = csv.reader(file, delimiter="\t")

        # csv data sorted by the identity
        csvData = sorted(
            csvData,
            key=operator.itemgetter(0),
            reverse=True,
        )
        # filter out the exclusions
        for line in csvData:
            query_id = line[4]
            if query_id in ex_list:
                continue
            else:
                # parse query-comment
                fields = line[5].split("|")
                segName = fields[2]
                if not segments[segName]:
                    segments[segName] = True
                    f.writelines("\t".join(line) + "\n")
    f.close()


if __name__ == "__main__":
    main()
