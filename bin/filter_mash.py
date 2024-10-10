#!/usr/bin/env python

import argparse
import csv, operator

def main():
    description = """
        program is used filter out the rows from the mash screen out:
            only keep the row with the best identity for each named segment
            for each type, for example, H1, H3, we keep one references which passed
            identity cutoff 
        """
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help=f"The screen file name from the mash screen output\n",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help=f"The filtered output file\n",
    )
    parser.add_argument(
        "-m",
        "--min_shared_hashes",
        type=int,
        default=100,
        help=f"The min shared hash value\n",
    )

    args = parser.parse_args()

    # [identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment]
    # 0.901759        114/1000        907     0       AB239126        Human|6|NA|H5N1|Viet_Nam|A/Hanoi/30408/2005|N1|na|na|na
    # 0.931235        224/1000        1046    0       AB434290        Pig|6|NA|H1N1|Thailand|A/swine/Ratchaburi/NIAH1481/2000|N1|na|na|na
    # open .tsv file

    f = open(f"{args.output}", "w", encoding="utf-8")
    #segments = {"PB1": False, "PB2": False, "PA": False, "HA": False, "NP": False, "NA": False, "M": False, "NS": False}
    segments = {}
    header = ['identity', 'shared-hashes', 'median-multiplicity', 'p-value', 'query-ID', 'query-comment']
    f.writelines("\t".join(header) + "\n")
    
    with open(args.input) as file:
        csvData = csv.reader(file, delimiter="\t")

        # csv data sorted by the identity
        csvData = sorted(
            csvData,
            key=operator.itemgetter(0),
            reverse=True,
        )
       
        for line in csvData:
            #shared-hashes
            hashes = line[1].split("/")
            shared_hashes = int(hashes[0])
            # parse query-comment
            fields = line[5].split("|")
            segName = fields[2]
            segType = fields[3]
            key = segName + "_" + segType
            # only keep the best hit for each type of the segment
            #if not segments[segName]:
            if shared_hashes >= args.min_shared_hashes and key not in segments:
                segments[key] = True
                f.writelines("\t".join(line) + "\n")
    f.close()


if __name__ == "__main__":
    main()
