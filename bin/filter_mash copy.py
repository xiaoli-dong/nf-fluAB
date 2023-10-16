#!/usr/bin/env python

import argparse
import csv, operator
import json


def main():

    description = "program is used filter out the rows from the mash screen out, whose access ids are in the exclusion list; then output the best record for each flu segments"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-i", "--input", required=True, help=f"The output file name from the mash screen command\n")
    parser.add_argument("-p", "--prefix", required=True, help=f"Prefix of the output file\n")
    parser.add_argument("-s", "--sid", required=True, help=f"sampleid\n")

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
    segments_ref = {
        "1_PB1_accession": "NA",
        "1_PB1_identity": "NA",
        "1_PB1_shared-hashes": "NA",
        "2_PB2_accession": "NA",
        "2_PB2_identity": "NA",
        "2_PB2_shared-hashes": "NA",
        "3_PA_accession": "NA",
        "3_PA_identity": "NA",
        "3_PA_shared-hashes": "NA",
        "4_HA_accession": "NA",
        "4_HA_identity": "NA",
        "4_HA_shared-hashes": "NA",
        "5_NP_accession": "NA",
        "5_NP_identity": "NA",
        "5_NP_shared-hashes": "NA",
        "6_NA_accession": "NA",
        "6_NA_identity": "NA",
        "6_NA_shared-hashes": "NA",
        "7_M_accession": "NA",
        "7_M_identity": "NA",
        "7_M_shared-hashes": "NA",
        "8_NS_accession": "NA",
        "8_NS_identity": "NA",
        "8_NS_shared-hashes": "NA",
    }

    # [identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment]
    # 0.901759        114/1000        907     0       AB239126        Human|6|NA|H5N1|Viet_Nam|A/Hanoi/30408/2005|N1|na|na|na
    # 0.931235        224/1000        1046    0       AB434290        Pig|6|NA|H1N1|Thailand|A/swine/Ratchaburi/NIAH1481/2000|N1|na|na|na
    # open .tsv file
    segments = {"PB1": False, "PB2": False, "PA": False, "HA": False, "NP": False, "NA": False, "M": False, "NS": False}

    f = open(f"{args.prefix}.screen", "w", encoding="utf-8")

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
            query_comment = line[5]
            if line[4] in ex_list:
                continue
            else:
                fields = query_comment.split("|")
                segName = fields[2]
                if not segments[segName]:
                    segments[segName] = True
                    # acc_list.append(query_id)
                    l = [fields[1], fields[2], line[4], line[0], line[1]]
                    f.writelines("\t".join(l) + "\n")

                    segments_ref[f"{fields[1]}_{segName}_accession"] = line[4]
                    segments_ref[f"{fields[1]}_{segName}_identity"] = line[0]
                    segments_ref[f"{fields[1]}_{segName}_shared-hashes"] = line[1]
    f.close()
    # Open the file in writing mode
    with open(f"{args.prefix}.ref.csv", "w", newline="") as f:
        # Create the writer with the dictionary keys as headers
        fieldnames = sorted(segments_ref.keys())
        fieldnames.insert(0, "sampleid")
        segments_ref["sampleid"] = args.sid

        writer = csv.DictWriter(f, fieldnames)
        # Write the header defined in the fieldnames argument
        writer.writeheader()
        # Write one or more rows
        writer.writerow(segments_ref)

    # jsonString = json.dumps(segments_ref, indent=4)
    # print(jsonString)

    f.close()


if __name__ == "__main__":
    main()
