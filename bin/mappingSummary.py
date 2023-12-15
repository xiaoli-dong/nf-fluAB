#!/usr/bin/env python
import argparse
import os
import csv


# [identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment]
# 0.901759        114/1000        907     0       AB239126        Human|6|NA|H5N1|Viet_Nam|A/Hanoi/30408/2005|N1|na|na|na
def parse_mash_screen_file(path_to_file):
    if os.path.exists(path_to_file) == False or os.path.isfile(path_to_file) == False:
        print(f"\nERROR: mash screen file {path_to_file} does not exist.\n")
        exit(1)
    if os.path.getsize(path_to_file) == 0:
        print(f"\nERROR: mash screen file {path_to_file} is empty.\n")
        exit(1)
    ref2segid = {}
    segments = {
        "1_PB1": {},
        "2_PB2": {},
        "3_PA": {},
        "4_HA": {},
        "5_NP": {},
        "6_NA": {},
        "7_M": {},
        "8_NS": {},
    }

    # print(path_to_file)
    with open(path_to_file, "r", encoding="utf8") as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter="\t")
        # print(tsv_reader)
        for tsv in tsv_reader:
            query_comment = tsv[5].split("|")
            segName = query_comment[2]
            segid = query_comment[1]
            ref2segid[tsv[4]] = f"{segid}_{segName}"
            segments[f"{segid}_{segName}"] = tsv

    return segments, ref2segid


def parse_samtools_coverage_file(path_to_file):
    if os.path.exists(path_to_file) == False or os.path.isfile(path_to_file) == False:
        print(f"\nERROR: coverage file {path_to_file} does not exist.\n")
        exit(1)
    if os.path.getsize(path_to_file) == 0:
        print(f"\nERROR: coverage file {path_to_file} is empty.\n")
        exit(1)

    refs = {}

    with open(path_to_file, "r", encoding="utf8") as coverage_file:
        tsv_reader = csv.DictReader(coverage_file, delimiter="\t")

        for tsv in tsv_reader:
            rname = tsv.pop("#rname")
            refs[rname] = tsv

    # print(refs)
    return refs


def main():
    description = "Parse mash scree output report and samtool coverage output report and merge the information together"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-s",
        "--sample-name",
        required=True,
        help=f"Sample name\n",
    )
    parser.add_argument(
        "-c",
        "--samtool-coverage-file",
        required=True,
        help=f"Path to samtools coveage output file with header\n",
    )
    parser.add_argument(
        "-m",
        "--mash-screen-file",
        required=True,
        help=f"Path to mash screen output file\n",
    )

    parser.add_argument(
        "-t",
        "--tsv-output",
        required=False,
        help=f"The output file name for tabular format\n",
    )

    args = parser.parse_args()
    results = parse_mash_screen_file(args.mash_screen_file)

    segments_dict = results[0]
    # jsonString = json.dumps(segments_dict, indent=4)
    # print(jsonString)
    ref2segid_dict = results[1]

    cov_dict = parse_samtools_coverage_file(args.samtool_coverage_file)

    sorted_ref2segid = dict(sorted(ref2segid_dict.items(), key=lambda x: x[1], reverse=False))
    fieldnames = [
        # "sname",
        "segmentid",
        # "segment_name",
        # "accession",
        "startpos",
        "endpos",
        "numreads",
        "covbases",
        "coverage",
        "meandepth",
        "meanbaseq",
        "meanmapq",
        "identity",
        "shared-hashes",
        "median-multiplicity",
        "p-value",
        "query-ID",
        "query-comment",
    ]
    # print("\t".join(fieldnames))
    # Open the file in writing mode
    with open(f"{args.tsv_output}", "w", newline="") as f:
        # using csv.writer method from CSV package
        write = csv.writer(f, delimiter="\t")
        write.writerow(fieldnames)
        for acc in sorted_ref2segid.keys():
            if acc in cov_dict.keys():
                # segname_segid
                x = ref2segid_dict[acc].split("_")
                print(x)
                line = [
                    # args.sample_name,
                    # args.sample_name + "_" + "segment_" + x[0],
                    args.sample_name + "_" + "segment" + x[0],
                    # x[1],
                    # acc,
                    cov_dict[acc]["startpos"],
                    cov_dict[acc]["endpos"],
                    cov_dict[acc]["numreads"],
                    cov_dict[acc]["covbases"],
                    cov_dict[acc]["coverage"],
                    cov_dict[acc]["meandepth"],
                    cov_dict[acc]["meanbaseq"],
                    cov_dict[acc]["meanmapq"],
                ]

                # get rid of accession
                screen_line = segments_dict[ref2segid_dict[acc]]
                # del screen_line[4]
                line.extend(screen_line)
                write.writerow(line)

    f.close()

    # jsonString = json.dumps(segments, indent=4)
    # print(jsonString)


if __name__ == "__main__":
    main()
