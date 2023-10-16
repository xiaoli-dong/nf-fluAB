#!/usr/bin/env python

import argparse
import csv, operator
import subprocess


def main():

    description = "program is used filter out the rows from the mash screen out, whose access ids are in the exclusion list; then output the best record for each flu segments"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-m", "--mash_screen_output", required=True, help=f"The output file name from the mash screen command\n"
    )

    parser.add_argument(
        "-r",
        "--reference_db_fasta",
        required=True,
        help=f"Fasta file of the reference db used to generated mash_screen_output\n",
    )
    parser.add_argument("-p", "--output_prefix", required=True, help=f"Prefix of the output file\n")

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
    segments = {"PB1": False, "PB2": False, "PA": False, "HA": False, "NP": False, "NA": False, "M": False, "NS": False}

    # 0.901759        114/1000        907     0       AB239126        Human|6|NA|H5N1|Viet_Nam|A/Hanoi/30408/2005|N1|na|na|na
    # 0.931235        224/1000        1046    0       AB434290        Pig|6|NA|H1N1|Thailand|A/swine/Ratchaburi/NIAH1481/2000|N1|na|na|na
    # open .tsv file
    acc_list = []

    f = open(f"{args.output_prefix}_segments.screen", "w", encoding="utf-8")

    with open(args.mash_screen_output) as file:

        # Passing the TSV file to
        # reader() function
        # with tab delimiter
        # This function will
        # read data from file
        csvData = csv.reader(file, delimiter="\t")
        csvData = sorted(
            csvData,
            key=operator.itemgetter(0),
            reverse=True,
        )
        # filter out the exclusions
        for line in csvData:
            if line[4] in ex_list:
                continue
            else:
                # print(line)
                fields = line[5].split("|")
                # print(fields[2])
                if not segments[fields[2]]:
                    segments[fields[2]] = True
                    acc_list.append(line[4])
                    l = [fields[1], fields[2], line[4], line[0], line[1]]
                    f.writelines("\t".join(l) + "\n")
                    # print(fields[1], fields[2], line[4], line[0], line[1], sep="\t")

    f.close()

    patterns = ",".join(acc_list)

    terminal_command = (
        f"seqkit grep -p {patterns}  {args.reference_db_fasta} -o  {args.output_prefix}_working_reference.fa"
    )

    completed_process = subprocess.run(
        terminal_command,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        shell=True,
        check=True,
    )

    if completed_process.returncode != 0:
        print(f"\nERROR: seqkit grep terminated with errors.\nseqkit grep error code: {completed_process.returncode}\n")
        exit(1)

    print()


if __name__ == "__main__":
    main()
