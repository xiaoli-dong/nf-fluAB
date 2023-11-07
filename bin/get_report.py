#!/usr/bin/env python

import itertools
import argparse
import sys
import os
from collections import defaultdict
import csv
import re
import json


def parse_seqkit_stats_file(path_to_file, paired_end):
    if path_to_file is None:
        seqstats = {"total_reads": "na", "total_bases": "na", "q20_rate": "na", "q30_rate": "na"}
        if paired_end:
            seqstats["read1_avg_len"] = "na"
            seqstats["read2_avg_len"] = "na"
        else:
            seqstats["read_avg_len"] = "na"

        return seqstats

    if os.path.exists(path_to_file) == False or os.path.isfile(path_to_file) == False:
        print(f"\nERROR: seqkit stats file {path_to_file} does not exist.\n")
        exit(1)
    if os.path.getsize(path_to_file) == 0:
        print(f"\nERROR: seqkit stats file {path_to_file} is empty.\n")
        exit(1)

    seqstats = {"total_reads": 0, "total_bases": 0, "q20_rate": 0, "q30_rate": 0}
    r1_pattern = re.compile(r".*_R*1[_.].*")
    r2_pattern = re.compile(r".*_R*2[_.].*")

    with open(path_to_file, "r", encoding="utf8") as stats_file:
        tsv_reader = csv.DictReader(stats_file, delimiter="\t")
        for stats in tsv_reader:
            file = stats["file"]
            seqstats["total_reads"] += int(stats["num_seqs"])
            seqstats["total_bases"] += int(stats["sum_len"])
            seqstats["q20_rate"] += float(stats["Q20(%)"])
            seqstats["q30_rate"] += float(stats["Q30(%)"])

            if paired_end:
                if r1_pattern.match(file):
                    seqstats["read1_avg_len"] = stats["avg_len"]

                elif r2_pattern.match(file):
                    seqstats["read2_avg_len"] = stats["avg_len"]
            else:
                seqstats["read_avg_len"] = stats["avg_len"]

            # print(f"{file} num_seqs: {stats} sum_len {sum_len}")

        if paired_end:
            seqstats["q20_rate"] = round(float(seqstats["q20_rate"]) / 2, 2)
            seqstats["q30_rate"] = round(float(seqstats["q30_rate"]) / 2, 2)
    # print(seqstats)
    return seqstats


def parse_map_summary_file(path_to_file):
    refs = {}

    if path_to_file is None:
        return refs

    if os.path.exists(path_to_file) == False or os.path.isfile(path_to_file) == False:
        print(f"\nERROR: mapping summary file {path_to_file} does not exist.\n")
        exit(1)
    if os.path.getsize(path_to_file) == 0:
        print(f"\nWarning: mapping summary file {path_to_file} is empty.\n")
        # exit(1)
        return refs

    with open(path_to_file, "r", encoding="utf8") as mapping_summary_file:
        reader = csv.DictReader(mapping_summary_file, delimiter=",")
        for row in reader:
            id = row.pop("#id")
            # Human|7|M|H1N1|Kenya|A/Kenya/035/2018|A|na|na|na
            query_id = row["query-ID"]
            x = query_id.split("|")
            row["host"] = x[0]
            row["segment_number"] = x[1]
            row["segment_name"] = x[2]
            row["serotype"] = x[3]
            row["virus_name"] = x[5]
            refs[id] = {}
            refs[id]["mapping"] = row
            # refs[id] = row
    sorted_refs = {i: refs[i] for i in sorted(refs.keys())}
    print(sorted_refs)
    return sorted_refs


def parse_seqkit_fx2tab_file(path_to_file):
    refs = {}

    if path_to_file is None:
        return refs

    if os.path.exists(path_to_file) == False or os.path.isfile(path_to_file) == False:
        print(f"\nERROR: coverage file {path_to_file} does not exist.\n")
        exit(1)
    if os.path.getsize(path_to_file) == 0:
        print(f"\nWarning: coverage file {path_to_file} is empty.\n")
        # exit(1)
        return refs

    with open(path_to_file, "r", encoding="utf8") as fx2tab_file:
        tsv_reader = csv.DictReader(fx2tab_file, delimiter="\t")

        for tsv in tsv_reader:
            id = tsv.pop("#id")
            tsv["percent_complete"] = round(100 * int(tsv["ATCG"]) / int(tsv["length"]), 2)
            tsv["percentN"] = round(100 * int(tsv["N"]) / int(tsv["length"]), 2)
            tsv["gene_length"] = tsv.pop("length")
            refs[id] = {}
            refs[id]["consensus_stats"] = tsv
    sorted_refs = {i: refs[i] for i in sorted(refs.keys())}
    print(sorted_refs)
    return sorted_refs


# blastn output: '6 std qlen slen qcovs'
def parse_typing_file(path_to_file):
    refs = {}

    if path_to_file is None:
        return refs

    if os.path.exists(path_to_file) == False or os.path.isfile(path_to_file) == False:
        print(f"\nERROR: coverage file {path_to_file} does not exist.\n")
        exit(1)
    if os.path.getsize(path_to_file) == 0:
        print(f"\nWarning: coverage file {path_to_file} is empty.\n")
        # exit(1)
        return refs

    # print(path_to_file)

    with open(path_to_file, "r", encoding="utf8") as blastn_typing_file:
        header = [
            # "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "qlen",
            "slen",
            "qcovs",
        ]
        tsv_reader = csv.reader(blastn_typing_file, delimiter="\t")
        # print(tsv_reader)
        for tsv in tsv_reader:
            record_dict = {}
            id = tsv.pop(0)
            # print(type(tsv))
            # print(tsv)
            for i in range(len(tsv)):
                if header[i] == "sseqid":
                    # CY163680~~4~~H3
                    x = tsv[i].split("~~")
                    record_dict["typing_accession"] = x[0]
                    record_dict["type"] = x[2]
                else:
                    record_dict[header[i]] = tsv[i]
            # refs[id] = record_dict
            refs[id] = {}
            refs[id]["typing"] = record_dict

    sorted_refs = {i: refs[i] for i in sorted(refs.keys())}
    print(sorted_refs)
    return sorted_refs


def parse_nextclade_csv_files(path_to_files):
    refs = {}

    if path_to_files is None:
        return refs

    print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    print(path_to_files)
    for path_to_file in path_to_files.split(sep=","):
        if os.path.exists(path_to_file) == False or os.path.isfile(path_to_file) == False:
            print(f"\nERROR: nextclade csv file {path_to_file} does not exist.\n")
            exit(1)
        if os.path.getsize(path_to_file) == 0:
            print(f"\nWarning: nextclade csv file {path_to_file} is empty.\n")
            # exit(1)
            continue

        # print(path_to_file)

        with open(path_to_file, "r", encoding="utf8") as nextclade_csv_file:
            reader = csv.DictReader(nextclade_csv_file, delimiter="\t")
            for row in reader:
                print(type(row))
                id = row.pop("seqName")
                # Human|7|M|H1N1|Kenya|A/Kenya/035/2018|A|na|na|na
                # refs[id] = slicedict(row, "clade")
                refs[id] = {}
                refs[id]["nextclade"] = slicedict(row, "clade")

    sorted_refs = {i: refs[i] for i in sorted(refs.keys())}
    print(sorted_refs)
    return sorted_refs


def slicedict(d, s):
    return {k: v for k, v in d.items() if k.startswith(s)}


import collections.abc


def dict_merge(*args, add_keys=True):
    assert len(args) >= 2, "dict_merge requires at least two dicts to merge"
    rtn_dct = args[0].copy()
    merge_dicts = args[1:]
    for merge_dct in merge_dicts:
        if add_keys is False:
            merge_dct = {key: merge_dct[key] for key in set(rtn_dct).intersection(set(merge_dct))}
        for k, v in merge_dct.items():
            if not rtn_dct.get(k):
                rtn_dct[k] = v
            elif k in rtn_dct and type(v) != type(rtn_dct[k]):
                raise TypeError(
                    f"Overlapping keys exist with different types: original is {type(rtn_dct[k])}, new value is {type(v)}"
                )
            elif isinstance(rtn_dct[k], dict) and isinstance(merge_dct[k], collections.abc.Mapping):
                rtn_dct[k] = dict_merge(rtn_dct[k], merge_dct[k], add_keys=add_keys)
            elif isinstance(v, list):
                for list_value in v:
                    if list_value not in rtn_dct[k]:
                        rtn_dct[k].append(list_value)
            else:
                rtn_dct[k] = v
    return rtn_dct


def main():
    description = "Generate report"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-s",
        "--sample-name",
        required=True,
        help=f"Sample name\n",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        required=True,
        help=f"prefix of the output files\n",
    )
    parser.add_argument(
        "-r",
        "--raw-stats-file",
        default=None,
        help=f"Path to the raw read seqkit stats output file in tabular format\n",
    )
    parser.add_argument(
        "-q",
        "--qc-stats-file",
        default=None,
        help=f"Path to the qc read seqkit stats output file in tabular format\n",
    )

    parser.add_argument(
        "-f",
        "--seqkit-fx2tab-file",
        default=None,
        help=f"Path to seqkt fx2table produced stats file\n",
    )
    parser.add_argument(
        "-b",
        "--blastn-typing-file",
        default=None,
        help=f"Path to blastn produced typing file\n",
    )
    parser.add_argument(
        "-j",
        "--json-output",
        default=None,
        help=f"The output file name for the json data\n",
    )

    parser.add_argument(
        "-t",
        "--tsv-output",
        default=None,
        help=f"The output file name for tabular format\n",
    )
    parser.add_argument(
        "-m",
        "--mapping-summary-file",
        default=None,
        help=f"Path to mapping_summary_file file\n",
    )
    parser.add_argument(
        "-n",
        "--nextclade-csv-files",
        default=None,
        help=f"Path to nextclade csv files\n",
    )

    # parser.add_argument("file", action="store", nargs=1)

    args = parser.parse_args()
    raw_read_stats = parse_seqkit_stats_file(args.raw_stats_file, True)
    qc_read_stats = parse_seqkit_stats_file(args.qc_stats_file, True)
    summary = {"sname": args.sample_name}
    summary["raw_reads"] = raw_read_stats
    summary["qc_reads"] = qc_read_stats

    consensus_stats = parse_seqkit_fx2tab_file(args.seqkit_fx2tab_file)
    mapping_summary = parse_map_summary_file(args.mapping_summary_file)
    types = parse_typing_file(args.blastn_typing_file)
    print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxttttttttttttttt")
    nextclade = parse_nextclade_csv_files(args.nextclade_csv_files)
    print("yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy")
    total_summary = dict_merge(summary, consensus_stats, mapping_summary, types, nextclade)
    jsonString = json.dumps(total_summary, indent=4)
    print(jsonString)
    json_summary_file = open(f"{args.prefix}.json", "w")
    json_summary_file.write(jsonString)
    json_summary_file.close()

    segment_summary_file = open(f"{args.prefix}.segments.csv", "w")
    header = [
        "contig",
        "contig_length",
        "total_ATCG",
        #'total_nonATCG'
        "totalN",
        "percent_complete",
        "percentN",
        "mapped_reads",
        "covbases",
        "coverage",
        "avg_depth",
        # typing information
        "type",
        "typing_pident",
        "typing_qcovs",
        "typing_accession",
        #
        "clade",
        # "clade_database",
        # "influenza_db_version",
        "ref_accession",
        "ref_segment_number",
        "ref_segment_name",
        "ref_serotype",
        "ref_virus_name",
        # "reference_typing_data",
        "ref_identity",  # "mash_distance",
        "ref_shared-hashes",  # "minhashes",
    ]
    # print(",".join(header))
    segment_summary_file.write(",".join(header))
    segment_summary_file.write("\n")
    for key, value_dict in total_summary.items():
        if key in ["sname", "raw_reads", "qc_reads"]:
            continue
        values = []
        values.append(key)
        values.append(value_dict["consensus_stats"]["gene_length"])
        values.append(value_dict["consensus_stats"]["ATCG"])
        values.append(value_dict["consensus_stats"]["N"])
        values.append(value_dict["consensus_stats"]["percent_complete"])
        values.append(value_dict["consensus_stats"]["percentN"])
        values.append(value_dict["mapping"]["numreads"])
        values.append(value_dict["mapping"]["covbases"])
        values.append(value_dict["mapping"]["coverage"])
        values.append(value_dict["mapping"]["meandepth"])
        if "typing" in value_dict.keys():
            values.append(value_dict["typing"]["type"])
            values.append(value_dict["typing"]["pident"])
            values.append(value_dict["typing"]["qcovs"])
            values.append(value_dict["typing"]["typing_accession"])
        else:
            values.extend(["na", "na", "na", "na"])

        if "clade" in value_dict.keys():
            values.append(value_dict["clade"]["clade"])
        else:
            values.append("na")
        values.append(value_dict["mapping"]["accession"])
        values.append(value_dict["mapping"]["segment_number"])
        values.append(value_dict["mapping"]["segment_name"])
        values.append(value_dict["mapping"]["serotype"])
        values.append(value_dict["mapping"]["virus_name"])
        values.append(value_dict["mapping"]["identity"])
        values.append(value_dict["mapping"]["shared-hashes"])

        # when list contain both number and str, it raises typeError for join
        # convert list to str
        all_strings = list(map(str, values))
        # print(",".join(all_strings))
        segment_summary_file.write(",".join(all_strings))
        segment_summary_file.write("\n")
    segment_summary_file.close()

    segment_summary_oneline_file = open(f"{args.prefix}.segments_oneline.csv", "w")
    header = [
        "contig",
        "contig_length",
        "percent_complete",
        "percentN",
        "mapped_reads",
        "covbases",
        "coverage",
        "avg_depth",
        "type",
        "typing_accession",
        "clade",
    ]

    # 8 segments
    one_liner_header = header * 8
    # print(",".join(one_liner_header))
    # print()
    segment_summary_oneline_file.write(",".join(one_liner_header))
    segment_summary_oneline_file.write("\n")
    values = []
    for key, value_dict in total_summary.items():
        if key in ["sname", "raw_reads", "qc_reads"]:
            continue

        values.append(key)
        values.append(value_dict["consensus_stats"]["gene_length"])
        values.append(value_dict["consensus_stats"]["percent_complete"])
        values.append(value_dict["consensus_stats"]["percentN"])
        values.append(value_dict["mapping"]["numreads"])
        values.append(value_dict["mapping"]["covbases"])
        values.append(value_dict["mapping"]["coverage"])
        values.append(value_dict["mapping"]["meandepth"])
        if "typing" in value_dict.keys():
            values.append(value_dict["typing"]["type"])
            values.append(value_dict["typing"]["typing_accession"])
        else:
            values.extend(["na", "na", "na", "na"])

        if "clade" in value_dict.keys():
            values.append(value_dict["clade"]["clade"])
        else:
            values.append("na")

    # when list contain both number and str, it raises typeError for join
    # convert list to str
    all_strings = list(map(str, values))
    # print(",".join(all_strings))
    # print()
    segment_summary_oneline_file.write(",".join(all_strings))
    segment_summary_oneline_file.write("\n")
    segment_summary_oneline_file.close()


if __name__ == "__main__":
    main()
