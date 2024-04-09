#!/usr/bin/env python

import argparse
import os
from collections import defaultdict
import csv
import json
import sys

"""
python ../../../../bin/consensus_summary.py --c-stats-file consensus/230721_S_I_314_87-bwa_bcftools.consensus_stats.txt --c-coverage-file mapping/230721_S_I_314_87-bwa.rmdup.coverage.txt --c-typing-file typing/230721_S_I_314_87-bwa_bcftools.typing.tsv  --c-nextclade-files ../230721_S_I_314_87-bwa_bcftools-segment_4-ref_accession_CY054670/nextclade/230721_S_I_314_87-bwa_bcftools-segment_4-ref_accession_CY054670.tsv,../230721_S_I_314_87-bwa_bcftools-segment_4-ref_accession_MZ198323/nextclade/230721_S_I_314_87-bwa_bcftools-segment_4-ref_accession_MZ198323.tsv
"""

def parse_consensus_stats(stats_file):
    """
    parse "SEQKIT_FX2TAB --only-id --name --length -C ATCG -C RYSWKMBDHV -C N -H" output
    #id     length  ATCG    RYSWKMBDHV      N
    230721_S_I_314_87-bwa_bcftools-segment_1-ref_accession_MK381163 2280    2238    42      0
    230721_S_I_314_87-bwa_bcftools-segment_1-ref_accession_MN560975 2280    2255    25      0
    """
    stats_dict = {}

    if stats_file is None:
        return stats_dict

    if os.path.exists(stats_file) == False or os.path.isfile(stats_file) == False:
        print(f"\nERROR: coverage file {stats_file} does not exist.\n", file=sys.stderr)
        exit(1)
    if os.path.getsize(stats_file) == 0:
        print(f"\nWarning: coverage file {stats_file} is empty.\n", file=sys.stderr)
        # exit(1)
        return stats_dict

    with open(stats_file, "r", encoding="utf8") as fx2tab_file:
        reader = csv.DictReader(fx2tab_file, delimiter="\t")
        for row in reader:
            seqid = row.pop("#id")
            stats_dict[seqid] = {}
            stats_dict[seqid]['length'] = row["length"]
            stats_dict[seqid]['pct_Ns'] = round(100 * int(row["N"]) / int(row["length"]), 2)
            stats_dict[seqid]['pct_completeness'] = 100 - stats_dict[seqid]['pct_Ns']

    sorted_stats_dict = {i: stats_dict[i] for i in sorted(stats_dict.keys())}        

    #print(sorted_stats_dict)
    for key, value_dict in sorted_stats_dict.items():
        print("key*********************" + key, file=sys.stderr)
        
    return sorted_stats_dict


def parse_consensus_coverage_file(coverage_file):
    
    """
    parse "samtools coverage" output
    #rname  startpos        endpos  numreads        covbases        coverage        meandepth       meanbaseq       meanmapq
    CY054670        1       1701    3358    1701    100     275.253 35.5    60
    KF805639        1       973     15684   973     100     2230.12 35.5    56.4

    """
    if os.path.exists(coverage_file) == False or os.path.isfile(coverage_file) == False:
        print(f"\nERROR: coverage file {coverage_file} does not exist.\n", file=sys.stderr)
        exit(1)
    if os.path.getsize(coverage_file) == 0:
        print(f"\nERROR: coverage file {coverage_file} is empty.\n", file=sys.stderr)
        exit(1)

    cov_dict = {}

    with open(coverage_file, "r", encoding="utf8") as csvfile:
            reader = csv.DictReader(csvfile, delimiter="\t")
            cols_to_keep = ["numreads", "covbases", "coverage", "meandepth"]
            for row in reader:
                #print(type(row))
                id = row.pop("#rname")
                cov_dict[id] = {}

                for col in cols_to_keep:
                    cov_dict[id][col] = row[col]

    sorted_cov_dict = {i: cov_dict[i] for i in sorted(cov_dict.keys())}
    #print(sorted_cov_dict)
    return sorted_cov_dict

def parse_typing_file(blast_tabular_output):

    """
    parse blastn output: '6 std qlen slen qcovs'
    qseqid  sseqid  pident  length  mismatch        gapopen qstart  qend    sstart  send    evalue  bitscore        qlen    slen    qcovs
    230721_S_I_314_87-bwa_bcftools-segment_4-ref_accession_CY054670 FJ966974~~4~~H1 94.709  1701    90      0       1       1701    1       1701    0.0     2643    1701    1701    100
    230721_S_I_314_87-bwa_bcftools-segment_4-ref_accession_MZ198323 CY163680~~4~~H3 93.831  1702    103     2       1       1701    1       1701    0.0     2560    1701    1701    100
    """

    typing_dict = {}

    if blast_tabular_output is None:
        return typing_dict

    if os.path.exists(blast_tabular_output) == False or os.path.isfile(blast_tabular_output) == False:
        print(f"\nERROR: coverage file {blast_tabular_output} does not exist.\n", file=sys.stderr)
        exit(1)
    if os.path.getsize(blast_tabular_output) == 0:
        print(f"\nWarning: coverage file {blast_tabular_output} is empty.\n", file=sys.stderr)
        # exit(1)
        return typing_dict

    # print(blast_tabular_output)
    with open(blast_tabular_output, "r", encoding="utf8") as csvfile:
            reader = csv.DictReader(csvfile, delimiter="\t")
            cols_to_keep = ["sseqid", "pident", "mismatch", "gapopen", "evalue", "qcovs"]
            for row in reader:
                #print(type(row))
                id = row.pop("qseqid")
                typing_dict[id] = {}
                for col in cols_to_keep:
                    typing_dict[id][col] = row[col]

    sorted_typing_dict = {i: typing_dict[i] for i in sorted(typing_dict.keys())}
    #print(sorted_typing_dict)
    return sorted_typing_dict


def parse_nextclade_files(nextclade_tsv_outputs, sep):
    
    """
    parse nextclade tsv output file
    
    """
    nextclade_dict = {}

    if nextclade_tsv_outputs is None:
        return nextclade_dict

    #print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    #print(nextclade_tsv_outputs)
    for nextclade_tsv_output in nextclade_tsv_outputs.split(sep=","):
        if os.path.exists(nextclade_tsv_output) == False or os.path.isfile(nextclade_tsv_output) == False:
            print(f"\nERROR: nextclade tsv file {nextclade_tsv_output} does not exist.\n", file=sys.stderr)
            exit(1)
        if os.path.getsize(nextclade_tsv_output) == 0:
            print(f"\nWarning: nextclade tsv file {nextclade_tsv_output} is empty.\n", file=sys.stderr)
            # exit(1)
            continue

        # print(nextclade_tsv_output)

        with open(nextclade_tsv_output, "r", encoding="utf8") as nextclade_tsv_file:
            reader = csv.DictReader(nextclade_tsv_file, delimiter=sep)
            for row in reader:
                #print(type(row))
                id = row.pop("seqName")
                nextclade_dict[id] = {}
                nextclade_dict[id]["clade"] = row["clade"]

    sorted_nextclade_dict = {i: nextclade_dict[i] for i in sorted(nextclade_dict.keys())}
    #print(sorted_nextclade_dict)
    return sorted_nextclade_dict



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
    description = "Consensus summary report"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-s",
        "--c-stats-file",
        default=None,
        help=f"Path to seqkt fx2table produced consensus stats file\n",
    )
    parser.add_argument(
        "-c",
        "--c-coverage-file",
        default=None,
        help=f"Path to samtool coverage produced tabular file\n",
    )
    parser.add_argument(
        "-t",
        "--c-typing-file",
        default=None,
        help=f"Path to blastn outfmt 6  produced output\n",
    )
    
    parser.add_argument(
        "-n",
        "--c-nextclade-files",
        default=None,
        help=f"Path to nextclade produced tsv files\n",
    )

    
    args = parser.parse_args()
    
    c_stats_dict = parse_consensus_stats(args.c_stats_file)
    c_cov_dict = parse_consensus_coverage_file(args.c_coverage_file)

    for cid in c_stats_dict.keys():
        refid = cid.split("-")[-1].split("_")[-1]
        if refid in c_cov_dict:
            c_stats_dict[cid].update(c_cov_dict[refid])
    


    c_typing_dict = parse_typing_file(args.c_typing_file)
    c_nextclade_dict = parse_nextclade_files(args.c_nextclade_files, "\t")

    total_summary = dict_merge(c_stats_dict, c_typing_dict, c_nextclade_dict)
    jsonString = json.dumps(total_summary, indent=4)
    print(jsonString, file=sys.stderr)
    """ json_summary_file = open(f"{args.prefix}.json", "w")
    json_summary_file.write(jsonString)
    json_summary_file.close() """
    
    field_names = ["clade", "length", "pct_Ns", "pct_completeness", "numreads", "covbases", "coverage", "meandepth", "sseqid", "pident", "mismatch", "gapopen", "evalue", "qcovs"]
    print("cid," + ",".join(field_names))

    for cid in total_summary.keys():
        values = [cid]
        for field in field_names:
            if field in total_summary[cid]:
                if isinstance(total_summary[cid][field], int) or  isinstance(total_summary[cid][field], float):
                    values.append(str(total_summary[cid][field]))
                else:
                    values.append(total_summary[cid][field])
            else:
                values.append('N/A')
        print(",".join(values))

if __name__ == "__main__":
    main()
