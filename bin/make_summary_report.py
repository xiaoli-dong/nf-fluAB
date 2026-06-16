#!/usr/bin/env python3

import argparse
import pandas as pd
import re
import os

def process_analysis(consensus_path, depth_path, nextclade_path, squirrel_path):
    """
    Processes mpox analysis. Returns contig-level DataFrame.
    """

    # 1. Consensus
    cons_df = pd.DataFrame()
    if os.path.exists(consensus_path) and os.path.getsize(consensus_path) > 0:
        cons = pd.read_csv(consensus_path, sep="\t")
        if not cons.empty:
            """ cons_df = cons.rename(columns={
                "#id": "contigid",
                "coverage": "consensus_coverage/10x",
                "completeness": "consensus_completeness/10x"
            })[["contigid", "consensus_coverage/10x", "consensus_completeness/10x"]].drop_duplicates() """
            cons_df = cons.rename(columns={
                "#id": "contigid",
                "coverage": "consensus_coverage",
                "completeness": "consensus_completeness"
            })[["contigid", "consensus_coverage", "consensus_completeness"]].drop_duplicates()

    # 2. Depth
    depth_df = pd.DataFrame()
    if os.path.exists(depth_path) and os.path.getsize(depth_path) > 0:
        depth = pd.read_csv(depth_path, sep="\t")
        if not depth.empty:
            #depth["contigid"] = depth["sample"].astype(str) + "|ref|" + depth["chrom"].astype(str)
            depth["contigid"] = depth["sample"].astype(str)
            depth_df = depth[["contigid", "mapped_reads", "mean_depth"]].drop_duplicates()

    # 3. Nextclade
    nextclade_df = pd.DataFrame()
    if os.path.exists(nextclade_path) and os.path.getsize(nextclade_path) > 0:
        nc = pd.read_csv(nextclade_path, sep="\t")
        if not nc.empty:
            nextclade_df = nc.rename(
                columns={"seqName": "contigid",
                        "clade": "nextclade-clade",
                        "lineage": "nextclade-lineage",
                        "qc.overallStatus": "nextclade-qc.overallStatus"
                        }
                )[["contigid", "nextclade-clade", "nextclade-lineage", "nextclade-qc.overallStatus"]].drop_duplicates()

    # 3.1 squirrels
    squirrel_df = pd.DataFrame()
    if os.path.exists(squirrel_path) and os.path.getsize(squirrel_path) > 0:
        sq = pd.read_csv(squirrel_path, sep=",")
        if not sq.empty:
            squirrel_df = sq.rename(
                columns={"taxon": "contigid",
                        "prediction": "squirrel-clade",
                        "score": "squirrel-score"
                        }
                )[["contigid", "squirrel-clade", "squirrel-score"]].drop_duplicates()

    squirrel_df["squirrel-score"] = squirrel_df["squirrel-score"].round(2)

    # 4. Merge
    if cons_df.empty and depth_df.empty and nextclade_df.empty:
        return pd.DataFrame()

    merged = pd.merge(cons_df, depth_df, on="contigid", how="outer")
    merged = pd.merge(merged, nextclade_df, on="contigid", how="outer")
    merged = pd.merge(merged, squirrel_df, on="contigid", how="outer")

    # Extract sample_id
    merged["sample_id"] = merged["contigid"].str.split("|").str[0]

    return merged


def extract_trailing_number(sample_id):
    if pd.isna(sample_id): return 999999
    match = re.search(r'(\d+)$', str(sample_id))
    return int(match.group(1)) if match else 999999


def main():
    parser = argparse.ArgumentParser(description="Consolidate Mpox Analysis with QC Stats.")
    parser.add_argument("--qc", required=True)
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--depth", required=True)
    parser.add_argument("--nextclade", required=True)
    parser.add_argument("--squirrel", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    # 1. QC
    qc = pd.read_csv(args.qc).rename(columns={"sample": "sample_id"})
    qc = qc[["sample_id", "input_num_seqs", "trimmed_num_seqs", "dehosted_num_seqs"]].drop_duplicates()

    # 2. Analysis
    analysis = process_analysis(args.consensus, args.depth, args.nextclade, args.squirrel)

    # 3. Merge QC + Analysis
    master = analysis.merge(qc, on="sample_id", how="outer")

    # Handle QC-only samples
    master["nextclade-clade"] = master["nextclade-clade"].fillna("No Call")

     # Handle QC-only samples
    master["squirrel-clade"] = master["squirrel-clade"].fillna("No Call")

    # 4. Sorting
    master["_num"] = master["sample_id"].apply(extract_trailing_number)
    master = master.sort_values(by="_num").drop(columns=["_num"])

    # Column order

    column_order = [
        "sample_id", "contigid", "nextclade-clade", "nextclade-lineage", "nextclade-qc.overallStatus","squirrel-clade", "squirrel-score",
        "input_num_seqs", "trimmed_num_seqs", "dehosted_num_seqs",
        "consensus_coverage", "consensus_completeness",
        "mapped_reads", "mean_depth"
    ]

    master = master[[c for c in column_order if c in master.columns]]

    # Split sample_id
    sample_split = master["sample_id"].str.split("-", n=1, expand=True)
    master.insert(0, "runid", sample_split[0])
    master.insert(1, "sample", sample_split[1])
    master = master.drop(columns=["sample_id"])

    # 5. Save
    master.to_csv(args.out, sep="\t", index=False, na_rep="NaN")

    print(f"✔ Final MPXV report written to: {args.out}")
    print(f"Total Rows: {len(master)}")


if __name__ == "__main__":
    main()
