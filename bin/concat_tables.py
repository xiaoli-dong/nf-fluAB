#!/usr/bin/env python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description="Concatenate delimited tables while preserving all columns."
)

parser.add_argument(
    "-o", "--output",
    required=True,
    help="Output file"
)

parser.add_argument(
    "--input-sep",
    default="\\t",
    help=r"Input field separator (default: '\t'). Examples: '\t', ','"
)

parser.add_argument(
    "--output-sep",
    default="\\t",
    help=r"Output field separator (default: '\t'). Examples: '\t', ','"
)

parser.add_argument(
    "--na",
    default="NA",
    help="Value used to fill missing columns (default: NA)"
)

parser.add_argument(
    "input_files",
    nargs="+",
    help="Input table files"
)

args = parser.parse_args()

# Convert escaped separators (e.g. '\t') into actual characters
input_sep = args.input_sep.encode().decode("unicode_escape")
output_sep = args.output_sep.encode().decode("unicode_escape")

dfs = []
master_cols = []

for f in args.input_files:
    print(f"Reading {f}")

    df = pd.read_csv(
        f,
        sep=input_sep,
        dtype=str,
        keep_default_na=False,
        comment="#"
    )

    dfs.append(df)

    # Use the file with the most columns as the master
    if len(df.columns) > len(master_cols):
        master_cols = list(df.columns)

# Append any columns not already present in the master
for df in dfs:
    for col in df.columns:
        if col not in master_cols:
            master_cols.append(col)

# Concatenate using the union of all columns
merged = pd.concat(
    dfs,
    ignore_index=True,
    sort=False
)

# Reorder columns to match the master column order
merged = merged.reindex(columns=master_cols)

# Fill missing values
merged = merged.fillna(args.na)

# Write output
merged.to_csv(
    args.output,
    sep=output_sep,
    index=False
)

print(f"\nMerged {len(args.input_files)} files")
print(f"Rows    : {len(merged)}")
print(f"Columns : {len(master_cols)}")
print(f"Output  : {args.output}")
