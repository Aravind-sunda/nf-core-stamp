#!/usr/bin/env python3
"""
Merge per-sample featureCounts outputs into a single count matrix.

Usage:
    python 06_merge_counts.py --indir 05_featurecounts --outfile counts_matrix.tsv

Each per-sample file has columns:
    Geneid | Chr | Start | End | Strand | Length | <sample.bam>

The merge keeps the annotation columns from the first file and appends
the count column from every subsequent file. Gene rows are identical
across all files (same GTF), so no join is needed — pure column-bind.

The BAM filename in the header is replaced with the clean sample name.
"""

import argparse
import glob
import os
import sys

import pandas as pd


# Everything that is not a per-sample count column. gene_id is emitted by
# featureCounts --extraAttributes; it must be listed here or load_featurecounts()
# mistakes it for a second count column and aborts.
ANNOTATION_COLS = ["Chr", "Start", "End", "Strand", "Length", "gene_id"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--indir",   required=True, help="Directory containing per-sample featureCounts .txt files")
    parser.add_argument("--outfile", required=True, help="Output TSV path for the merged count matrix")
    parser.add_argument("--suffix",  default=".featurecounts.txt", help="File suffix to glob for (default: .featurecounts.txt)")
    # Deliberately not exposed by the pipeline: MERGE_COUNTS never passes this.
    # helper_normalize_edits_bulk.py reads Length from this same matrix to compute
    # EPKM/EPKMR, and its process_feature_counts() drops the annotation columns
    # unconditionally — so a matrix built with this flag breaks NORMALIZE_EDITS_BULK.
    # Kept for standalone use only.
    parser.add_argument("--drop-annotation", action="store_true",
                        help="Drop Chr/Start/End/Strand/Length columns from output. "
                             "Standalone use only — the pipeline never sets this, because "
                             "NORMALIZE_EDITS_BULK requires the Length column.")
    return parser.parse_args()


def sample_name(filepath, suffix):
    return os.path.basename(filepath).replace(suffix, "")


def load_featurecounts(filepath, suffix):
    df = pd.read_csv(filepath, sep="\t", comment="#", index_col="Geneid")
    # rename the BAM count column to the clean sample name
    bam_col = [c for c in df.columns if c not in ANNOTATION_COLS]
    if len(bam_col) != 1:
        sys.exit(f"[ERROR] Expected 1 count column in {filepath}, found: {bam_col}")
    df = df.rename(columns={bam_col[0]: sample_name(filepath, suffix)})
    return df


def main():
    args = parse_args()

    files = sorted(glob.glob(os.path.join(args.indir, f"*{args.suffix}")))
    if not files:
        sys.exit(f"[ERROR] No files matching *{args.suffix} found in {args.indir}")

    print(f"[LOG] Found {len(files)} sample file(s)")

    # load first file with annotation columns
    merged = load_featurecounts(files[0], args.suffix)
    print(f"[LOG] Loaded {sample_name(files[0], args.suffix)}")

    # append count column from each remaining file
    for f in files[1:]:
        name = sample_name(f, args.suffix)
        df = load_featurecounts(f, args.suffix)
        # sanity check: gene order must match
        if not merged.index.equals(df.index):
            sys.exit(f"[ERROR] Gene index mismatch between first file and {f} — check that all files used the same GTF")
        merged[name] = df[name]
        print(f"[LOG] Loaded {name}")

    # errors="ignore" so matrices produced before --extraAttributes was added
    # (no gene_id column) still work.
    n_annotation = sum(c in merged.columns for c in ANNOTATION_COLS)
    if args.drop_annotation:
        merged = merged.drop(columns=ANNOTATION_COLS, errors="ignore")
        n_annotation = 0

    os.makedirs(os.path.dirname(os.path.abspath(args.outfile)), exist_ok=True)
    merged.to_csv(args.outfile, sep="\t")
    print(f"[DONE] Count matrix written to {args.outfile}")
    print(f"[DONE] Shape: {merged.shape[0]} genes x {merged.shape[1] - n_annotation} samples")


if __name__ == "__main__":
    main()
