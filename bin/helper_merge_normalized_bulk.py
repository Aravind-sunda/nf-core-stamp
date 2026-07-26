#!/usr/bin/env python3
"""
Merge per-sample normalized bulk edit tables into gene x sample matrices.

Usage:
    helper_merge_normalized_bulk.py --indir . --outdir .

Each input is one sample's *.EPR_EPKM_normalized.tsv from
helper_normalize_edits_bulk.py, in long format with one row per gene:

    sample_name | contig | feature_name | feature_strand | feature_type |
    Length | strand | strand_conversion | total_edits | total_coverage |
    featureCount_count | EPR | EPKM | EPKMR | EPM | EPMR

This writes one wide matrix per metric (genes as rows, samples as columns),
plus a gene annotation sidecar carrying the per-gene columns that the wide
matrices drop.

Missing cells: a gene absent from a sample's table is written as empty (NaN)
rather than 0. Absence means either "no edits called" (true 0) or "dropped
because the gene had zero featureCounts, making the metric undefined", and the
long table does not record which. Use --fill 0 if downstream analysis needs
zeros and that conflation is acceptable.
"""

import argparse
import glob
import os
import sys

import pandas as pd


# Written as one wide matrix each. total_edits and featureCount_count are the
# inputs to every ratio above them, so they are emitted too — comparing a
# suspicious metric against its numerator/denominator is the usual first check.
METRIC_COLS = ["EPR", "EPKM", "EPKMR", "EPM", "EPMR", "total_edits", "featureCount_count"]

# Per-gene, sample-invariant columns; collapsed into a single annotation table.
ANNOTATION_COLS = ["contig", "feature_strand", "feature_type", "Length"]

GENE_COL   = "feature_name"
SAMPLE_COL = "sample_name"


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--indir", required=True,
                        help="Directory containing per-sample *.EPR_EPKM_normalized.tsv files")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for the merged matrices")
    parser.add_argument("--suffix", default=".EPR_EPKM_normalized.tsv",
                        help="File suffix to glob for (default: .EPR_EPKM_normalized.tsv)")
    parser.add_argument("--fill", default=None,
                        help="Value for genes absent from a sample (default: leave empty/NaN)")
    return parser.parse_args()


def load_one(filepath):
    df = pd.read_csv(filepath, sep="\t")
    if df.empty:
        print(f"[WARN] {os.path.basename(filepath)} has no rows — skipping")
        return None

    missing = {GENE_COL, SAMPLE_COL}.difference(df.columns)
    if missing:
        sys.exit(f"[ERROR] {filepath} is missing required column(s): {sorted(missing)}")

    # One row per gene is what normalizing_marine_edits() produces (it groups by
    # contig/feature_name). Duplicates would make the pivot silently pick one value,
    # so fail instead of guessing which is correct.
    dups = df[GENE_COL].duplicated().sum()
    if dups:
        offenders = df.loc[df[GENE_COL].duplicated(keep=False), GENE_COL].unique()[:5]
        sys.exit(
            f"[ERROR] {filepath} has {dups} duplicated {GENE_COL} rows "
            f"(e.g. {offenders.tolist()}). Cannot pivot unambiguously."
        )
    return df


def main():
    args = parse_args()

    files = sorted(glob.glob(os.path.join(args.indir, f"*{args.suffix}")))
    if not files:
        sys.exit(f"[ERROR] No files matching *{args.suffix} found in {args.indir}")
    print(f"[LOG] Found {len(files)} sample file(s)")

    frames = [df for df in (load_one(f) for f in files) if df is not None]
    if not frames:
        sys.exit("[ERROR] Every input file was empty — nothing to merge")

    combined = pd.concat(frames, ignore_index=True)
    samples = sorted(combined[SAMPLE_COL].unique())
    genes   = combined[GENE_COL].nunique()
    print(f"[LOG] {len(samples)} sample(s), {genes:,} distinct genes")

    os.makedirs(args.outdir, exist_ok=True)

    for metric in METRIC_COLS:
        if metric not in combined.columns:
            print(f"[WARN] Column '{metric}' absent from inputs — skipping")
            continue

        wide = combined.pivot(index=GENE_COL, columns=SAMPLE_COL, values=metric)
        wide = wide.reindex(columns=samples)          # stable, sorted sample order
        if args.fill is not None:
            wide = wide.fillna(float(args.fill))
        wide.columns.name = None

        outfile = os.path.join(args.outdir, f"normalized_matrix_{metric}.tsv")
        wide.to_csv(outfile, sep="\t")
        filled = int(wide.notna().sum().sum())
        total  = wide.shape[0] * wide.shape[1]
        print(f"[DONE] {metric:19s} -> {os.path.basename(outfile)}  "
              f"({wide.shape[0]:,} genes x {wide.shape[1]} samples, "
              f"{100 * filled / total:.1f}% populated)")

    # Annotation is per-gene, so collapse to the first observation of each gene.
    annot_cols = [c for c in ANNOTATION_COLS if c in combined.columns]
    annotation = (
        combined[[GENE_COL] + annot_cols]
        .drop_duplicates(subset=[GENE_COL])
        .set_index(GENE_COL)
        .sort_index()
    )
    annot_file = os.path.join(args.outdir, "normalized_matrix_gene_annotation.tsv")
    annotation.to_csv(annot_file, sep="\t")
    print(f"[DONE] gene annotation      -> {os.path.basename(annot_file)} "
          f"({annotation.shape[0]:,} genes)")


if __name__ == "__main__":
    main()
