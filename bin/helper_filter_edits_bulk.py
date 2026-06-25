#!/usr/bin/env python3
"""
filter_edits_bulk.py  —  Filter MARINE RNA editing sites for bulk RNA-seq.

Author: Aravind Sundaravadivelu
Date: 2026-06-08

Description
-----------
Filters MARINE output for bulk RNA-seq data. Optionally annotates unannotated
sites (strandedness=0 only) before applying filters. Output is passed downstream
to EPR normalization.

Steps
-----
1. Load MARINE TSV and compute edit_fraction = count / coverage
2. Annotate with BED file if input is unannotated and strandedness == 0
   (strandedness 1/2 requires a pre-annotated MARINE output)
3. Filter 1 — Remove multi-allelic sites (>1 alt base at same contig/position/ref)
4. Filter 2 — Remove sites overlapping dbSNP
5. Filter 3 — Remove sites where edit_fraction > max_frac
6. Filter 4 — Remove unannotated sites (feature_type == -1)

Required inputs
---------------
--marine-results   Path to MARINE TSV (annotated or unannotated)
--strandedness     0 = unstranded, 1 = forward, 2 = reverse
--dbsnp-bed        Path to dbSNP BED file (3-column BED)

Conditional inputs
------------------
--annotation-bed   BED6 gene annotation — required only when --strandedness 0

Optional inputs
---------------
--max-frac         Max editing fraction cutoff (default: 0.10)
--output-dir / -o  Directory to write TSV outputs and plots

Example
-------
python filter_edits_bulk.py \\
    --marine-results results/final_filtered_site_info_annotated.tsv \\
    --strandedness 2 \\
    --dbsnp-bed reference/hg38_dbsnp_combined.bed3 \\
    --max-frac 0.10 \\
    --output-dir filtered_output/
"""

import argparse
import math
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt
import pybedtools


# ---------------------------------------------------------------------------
# Annotation
# ---------------------------------------------------------------------------

def annotate_edit_sites(df, annotation_bed):
    """Annotate edit sites with gene info — strand-unaware (strandedness=0 only)."""
    print(f"Annotating sites with {annotation_bed} ...")

    sites = df[["contig", "position"]].drop_duplicates().copy()
    sites["start"] = sites["position"] - 1   # 1-based → 0-based BED start
    sites["end"]   = sites["position"]
    sites["name"]  = sites["contig"] + "_" + sites["position"].astype(str)

    edits_bt = pybedtools.BedTool.from_dataframe(sites[["contig", "start", "end", "name"]])
    annot_bt = pybedtools.BedTool(annotation_bed)

    intersect = edits_bt.intersect(annot_bt, wb=True, loj=True).to_dataframe()
    intersect.columns = [
        "contig", "start", "end", "name",
        "feature_chrom", "feature_start", "feature_end",
        "feature_name", "feature_type", "feature_strand",
    ]
    intersect = intersect[["name", "feature_name", "feature_type", "feature_strand"]]

    df = df.copy()
    df["_site_key"] = df["contig"] + "_" + df["position"].astype(str)
    intersect = intersect.rename(columns={"name": "_site_key"})

    df = df.merge(intersect, on="_site_key").drop(columns=["_site_key"])
    print(f"  Sites after annotation: {len(df):,}")
    return df


# ---------------------------------------------------------------------------
# Filter functions
# ---------------------------------------------------------------------------

def filter_alternate_edits(df):
    """Remove sites where more than one alt base appears at the same contig/position/ref."""
    mask = df.groupby(["contig", "position", "ref"])["alt"].transform("nunique") == 1
    return df[mask].copy()


def filter_dbsnp(df, dbsnp_bed_path):
    """Remove edit sites whose position overlaps a dbSNP entry."""
    sites = df[["contig", "position"]].drop_duplicates().copy()
    sites["start"] = sites["position"] - 1
    sites["end"]   = sites["position"]
    sites["name"]  = sites["contig"] + "_" + sites["position"].astype(str)

    edits_bt = pybedtools.BedTool.from_dataframe(sites[["contig", "start", "end", "name"]])
    dbsnp_bt = pybedtools.BedTool(dbsnp_bed_path)

    non_overlapping = edits_bt.intersect(dbsnp_bt, v=True)
    if non_overlapping.count() == 0:
        print("  Warning: all sites overlap dbSNP — returning empty DataFrame.")
        return df.iloc[0:0].copy()
    keep_keys = set(non_overlapping.to_dataframe()["name"])

    site_key = df["contig"] + "_" + df["position"].astype(str)
    return df[site_key.isin(keep_keys)].copy()


def filter_max_editing_fraction(df, max_frac=0.10):
    """Remove rows where edit_fraction > max_frac."""
    return df[df["edit_fraction"] <= max_frac].copy()


def filter_unannotated(df):
    """Remove sites with no gene annotation (feature_type == -1)."""
    return df[df["feature_type"].astype(str) != "-1"].copy()


# ---------------------------------------------------------------------------
# Logging helper
# ---------------------------------------------------------------------------

def _report(df_before, df_after):
    n_before  = len(df_before)
    n_after   = len(df_after)
    n_removed = n_before - n_after
    pct       = n_removed / n_before * 100 if n_before > 0 else 0
    print(f"  Before: {n_before:,}  |  After: {n_after:,}  |  Removed: {n_removed:,} ({pct:.1f}%)")


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _pie_grid(steps, title, output_path):
    """Render one pie chart per (label, DataFrame) pair as a grid image."""
    n     = len(steps)
    ncols = min(n, 3)
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
    axes_flat = axes.flatten() if n > 1 else [axes]

    for ax, (label, d) in zip(axes_flat, steps):
        counts = d["strand_conversion"].value_counts()
        ax.pie(counts, labels=counts.index, autopct="%1.1f%%", startangle=90)
        ax.set_title(f"{label}\n({len(d):,} edits)", fontsize=10)

    for ax in axes_flat[n:]:
        ax.set_visible(False)

    fig.suptitle(title, fontsize=13, y=1.01)
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, bbox_inches="tight", dpi=150)
        print(f"  Saved: {output_path}")
    else:
        plt.show()
    plt.close()


def _hist_grid(steps, output_path):
    """Render edit_fraction histograms for every step after 'Raw input'."""
    plot_steps = [(lbl, d) for lbl, d in steps if lbl != "Raw input"]
    n     = len(plot_steps)
    ncols = min(n, 2)
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
    axes_flat = axes.flatten() if n > 1 else [axes]

    for ax, (label, d) in zip(axes_flat, plot_steps):
        ax.hist(d["edit_fraction"], bins=50, range=(0, 0.2),
                color="skyblue", edgecolor="black")
        ax.set_title(label, fontsize=10)
        ax.set_xlabel("Edit Fraction")
        ax.set_ylabel("Frequency")

    for ax in axes_flat[n:]:
        ax.set_visible(False)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, bbox_inches="tight", dpi=150)
        print(f"  Saved: {output_path}")
    else:
        plt.show()
    plt.close()


def _plot_postfilter_distributions(df, output_dir):
    """Save post-filter edit_fraction and coverage distribution plots."""
    plots_dir = os.path.join(output_dir, "plots")

    plt.figure(figsize=(10, 6))
    plt.hist(df["edit_fraction"], bins=50, range=(0, 0.1), color="skyblue", edgecolor="black")
    plt.title("Distribution of Edit Fractions After All Filters")
    plt.xlabel("Edit Fraction")
    plt.ylabel("Frequency")
    path = os.path.join(plots_dir, "edit_fraction_distribution_after_all_filters.png")
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Saved: {path}")

    plt.figure(figsize=(10, 6))
    plt.hist(df["coverage"], bins=50, color="salmon", edgecolor="black")
    plt.title("Distribution of Coverage After All Filters")
    plt.xlabel("Coverage")
    plt.ylabel("Frequency")
    path = os.path.join(plots_dir, "coverage_distribution_after_all_filters.png")
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Filter MARINE RNA editing sites for bulk RNA-seq.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--marine-results", required=True,
                        help="Path to MARINE TSV (annotated or unannotated)")
    parser.add_argument("--strandedness", type=int, choices=[0, 1, 2], required=True,
                        help="Library strandedness: 0=unstranded, 1=forward, 2=reverse")
    parser.add_argument("--annotation-bed", default=None,
                        help="BED6 gene annotation file — required when --strandedness 0 and input is unannotated")
    parser.add_argument("--dbsnp-bed", required=True,
                        help="Path to dbSNP BED file")
    parser.add_argument("--max-frac", type=float, default=0.10,
                        help="Max edit fraction cutoff (Filter 3)")
    parser.add_argument("--output-dir", "-o", default=None,
                        help="Directory to save all outputs (TSVs + plots)")
    # ── Per-filter on/off controls (set flag to skip that filter entirely) ────
    parser.add_argument("--no-filter-multiallelic", action="store_true", default=False,
                        help="Skip Filter 1 — multi-allelic site removal")
    parser.add_argument("--no-filter-dbsnp",        action="store_true", default=False,
                        help="Skip Filter 2 — dbSNP overlap removal")
    parser.add_argument("--no-filter-max-frac",     action="store_true", default=False,
                        help="Skip Filter 3 — max editing-fraction threshold")
    parser.add_argument("--no-filter-unannotated",  action="store_true", default=False,
                        help="Skip Filter 4 — unannotated site removal")
    args = parser.parse_args()

    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
        os.makedirs(os.path.join(args.output_dir, "plots"), exist_ok=True)
        print(f"Output directory: {args.output_dir}")

    # ---- Load ---------------------------------------------------------------
    print(f"\nLoading {args.marine_results} ...")
    df = pd.read_csv(args.marine_results, sep="\t")
    zero_cov = (df["coverage"] == 0).sum()
    if zero_cov > 0:
        print(f"  Warning: {zero_cov:,} rows have coverage == 0 and will produce NaN edit_fraction. Consider removing them upstream.")
    df["edit_fraction"] = df["count"] / df["coverage"]
    print(f"Loaded {len(df):,} edit entries")

    # ---- Annotate if needed -------------------------------------------------
    required_annot_cols = {"feature_name", "feature_type", "feature_strand"}
    if not required_annot_cols.issubset(df.columns):
        if args.strandedness in (1, 2):
            print("Error: Input is unannotated but --strandedness is 1 or 2.")
            print("       For stranded libraries, provide a pre-annotated MARINE output file.")
            sys.exit(1)
        if not args.annotation_bed:
            print("Error: Input is unannotated and --annotation-bed was not provided.")
            sys.exit(1)
        print(f"\nAnnotating sites (strandedness=0, strand-unaware) ...")
        df = annotate_edit_sites(df, args.annotation_bed)
    else:
        print("Annotation columns present — skipping annotation step.")

    steps = [("Raw input", df)]

    # ---- Filter 1: Multi-allelic --------------------------------------------
    if args.no_filter_multiallelic:
        print("\nFilter 1: SKIPPED (--no-filter-multiallelic)")
        df_01 = df
    else:
        print("\nFilter 1: Remove multi-allelic sites")
        df_01 = filter_alternate_edits(df)
        _report(df, df_01)
    steps.append(("After F1\n(multi-allelic)", df_01))

    # ---- Filter 2: dbSNP ----------------------------------------------------
    if args.no_filter_dbsnp:
        print("\nFilter 2: SKIPPED (--no-filter-dbsnp)")
        df_02 = df_01
    else:
        print("\nFilter 2: Remove sites overlapping dbSNP")
        df_02 = filter_dbsnp(df_01, args.dbsnp_bed)
        _report(df_01, df_02)
    steps.append(("After F2\n(dbSNP overlap)", df_02))

    # ---- Filter 3: Max editing fraction -------------------------------------
    if args.no_filter_max_frac:
        print("\nFilter 3: SKIPPED (--no-filter-max-frac)")
        df_03 = df_02
    else:
        print(f"\nFilter 3: Remove sites with edit_fraction > {args.max_frac}")
        df_03 = filter_max_editing_fraction(df_02, max_frac=args.max_frac)
        _report(df_02, df_03)
    steps.append((f"After F3\n(frac > {args.max_frac})", df_03))

    # ---- Filter 4: Unannotated sites ----------------------------------------
    if args.no_filter_unannotated:
        print("\nFilter 4: SKIPPED (--no-filter-unannotated)")
        df_04 = df_03
    else:
        print("\nFilter 4: Remove unannotated sites (feature_type == -1)")
        df_04 = filter_unannotated(df_03)
        _report(df_03, df_04)
    steps.append(("After F4\n(unannotated)", df_04))

    df_final = df_04

    # ---- Summary table ------------------------------------------------------
    print("\n--- Summary ---")
    summary = pd.DataFrame([
        {
            "Step":             label.replace("\n", " "),
            "Edit entries":     len(d),
            "Removed from raw": len(df) - len(d),
            "% remaining":      f"{len(d) / len(df) * 100:.1f}%",
        }
        for label, d in steps
    ])
    print(summary.to_string(index=False))

    if args.output_dir:
        summary_path = os.path.join(args.output_dir, "filter_summary.tsv")
        summary.to_csv(summary_path, sep="\t", index=False)
        print(f"\nSaved summary  → {summary_path}")

        final_path = os.path.join(args.output_dir, "filtered_edits.tsv")
        df_final.to_csv(final_path, sep="\t", index=False)
        print(f"Saved filtered → {final_path}")

    # ---- Plots --------------------------------------------------------------
    print("\nGenerating plots ...")

    pie_path  = os.path.join(args.output_dir, "plots", "piecharts_all_steps.png")      if args.output_dir else None
    hist_path = os.path.join(args.output_dir, "plots", "edit_fraction_histograms.png") if args.output_dir else None

    _pie_grid(steps, "Strand-conversion distribution at each filtering step", pie_path)
    _hist_grid(steps, hist_path)

    if args.output_dir:
        _plot_postfilter_distributions(df_final, args.output_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
