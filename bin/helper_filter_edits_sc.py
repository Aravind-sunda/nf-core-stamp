#!/usr/bin/env python3
"""
filter_edits.py  —  Filter MARINE RNA editing sites step by step.

Required inputs
---------------
--marine-results   Path to MARINE TSV (final_filtered_site_info_annotated.tsv)
--dbsnp-bed        Path to dbSNP BED file (3-column BED)

Optional inputs
---------------
--min-count        Minimum total edited reads per site across all cells (default: 3)
--max-frac         Max editing fraction for the primary fraction filter (default: 0.10)
--output-dir / -o  If given, all TSV outputs and images are written to this directory

Example
-------
python filter_edits.py \\
    --marine-results results/final_filtered_site_info_annotated.tsv \\
    --dbsnp-bed      reference/mm10_dbsnp_combined.bed3 \\
    --min-count 3    \\
    --max-frac 0.10  \\
    --output-dir     filtered_output/
"""

import argparse
import math
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt
import pybedtools


# ---------------------------------------------------------------------------
# Filter functions
# ---------------------------------------------------------------------------

def filter_multi_conversion(df):
    """Remove rows where a barcode has >1 strand_conversion at the same site."""
    mask = (
        df.groupby(["barcode", "contig", "position"])["strand_conversion"]
        .transform("nunique") == 1
    )
    return df[mask].copy()


def filter_dbsnp(df, dbsnp_bed_path):
    """Remove edit sites whose position overlaps a dbSNP entry."""
    sites = df[["contig", "position"]].drop_duplicates().copy()
    sites["start"] = sites["position"] - 1          # 1-based → 0-based BED start
    sites["end"]   = sites["position"]
    sites["name"]  = sites["contig"] + "_" + sites["position"].astype(str)

    edits_bt = pybedtools.BedTool.from_dataframe(sites[["contig", "start", "end", "name"]])
    dbsnp_bt = pybedtools.BedTool(dbsnp_bed_path)

    non_overlapping = edits_bt.intersect(dbsnp_bt, v=True)
    keep_keys = set(non_overlapping.to_dataframe()["name"])

    site_key = df["contig"] + "_" + df["position"].astype(str)
    return df[site_key.isin(keep_keys)].copy()


def filter_min_total_edits(df, min_count=3):
    """Drop sites whose total edited-read count across all cells is below min_count."""
    site_totals = df.groupby(["contig", "position"])["count"].transform("sum")
    return df[site_totals >= min_count].copy()


def filter_max_editing_fraction(df, max_frac=0.05):
    """Remove rows where edit_fraction > max_frac."""
    return df[df["edit_fraction"] <= max_frac].copy()


def filter_unannotated(df):
    """Remove sites with no gene annotation (feature_type == -1)."""
    return df[df["feature_type"].astype(str) != "-1"].copy()


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _pie_grid(steps, title, output_path):
    """Render one pie chart per (label, DataFrame) pair in *steps* as a grid image."""
    n = len(steps)
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
    """Render editing-fraction histograms for every step after 'Raw input'."""
    plot_steps = [(lbl, d) for lbl, d in steps if lbl != "Raw input"]
    n = len(plot_steps)
    ncols = min(n, 2)
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
    axes_flat = axes.flatten() if n > 1 else [axes]

    for ax, (label, d) in zip(axes_flat, plot_steps):
        ax.hist(d["edit_fraction"], bins=50, range=(0, 0.2),
                color="skyblue", edgecolor="black")
        ax.set_title(label, fontsize=10)
        ax.set_xlabel("Editing Fraction")
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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Filter MARINE RNA editing sites step by step.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--marine-results", required=True,
                        help="Path to MARINE annotated TSV file")
    parser.add_argument("--dbsnp-bed", required=True,
                        help="Path to dbSNP BED file")
    parser.add_argument("--min-count", type=int, default=3,
                        help="Min total edited reads per site (Filter 3)")
    parser.add_argument("--max-frac", type=float, default=0.10,
                        help="Max editing fraction (Filter 4)")
    parser.add_argument("--output-dir", "-o", default=None,
                        help="Directory to save all outputs (TSVs + images)")
    # ── Per-filter on/off controls (set flag to skip that filter entirely) ────
    parser.add_argument("--no-filter-multi-conversion", action="store_true", default=False,
                        help="Skip Filter 1 — multiple conversion types per barcode+site")
    parser.add_argument("--no-filter-dbsnp",            action="store_true", default=False,
                        help="Skip Filter 2 — dbSNP overlap removal")
    parser.add_argument("--no-filter-min-count",        action="store_true", default=False,
                        help="Skip Filter 3 — minimum total edited reads per site")
    parser.add_argument("--no-filter-max-frac",         action="store_true", default=False,
                        help="Skip Filter 4 — max editing-fraction threshold")
    parser.add_argument("--no-filter-unannotated",      action="store_true", default=False,
                        help="Skip Filter 5 — unannotated site removal")
    args = parser.parse_args()

    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
        print(f"Output directory: {args.output_dir}")

    # ---- Load ---------------------------------------------------------------
    print(f"\nLoading {args.marine_results} ...")
    df = pd.read_csv(args.marine_results, sep="\t")
    df["edit_fraction"] = df["count"] / df["coverage"]
    print(f"Loaded {len(df):,} edit entries across {df['barcode'].nunique():,} cells")

    steps = [("Raw input", df)]

    # ---- Filter 1 -----------------------------------------------------------
    if args.no_filter_multi_conversion:
        print("\nFilter 1: SKIPPED (--no-filter-multi-conversion)")
        df_01 = df
    else:
        print("\nFilter 1: Remove multiple conversion types per barcode + site")
        df_01 = filter_multi_conversion(df)
        _report(df, df_01)
    steps.append(("After F1\n(multi-conversion)", df_01))

    # ---- Filter 2 -----------------------------------------------------------
    if args.no_filter_dbsnp:
        print("\nFilter 2: SKIPPED (--no-filter-dbsnp)")
        df_02 = df_01
    else:
        print("\nFilter 2: Remove sites overlapping dbSNP")
        df_02 = filter_dbsnp(df_01, args.dbsnp_bed)
        _report(df_01, df_02)
    steps.append(("After F2\n(dbSNP overlap)", df_02))

    # ---- Filter 3 -----------------------------------------------------------
    if args.no_filter_min_count:
        print("\nFilter 3: SKIPPED (--no-filter-min-count)")
        df_03 = df_02
    else:
        print(f"\nFilter 3: Remove sites with < {args.min_count} total edits")
        df_03 = filter_min_total_edits(df_02, min_count=args.min_count)
        _report(df_02, df_03)
    steps.append((f"After F3\n(<{args.min_count} total edits)", df_03))

    # ---- Filter 4 -----------------------------------------------------------
    if args.no_filter_max_frac:
        print("\nFilter 4: SKIPPED (--no-filter-max-frac)")
        df_04 = df_03
    else:
        print(f"\nFilter 4: Remove editing fraction > {args.max_frac}")
        df_04 = filter_max_editing_fraction(df_03, max_frac=args.max_frac)
        _report(df_03, df_04)
    steps.append((f"After F4\n(frac > {args.max_frac})", df_04))

    # ---- Filter 5 -----------------------------------------------------------
    if args.no_filter_unannotated:
        print("\nFilter 5: SKIPPED (--no-filter-unannotated)")
        df_05 = df_04
    else:
        print("\nFilter 5: Remove unannotated sites (feature_type == -1)")
        df_05 = filter_unannotated(df_04)
        _report(df_04, df_05)
    steps.append(("After F5\n(unannotated)", df_05))

    df_final = df_05

    # ---- Summary table ------------------------------------------------------
    print("\n--- Summary ---")
    summary = pd.DataFrame([
        {
            "Step":            label.replace("\n", " "),
            "Edit entries":    len(d),
            "Removed from raw": len(df) - len(d),
            "% remaining":     f"{len(d) / len(df) * 100:.1f}%",
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
        print(f"Saved filtered data → {final_path}")

    # ---- Plots --------------------------------------------------------------
    print("\nGenerating plots ...")

    pie_path  = os.path.join(args.output_dir, "piecharts.png")  if args.output_dir else None
    hist_path = os.path.join(args.output_dir, "edit_fraction_histograms.png") if args.output_dir else None

    # add sample name to the title of the pie charts
    sample_name = os.path.basename(os.path.dirname(args.marine_results))
    _pie_grid(steps, f"Strand-conversion distribution at each filtering step ({sample_name})", pie_path)
    _hist_grid(steps, hist_path)

    print("\nDone.")


def _report(df_before, df_after):
    n_before  = len(df_before)
    n_after   = len(df_after)
    n_removed = n_before - n_after
    pct       = n_removed / n_before * 100 if n_before > 0 else 0
    print(f"  Before: {n_before:,}  |  After: {n_after:,}  |  Removed: {n_removed:,} ({pct:.1f}%)")


if __name__ == "__main__":
    main()
