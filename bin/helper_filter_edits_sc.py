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
--max-frac         Max per-cell editing fraction at a site (default: 0.10)
--site-max-frac    Max per-site editing fraction across all cells, Filter 6
                   (default: 0.05; only applied with --filter-site-max-frac)
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

# Previous implementation -- grouped by barcode as well, so a mixed site only
# lost the rows of the cells that showed the mixture and survived through the
# rest of the population. Sequencing-error hotspots in deeply covered genes are
# exactly that.
#
# def filter_multi_conversion(df):
#     """Remove rows where a barcode has >1 strand_conversion at the same site."""
#     mask = (
#         df.groupby(["barcode", "contig", "position"])["strand_conversion"]
#         .transform("nunique") == 1
#     )
#     return df[mask].copy()


def filter_multi_conversion(df):
    """Drop sites showing more than one strand_conversion across all cells.

    Per site, not per (cell, site): Methods say "Edit sites (sites with more than
    one edit type) [...] were filtered out". A position with several conversion
    types across the population is a sequencing-error hotspot rather than an edit
    site.
    """
    mask = (
        df.groupby(["contig", "position"])["strand_conversion"]
        .transform("nunique") == 1
    )
    return df[mask].copy()


def filter_dbsnp(df, dbsnp_bed_path):
    """Remove edit sites whose position overlaps a dbSNP entry."""
    sites = df[["contig", "position"]].drop_duplicates().copy()
    sites["start"] = sites["position"] - 1          # 1-based → 0-based BED start
    sites["end"]   = sites["position"]
    sites["name"]  = sites["contig"] + "_" + sites["position"].astype(str)

    # bedtools `-sorted` streams both files in lockstep instead of building an
    # in-memory interval tree of the whole dbSNP (>100 GB for hg38). It requires
    # both sides in the same order, and PREPARE_DBSNP sorts the dbSNP with
    # `sort -k1,1 -k2,2n` — lexicographic contig, numeric start — so match that here.
    sites = sites.sort_values(["contig", "start"], kind="mergesort")

    edits_bt = pybedtools.BedTool.from_dataframe(sites[["contig", "start", "end", "name"]])
    dbsnp_bt = pybedtools.BedTool(dbsnp_bed_path)

    non_overlapping = edits_bt.intersect(dbsnp_bt, v=True, sorted=True)
    if non_overlapping.count() == 0:
        # Never a legitimate result on real data: every candidate edit being a known
        # SNP means the intersect failed (historically an OOM kill) rather than
        # genuinely matching every site. Without this guard the empty result reaches
        # .to_dataframe() below and surfaces as an unrelated-looking KeyError.
        raise RuntimeError(
            f"dbSNP filter removed all {len(sites):,} candidate sites, which is not a "
            f"plausible result. The bedtools intersect against {dbsnp_bed_path} most "
            f"likely failed (out of memory, or a truncated/malformed dbSNP BED) rather "
            f"than genuinely matching every site. Refusing to emit an empty edit set."
        )
    keep_keys = set(non_overlapping.to_dataframe()["name"])

    site_key = df["contig"] + "_" + df["position"].astype(str)
    return df[site_key.isin(keep_keys)].copy()


def filter_min_total_edits(df, min_count=3):
    """Drop sites whose total edited-read count across all cells is below min_count."""
    site_totals = df.groupby(["contig", "position"])["count"].transform("sum")
    return df[site_totals >= min_count].copy()


def filter_max_editing_fraction(df, max_frac=0.10):
    """Remove rows where a cell's edit_fraction at a site exceeds max_frac.

    Per (cell, site), not per site. Doubles as a germline-variant filter: a SNP
    sits at edit_fraction ~1.0 in the cells carrying it, so this drops it even
    when dbSNP does not list it. Costs sensitivity in return -- a row needs
    coverage >= 1/max_frac in that one cell to survive with a single edit.
    """
    return df[df["edit_fraction"] <= max_frac].copy()


def filter_max_site_editing_fraction(df, max_frac=0.05):
    """Drop sites whose edited-read fraction across all cells exceeds max_frac.

    Off by default. MARINE emits rows only for cells that HAVE an edit at a
    site, so summing `coverage` here gives coverage in edited cells only, not
    total coverage at the site -- the resulting fraction is biased low for genes
    with many deeply covered edited cells. Measured on rep1 chr19 it retains 7
    genes (Malat1 dominating) and halves C>T purity, so it is exposed as an
    opt-in experiment rather than part of the default chain. Computing this
    correctly needs MARINE's full coverage matrix, including unedited cells.
    """
    grouped = df.groupby(["contig", "position"])
    site_frac = (
        grouped["count"].transform("sum") / grouped["coverage"].transform("sum")
    )
    return df[site_frac <= max_frac].copy()


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
                        help="Max per-cell editing fraction at a site (Filter 4)")
    parser.add_argument("--site-max-frac", type=float, default=0.05,
                        help="Max per-site editing fraction across all cells (Filter 6)")
    parser.add_argument("--output-dir", "-o", default=None,
                        help="Directory to save all outputs (TSVs + images)")
    # ── Per-filter on/off controls (set flag to skip that filter entirely) ────
    parser.add_argument("--no-filter-multi-conversion", action="store_true", default=False,
                        help="Skip Filter 1 — sites showing multiple conversion types")
    parser.add_argument("--no-filter-dbsnp",            action="store_true", default=False,
                        help="Skip Filter 2 — dbSNP overlap removal")
    parser.add_argument("--no-filter-min-count",        action="store_true", default=False,
                        help="Skip Filter 3 — minimum total edited reads per site")
    parser.add_argument("--no-filter-max-frac",         action="store_true", default=False,
                        help="Skip Filter 4 — max editing-fraction threshold")
    parser.add_argument("--no-filter-unannotated",      action="store_true", default=False,
                        help="Skip Filter 5 — unannotated site removal")
    # Filter 6 is opt-in rather than opt-out: see filter_max_site_editing_fraction.
    parser.add_argument("--filter-site-max-frac",       action="store_true", default=False,
                        help="Enable Filter 6 — per-site editing-fraction threshold")
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
        print("\nFilter 1: Remove sites with multiple conversion types across all cells")
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
        print(f"\nFilter 4: Remove cell-site entries with edit fraction > {args.max_frac}")
        df_04 = filter_max_editing_fraction(df_03, max_frac=args.max_frac)
        _report(df_03, df_04)
    steps.append((f"After F4\n(cell frac > {args.max_frac})", df_04))

    # ---- Filter 5 -----------------------------------------------------------
    if args.no_filter_unannotated:
        print("\nFilter 5: SKIPPED (--no-filter-unannotated)")
        df_05 = df_04
    else:
        print("\nFilter 5: Remove unannotated sites (feature_type == -1)")
        df_05 = filter_unannotated(df_04)
        _report(df_04, df_05)
    steps.append(("After F5\n(unannotated)", df_05))

    # ---- Filter 6 (opt-in) --------------------------------------------------
    if not args.filter_site_max_frac:
        print("\nFilter 6: SKIPPED (enable with --filter-site-max-frac)")
        df_06 = df_05
    else:
        print(f"\nFilter 6: Remove sites with per-site edit fraction > {args.site_max_frac}")
        df_06 = filter_max_site_editing_fraction(df_05, max_frac=args.site_max_frac)
        _report(df_05, df_06)
    # Appended unconditionally, like every other filter, so the summary table and
    # the pie/histogram grids always show all six steps -- a skipped filter shows
    # as a no-op row rather than vanishing from the figures.
    steps.append((f"After F6\n(site frac > {args.site_max_frac})", df_06))

    df_final = df_06

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
