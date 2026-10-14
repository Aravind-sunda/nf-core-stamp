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
--output-dir / -o  If given, all TSV outputs are written to this directory

There is no editing-fraction filter here. A per-cell fraction is meaningless at
1-3 reads per cell per site (one real edit reads as 0.33-1.0), and a per-site
fraction needs true read depth, which MARINE's edit table cannot supply. That
filter (F5) runs downstream in helper_filter_site_frac_sc.py on the sites.bed
written here.

Example
-------
python filter_edits.py \\
    --marine-results results/final_filtered_site_info_annotated.tsv \\
    --dbsnp-bed      reference/mm10_dbsnp_combined.bed3 \\
    --min-count 3    \\
    --output-dir     filtered_output/
"""

import argparse
import os
import sys

import pandas as pd
import pybedtools

# Same directory (bin/); shared so F1–F4 and F5 count their steps identically.
from helper_filter_site_frac_sc import step_stats


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


def filter_unannotated(df):
    """Remove sites with no gene annotation (feature_type == -1)."""
    return df[df["feature_type"].astype(str) != "-1"].copy()


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
    parser.add_argument("--output-dir", "-o", default=None,
                        help="Directory to save all outputs")
    # ── Per-filter on/off controls (set flag to skip that filter entirely) ────
    parser.add_argument("--no-filter-multi-conversion", action="store_true", default=False,
                        help="Skip Filter 1 — sites showing multiple conversion types")
    parser.add_argument("--no-filter-dbsnp",            action="store_true", default=False,
                        help="Skip Filter 2 — dbSNP overlap removal")
    parser.add_argument("--no-filter-min-count",        action="store_true", default=False,
                        help="Skip Filter 3 — minimum total edited reads per site")
    parser.add_argument("--no-filter-unannotated",      action="store_true", default=False,
                        help="Skip Filter 4 — unannotated site removal")
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

        sites = df_final[["contig", "position"]].drop_duplicates()
        sites = sites.sort_values(["contig", "position"], kind="mergesort")
        bed_path = os.path.join(args.output_dir, "sites.bed")
        pd.DataFrame({"contig": sites["contig"],
                      "start":  sites["position"] - 1,     # 1-based → 0-based BED start
                      "end":    sites["position"]}).to_csv(
            bed_path, sep="\t", header=False, index=False)
        print(f"Saved {len(sites):,} surviving sites → {bed_path}")

        # Plots are drawn by helper_filter_site_frac_sc.py (FILTER_SITE_FRAC_SC),
        # which runs whether or not F5 is on, so F1–F5 share one set of figures.
        # These per-step counts are all it needs from here.
        stats_path = os.path.join(args.output_dir, "filter_step_stats_f1_f4.tsv")
        step_stats(steps).to_csv(stats_path, sep="\t", index=False)
        print(f"Saved per-step counts for plotting → {stats_path}")

    print("\nDone.")


def _report(df_before, df_after):
    n_before  = len(df_before)
    n_after   = len(df_after)
    n_removed = n_before - n_after
    pct       = n_removed / n_before * 100 if n_before > 0 else 0
    print(f"  Before: {n_before:,}  |  After: {n_after:,}  |  Removed: {n_removed:,} ({pct:.1f}%)")


if __name__ == "__main__":
    main()
