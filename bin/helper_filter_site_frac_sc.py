#!/usr/bin/env python3
"""
helper_filter_site_frac_sc.py  —  SC Filter 5: per-site editing fraction against true depth.

Drops sites where more than --site-max-frac of the reads covering the site are
edited, summed across all cells in --barcodes. The denominator is read depth
from `samtools depth` over the BAM (SITE_DEPTH_SC), not MARINE's `coverage`
column: MARINE emits rows only for cells that HAVE an edit at a site, so summing
that column conditions the denominator on the outcome and biases the fraction low.

The edit rows are restricted to --barcodes too, so edited reads and depth are
counted over the same cells.

F3 (>= min-count edited reads per site) ran upstream across all CellRanger cells.
With a narrower barcode list, a site can pass F3 on reads from cells outside the
list and be left with 1-2 edited reads here, so the threshold is re-applied
within the list unless --no-recount-min-count is given.

Inputs
------
--filtered-edits        F1–F4 output of helper_filter_edits_sc.py (filtered_edits.tsv)
--site-depth            contig / position / depth TSV, header included (1-based positions)
--barcodes              One barcode per line, plain or .gz; the cells depth was measured over
--site-max-frac         Max edited-read fraction per site (default: 0.05)
--min-count             Min edited reads per site within --barcodes (default: 3)
--no-recount-min-count  Keep F3's all-cell count; do not re-apply --min-count here
--output-dir / -o       Output directory
"""

import argparse
import os

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Filter SC edit sites by editing fraction against true read depth.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--filtered-edits", required=True)
    parser.add_argument("--site-depth",     required=True)
    parser.add_argument("--barcodes",       required=True)
    parser.add_argument("--site-max-frac",  type=float, default=0.05)
    parser.add_argument("--min-count",      type=int,   default=3)
    parser.add_argument("--no-recount-min-count", action="store_true", default=False)
    parser.add_argument("--output-dir", "-o", default=".")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    edits = pd.read_csv(args.filtered_edits, sep="\t")
    barcodes = set(pd.read_csv(args.barcodes, header=None, usecols=[0], sep="\t")[0])
    depth = pd.read_csv(args.site_depth, sep="\t")

    n_raw = len(edits)
    sites_raw = edits.drop_duplicates(["contig", "position"])
    edits = edits[edits["barcode"].isin(barcodes)]
    print(f"Barcodes: {len(barcodes):,} in list; {edits['barcode'].nunique():,} with edits; "
          f"{n_raw - len(edits):,} of {n_raw:,} edit rows outside the list removed")

    sites = (edits.groupby(["contig", "position"])
                  .agg(strand_conversion=("strand_conversion", "first"),
                       feature_name=("feature_name", "first"),
                       n_cells=("barcode", "nunique"),
                       edited_reads=("count", "sum"))
                  .reset_index()
                  .merge(depth, on=["contig", "position"], how="left"))
    # Positions absent from the depth table have no qualifying reads at all.
    sites["depth"] = sites["depth"].fillna(0).astype(int)
    sites["site_frac"] = sites["edited_reads"] / sites["depth"]     # depth 0 → inf, dropped
    sites["pass_min_count"] = (args.no_recount_min_count
                               | (sites["edited_reads"] >= args.min_count))
    sites["pass_site_frac"] = sites["site_frac"] <= args.site_max_frac
    sites["pass"] = sites["pass_min_count"] & sites["pass_site_frac"]

    # Consistency check: MARINE's edited reads are a subset of the reads samtools
    # counts, so edited_reads > depth means the two read selections have diverged.
    n_over = int((sites["edited_reads"] > sites["depth"]).sum())
    if n_over:
        print(f"WARNING: {n_over:,} of {len(sites):,} sites have more edited reads than depth "
              f"(these fail F5). A few are expected; more than ~1% means the samtools read "
              f"selection no longer mirrors MARINE's.")

    keys = ["contig", "position"]
    counted = sites[sites["pass_min_count"]]
    keep = sites[sites["pass"]]
    out = edits.merge(keep[keys], on=keys)

    recount_label = ("Min-count not re-applied" if args.no_recount_min_count
                     else f"< {args.min_count} edited reads within barcodes")
    summary = pd.DataFrame([
        {"Step": "F1–F4 input",                          "Sites": len(sites_raw),
         "Edit entries": n_raw},
        {"Step": "Restricted to barcodes",               "Sites": len(sites),
         "Edit entries": len(edits)},
        {"Step": f"After {recount_label}",               "Sites": len(counted),
         "Edit entries": len(edits.merge(counted[keys], on=keys))},
        {"Step": f"After F5 (site frac > {args.site_max_frac})", "Sites": len(keep),
         "Edit entries": len(out)},
    ])
    ct = lambda d: 100 * (d["strand_conversion"] == "C>T").mean() if len(d) else 0.0
    summary["C>T % of sites"] = [ct(sites_raw), ct(sites), ct(counted), ct(keep)]
    print(summary.to_string(index=False))

    sites.to_csv(os.path.join(args.output_dir, "site_frac.tsv"), sep="\t", index=False)
    summary.to_csv(os.path.join(args.output_dir, "site_frac_summary.tsv"), sep="\t", index=False)
    out.to_csv(os.path.join(args.output_dir, "filtered_edits_site_frac.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
