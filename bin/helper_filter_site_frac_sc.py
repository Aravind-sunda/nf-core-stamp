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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _pie_grid(steps, title, output_path):
    """One pie of conversion types per step, counted by distinct site (the paper's C>T purity)."""
    fig, axes = plt.subplots(1, len(steps), figsize=(6 * len(steps), 5))
    for ax, (label, d) in zip(axes, steps):
        counts = d["strand_conversion"].value_counts()
        small = counts / counts.sum() < 0.02     # unreadable as separate slices
        if small.sum() > 1:
            counts = pd.concat([counts[~small], pd.Series({"other": counts[small].sum()})])
        ax.pie(counts, labels=counts.index, autopct="%1.1f%%", startangle=90)
        ax.set_title(f"{label}\n({len(d):,} unique edit sites)", fontsize=10)
    fig.suptitle(title, fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", dpi=150)
    plt.close()


def _hist_grid(steps, max_frac, output_path):
    """Per-site edited fraction (edited reads / depth), 0 to 1, with the F5 cutoff marked."""
    fig, axes = plt.subplots(1, len(steps), figsize=(6 * len(steps), 5))
    for ax, (label, d) in zip(axes, steps):
        finite = d["site_frac"][np.isfinite(d["site_frac"])]
        ax.hist(finite.clip(upper=1), bins=50, range=(0, 1), color="skyblue", edgecolor="black")
        ax.axvline(max_frac, color="crimson", linestyle="--", label=f"cutoff {max_frac}")
        ax.set_yscale("log")
        n_zero = len(d) - len(finite)
        ax.set_title(f"{label}\n({len(d):,} sites"
                     f"{f'; {n_zero:,} with depth 0 not shown' if n_zero else ''})", fontsize=10)
        ax.set_xlabel("Site editing fraction (edited reads / depth)")
        ax.set_ylabel("Sites (log scale)")
        ax.legend()
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", dpi=150)
    plt.close()


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
    parser.add_argument("--sample", default="", help="Sample name for plot titles")
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

    # Named site_frac_* so they do not overwrite helper_filter_edits_sc.py's plots,
    # which publish to the same directory.
    frac_steps = [("Restricted to barcodes", sites),
                  (f"After {recount_label}", counted),
                  (f"After F5 (site frac > {args.site_max_frac})", keep)]
    _pie_grid([("F1–F4 input", sites_raw)] + frac_steps,
              f"Conversion types by unique edit site at each F5 step{f' ({args.sample})' if args.sample else ''}",
              os.path.join(args.output_dir, "site_frac_piecharts.png"))
    _hist_grid(frac_steps, args.site_max_frac,
               os.path.join(args.output_dir, "site_frac_histograms.png"))


if __name__ == "__main__":
    main()
