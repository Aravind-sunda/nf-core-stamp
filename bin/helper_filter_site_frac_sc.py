#!/usr/bin/env python3
"""
helper_filter_site_frac_sc.py  —  SC Filter 5 (per-site editing fraction against true
depth) and the plots for every SC filter step.

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

Plots are drawn here for F1–F5 alike, from the per-step counts that
helper_filter_edits_sc.py writes (--step-stats) plus this script's own steps.
Without --site-depth, F5 is off: nothing is filtered and only the F1–F4 plots are drawn.

Inputs
------
--step-stats            filter_step_stats_f1_f4.tsv from helper_filter_edits_sc.py; with
                        F5 on, written back out with F5's steps as filter_step_stats.tsv
--filtered-edits        F1–F4 output of helper_filter_edits_sc.py (filtered_edits.tsv)
--site-depth            contig / position / depth TSV, header included (1-based positions)
--barcodes              One barcode per line, plain or .gz; the cells depth was measured over
--site-max-frac         Max edited-read fraction per site (default: 0.05)
--min-count             Min edited reads per site within --barcodes (default: 3)
--no-recount-min-count  Keep F3's all-cell count; do not re-apply --min-count here
--sample                Sample name for plot titles
--output-dir / -o       Output directory
"""

import argparse
import math
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# One pie file per measure.
PIE_MEASURES = {                  # measure: (panel unit, figure-title unit)
    "unique_sites": ("unique edit sites", "unique edit site"),
    "edit_rows":    ("edit rows",         "edit row (one per cell per site)"),
}
CELL_FRAC_BINS = np.linspace(0, 1, 51)

# Fixed colour per conversion type so a slice keeps its colour across steps.
_TYPES = ["C>T", "G>A", "A>G", "T>C", "T>A", "A>T", "C>A", "G>T", "A>C", "T>G", "C>G", "G>C"]
_DISTINCT = [c for i, c in enumerate(plt.get_cmap("tab10").colors) if i != 7]  # grey kept for "other"
PIE_COLORS = dict(zip(_TYPES, _DISTINCT + ["#393b79", "#8c6d31", "#637939"]))
PIE_COLORS["other"] = "#c7c7c7"


# ---------------------------------------------------------------------------
# Per-step counts (also imported by helper_filter_edits_sc.py)
# ---------------------------------------------------------------------------

def step_stats(steps):
    """Long table of what the plots need from each (label, edits DataFrame) step.

    measure = unique_sites / edit_rows, key = conversion type;
    measure = cell_fraction_bin, key = bin start (per-row edit_fraction, 0-1).
    Unique sites are deduplicated on the conversion type too: before F1 a site
    can carry several types and then counts once per type.
    """
    rows = []
    for label, d in steps:
        label = label.replace("\n", " ")
        sites = d.drop_duplicates(["contig", "position", "strand_conversion"])
        for measure, counts in [
            ("unique_sites", sites["strand_conversion"].value_counts()),
            ("edit_rows",    d["strand_conversion"].value_counts()),
        ]:
            rows += [(label, measure, k, int(v)) for k, v in counts.items()]
        hist, _ = np.histogram(d["edit_fraction"].dropna().clip(upper=1), bins=CELL_FRAC_BINS)
        rows += [(label, "cell_fraction_bin", f"{e:.2f}", int(v))
                 for e, v in zip(CELL_FRAC_BINS[:-1], hist)]
    return pd.DataFrame(rows, columns=["step", "measure", "key", "value"])


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def _grid(n):
    """Two rows once there are more than four panels."""
    nrows = 2 if n > 4 else 1
    ncols = math.ceil(n / nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5.5 * nrows), squeeze=False)
    axes = axes.flatten()
    for ax in axes[n:]:
        ax.set_visible(False)
    return fig, axes


def _pie_grid(stats, measure, unit, title, output_path):
    """One pie of conversion types per filter step, for one measure."""
    d = stats[stats.measure == measure]
    step_order = list(dict.fromkeys(d.step))
    fig, axes = _grid(len(step_order))
    for ax, step in zip(axes, step_order):
        counts = d[d.step == step].set_index("key")["value"].sort_values(ascending=False)
        small = counts / counts.sum() < 0.05     # unreadable as separate slices
        if small.sum() > 1:
            counts = pd.concat([counts[~small], pd.Series({"other": counts[small].sum()})])
        ax.pie(counts, labels=counts.index, autopct="%1.1f%%", startangle=90,
               colors=[PIE_COLORS.get(k, "#7f7f7f") for k in counts.index])
        ax.set_title(f"{step}\n({counts.sum():,} {unit})", fontsize=10)
    fig.suptitle(title, fontsize=13, y=1.01)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", dpi=150)
    plt.close()


def _cell_hist_grid(stats, title, output_path):
    """Per-row (cell, site) edit_fraction from MARINE's coverage, 0 to 1, per filter step."""
    d = stats[stats.measure == "cell_fraction_bin"]
    step_order = list(dict.fromkeys(d.step))
    width = CELL_FRAC_BINS[1] - CELL_FRAC_BINS[0]
    fig, axes = _grid(len(step_order))
    for ax, step in zip(axes, step_order):
        h = d[d.step == step]
        ax.bar(h.key.astype(float), h.value, width=width, align="edge",
               color="skyblue", edgecolor="black")
        ax.set_yscale("log")
        ax.set_xlim(0, 1)
        ax.set_title(f"{step}\n({h.value.sum():,} edit rows)", fontsize=10)
        ax.set_xlabel("Cell editing fraction (edited reads / coverage in that cell)")
        ax.set_ylabel("Edit rows (log scale)")
    fig.suptitle(title, fontsize=13, y=1.01)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", dpi=150)
    plt.close()


def _site_hist_grid(steps, max_frac, title, output_path):
    """Per-site edited fraction (edited reads / depth), 0 to 1, with the F5 cutoff marked."""
    fig, axes = _grid(len(steps))
    for ax, (label, d) in zip(axes, steps):
        finite = d["site_frac"][np.isfinite(d["site_frac"])]
        ax.hist(finite.clip(upper=1), bins=50, range=(0, 1), color="skyblue", edgecolor="black")
        ax.axvline(max_frac, color="crimson", linestyle="--", label=f"cutoff {max_frac}")
        ax.set_yscale("log")
        n_zero = len(d) - len(finite)
        ax.set_title(f"{label}\n({len(d):,} unique edit sites"
                     f"{f'; {n_zero:,} with depth 0 not shown' if n_zero else ''})", fontsize=10)
        ax.set_xlabel("Site editing fraction (edited reads / depth)")
        ax.set_ylabel("Sites (log scale)")
        ax.legend()
    fig.suptitle(title, fontsize=13, y=1.01)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", dpi=150)
    plt.close()


def plot_all(stats, sample, output_dir):
    """Pies and the cell-fraction histogram, over every step in *stats*."""
    tag = f" ({sample})" if sample else ""
    for measure, (unit, title_unit) in PIE_MEASURES.items():
        note = (" — before F1, a site with several conversion types counts once per type"
                if measure == "unique_sites" else "")
        _pie_grid(stats, measure, unit, f"Conversion types by {title_unit} at each filter step{tag}{note}",
                  os.path.join(output_dir, f"piecharts_{measure}.png"))
    _cell_hist_grid(stats, f"Cell editing fraction at each filter step{tag}",
                    os.path.join(output_dir, "histograms_cell_fraction.png"))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="SC Filter 5 (site editing fraction against true depth) and the SC filter plots.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--step-stats",     required=True)
    parser.add_argument("--filtered-edits")
    parser.add_argument("--site-depth")
    parser.add_argument("--barcodes")
    parser.add_argument("--site-max-frac",  type=float, default=0.05)
    parser.add_argument("--min-count",      type=int,   default=3)
    parser.add_argument("--no-recount-min-count", action="store_true", default=False)
    parser.add_argument("--sample", default="", help="Sample name for plot titles")
    parser.add_argument("--output-dir", "-o", default=".")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    stats = pd.read_csv(args.step_stats, sep="\t", dtype={"key": str})

    if not args.site_depth:
        print("F5 off (no --site-depth): plotting F1–F4 only")
        plot_all(stats, args.sample, args.output_dir)
        return
    if not (args.filtered_edits and args.barcodes):
        parser.error("--site-depth needs --filtered-edits and --barcodes")

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
    edits_counted = edits.merge(counted[keys], on=keys)
    out = edits.merge(keep[keys], on=keys)

    recount_label = ("Min-count not re-applied" if args.no_recount_min_count
                     else f"< {args.min_count} edited reads within barcodes")
    summary = pd.DataFrame([
        {"Step": "F1–F4 input",                          "Sites": len(sites_raw),
         "Edit entries": n_raw},
        {"Step": "Restricted to barcodes",               "Sites": len(sites),
         "Edit entries": len(edits)},
        {"Step": f"After {recount_label}",               "Sites": len(counted),
         "Edit entries": len(edits_counted)},
        {"Step": f"After F5 (site frac > {args.site_max_frac})", "Sites": len(keep),
         "Edit entries": len(out)},
    ])
    ct = lambda d: 100 * (d["strand_conversion"] == "C>T").mean() if len(d) else 0.0
    summary["C>T % of sites"] = [ct(sites_raw), ct(sites), ct(counted), ct(keep)]
    print(summary.to_string(index=False))

    sites.to_csv(os.path.join(args.output_dir, "site_frac.tsv"), sep="\t", index=False)
    summary.to_csv(os.path.join(args.output_dir, "site_frac_summary.tsv"), sep="\t", index=False)
    out.to_csv(os.path.join(args.output_dir, "filtered_edits_site_frac.tsv"), sep="\t", index=False)

    # F5's steps appended to F1–F4's, so every plot runs raw → F5 in one figure.
    f5_steps = [("F5: restricted to barcodes", edits),
                (f"F5: after {recount_label}", edits_counted),
                (f"F5: after site frac > {args.site_max_frac}", out)]
    all_stats = pd.concat([stats, step_stats(f5_steps)], ignore_index=True)
    all_stats.to_csv(os.path.join(args.output_dir, "filter_step_stats.tsv"), sep="\t", index=False)
    plot_all(all_stats, args.sample, args.output_dir)

    tag = f" ({args.sample})" if args.sample else ""
    _site_hist_grid([("F5: restricted to barcodes", sites),
                     (f"F5: after {recount_label}", counted),
                     (f"F5: after site frac > {args.site_max_frac}", keep)],
                    args.site_max_frac, f"Site editing fraction at each F5 step{tag}",
                    os.path.join(args.output_dir, "histograms_site_fraction.png"))


if __name__ == "__main__":
    main()
