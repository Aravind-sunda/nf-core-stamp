#!/usr/bin/env python3
"""
metagene_plot.py — overlay metagene distributions from one or more
metaPlotR-style "distance measures" files (the output of metaplot_dist.py).

Each input file must already be collapsed to the longest isoform per gene
(i.e. each gene_name maps to exactly one refseqID). This is enforced before
plotting; regenerate the input with metaplot_dist.py's default (longest
isoform) if the check fails.

The x-axis is the metagene coordinate rel_location: 0-1 = 5'UTR, 1-2 = CDS,
2-3 = 3'UTR. With --rescale, the 5'UTR and 3'UTR widths are shrunk relative to
the CDS (set to width 1) using the median region lengths, reproducing the
metaPlotR rescaling:
    utr5.SF = median(utr5_size) / median(cds_size)
    utr3.SF = median(utr3_size) / median(cds_size)
    5'UTR coords [0,1] -> [1-utr5.SF, 1];  3'UTR coords [2,3] -> [2, 2+utr3.SF]

Default plot is a per-library KDE density (each library normalised to area 1),
so libraries with different numbers of sites are compared by shape, not depth.
Use --counts to scale each KDE by its site count (area under curve = N sites),
making between-sample editing depth differences visible while keeping smooth
lines. Add --histogram to switch from smooth KDE lines to binned bars instead.

Requires: pandas, numpy, scipy, matplotlib.
"""

import argparse
import math
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

REQUIRED_COLS = ["gene_name", "refseqID", "rel_location",
                 "utr5_size", "cds_size", "utr3_size"]


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def load_and_check(path, allow_multi_isoform=False):
    """Read a dist-measures file; verify one isoform per gene; return df."""
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        log(f"ERROR: {path} is missing required column(s): {missing}")
        log(f"       found columns: {list(df.columns)}")
        sys.exit(1)

    # rel_location is NA for sites outside annotated regions (e.g. --keep-incomplete)
    n_before = len(df)
    df = df[df["rel_location"].notna()].copy()
    dropped = n_before - len(df)
    if dropped:
        log(f"  {os.path.basename(path)}: dropped {dropped} rows with NA rel_location")

    # Longest-isoform enforcement: each gene_name -> exactly one refseqID
    per_gene = df.groupby("gene_name")["refseqID"].nunique()
    bad = per_gene[per_gene > 1]
    if len(bad) > 0:
        log(f"ERROR: {path} is NOT collapsed to one isoform per gene.")
        log(f"       {len(bad)} gene(s) have multiple transcripts, e.g.: "
            f"{list(bad.index[:5])}")
        log("       Regenerate with metaplot_dist.py using its default "
            "(longest isoform), i.e. WITHOUT --all-isoforms.")
        if not allow_multi_isoform:
            sys.exit(1)
        log("       --allow-multi-isoform set: continuing anyway.")
    return df


def compute_scale_factors(dfs):
    """Pool unique transcripts across all libraries and derive UTR scale
    factors from median region lengths (one row per transcript)."""
    uniq = pd.concat(
        [d[["refseqID", "utr5_size", "cds_size", "utr3_size"]] for d in dfs],
        ignore_index=True,
    ).drop_duplicates(subset="refseqID")
    m5 = np.nanmedian(uniq["utr5_size"])
    mc = np.nanmedian(uniq["cds_size"])
    m3 = np.nanmedian(uniq["utr3_size"])
    utr5_sf = m5 / mc
    utr3_sf = m3 / mc
    log(f"  median region lengths (per transcript): 5'UTR={m5:.0f}, "
        f"CDS={mc:.0f}, 3'UTR={m3:.0f}")
    log(f"  scale factors: utr5.SF={utr5_sf:.3f}, utr3.SF={utr3_sf:.3f}")
    return utr5_sf, utr3_sf


def rescale_coords(r, utr5_sf, utr3_sf):
    """Map rel_location into rescaled space (CDS fixed to [1,2])."""
    r = np.asarray(r, dtype=float)
    out = r.copy()
    in5 = r < 1
    in3 = r >= 2
    out[in5] = (1 - utr5_sf) + r[in5] * utr5_sf       # [0,1] -> [1-sf, 1]
    out[in3] = 2 + (r[in3] - 2) * utr3_sf             # [2,3] -> [2, 2+sf]
    return out


def draw_series(ax, c, color, grid, xlo, xhi, histogram, counts, bins,
                bw_adjust, label=None):
    """Draw one library's distribution on ax. Returns False if a KDE could
    not be computed (too few/degenerate points).

    With counts=True and histogram=False: KDE is scaled by N (number of sites)
    so the area under the curve equals N. This makes between-sample differences
    in editing depth visible while keeping smooth lines.
    """
    if histogram:
        ax.hist(c, bins=bins, range=(xlo, xhi), density=not counts,
                histtype="step", linewidth=1.6, color=color, label=label)
        return True
    if len(c) < 2 or np.ptp(c) == 0:
        return False
    kde = gaussian_kde(c)
    kde.set_bandwidth(kde.factor * bw_adjust)
    y = kde(grid)
    if counts:
        # Scale density by N: area under curve = N sites.
        # y-axis becomes "sites per unit metagene coordinate" — between-sample
        # depth differences are visible while the curve remains smooth.
        y = y * len(c)
    ax.plot(grid, y, color=color, linewidth=2.0, label=label)
    ax.fill_between(grid, y, color=color, alpha=0.08)
    return True


def decorate_axis(ax, xlo, xhi, center5, center3, ytop, label_regions=True):
    """Shade UTR regions, draw CDS boundaries, and (optionally) label regions."""
    ax.axvspan(xlo, 1, color="0.5", alpha=0.05, zorder=0)
    ax.axvspan(2, xhi, color="0.5", alpha=0.05, zorder=0)
    for vx in (1, 2):
        ax.axvline(vx, color="grey", linestyle="-", linewidth=0.9, alpha=0.7)
    if label_regions:
        for cx, name in ((center5, "5'UTR"), (1.5, "CDS"), (center3, "3'UTR")):
            ax.text(cx, ytop * 0.97, name, ha="center", va="top",
                    fontsize=9, color="0.25")
    ax.set_xlim(xlo, xhi)
    ax.spines[["top", "right"]].set_visible(False)


def main():
    ap = argparse.ArgumentParser(
        description="Overlay metagene distributions from one or more "
                    "dist-measures files.")
    ap.add_argument("files", nargs="+",
                    help="dist-measures file(s) from metaplot_dist.py.")
    ap.add_argument("-o", "--out", default="metagene.png",
                    help="Output image (extension sets format: .png/.pdf/.svg).")
    ap.add_argument("--labels", default=None,
                    help="Comma-separated legend labels, one per file "
                         "(default: file basenames).")
    ap.add_argument("--rescale", action="store_true",
                    help="Rescale 5'UTR/3'UTR widths by median length "
                         "relative to CDS (metaPlotR-style).")
    ap.add_argument("--histogram", action="store_true",
                    help="Draw binned histograms instead of KDE density lines.")
    ap.add_argument("--counts", action="store_true",
                    help="Scale each curve/bar by its site count (area = N sites). "
                         "Makes between-sample editing depth differences visible. "
                         "Works for both KDE (default smooth lines) and --histogram.")
    ap.add_argument("--bins", type=int, default=100,
                    help="Histogram bin count (default 100).")
    ap.add_argument("--bw-adjust", type=float, default=1.0,
                    help="KDE bandwidth multiplier (default 1.0; smaller = "
                         "sharper).")
    ap.add_argument("--facet", action="store_true",
                    help="Draw one panel per library (small multiples) instead "
                         "of overlaying all curves on one axis.")
    ap.add_argument("--ncols", type=int, default=None,
                    help="Columns in the facet grid (default: min(n_files, 3)).")
    ap.add_argument("--title", default="Metagene distribution of sites")
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--figsize", default=None,
                    help="Figure size W,H in inches. Overlay default 8,4.5; "
                         "facet default scales with the grid.")
    ap.add_argument("--allow-multi-isoform", action="store_true",
                    help="Override the one-isoform-per-gene requirement "
                         "(not recommended).")
    args = ap.parse_args()

    # Labels
    if args.labels:
        labels = [x.strip() for x in args.labels.split(",")]
        if len(labels) != len(args.files):
            log(f"ERROR: got {len(labels)} labels for {len(args.files)} files.")
            sys.exit(1)
    else:
        labels = [os.path.splitext(os.path.basename(f))[0] for f in args.files]

    # Load + check every file
    log("Loading and checking input files...")
    dfs = [load_and_check(f, args.allow_multi_isoform) for f in args.files]
    for lab, d in zip(labels, dfs):
        log(f"  {lab}: {len(d)} sites")

    # Optional rescaling
    if args.rescale:
        utr5_sf, utr3_sf = compute_scale_factors(dfs)
        coords = [rescale_coords(d["rel_location"].values, utr5_sf, utr3_sf)
                  for d in dfs]
        xlo, xhi = 1 - utr5_sf, 2 + utr3_sf
        center5 = (xlo + 1) / 2
        center3 = (2 + xhi) / 2
    else:
        coords = [d["rel_location"].values for d in dfs]
        xlo, xhi = 0.0, 3.0
        center5, center3 = 0.5, 2.5

    # Common plotting setup
    cmap = plt.get_cmap("tab10" if len(dfs) <= 10 else "tab20")
    grid = np.linspace(xlo, xhi, 512)
    xlabel = ("Metagene coordinate"
              + (" (rescaled by median region length)" if args.rescale
                 else "  (0-1 5'UTR | 1-2 CDS | 2-3 3'UTR)"))
    ylabel = "Site density (count-scaled)" if args.counts else "Density"
    fw = fh = None
    if args.figsize:
        fw, fh = (float(x) for x in args.figsize.split(","))

    if args.facet:
        n = len(dfs)
        ncols = args.ncols if args.ncols else min(n, 3)
        nrows = math.ceil(n / ncols)
        figsize = (fw, fh) if fw else (ncols * 4.2, nrows * 3.0)
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize,
                                 sharex=True, sharey=True, squeeze=False)
        flat = axes.flatten()
        for i, (lab, c) in enumerate(zip(labels, coords)):
            ok = draw_series(flat[i], c, cmap(i % cmap.N), grid, xlo, xhi,
                             args.histogram, args.counts, args.bins,
                             args.bw_adjust)
            flat[i].set_title(f"{lab} (n={len(c)})", fontsize=10)
            if not ok:
                flat[i].text(0.5, 0.5, "insufficient data for KDE",
                             transform=flat[i].transAxes, ha="center",
                             va="center", fontsize=9, color="0.4")
        for j in range(n, len(flat)):
            flat[j].axis("off")
        ytop = max(flat[i].get_ylim()[1] for i in range(n)) * 1.08
        for i in range(n):
            flat[i].set_ylim(top=ytop)
            decorate_axis(flat[i], xlo, xhi, center5, center3, ytop,
                          label_regions=(i // ncols == 0))
        fig.supxlabel(xlabel, fontsize=10)
        fig.supylabel(ylabel, fontsize=10)
        fig.suptitle(args.title, fontsize=12)
        fig.tight_layout()
    else:
        fig, ax = plt.subplots(figsize=(fw, fh) if fw else (8, 4.5))
        for i, (lab, c) in enumerate(zip(labels, coords)):
            drawn = draw_series(ax, c, cmap(i % cmap.N), grid, xlo, xhi,
                                args.histogram, args.counts, args.bins,
                                args.bw_adjust, label=f"{lab} (n={len(c)})")
            if not drawn:
                log(f"WARN: {labels[i]} has too few/invariant points for KDE; "
                    "skipping.")
        ytop = ax.get_ylim()[1] * 1.08
        ax.set_ylim(top=ytop)
        decorate_axis(ax, xlo, xhi, center5, center3, ytop, label_regions=True)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(args.title)
        ax.legend(frameon=False, fontsize=9, loc="center left",
                  bbox_to_anchor=(1.01, 0.5))
        fig.tight_layout()

    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    log(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
