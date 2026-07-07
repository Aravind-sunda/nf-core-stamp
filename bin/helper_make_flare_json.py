#!/usr/bin/env python3
"""
Generate a FLARE JSON config for a single sample.

FLARE runs downstream of SAILOR and consumes that sample's SAILOR outputs
(ranked BED, forward/reverse bigwigs, filtered-merged BAM). This helper writes
the per-sample config the bundled FLARE Snakefile expects, for either mode:

  * Cluster-identification mode  (--regions <folder>):  discovers and scores new
    edited clusters. Sets "keep_all": false and the "regions" key.
  * Edit-fraction mode           (--regions_file <tsv>): quantifies editing within
    pre-defined regions. Sets "keep_all": true and the "regions_file" key.

Exactly one of --regions / --regions_file must be given; the mode (and keep_all)
is selected automatically from which one is supplied.

Example (cluster identification):
    helper_make_flare_json.py \
        --label       sample1 \
        --output_json sample1_flare.json \
        --output_folder ./flare_output \
        --stamp_sites_file sample1.combined.....ranked.bed \
        --forward_bw  sample1.fwd.sorted.bw \
        --reverse_bw  sample1.rev.sorted.bw \
        --bam         sample1_filtered_merged.sorted.bam \
        --fasta       genome.fa \
        --regions     hg38_regions \
        --edit_type   CT \
        --fdr_threshold 0.1 \
        --max_merge_dist 15

Example (edit fraction): as above but replace --regions with
        --regions_file my_regions_of_interest.tsv
"""

import argparse
import json
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("--label",           required=True,
                   help="Sample ID; used as the FLARE output label.")
    p.add_argument("--output_json",     required=True,
                   help="Path to write the generated JSON config.")
    p.add_argument("--output_folder",   required=True,
                   help="Directory where FLARE will write its outputs.")

    p.add_argument("--stamp_sites_file", required=True,
                   help="SAILOR ranked BED for this sample (the STAMP/edit sites).")
    p.add_argument("--forward_bw",      required=True,
                   help="SAILOR forward-strand bigwig (.fwd.sorted.bw).")
    p.add_argument("--reverse_bw",      required=True,
                   help="SAILOR reverse-strand bigwig (.rev.sorted.bw).")
    p.add_argument("--bam",             required=True,
                   help="SAILOR filtered-merged sorted BAM (_filtered_merged.sorted.bam).")
    p.add_argument("--fasta",           required=True,
                   help="Reference genome FASTA.")

    # Mode selector: exactly one of these. --regions → cluster identification,
    # --regions_file → edit fraction.
    p.add_argument("--regions",         default=None,
                   help="FLARE regions FOLDER from generate_regions.py → cluster-identification mode (keep_all=false).")
    p.add_argument("--regions_file",    default=None,
                   help="Pre-defined regions FILE (.tsv) → edit-fraction mode (keep_all=true).")

    p.add_argument("--edit_type",       required=True,
                   help="Edit type in 2-char form, e.g. CT (C>T) or AG (A>G).")
    p.add_argument("--fdr_threshold",   default=0.1, type=float,
                   help="FDR threshold for cluster/peak filtering (default: 0.1).")
    p.add_argument("--max_merge_dist",  default=15, type=int,
                   help="Max distance (bp) for merging adjacent peaks (default: 15).")
    return p.parse_args()


def main():
    args = parse_args()

    # Select the mode from the regions argument supplied. Exactly one is allowed:
    # a folder → cluster identification; a file → edit fraction (keep_all=true).
    # Neither given is an error, and both given at once is also an error.
    if not args.regions and not args.regions_file:
        sys.exit(
            "[ERROR] Missing regions: provide --regions <folder> (cluster-identification "
            "mode) or --regions_file <file> (edit-fraction mode). Neither was given."
        )
    if args.regions and args.regions_file:
        sys.exit(
            "[ERROR] Provide only ONE of --regions (cluster-identification mode) or "
            "--regions_file (edit-fraction mode), not both."
        )

    config = {
        "label":          args.label,
        "output_folder":  args.output_folder,
        "stamp_sites_file": args.stamp_sites_file,
        "forward_bw":     args.forward_bw,
        "reverse_bw":     args.reverse_bw,
        "bam":            args.bam,
        "fasta":          args.fasta,
        "edit_type":      args.edit_type.upper(),
        "fdr_threshold":  args.fdr_threshold,
        "max_merge_dist": args.max_merge_dist,
    }

    if args.regions:
        # Cluster-identification mode: a regions *folder*, keep_all disabled.
        config["regions"]  = args.regions
        config["keep_all"] = False
    else:
        # Edit-fraction mode: a regions *file*, keep_all enabled.
        config["regions_file"] = args.regions_file
        config["keep_all"]     = True

    out_path = Path(args.output_json)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(config, f, indent=2)

    print(f"[DONE] FLARE config written to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
