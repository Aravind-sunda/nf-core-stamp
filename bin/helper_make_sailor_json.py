#!/usr/bin/env python3
"""
Generate a SAILOR JSON config file.

Auto-discovers all .bam files in --samples_path and writes a ready-to-use
config. Override any default with the corresponding flag.

Example:
    python make_sailor_json.py \
        --samples_path /path/to/03_star \
        --output_dir   /path/to/04_sailor \
        --output_json  /path/to/scripts/sailor_config.json
"""

import argparse
import json
import sys
from pathlib import Path


DEFAULTS = {
    # These are now required CLI arguments — no defaults applied:
    # "edit_type":            "CT",
    # "reverse_stranded":     True,
    # "library":              "single",
    # "reference_fasta":      "/home/tmhaxs421/brannanlab/tmhaxs421/reference/cellranger_reference/custom_cellranger/GRCh38-2024-A_mruby_gfp/genome.fa",
    # "known_snps":           "/home/tmhaxs421/brannanlab/Vrutant/Genomes/hg38_bed/hg38_dbsnp_combined.bed3",
    "remove_duplicates":    False,
    "min_variant_coverage": 5,
    "edit_fraction":        0.01,
    "alpha":                0,
    "beta":                 0,
    "junction_overhang":    10,
    "edge_mutation":        5,
    "mm_tolerance":         1,
    "dp":                   "DP4",
    "keep_all_edited":      False,
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("--samples_path",  required=True,
                   help="Directory containing input BAM files.")
    p.add_argument("--output_dir",    required=True,
                   help="Directory where SAILOR will write its outputs.")
    p.add_argument("--output_json",   required=True,
                   help="Path to write the generated JSON config.")
    p.add_argument("--bam_pattern",   default="*.bam",
                   help="Glob pattern for BAM files inside samples_path (default: *.bam).")

    p.add_argument("--edit_type",            required=True,
                   help="Edit type, e.g. CT (C>T / A-to-I) or AG (A>G).")
    p.add_argument("--reverse_stranded",     required=True,
                   type=lambda x: x.lower() not in ("false", "0", "no"),
                   help="Reverse-stranded library: true | false.")
    p.add_argument("--library",              required=True,
                   choices=["single", "paired"],
                   help="Library type: single | paired.")
    p.add_argument("--reference_fasta",      required=True,
                   help="Path to reference genome FASTA.")
    p.add_argument("--known_snps",           required=True,
                   help="Path to known SNPs BED3 file.")
    
    p.add_argument("--remove_duplicates",    default=DEFAULTS["remove_duplicates"],
                   type=lambda x: x.lower() not in ("false", "0", "no"),
                   help=f"Remove PCR duplicates (default: {DEFAULTS['remove_duplicates']}).")
    p.add_argument("--min_variant_coverage", default=DEFAULTS["min_variant_coverage"],
                   type=int,
                   help=f"Minimum read coverage to call a variant (default: {DEFAULTS['min_variant_coverage']}).")
    p.add_argument("--edit_fraction",        default=DEFAULTS["edit_fraction"],
                   type=float,
                   help=f"Minimum edit fraction threshold (default: {DEFAULTS['edit_fraction']}).")
    p.add_argument("--alpha",                default=DEFAULTS["alpha"],
                   type=int,
                   help=f"Alpha pseudocount for Bayesian scoring (default: {DEFAULTS['alpha']}).")
    p.add_argument("--beta",                 default=DEFAULTS["beta"],
                   type=int,
                   help=f"Beta pseudocount for Bayesian scoring (default: {DEFAULTS['beta']}).")
    p.add_argument("--junction_overhang",    default=DEFAULTS["junction_overhang"],
                   type=int,
                   help=f"Minimum junction overhang in nt (default: {DEFAULTS['junction_overhang']}).")
    p.add_argument("--edge_mutation",        default=DEFAULTS["edge_mutation"],
                   type=int,
                   help=f"Minimum nt from read end to call a mismatch (default: {DEFAULTS['edge_mutation']}).")
    p.add_argument("--mm_tolerance",         default=DEFAULTS["mm_tolerance"],
                   type=int,
                   help=f"Max non-target mismatches allowed per read (default: {DEFAULTS['mm_tolerance']}).")
    p.add_argument("--dp",                   default=DEFAULTS["dp"],
                   choices=["DP", "DP4"],
                   help=f"Coverage metric for variant filtering (default: {DEFAULTS['dp']}).")
    p.add_argument("--keep_all_edited",      default=DEFAULTS["keep_all_edited"],
                   type=lambda x: x.lower() not in ("false", "0", "no"),
                   help=f"Keep 100%%-edited sites instead of flagging as SNPs (default: {DEFAULTS['keep_all_edited']}).")
    return p.parse_args()


def discover_bams(samples_path: str, pattern: str) -> list:
    bam_dir = Path(samples_path)
    if not bam_dir.is_dir():
        print(f"[ERROR] samples_path does not exist: {samples_path}", file=sys.stderr)
        sys.exit(1)
    bams = sorted(p.name for p in bam_dir.glob(pattern))
    if not bams:
        print(f"[ERROR] No files matching '{pattern}' found in {samples_path}", file=sys.stderr)
        sys.exit(1)
    return bams


def main():
    args = parse_args()

    bams = discover_bams(args.samples_path, args.bam_pattern)
    print(f"[INFO] Found {len(bams)} BAM file(s) in {args.samples_path}:")
    for b in bams:
        print(f"       {b}")

    config = {
        "samples_path":         args.samples_path,
        "samples":              bams,
        "reverse_stranded":     args.reverse_stranded,
        "library":              args.library,
        "edit_type":            args.edit_type.upper(),
        "reference_fasta":      args.reference_fasta,
        "known_snps":           args.known_snps,
        "output_dir":           args.output_dir,
        "remove_duplicates":    args.remove_duplicates,
        "min_variant_coverage": args.min_variant_coverage,
        "edit_fraction":        args.edit_fraction,
        "alpha":                args.alpha,
        "beta":                 args.beta,
        "junction_overhang":    args.junction_overhang,
        "edge_mutation":        args.edge_mutation,
        "mm_tolerance":         args.mm_tolerance,
        "dp":                   args.dp,
        "keep_all_edited":      args.keep_all_edited,
    }

    out_path = Path(args.output_json)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(config, f, indent=2)

    print(f"[DONE] SAILOR config written to {out_path}")


if __name__ == "__main__":
    main()
