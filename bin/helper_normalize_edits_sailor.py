#!/usr/bin/env python3

"""
Normalize SAILOR ranked edit sites to EPR, EPKM, EPKMR, EPM, EPMR metrics.

The SAILOR counterpart of helper_normalize_edits_bulk.py. It takes SAILOR's
combined ranked BED instead of MARINE's site table, annotates each site to a gene,
and then applies the identical normalisation formulas so the two callers' outputs
are directly comparable.

Steps:
1) Read the combined ranked BED and parse per-site edited/total read counts
2) Drop sites below the confidence threshold
3) Annotate sites with genes from the BED6 gene model
4) Attach gene lengths and per-sample library counts from the featureCounts matrix
5) Compute EPR, EPKM, EPKMR, EPM, EPMR per gene

Input BED layout (produced by SAILOR's join_edits rule):

    chrom  start  end  confidence  edited,total  strand
    chr1   629183 629184  0.860058355   1,16      +

Note this differs from the per-strand BEDs under 7_scored_outputs/, where column 4
is 'coverage|conversion|edit_fraction' and column 5 is the confidence. Only the
combined BED is consumed here.
"""

import sys
import argparse
import logging
import os

import pandas as pd
import pybedtools


logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)


BED_COLS = ["contig", "start", "position", "confidence", "reads", "strand"]


def read_ranked_bed(path):
    """Read a SAILOR combined ranked BED and split the 'edited,total' field."""
    df = pd.read_csv(path, sep="\t", header=None, names=BED_COLS,
                     dtype={"contig": str})
    if df.empty:
        return df.assign(count=pd.Series(dtype=int), coverage=pd.Series(dtype=int))

    reads = df["reads"].astype(str).str.split(",", expand=True)
    if reads.shape[1] != 2:
        raise ValueError(
            f"{path}: column 5 should look like 'edited,total' but got "
            f"{df['reads'].iloc[0]!r}. This script expects SAILOR's *combined* "
            f"ranked BED, not the per-strand files in 7_scored_outputs/."
        )
    df["count"]    = pd.to_numeric(reads[0], errors="coerce")
    df["coverage"] = pd.to_numeric(reads[1], errors="coerce")

    bad = df[["count", "coverage"]].isna().any(axis=1).sum()
    if bad:
        logger.warning(f"{bad:,} rows had unparseable read counts — dropping.")
        df = df.dropna(subset=["count", "coverage"])

    df["count"]    = df["count"].astype(int)
    df["coverage"] = df["coverage"].astype(int)
    return df


def annotate_sites(df, annotation_bed, strand_aware):
    """Annotate sites with gene info by intersecting the BED6 gene model."""
    sites = df[["contig", "start", "position", "strand"]].drop_duplicates().copy()
    sites["name"] = sites["contig"] + "_" + sites["position"].astype(str)

    # bedtools `-sorted` streams both files instead of building an interval tree of
    # the whole annotation. It needs both sides in the same order, and PREPARE_DBSNP
    # sorts references with `sort -k1,1 -k2,2n` (lexicographic contig, numeric start),
    # which is what pandas' string ordering gives us here.
    sites = sites.sort_values(["contig", "start"], kind="mergesort")

    cols = ["contig", "start", "position", "name", "score", "strand"]
    sites["score"] = 0
    edits_bt = pybedtools.BedTool.from_dataframe(sites[cols])
    annot_bt = pybedtools.BedTool(annotation_bed)

    # loj keeps unannotated sites (feature_name '.') so the caller can report them;
    # -s restricts matches to the same strand for stranded libraries.
    intersect = edits_bt.intersect(annot_bt, wb=True, loj=True, s=strand_aware).to_dataframe(
        header=None,
        names=[
            "contig", "start", "position", "name", "score", "strand",
            "feature_chrom", "feature_start", "feature_end",
            "feature_name", "feature_type", "feature_strand",
        ],
    )
    intersect = intersect[["name", "feature_name", "feature_type", "feature_strand"]]

    df = df.copy()
    df["_site_key"] = df["contig"] + "_" + df["position"].astype(str)
    intersect = intersect.rename(columns={"name": "_site_key"})

    df = df.merge(intersect, on="_site_key").drop(columns=["_site_key"])

    unannotated = (df["feature_name"] == ".").sum()
    if unannotated:
        logger.info(f"  {unannotated:,} site-rows fell outside any gene — dropping.")
        df = df[df["feature_name"] != "."]
    logger.info(f"  Sites after annotation: {len(df):,}")
    return df


def extract_gene_lengths_from_featurecounts(feature_counts_path):
    """Extract Geneid and Length columns from a featureCounts matrix."""
    fc = pd.read_csv(feature_counts_path, sep="\t", comment="#")
    return fc[["Geneid", "Length"]].copy()


def process_feature_counts(feature_counts_path):
    """Read the featureCounts matrix and return Geneid + sample count columns."""
    feature_counts = pd.read_csv(feature_counts_path, sep="\t", comment="#")
    # gene_id comes from featureCounts --extraAttributes and is annotation, not a
    # sample; leaving it in would make it look like a count column downstream.
    # errors="ignore" keeps matrices built before that flag was added working.
    return feature_counts.drop(
        columns=["Chr", "Start", "End", "Strand", "Length", "gene_id"], errors="ignore"
    )


def normalizing_sailor_edits(sites, feature_counts, sample_name):
    """Aggregate per-gene and compute the same metrics as the MARINE bulk path."""
    if sample_name not in feature_counts.columns:
        logger.error(f"Sample '{sample_name}' not found in featureCounts.")
        logger.error(f"Available samples: {feature_counts.columns.tolist()}")
        raise ValueError(f"Sample {sample_name} not found in feature counts.")

    # Library size comes from all genes before any filtering so the EPKM/EPM
    # denominators reflect true depth, matching helper_normalize_edits_bulk.py.
    total_reads = feature_counts[sample_name].sum()
    if total_reads == 0:
        raise ValueError(f"Total read count for sample '{sample_name}' is 0 — cannot normalize.")
    logger.info(f"Total library size for '{sample_name}': {total_reads:,} reads")

    feature_counts = feature_counts[["Geneid", sample_name]]

    normalized_matrix = sites.groupby(["contig", "feature_name"]).agg(
        total_edits=("count", "sum"),
        total_coverage=("coverage", "sum"),
        feature_type=("feature_type", "first"),
        feature_strand=("feature_strand", "first"),
        strand=("strand", lambda x: ",".join(x.unique()) if x.nunique() > 1 else x.iloc[0]),
        mean_confidence=("confidence", "mean"),
        n_sites=("count", "size"),
        Length=("Length", "first"),
    ).reset_index()

    normalized_matrix = normalized_matrix.merge(
        feature_counts, left_on="feature_name", right_on="Geneid", how="left"
    ).drop(columns=["Geneid"])

    missing_counts = normalized_matrix[sample_name].isna().sum()
    if missing_counts > 0:
        logger.warning(f"{missing_counts:,} genes have no featureCounts entry — dropping.")
    normalized_matrix = normalized_matrix.dropna(subset=[sample_name])

    zero_counts = (normalized_matrix[sample_name] == 0).sum()
    if zero_counts > 0:
        logger.warning(f"{zero_counts:,} genes have zero featureCounts — dropping to avoid division by zero.")
    normalized_matrix = normalized_matrix[normalized_matrix[sample_name] != 0]

    normalized_matrix["EPR"]   = normalized_matrix["total_edits"] / normalized_matrix[sample_name]
    normalized_matrix["EPKM"]  = normalized_matrix["total_edits"] / (total_reads / 1e6 * (normalized_matrix["Length"] / 1e3))
    normalized_matrix["EPKMR"] = normalized_matrix["EPKM"] / normalized_matrix[sample_name]
    normalized_matrix["EPM"]   = normalized_matrix["total_edits"] / (total_reads / 1e6)
    normalized_matrix["EPMR"]  = normalized_matrix["EPM"] / normalized_matrix[sample_name]

    normalized_matrix["sample_name"] = sample_name
    normalized_matrix = normalized_matrix.rename(columns={sample_name: "featureCount_count"})

    logger.info(f"Final output: {len(normalized_matrix):,} rows across "
                f"{normalized_matrix['feature_name'].nunique()} genes")

    # Column order matches helper_normalize_edits_bulk.py, plus the two SAILOR-only
    # columns, so the tables can be concatenated or compared directly.
    return normalized_matrix[[
        "sample_name", "contig", "feature_name", "feature_strand", "feature_type",
        "Length", "strand", "total_edits", "total_coverage",
        "featureCount_count", "EPR", "EPKM", "EPKMR", "EPM", "EPMR",
        "n_sites", "mean_confidence",
    ]]


def main():
    parser = argparse.ArgumentParser(
        description="Normalize SAILOR ranked edits to EPR/EPKM/EPKMR/EPM/EPMR metrics."
    )
    parser.add_argument("-i", "--input_bed", required=True,
                        help="SAILOR combined ranked BED", dest="input_bed")
    parser.add_argument("-c", "--feature_counts_matrix", required=True,
                        help="Path to the featureCounts matrix", dest="feature_counts_matrix")
    parser.add_argument("-a", "--annotation_bed", required=True,
                        help="BED6 gene model used to assign sites to genes", dest="annotation_bed")
    parser.add_argument("-s", "--sample_name", required=True,
                        help="Sample name matching a column in the featureCounts matrix", dest="sample_name")
    parser.add_argument("-str", "--strandedness", type=int, choices=[0, 1, 2], required=True,
                        help="Strandedness mode: 0=unstranded, 1=forward, 2=reverse", dest="strandedness")
    parser.add_argument("--min-confidence", type=float, default=0.5,
                        help="Drop sites with SAILOR confidence below this (default: 0.5)",
                        dest="min_confidence")
    parser.add_argument("-d", "--output_directory", default=".",
                        help="Output directory (default: current directory)", dest="output_directory")

    args = parser.parse_args()
    os.makedirs(args.output_directory, exist_ok=True)

    # ---- Load ---------------------------------------------------------------
    logger.info(f"Loading {args.input_bed} ...")
    df = read_ranked_bed(args.input_bed)
    logger.info(f"Loaded {len(df):,} ranked sites")

    # ---- Confidence filter --------------------------------------------------
    before = len(df)
    df = df[df["confidence"] >= args.min_confidence]
    logger.info(f"Confidence >= {args.min_confidence}: {len(df):,} / {before:,} sites retained")

    out_name = f"{args.sample_name}.sailor.EPR_EPKM_normalized.tsv"
    out_path = os.path.join(args.output_directory, out_name)

    if df.empty:
        logger.warning("No sites passed the confidence filter — writing an empty table.")
        pd.DataFrame(columns=[
            "sample_name", "contig", "feature_name", "feature_strand", "feature_type",
            "Length", "strand", "total_edits", "total_coverage",
            "featureCount_count", "EPR", "EPKM", "EPKMR", "EPM", "EPMR",
            "n_sites", "mean_confidence",
        ]).to_csv(out_path, sep="\t", index=False)
        logger.info(f"Saved normalized edits → {out_name}")
        return

    # ---- Annotate sites to genes -------------------------------------------
    # Unstranded libraries cannot resolve which strand a site's gene is on, so the
    # intersect is left strand-unaware there, mirroring the MARINE bulk behaviour.
    strand_aware = args.strandedness != 0
    logger.info(f"Annotating sites with {args.annotation_bed} "
                f"({'strand-aware' if strand_aware else 'strand-unaware'}) ...")
    df = annotate_sites(df, args.annotation_bed, strand_aware)

    # ---- Gene lengths -------------------------------------------------------
    logger.info("Extracting gene lengths from featureCounts Length column ...")
    gene_lengths = extract_gene_lengths_from_featurecounts(args.feature_counts_matrix)
    df = df.merge(gene_lengths, left_on="feature_name", right_on="Geneid", how="left").drop(columns=["Geneid"])

    missing = df["Length"].isna().sum()
    if missing:
        logger.warning(f"{missing:,} rows have no gene length (not in featureCounts) — dropping.")
        df = df.dropna(subset=["Length"])

    # ---- Normalize ----------------------------------------------------------
    logger.info("Processing featureCounts matrix ...")
    feature_counts = process_feature_counts(args.feature_counts_matrix)

    logger.info("Normalizing edits ...")
    normalized = normalizing_sailor_edits(df, feature_counts, args.sample_name)

    normalized.to_csv(out_path, sep="\t", index=False)
    logger.info(f"Saved normalized edits → {out_name}")


if __name__ == "__main__":
    if sys.version_info < (3, 8, 0):
        sys.stderr.write("You need Python 3.8 or later to run this script\n")
        sys.exit(1)
    main()
