#!/usr/bin/env python3

"""
Author: Aravind Sundaravadivelu
Date: 2026-06-08
Description: Normalize MARINE bulk RNA-seq editing output to EPR, EPKM, EPKMR, EPM, EPMR metrics.

The script follows these steps:
1) Read the filtered MARINE input file
2) Filter rows to the requested strand_conversion type (strandedness-aware)
3) Extract gene lengths directly from the featureCounts Length column
4) Annotate the marine edits with gene lengths
5) Read the featureCounts matrix and annotate edits with sample read counts
6) Compute EPR, EPKM, EPKMR, EPM, EPMR per gene per sample

NOTE: The working of this script can be tested in
      /home/tmhaxs421/brannanlab/tmhaxs421/scripts/Marine_scripts/bulk/02_EPR_normalization_trial.ipynb
"""

import sys
import argparse
import logging
import os

import pandas as pd


logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)


def complement_edit(edit: str) -> str:
    """Return the strand-complement of an edit string e.g. 'C>T' -> 'G>A'."""
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    ref, alt = edit.split(">")
    return f"{comp[ref]}>{comp[alt]}"


def _join_unique(series):
    """Collapse a group's values to a single cell, comma-joining genuine variation."""
    uniq = series.unique()
    return ",".join(str(v) for v in uniq) if len(uniq) > 1 else series.iloc[0]


def _warn_multi_locus_genes(df):
    """Log gene symbols whose edits span more than one contig.

    A symbol can map to several distinct gene_ids — 21 do in GRCh38-2024-A, e.g.
    LSP1P5 is ENSG00000288905 on chr1 and ENSG00000291150 on chr2. featureCounts
    (-g gene_name) already merges those into one count, so this script has to
    merge them too; the alternative divides one locus' edits by every locus'
    reads. The merge is therefore correct given the counting scheme, but it does
    conflate distinct genes, so name them rather than let it pass silently.
    Keying the whole pipeline on gene_id would be the way to remove the ambiguity.
    """
    spans = df.groupby("feature_name")["contig"].nunique()
    multi = spans[spans > 1]
    if len(multi):
        logger.warning(
            f"{len(multi)} gene symbol(s) have edits on more than one contig and will be "
            f"merged into a single row, matching how featureCounts counted them: "
            f"{', '.join(multi.index.astype(str)[:10])}"
            f"{' ...' if len(multi) > 10 else ''}"
        )


def extract_gene_lengths_from_featurecounts(feature_counts_path):
    """Extract Geneid and Length columns from a featureCounts output file."""
    fc = pd.read_csv(feature_counts_path, sep="\t", comment="#")
    return fc[["Geneid", "Length"]].copy()


def annotate_gene_length(df, gene_lengths_df):
    """Merge featureCounts gene lengths into marine edits on feature_name. Drops rows with no match."""
    annotated = df.merge(
        gene_lengths_df, left_on="feature_name", right_on="Geneid", how="left"
    ).drop(columns=["Geneid"])

    missing = annotated["Length"].isna().sum()
    if missing > 0:
        missing_genes = annotated.loc[annotated["Length"].isna(), "feature_name"].unique()
        logger.warning(f"{missing:,} rows have no gene length (feature_name not in featureCounts) — dropping.")
        logger.warning(f"  Sample genes not found: {missing_genes[:10].tolist()}")
    annotated = annotated.dropna(subset=["Length"])
    return annotated


def _dedupe_annotation_lists(names, types, strands):
    """Drop repeats of the same gene from one site's annotation lists.

    A site inside two overlapping loci that share a symbol is annotated
    'MKKS,MKKS,MKKS'. Exploding that as-is yields three rows for one site, and the
    per-gene sum then counts its reads three times. Collapsing to the first entry
    per name keeps one row per site-gene pair.

    Deliberately keyed on the name alone, not the (name, type, strand) triple:
    four ambiguous symbols have loci with different biotypes that overlap (e.g.
    ELFN2 is lncRNA and protein_coding on chr22), so a triple would stay distinct
    and the double count would survive. The discarded biotype costs nothing —
    feature_type is only ever tested for '-1' (unannotated) and otherwise carried
    through as a label; it never enters a metric.
    """
    seen = set()
    kept = ([], [], [])
    for name, ftype, strand in zip(names, types, strands):
        if name in seen:
            continue
        seen.add(name)
        for bucket, value in zip(kept, (name, ftype, strand)):
            bucket.append(value)
    return kept


def split_multi_annotated_sites(df):
    """Explode rows where feature_name/feature_type/feature_strand hold comma-separated
    values (site overlapping multiple genes) into one row per annotation.
    All other columns (count, coverage, etc.) are duplicated across the new rows.
    Repeats of the same gene within a site are collapsed first, so a site is only
    ever counted once per gene.
    """
    multi = df["feature_name"].str.contains(",", na=False)
    n_multi = multi.sum()
    if n_multi == 0:
        return df
    logger.info(f"Splitting {n_multi:,} multi-annotated rows into per-gene entries ...")
    df = df.copy()
    cols = ["feature_name", "feature_type", "feature_strand"]
    for col in cols:
        df[col] = df[col].str.split(",")

    # Ragged annotations would silently truncate under zip(), so refuse them.
    lengths = pd.DataFrame({c: df[c].str.len() for c in cols})
    ragged = (lengths.nunique(axis=1) > 1).sum()
    if ragged:
        raise ValueError(
            f"{ragged:,} row(s) have differing counts of comma-separated values across "
            f"{cols}. Cannot align annotations to genes."
        )

    n_before = df[cols[0]].str.len().sum()
    deduped = [
        _dedupe_annotation_lists(n, t, s)
        for n, t, s in zip(df["feature_name"], df["feature_type"], df["feature_strand"])
    ]
    for i, col in enumerate(cols):
        df[col] = [d[i] for d in deduped]
    n_after = df[cols[0]].str.len().sum()
    if n_after != n_before:
        logger.info(
            f"  Collapsed {n_before - n_after:,} repeated same-gene annotation(s) "
            f"so each site counts once per gene"
        )

    df = df.explode(cols).reset_index(drop=True)
    logger.info(f"  Rows after split: {len(df):,}")
    return df


def process_feature_counts(feature_counts_path):
    """Read featureCounts matrix and return only Geneid + sample count columns."""
    feature_counts = pd.read_csv(feature_counts_path, sep="\t", comment="#")
    # gene_id comes from featureCounts --extraAttributes and is annotation, not a
    # sample; leaving it in would make it look like a count column downstream.
    # errors="ignore" keeps matrices built before that flag was added working.
    feature_counts = feature_counts.drop(
        columns=["Chr", "Start", "End", "Strand", "Length", "gene_id"], errors="ignore"
    )
    return feature_counts


def normalizing_marine_edits(marine_counts, feature_counts, sample_name):
    marine_counts = marine_counts[[
        "contig", "feature_name", "feature_strand", "feature_type",
        "Length", "count", "coverage", "strand_conversion", "strand",
    ]]

    if sample_name not in feature_counts.columns:
        logger.error(f"Sample '{sample_name}' not found in featureCounts.")
        logger.error(f"Available samples: {feature_counts.columns.tolist()}")
        raise ValueError(f"Sample {sample_name} not found in feature counts.")

    # Compute library size from all genes before any filtering so that EPKM/EPM
    # denominators reflect the true total read depth, not the filtered subset.
    total_reads = feature_counts[sample_name].sum()
    
    if total_reads == 0:
        raise ValueError(f"Total read count for sample '{sample_name}' is 0 — cannot normalize.")
    logger.info(f"Total library size for '{sample_name}': {total_reads:,} reads")

    feature_counts = feature_counts[["Geneid", sample_name]]

    _warn_multi_locus_genes(marine_counts)

    # Grouped by feature_name alone, deliberately NOT by (contig, feature_name).
    # featureCounts runs with -g gene_name, so it collapses every locus sharing a
    # symbol into ONE meta-feature — one count, one Length. Splitting the numerator
    # per contig while the denominator stays merged divides one locus' edits by all
    # loci's reads, which is wrong in every resulting row. Grouping by symbol keeps
    # numerator and denominator on the same key.
    normalized_matrix = marine_counts.groupby("feature_name").agg(
        contig=("contig", _join_unique),
        total_edits=("count", "sum"),
        total_coverage=("coverage", "sum"),
        feature_type=("feature_type", "first"),
        feature_strand=("feature_strand", _join_unique),
        strand_conversion=("strand_conversion", _join_unique),
        strand=("strand", _join_unique),
        Length=("Length", "first"),
    ).reset_index()

    # Annotate with feature counts for this sample
    normalized_matrix = normalized_matrix.merge(
        feature_counts, left_on="feature_name", right_on="Geneid", how="left"
    ).drop(columns=["Geneid"])

    # Log and drop missing sample counts
    missing_counts = normalized_matrix[sample_name].isna().sum()
    if missing_counts > 0:
        logger.warning(f"{missing_counts:,} genes have no featureCounts entry — dropping.")
    normalized_matrix = normalized_matrix.dropna(subset=[sample_name])

    # Log and drop zero sample counts to avoid division by zero
    zero_counts = (normalized_matrix[sample_name] == 0).sum()
    
    if zero_counts > 0:
        top_genes = normalized_matrix.loc[normalized_matrix[sample_name] == 0, "feature_name"].value_counts().head(5).to_dict()
        logger.warning(f"{zero_counts:,} genes have zero featureCounts — dropping to avoid division by zero.")
        logger.warning(f"  Most common: {top_genes}")
        
    normalized_matrix = normalized_matrix[normalized_matrix[sample_name] != 0]

    # total_reads = normalized_matrix[sample_name].sum()

    normalized_matrix["EPR"]   = normalized_matrix["total_edits"] / normalized_matrix[sample_name]
    normalized_matrix["EPKM"]  = normalized_matrix["total_edits"] / (total_reads / 1e6 * (normalized_matrix["Length"] / 1e3))
    normalized_matrix["EPKMR"] = normalized_matrix["EPKM"] / normalized_matrix[sample_name]
    normalized_matrix["EPM"]   = normalized_matrix["total_edits"] / (total_reads / 1e6)
    normalized_matrix["EPMR"]  = normalized_matrix["EPM"] / normalized_matrix[sample_name]

    normalized_matrix["sample_name"] = sample_name
    normalized_matrix = normalized_matrix.rename(columns={sample_name: "featureCount_count"})

    logger.info(f"Final output: {len(normalized_matrix):,} rows across "
                f"{normalized_matrix['feature_name'].nunique()} genes")

    normalized_matrix = normalized_matrix[[
        "sample_name", "contig", "feature_name", "feature_strand", "feature_type",
        "Length", "strand", "strand_conversion", "total_edits", "total_coverage",
        "featureCount_count", "EPR", "EPKM", "EPKMR", "EPM", "EPMR",
    ]]
    return normalized_matrix


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(
        description="Normalize MARINE bulk RNA-seq edits to EPR/EPKM/EPKMR metrics."
    )

    parser.add_argument("-i", "--input_file", required=True,
                        help="Path to the filtered MARINE input file", dest="input_file")
    parser.add_argument("-d", "--output_directory", default=script_dir,
                        help=f"Output directory (default: {script_dir})", dest="output_directory")
    parser.add_argument("-e", "--edit_type", default="C>T",
                        help="strand_conversion type to keep (default: 'C>T')", dest="edit_type")
    parser.add_argument("-c", "--feature_counts_matrix", required=True,
                        help="Path to the featureCounts matrix file", dest="feature_counts_matrix")
    parser.add_argument("-s", "--sample_name", required=True,
                        help="Sample name matching a column in the featureCounts matrix", dest="sample_name")
    parser.add_argument("-str", "--strandedness", type=int, choices=[0, 1, 2], required=True,
                        help="Strandedness mode: 0=unstranded, 1=forward, 2=reverse", dest="strandedness")
    parser.add_argument("--no-bedgraph", action="store_true", default=False,
                        help="Skip saving the edit_fraction bedgraph (default: bedgraph is saved)", dest="no_bedgraph")

    args = parser.parse_args()
    os.makedirs(args.output_directory, exist_ok=True)

    edit_type = args.edit_type
    if ">" not in edit_type or len(edit_type) != 3:
        raise ValueError(f"edit_type should look like 'C>T', got: {edit_type}")

    # ---- Load ---------------------------------------------------------------
    logger.info(f"Loading {args.input_file} ...")
    df = pd.read_csv(args.input_file, sep="\t")
    logger.info(f"Loaded {len(df):,} edit entries")

    # ---- Filter by strand_conversion ----------------------------------------
    if args.strandedness == 0:
        keep = {edit_type, complement_edit(edit_type)}
        df = df[df["strand_conversion"].isin(keep)]
        logger.info(f"Filtered to strand_conversion in {keep}: {len(df):,} entries remaining")
    else:
        df = df[df["strand_conversion"] == edit_type]
        logger.info(f"Filtered to strand_conversion == '{edit_type}': {len(df):,} entries remaining")

    # ---- Save bedgraph + BED6 (position-level, before gene-level split) -----
    if not args.no_bedgraph:
        # MARINE positions are 1-based; BED/bedgraph is 0-based half-open: start = position-1, end = position
        bedgraph_dir = os.path.join(args.output_directory, "bedgraphs")
        os.makedirs(bedgraph_dir, exist_ok=True)

        pos_df = df[["contig", "position", "feature_name", "edit_fraction", "strand", "strand_conversion"]].drop_duplicates(subset=["contig", "position"]).copy()
        pos_df["start"] = pos_df["position"] - 1

        # For strandedness=0, MARINE assigns strand="+" to all paired-end reads
        # so the strand column is not meaningful. Derive BED strand from
        # strand_conversion instead — it encodes the edit's genomic strand
        # directly via the reference base (C>T → + strand, G>A → - strand).
        # For stranded libraries (1/2), the MARINE strand column is correct.
        if args.strandedness == 0:
            comp = complement_edit(edit_type)
            pos_df["bed_strand"] = pos_df["strand_conversion"].map({edit_type: "+", comp: "-"})
        else:
            pos_df["bed_strand"] = pos_df["strand"]

        # 4-column bedgraph: chr, start, end, edit_fraction
        bedgraph = pos_df[["contig", "start", "position", "edit_fraction"]].rename(columns={"position": "end"})
        bedgraph = bedgraph.sort_values(["contig", "start"])
        bedgraph_name = f"{args.sample_name}.{edit_type.replace('>', '_')}.edit_fraction.bedgraph"
        bedgraph.to_csv(os.path.join(bedgraph_dir, bedgraph_name), sep="\t", index=False, header=False)
        logger.info(f"Saved bedgraph ({len(bedgraph):,} sites) → bedgraphs/{bedgraph_name}")

        # 6-column BED: chr, start, end, feature_name, edit_fraction, strand
        # Used as input to metaPlotR (strand-aware intersection via helper_metaplotdist.py)
        bed6 = pos_df[["contig", "start", "position", "feature_name", "edit_fraction", "bed_strand"]].rename(columns={"position": "end", "bed_strand": "strand"})
        bed6 = bed6.sort_values(["contig", "start"])
        bed6_name = f"{args.sample_name}.{edit_type.replace('>', '_')}.edit_fraction.bed"
        bed6.to_csv(os.path.join(bedgraph_dir, bed6_name), sep="\t", index=False, header=False)
        logger.info(f"Saved BED6 ({len(bed6):,} sites) → bedgraphs/{bed6_name}")
    else:
        logger.info("Skipping bedgraph/BED6 output (--no-bedgraph set)")

    # ---- Split multi-annotated sites ----------------------------------------
    df = split_multi_annotated_sites(df)

    # ---- For unstranded data: drop edits assigned to the wrong-strand gene --
    # LOJ annotation assigns each site to ALL overlapping genes regardless of
    # strand. For paired-end strandedness=0 MARINE sets strand="+" for every
    # read, so the strand column cannot be used to resolve this. Instead use
    # strand_conversion: C>T edits can only come from + strand RNA (ref=C on +
    # strand), G>A edits can only come from - strand RNA (ref=G on + strand =
    # C on - strand). Keep only rows where the gene strand matches the edit.
    if args.strandedness == 0:
        comp = complement_edit(edit_type)
        before = len(df)
        mask = (
            ((df["strand_conversion"] == edit_type)  & (df["feature_strand"] == "+")) |
            ((df["strand_conversion"] == comp)        & (df["feature_strand"] == "-"))
        )
        df = df[mask]
        logger.info(f"Strandedness=0: strand_conversion→feature_strand filter kept {len(df):,} / {before:,} rows "
                    f"(removed {before - len(df):,} cross-strand assignments)")

    # ---- Extract gene lengths from featureCounts ----------------------------
    logger.info("Extracting gene lengths from featureCounts Length column ...")
    gene_lengths_df = extract_gene_lengths_from_featurecounts(args.feature_counts_matrix)

    # ---- Annotate gene lengths ----------------------------------------------
    logger.info("Annotating edits with gene lengths ...")
    df = annotate_gene_length(df, gene_lengths_df)
    logger.info(f"After gene length annotation: {len(df):,} entries")

    # ---- Save length-annotated intermediate ---------------------------------
    df["sample_name"] = args.sample_name
    file_name = args.sample_name + ".edits.length_annotated.tsv"
    df.to_csv(os.path.join(args.output_directory, file_name), sep="\t", index=False)
    logger.info(f"Saved length-annotated edits → {file_name}")

    # ---- Process feature counts ---------------------------------------------
    logger.info("Processing featureCounts matrix ...")
    feature_counts = process_feature_counts(feature_counts_path=args.feature_counts_matrix)

    # ---- Normalize ----------------------------------------------------------
    logger.info("Normalizing edits ...")
    normalized_marine = normalizing_marine_edits(
        marine_counts=df,
        feature_counts=feature_counts,
        sample_name=args.sample_name,
    )

    file_name = args.sample_name + ".EPR_EPKM_normalized.tsv"
    normalized_marine.to_csv(os.path.join(args.output_directory, file_name), sep="\t", index=False)
    logger.info(f"Saved normalized edits → {file_name}")


if __name__ == "__main__":
    if sys.version_info < (3, 8, 0):
        sys.stderr.write("You need Python 3.8 or later to run this script\n")
        sys.exit(1)
    main()
