#!/usr/bin/env python3
import pandas as pd
import scanpy as sc
import numpy as np
import argparse
import logging
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)


def process_counts_matrix(counts_matrix_path):
    """
    Reads a 10x mtx directory and returns a counts DataFrame with total reads per cell.
    """
    adata = sc.read_10x_mtx(counts_matrix_path, var_names='gene_symbols')
    count_matrix = adata.to_df()
    count_matrix['total_reads'] = count_matrix.sum(axis=1)
    return count_matrix

def calculate_length_from_bed(bed):
    """
    Reads a BED file, calculates gene lengths, and returns a DataFrame with feature_name and length.
    """
    bed_df = pd.read_csv(bed, sep="\t", header=None)
    bed_df.columns = ['contig', 'start', 'end', 'feature_name', 'region', 'strand']
    bed_df['length'] = bed_df['end'] - bed_df['start']
    return bed_df[['feature_name', 'length']]

def annotate_edits_with_gene_length(edits, bed):
    """
    Merges the edits DataFrame with gene lengths based on the feature_name.
    """
    edits_annotated = pd.merge(edits, bed, on='feature_name', how='left')
    return edits_annotated

def annotate_and_normalize_edits(filtered_edits, counts_matrix):
    """
    Annotates the filtered edits with total cell counts and computes:
      - EPR: total_edits / cellranger_count
      - EPKM: total_edits normalized by total_reads and gene length
      - EPKMR: EPKM normalized by cellranger_count
      - EPM: total_edits normalized by total_reads (in millions)
      - EPMR: EPM normalized by cellranger_count
    """
    # Group the filtered edits by barcode, contig, and feature_name
    grouped_filtered_edits = filtered_edits.groupby(['barcode', 'contig', 'feature_name']).agg(
        total_edits=('count', 'sum'),
        total_coverage=('coverage', 'sum'),
        feature_type=('feature_type', 'first'),
        length=('length', 'first')
    ).reset_index()

    # Calculate the total counts for each cell by summing over all gene columns
    # double counting since this is done in the process_counts_matrix function
    # counts_matrix['total_reads'] = counts_matrix.sum(axis=1)

    # Map total_reads to the grouped edits
    grouped_filtered_edits['total_reads'] = grouped_filtered_edits['barcode'].map(counts_matrix['total_reads'])

    # Fetch the count for the specific gene in each cell (barcode)
    grouped_filtered_edits['cellranger_count'] = grouped_filtered_edits.apply(
        lambda row: counts_matrix.loc[row['barcode'], row['feature_name']]
        if (row['barcode'] in counts_matrix.index and row['feature_name'] in counts_matrix.columns)
        else None,
        axis=1
    )
    # before dropna on cellranger_count
    missing = grouped_filtered_edits[grouped_filtered_edits['cellranger_count'].isna()]
    logger.info(f"Missing cellranger_count: {len(missing)} rows removed")
    logger.info(f"  Sample barcodes: {missing['barcode'].unique()[:5].tolist()}")
    logger.info(f"  Sample genes: {missing['feature_name'].unique()[:5].tolist()}")
    
    # Remove rows with missing cellranger_count
    grouped_filtered_edits = grouped_filtered_edits.dropna(subset=['cellranger_count'])

    # Log zero-count rows BEFORE removing them so the count is accurate
    zero = grouped_filtered_edits[grouped_filtered_edits['cellranger_count'] == 0]
    logger.info(f"Zero cellranger_count: {len(zero)} rows removed")
    logger.info(f"  Most common genes: {zero['feature_name'].value_counts().head(5).to_dict()}")
    grouped_filtered_edits = grouped_filtered_edits[grouped_filtered_edits['cellranger_count'] != 0]

    grouped_filtered_edits = grouped_filtered_edits.dropna(subset=['total_reads'])
    grouped_filtered_edits = grouped_filtered_edits[grouped_filtered_edits['total_reads'] != 0]

    no_length = grouped_filtered_edits[grouped_filtered_edits['length'].isna()]
    logger.info(f"Missing gene length (not in BED): {len(no_length)} rows removed")
    logger.info(f"  Genes not in BED: {no_length['feature_name'].unique()[:10].tolist()}")
    grouped_filtered_edits = grouped_filtered_edits.dropna(subset=['length'])

    # Calculate EPR as the ratio of total_edits to cellranger_count
    grouped_filtered_edits['EPR'] = (
        grouped_filtered_edits['total_edits'].astype(int) /
        grouped_filtered_edits['cellranger_count'].astype(int)
    )

    # Calculate EPKM using total_reads (in millions) and length (in kilobases)
    grouped_filtered_edits['EPKM'] = (
        grouped_filtered_edits['total_edits'].astype(int) /
        ((grouped_filtered_edits['total_reads'].astype(int) / 10**6) *
         (grouped_filtered_edits['length'].astype(int) / 10**3))
    )

    # Calculate EPKMR as EPKM normalized by cellranger_count
    grouped_filtered_edits['EPKMR'] = (
        grouped_filtered_edits['EPKM'].astype(float) /
        grouped_filtered_edits['cellranger_count'].astype(int)
    )

    # Calculate EPM: total_edits normalized by cell total reads (in millions)
    grouped_filtered_edits['EPM'] = (
        grouped_filtered_edits['total_edits'].astype(int) /
        (grouped_filtered_edits['total_reads'].astype(int) / 10**6)
    )

    # Calculate EPMR: EPM normalized by cellranger_count
    grouped_filtered_edits['EPMR'] = (
        grouped_filtered_edits['EPM'] /
        grouped_filtered_edits['cellranger_count'].astype(int)
    )
    logger.info(f"Final output: {len(grouped_filtered_edits)} rows across "
            f"{grouped_filtered_edits['barcode'].nunique()} cells and "
            f"{grouped_filtered_edits['feature_name'].nunique()} genes")

    return grouped_filtered_edits

def expand_multi_annotations(df):
    """Expand rows with comma-separated multi-gene annotations into one row per gene."""
    for col in ("feature_name", "feature_type", "feature_strand"):
        df[col] = df[col].astype(str).str.split(",")
    df = df.explode(["feature_name", "feature_type", "feature_strand"])
    return df.reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(
        description="Annotate and normalize filtered edits using a counts matrix and BED file."
    )
    parser.add_argument("--filtered_edits", required=True,
                        help="Path to the filtered_edits TSV file")
    parser.add_argument("--counts_matrix", required=True,
                        help="Path to the 10x mtx directory for the counts matrix")
    parser.add_argument("--bed", required=True,
                        help="Path to the BED file for gene lengths")
    parser.add_argument("--output_dir", required=True,
                        help="Path to save the annotated edits TSV file")
    parser.add_argument("--edit-type", default="C>T",
                        help="Strand conversion type to retain (e.g. 'C>T' for A-to-I on reverse strand, "
                             "'A>G' for sense strand). Default: C>T")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    # Read and filter the edits file
    filtered_edits = pd.read_csv(args.filtered_edits, sep="\t")

    n_before = len(filtered_edits)
    filtered_edits = filtered_edits.query("strand_conversion == @args.edit_type").copy()
    logger.info(f"Strand filter: {n_before - len(filtered_edits)} rows removed, {len(filtered_edits)} remaining")

    n_before = len(filtered_edits)
    filtered_edits = expand_multi_annotations(filtered_edits)
    logger.info(f"Multi-annotation expansion: {n_before} → {len(filtered_edits)} rows")

    # Process the counts matrix
    counts_matrix = process_counts_matrix(args.counts_matrix)

    # Calculate gene lengths from the BED file and annotate the edits
    bed_length = calculate_length_from_bed(args.bed)
    filtered_edits = annotate_edits_with_gene_length(filtered_edits, bed_length)

    # Annotate and normalize the edits
    normalized_edits = annotate_and_normalize_edits(filtered_edits, counts_matrix)

    # Save the annotated edits to a file
    normalized_edits.to_csv(args.output_dir + "/normalized_edits.tsv", sep="\t", index=False)
    print(f"Normalized edits saved to {args.output_dir}")

if __name__ == '__main__':
    main()
