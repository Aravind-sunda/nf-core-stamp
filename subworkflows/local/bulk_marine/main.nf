/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BULK MARINE SUBWORKFLOW
    Steps: MARINE_BULK → filter_edits → normalize_edits → merge_normalized → metaPlotR

    Alignment, strandedness inference and gene quantification happen upstream in
    BULK_PREPROCESS, which is shared with BULK_SAILOR.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MARINE_BULK           } from '../../../modules/local/marine_bulk/main'
include { FILTER_EDITS_BULK     } from '../../../modules/local/filter_edits_bulk/main'
include { NORMALIZE_EDITS_BULK  } from '../../../modules/local/normalize_edits_bulk/main'
include { MERGE_NORMALIZED_BULK } from '../../../modules/local/merge_normalized_bulk/main'
include { METAPLOTR_BULK        } from '../../../modules/local/metaplotr_bulk/main'

workflow BULK_MARINE {

    take:
    ch_bams_with_strand  // channel: [ meta, bam, bai, strandedness ] from BULK_PREPROCESS
    counts_matrix        // path:    merged featureCounts matrix from BULK_PREPROCESS
    gene_bed             // path:    BED6 gene model file
    dbsnp_bed            // path:    sorted dbSNP BED file
    genepred             // path:    genePred file for metaPlotR
    fasta                // path:    reference FASTA (for samtools calmd when BAMs lack MD tags)

    main:

    def ch_versions = channel.empty()

    // ── MARINE edit calling ───────────────────────────────────────────────────
    MARINE_BULK(ch_bams_with_strand, fasta, gene_bed)
    ch_versions = ch_versions.mix(MARINE_BULK.out.versions)

    // ── Filter MARINE edits ───────────────────────────────────────────────────
    // Join MARINE results with their strand code
    def ch_filter_input = MARINE_BULK.out.results
        .map { meta, strand, dir -> [ meta, strand, dir ] }

    FILTER_EDITS_BULK(ch_filter_input, dbsnp_bed, gene_bed)
    ch_versions = ch_versions.mix(FILTER_EDITS_BULK.out.versions)

    // ── Normalize edits by gene expression ───────────────────────────────────
    // combine() broadcasts the single shared matrix to every sample channel item,
    // producing [ meta, strandedness, filtered_edits, counts_matrix ] per sample.
    def ch_normalize_input = FILTER_EDITS_BULK.out.filtered
        .combine(counts_matrix)

    NORMALIZE_EDITS_BULK(ch_normalize_input)
    ch_versions = ch_versions.mix(NORMALIZE_EDITS_BULK.out.versions)

    // ── Combine per-sample normalized tables into gene x sample matrices ──────
    MERGE_NORMALIZED_BULK(NORMALIZE_EDITS_BULK.out.normalized.map { _m, f -> f }.collect())
    ch_versions = ch_versions.mix(MERGE_NORMALIZED_BULK.out.versions)

    // ── metaPlotR metagene distances ──────────────────────────────────────────
    METAPLOTR_BULK(NORMALIZE_EDITS_BULK.out.bedgraph_dir, genepred)
    ch_versions = ch_versions.mix(METAPLOTR_BULK.out.versions)

    emit:
    versions = ch_versions
}
