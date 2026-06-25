/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SC MARINE SUBWORKFLOW
    Steps (FASTQ-start): CellRanger → MARINE_SC → FILTER_EDITS_SC → NORMALIZE_EDITS_SC
    Steps (BAM-start):              → MARINE_SC → FILTER_EDITS_SC → NORMALIZE_EDITS_SC

    Strandedness is hardcoded to 2 for 10x STAMP data; it is not inferred here.
    matrix_dir is carried as a keyed side-channel (joined by meta.id) so it is
    available for both MARINE_SC (barcodes) and NORMALIZE_EDITS_SC (UMI counts).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CELLRANGER          } from '../../../modules/local/cellranger/main'
include { MARINE_SC           } from '../../../modules/local/marine_sc/main'
include { FILTER_EDITS_SC     } from '../../../modules/local/filter_edits_sc/main'
include { NORMALIZE_EDITS_SC  } from '../../../modules/local/normalize_edits_sc/main'

workflow SC_MARINE {

    take:
    ch_input        // channel: items from samplesheet classification
                    //   FASTQ-start: [ meta{mode:'sc',start:'fastq'}, fastq_dir ]
                    //   BAM-start:   [ meta{mode:'sc',start:'bam'},   bam, matrix_dir ]
    cellranger_ref  // path: pre-built CellRanger reference (FASTQ-start only)
    fasta           // path: reference FASTA (for samtools calmd in MARINE_SC)
    gene_bed        // path: BED6 gene model
    dbsnp_bed       // path: dbSNP BED for SNP filtering

    main:

    def ch_versions      = channel.empty()
    def ch_multiqc_files = channel.empty()

    // ── Split input by start point ────────────────────────────────────────────
    ch_input.branch { entry ->
        fastq: entry[0].start == 'fastq'
        bam:   entry[0].start == 'bam'
    }.set { ch_by_start }

    // ── FASTQ START: CellRanger ───────────────────────────────────────────────
    // ch_by_start.fastq: [ meta, fastq_dir ]
    CELLRANGER(ch_by_start.fastq, cellranger_ref)
    ch_versions      = ch_versions.mix(CELLRANGER.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(CELLRANGER.out.outs_dir.map { _m, d -> d })

    // CellRanger BAM output has the BAI embedded; emit as [ meta, bam, bai, matrix_dir ]
    def ch_bams_from_fastq = CELLRANGER.out.bam
        .join(CELLRANGER.out.matrix_dir)
        .map { meta, bam, bai, matrix_dir -> [ meta, bam, bai, matrix_dir ] }

    // ── BAM START: validate BAI exists ────────────────────────────────────────
    // ch_by_start.bam: [ meta, bam, matrix_dir ] (3-element tuple from samplesheet)
    def ch_bams_from_bam = ch_by_start.bam
        .map { meta, bam, matrix_dir ->
            def bai = file("${bam}.bai")
            if (!bai.exists()) {
                bai = file("${bam.toString().replace('.bam', '.bai')}")
            }
            if (!bai.exists()) {
                error("Sample '${meta.id}' [SC BAM]: BAI index not found. " +
                      "Expected '${bam}.bai' or '${bam.toString().replace('.bam', '.bai')}'.")
            }
            [ meta, bam, bai, matrix_dir ]
        }

    // ── Merge both start points: [ meta, bam, bai, matrix_dir ] ─────────────
    def ch_all_sc_bams = ch_bams_from_fastq.mix(ch_bams_from_bam)

    // Carry matrix_dir separately by key so it is available after MARINE_SC
    def ch_matrix_dirs = ch_all_sc_bams
        .map { meta, _bam, _bai, matrix_dir -> [ meta.id, matrix_dir ] }

    // ── MARINE edit calling ───────────────────────────────────────────────────
    MARINE_SC(
        ch_all_sc_bams.map { meta, bam, bai, matrix_dir -> [ meta, bam, bai, matrix_dir ] },
        fasta,
        gene_bed
    )
    ch_versions = ch_versions.mix(MARINE_SC.out.versions)

    // ── Filter MARINE SC edits ────────────────────────────────────────────────
    FILTER_EDITS_SC(MARINE_SC.out.results, dbsnp_bed)
    ch_versions = ch_versions.mix(FILTER_EDITS_SC.out.versions)

    // ── Normalize by per-cell UMI counts ─────────────────────────────────────
    // Join filtered edits with the matrix_dir channel keyed by meta.id
    def ch_normalize_input = FILTER_EDITS_SC.out.filtered
        .map { meta, filtered -> [ meta.id, meta, filtered ] }
        .join(ch_matrix_dirs)
        .map { _id, meta, filtered, matrix_dir -> [ meta, filtered, matrix_dir ] }

    NORMALIZE_EDITS_SC(ch_normalize_input, gene_bed)
    ch_versions = ch_versions.mix(NORMALIZE_EDITS_SC.out.versions)

    emit:
    normalized   = NORMALIZE_EDITS_SC.out.normalized  // [ meta, normalized_edits.tsv ]
    multiqc      = ch_multiqc_files
    versions     = ch_versions
}
