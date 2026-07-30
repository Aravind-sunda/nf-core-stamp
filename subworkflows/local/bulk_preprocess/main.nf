/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BULK PREPROCESSING SUBWORKFLOW
    Steps: FastQC → fastp → STAR → infer_strandedness → featureCounts → merge_counts
    Supports both FASTQ-start and BAM-start samples (detected via meta.start).

    Shared by BULK_MARINE and BULK_SAILOR so that alignment, strandedness and gene
    quantification happen exactly once regardless of which edit callers are enabled.
    Keeping featureCounts here (rather than inside BULK_MARINE) is what lets SAILOR
    normalise against expression without requiring --run_marine.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC             } from '../../../modules/nf-core/fastqc/main'
include { FASTP              } from '../../../modules/local/fastp/main'
include { STAR_ALIGN         } from '../../../modules/local/star_align/main'
include { INFER_STRANDEDNESS } from '../../../modules/local/infer_strandedness/main'
include { FEATURECOUNTS      } from '../../../modules/local/featurecounts/main'
include { MERGE_COUNTS       } from '../../../modules/local/merge_counts/main'

workflow BULK_PREPROCESS {

    take:
    ch_input     // channel: [ meta, reads/bam ] from samplesheet validation
                 //   FASTQ-start: meta.start == 'fastq', data = [ fastq_1 ] or [ fastq_1, fastq_2 ]
                 //   BAM-start:   meta.start == 'bam',   data = bam_file
    star_index   // path: STAR genome index directory
    gtf          // path: GTF annotation file
    gene_bed     // path: BED6 gene model file

    main:

    def ch_versions      = channel.empty()
    def ch_multiqc_files = channel.empty()

    // ── Split input by start point ────────────────────────────────────────────
    ch_input.branch { entry ->
        fastq: entry[0].start == 'fastq'
        bam:   entry[0].start == 'bam'
    }.set { ch_by_start }

    // ── FASTQ START: QC → trim → align ───────────────────────────────────────
    FASTQC(ch_by_start.fastq)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { _m, f -> f })
    // FASTQC versions are collected via channel.topic('versions') in the main workflow

    FASTP(ch_by_start.fastq)
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json)
    ch_versions      = ch_versions.mix(FASTP.out.versions)

    STAR_ALIGN(FASTP.out.reads, star_index, gtf)
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final)
    ch_versions      = ch_versions.mix(STAR_ALIGN.out.versions)

    // Merge BAMs from FASTQ path and direct BAM input into one channel.
    // BAM-start samples need BAI files staged alongside — check and index if missing.
    def ch_bams_from_fastq = STAR_ALIGN.out.bam
        .join(STAR_ALIGN.out.bai)
        .map { meta, bam, bai -> [ meta, bam, bai ] }

    // For BAM-start: bai is expected to exist beside the bam; validate at runtime
    def ch_bams_from_bam = ch_by_start.bam
        .map { meta, bam ->
            def bai = file("${bam}.bai")
            if (!bai.exists()) {
                bai = file("${bam.toString().replace('.bam', '.bai')}")
            }
            if (!bai.exists()) {
                error("Sample '${meta.id}' [bulk BAM]: BAI index not found. " +
                      "Expected '${bam}.bai' or '${bam.toString().replace('.bam', '.bai')}'.")
            }
            [ meta, bam, bai ]
        }

    def ch_all_bams = ch_bams_from_fastq.mix(ch_bams_from_bam)

    // ── Infer strandedness per sample ─────────────────────────────────────────
    INFER_STRANDEDNESS(ch_all_bams, gene_bed)
    ch_versions      = ch_versions.mix(INFER_STRANDEDNESS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(INFER_STRANDEDNESS.out.strandedness_txt.map { _m, f -> f })

    // Parse strand_code_txt to an integer value and join back to BAMs
    def ch_strand_codes = INFER_STRANDEDNESS.out.strand_code_txt
        .map { meta, txt ->
            def raw = txt.text.trim()
            if (!raw.isInteger()) {
                error("INFER_STRANDEDNESS produced unexpected output for sample '${meta.id}': '${raw}'. Expected 0, 1, or 2.")
            }
            [ meta.id, raw.toInteger() ]
        }

    def ch_bams_with_strand = ch_all_bams
        .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
        .join(ch_strand_codes)
        .map { _id, meta, bam, bai, strand -> [ meta, bam, bai, strand ] }

    // ── featureCounts (gene expression) ──────────────────────────────────────
    // Needed by both edit callers: MARINE and SAILOR normalise edits to expression
    // using the same EPR/EPKM/EPM metrics, so the matrix is built once here.
    def ch_fc_input = ch_bams_with_strand
        .map { meta, bam, _bai, strand -> [ meta, bam, strand ] }

    FEATURECOUNTS(ch_fc_input, gtf)
    ch_multiqc_files = ch_multiqc_files.mix(FEATURECOUNTS.out.summary.map { _m, f -> f })
    ch_versions      = ch_versions.mix(FEATURECOUNTS.out.versions)

    // ── Merge all per-sample counts into one matrix ───────────────────────────
    MERGE_COUNTS(FEATURECOUNTS.out.counts.map { _m, f -> f }.collect())
    ch_versions = ch_versions.mix(MERGE_COUNTS.out.versions)

    emit:
    bams             = ch_all_bams              // [ meta, bam, bai ]
    strand_codes     = ch_strand_codes          // [ sample_id, strand_int ]
    bams_with_strand = ch_bams_with_strand      // [ meta, bam, bai, strand ]
    counts_matrix    = MERGE_COUNTS.out.matrix  // counts_matrix_combined.tsv (value channel)
    multiqc          = ch_multiqc_files
    versions         = ch_versions
}
