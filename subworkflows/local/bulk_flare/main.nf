/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BULK FLARE SUBWORKFLOW (RBP-STAMP)
    Steps: (per-sample SAILOR outputs) + regions folder → FLARE (Snakemake, per sample)

    FLARE runs cluster-identification mode downstream of SAILOR. For each sample it
    pairs the SAILOR ranked BED (STAMP sites) with that sample's bigwigs and
    filtered-merged BAM (resolved from the SAILOR output folder) plus the shared
    genome-level regions folder.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FLARE } from '../../../modules/local/flare/main'

workflow BULK_FLARE {

    take:
    ch_ranked_beds  // channel: [ meta, ranked_bed ] per sample (from BULK_SAILOR)
    ch_sailor_dir   // channel: sailor_output/ folder (single value — SAILOR runs as one batch)
    regions         // path:    FLARE regions folder (cluster-identification mode)
    fasta           // path:    reference FASTA
    fasta_fai       // path:    FASTA index (.fai)
    snakefile       // path:    FLARE Snakefile

    main:

    def ch_versions = channel.empty()

    // Pair each sample's ranked BED with its bigwigs/BAM, resolved from the
    // SAILOR output folder's 8_bw_and_bam/ subdirectory by sample ID.
    def ch_flare_inputs = ch_ranked_beds
        .combine(ch_sailor_dir)
        .map { meta, bed, sdir ->
            def bw_dir = "${sdir}/8_bw_and_bam"
            def fwd_bw = file("${bw_dir}/${meta.id}.fwd.sorted.bw")
            def rev_bw = file("${bw_dir}/${meta.id}.rev.sorted.bw")
            def bam    = file("${bw_dir}/${meta.id}_filtered_merged.sorted.bam")
            def bai    = file("${bw_dir}/${meta.id}_filtered_merged.sorted.bam.bai")
            [ meta, bed, fwd_bw, rev_bw, bam, bai ]
        }

    FLARE(ch_flare_inputs, fasta, fasta_fai, regions, snakefile)
    ch_versions = ch_versions.mix(FLARE.out.versions)

    emit:
    scored_peaks = FLARE.out.scored_peaks   // [ meta, scored.tsv ] per sample
    versions     = ch_versions
}
