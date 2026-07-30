/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BULK SAILOR SUBWORKFLOW
    Steps: (collected BAMs) → SAILOR (Snakemake) → METAPLOTR_SAILOR (3 confidence tiers)

    Receives pre-resolved strandedness (int) from the calling workflow.
    When MARINE ran first: caller derives consensus from BULK_MARINE.out.strand_codes.
    When running SAILOR only from BAMs: caller uses params.strandedness directly.

    All samples MUST share the same strandedness and library type (SE or PE).
    SAILOR processes all BAMs as a single Snakemake batch job.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAILOR                  } from '../../../modules/local/sailor/main'
include { METAPLOTR_SAILOR        } from '../../../modules/local/metaplotr_sailor/main'
include { NORMALIZE_EDITS_SAILOR  } from '../../../modules/local/normalize_edits_sailor/main'

workflow BULK_SAILOR {

    take:
    ch_all_bams   // channel: [ meta, bam, bai ] — all bulk BAMs for this run
    strandedness  // val: consensus strandedness integer (0/1/2), pre-resolved by caller
    counts_matrix // path: merged featureCounts matrix from BULK_PREPROCESS
    fasta         // path: reference FASTA (for SAILOR Snakemake rules)
    fasta_fai     // path: FASTA index (.fai) — staged alongside fasta so bam_to_bw.sh finds it
    dbsnp_bed     // path: dbSNP BED for SNP filtering
    snakefile     // path: SAILOR Snakefile
    gene_bed      // path: BED6 gene model, for annotating ranked sites to genes
    genepred      // path: genePred file for metaPlotR

    main:

    def ch_versions = channel.empty()

    // ── Collect all BAMs and BAIs into single lists for the SAILOR batch ─────
    // .collect() produces a value channel (singleton list) — SAILOR runs once.
    def ch_collected_bams = ch_all_bams.map { _m, bam, _bai -> bam }.collect()
    def ch_collected_bais = ch_all_bams.map { _m, _bam, bai -> bai }.collect()

    // ── Derive library type — must be uniform across all samples ─────────────
    def ch_library_type = ch_all_bams
        .map { meta, _bam, _bai -> meta.single_end ? 'single' : 'paired' }
        .collect()
        .map { types ->
            def unique_types = types.unique()
            if (unique_types.size() > 1) {
                error("SAILOR: all samples must share library type (SE or PE), " +
                      "but found mixed types: ${types}. " +
                      "Run SAILOR separately for SE and PE cohorts.")
            }
            unique_types[0]
        }

    // ── Run SAILOR ───────────────────────────────────────────────────────────
    // All inputs are value channels → SAILOR executes once for the full batch.
    SAILOR(
        ch_collected_bams,
        ch_collected_bais,
        strandedness,
        ch_library_type,
        fasta,
        fasta_fai,
        dbsnp_bed,
        snakefile
    )
    ch_versions = ch_versions.mix(SAILOR.out.versions)

    // ── Flatten per-sample ranked BED files for metaPlotR ───────────────────
    // SAILOR emits all BEDs as a list; flatten emits them one at a time.
    // Sample ID is reconstructed from the filename pattern:
    //   {id}.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed
    def ch_ranked_beds = SAILOR.out.ranked_beds
        .flatten()
        .map { bed ->
            def sample_id = bed.name.replaceAll(/\.combined\..*$/, '')
            [ [id: sample_id], bed ]
        }

    // strandedness is a val (broadcast to each queue channel item automatically)
    METAPLOTR_SAILOR(ch_ranked_beds, strandedness, genepred)
    ch_versions = ch_versions.mix(METAPLOTR_SAILOR.out.versions)

    // ── Annotate ranked sites to genes and normalise to expression ───────────
    // ch_ranked_beds carries SAILOR's own sample name (the BAM filename), which
    // BULK_FLARE relies on to locate matching bigwigs — so it must stay as-is.
    // featureCounts columns use the samplesheet id instead, and the two only
    // coincide for FASTQ-start samples: SAILOR renames STAR BAMs to '{id}.bam',
    // while BAM-start samples keep their own filename (sample 'exp' may arrive as
    // 'exp_chr1.bam'). Recover the real meta by joining on that filename, applying
    // the same rename rule as modules/local/sailor/main.nf.
    def ch_bam_key = ch_all_bams.map { meta, bam, _bai ->
        def sailor_name = bam.name.endsWith('.Aligned.sortedByCoord.out.bam')
            ? bam.name.replace('.Aligned.sortedByCoord.out.bam', '.bam')
            : bam.name
        [ sailor_name, meta ]
    }

    def ch_normalize_input = ch_ranked_beds
        .map { meta, bed -> [ meta.id, bed ] }
        .join(ch_bam_key)
        .map { _sailor_name, bed, meta -> [ meta, bed ] }
        .combine(counts_matrix)

    NORMALIZE_EDITS_SAILOR(ch_normalize_input, strandedness, gene_bed)
    ch_versions = ch_versions.mix(NORMALIZE_EDITS_SAILOR.out.versions)

    emit:
    ranked_beds = ch_ranked_beds               // [ meta, bed ] per sample
    output_dir  = SAILOR.out.output_dir        // sailor_output/ folder (for downstream FLARE)
    versions    = ch_versions
}
