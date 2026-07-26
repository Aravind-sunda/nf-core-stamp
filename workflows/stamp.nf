/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BULK_MARINE            } from '../subworkflows/local/bulk_marine/main'
include { BULK_SAILOR            } from '../subworkflows/local/bulk_sailor/main'
include { BULK_FLARE             } from '../subworkflows/local/bulk_flare/main'
include { SC_MARINE              } from '../subworkflows/local/sc_marine/main'
include { FLARE_GENERATE_REGIONS } from '../modules/local/flare_generate_regions/main'
include { PREPARE_DBSNP          } from '../modules/local/prepare_dbsnp/main'
include { SAMTOOLS_FAIDX         } from '../modules/local/samtools_faidx/main'
include { GUNZIP as GUNZIP_GTF      } from '../modules/local/gunzip/main'
include { GUNZIP as GUNZIP_GENE_BED } from '../modules/local/gunzip/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_stamp_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENEPRED ASSET LOOKUP
    Built-in genePred files for supported genomes (in assets/genepred/).
    Override with --genepred to supply your own file; select built-in with --genome.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def resolveGenePred(projectDir) {
    def assets = [
        hg19:     "${projectDir}/assets/genepred/hg19_ncbiRefSeqCurated.txt.gz",
        hg38:     "${projectDir}/assets/genepred/hg38_ncbiRefSeqCurated.txt.gz",
        hg38_V44: "${projectDir}/assets/genepred/hg38_V44_E110_basic_knownGene_genePred.txt.gz",
        mm10:     "${projectDir}/assets/genepred/mm10_ncbiRefSeqCurated.txt.gz",
        mm39:     "${projectDir}/assets/genepred/mm39_ncbiRefSeqCurated.txt.gz",
    ]
    if (params.genepred) {
        return params.genepred
    } else if (params.genome && assets.containsKey(params.genome)) {
        return assets[params.genome]
    }
    return null
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STAMP {

    take:
    ch_samplesheet               // channel: typed tuples from classifyAndValidateRow()
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions      = channel.empty()
    def ch_multiqc_files = channel.empty()

    // ── Reference file channels (validated in validateInputParameters()) ──────
    // channel.value(file()) creates a value channel directly — channel.fromPath().first()
    // triggers "first is useless on a value channel" in Nextflow >= 24.x.
    def ch_fasta = params.fasta
        ? channel.value(file(params.fasta, checkIfExists: true))
        : channel.empty()

    // Use pre-existing .fai if present; otherwise generate it with samtools faidx.
    def ch_fasta_fai = channel.empty()
    if (params.fasta) {
        def fai_path = file(params.fasta + '.fai')
        if (fai_path.exists()) {
            ch_fasta_fai = channel.value(fai_path)
        } else {
            SAMTOOLS_FAIDX(ch_fasta)
            ch_fasta_fai = SAMTOOLS_FAIDX.out.fai
            ch_versions  = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        }
    }

    def ch_star_index = params.star_index
        ? channel.value(file(params.star_index, checkIfExists: true))
        : channel.empty()

    // GTF and BED: decompress on-the-fly if .gz so that tools that don't
    // support gzipped inputs (RSeQC infer_experiment.py, MARINE, etc.) receive
    // a plain file. .first() converts the single-element process output queue
    // channel back into a value channel for broadcasting to all consumers.
    def ch_gtf = channel.empty()
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            GUNZIP_GTF(channel.value(file(params.gtf, checkIfExists: true)))
            ch_gtf      = GUNZIP_GTF.out.file
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = channel.value(file(params.gtf, checkIfExists: true))
        }
    }

    def ch_gene_bed = channel.empty()
    if (params.gene_bed) {
        if (params.gene_bed.endsWith('.gz')) {
            GUNZIP_GENE_BED(channel.value(file(params.gene_bed, checkIfExists: true)))
            ch_gene_bed = GUNZIP_GENE_BED.out.file
            ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_gene_bed = channel.value(file(params.gene_bed, checkIfExists: true))
        }
    }

    // Sorted once here rather than per sample: every consumer (bulk/sc filters and
    // SAILOR) shares the same prepared file. .first() converts the process output
    // back into a value channel so it can broadcast to all of them.
    def ch_dbsnp_bed = channel.empty()
    if (params.dbsnp_bed) {
        PREPARE_DBSNP(channel.value(file(params.dbsnp_bed, checkIfExists: true)))
        ch_dbsnp_bed = PREPARE_DBSNP.out.bed.first()
        ch_versions  = ch_versions.mix(PREPARE_DBSNP.out.versions)
    }

    def ch_snakefile = params.sailor_snakefile
        ? channel.value(file(params.sailor_snakefile, checkIfExists: true))
        : channel.empty()

    def ch_flare_snakefile = params.flare_snakefile
        ? channel.value(file(params.flare_snakefile, checkIfExists: true))
        : channel.empty()

    def ch_cellranger_ref = params.cellranger_ref
        ? channel.value(file(params.cellranger_ref, checkIfExists: true))
        : channel.empty()

    // ── GenePred (required for metaPlotR in bulk mode) ───────────────────────
    def genepred_path = resolveGenePred(projectDir)
    def ch_genepred = genepred_path
        ? channel.value(file(genepred_path, checkIfExists: true))
        : channel.empty()

    // ── Split samplesheet by pipeline mode ───────────────────────────────────
    ch_samplesheet.branch { entry ->
        bulk: entry[0].mode == 'bulk'
        sc:   entry[0].mode == 'sc'
    }.set { ch_by_mode }

    // ════════════════════════════════════════════════════════════════════════════
    //  BULK MODE
    // ════════════════════════════════════════════════════════════════════════════

    if (params.mode == 'bulk') {

        if (!params.run_marine && !params.run_sailor) {
            log.warn("Both --run_marine and --run_sailor are false in bulk mode. Nothing to run.")
        }

        // ── MARINE: FASTQ/BAM → edits → expression → filter → normalize → metaPlotR
        if (params.run_marine) {
            BULK_MARINE(
                ch_by_mode.bulk,
                ch_star_index,
                ch_gtf,
                ch_gene_bed,
                ch_dbsnp_bed,
                ch_genepred,
                ch_fasta
            )
            ch_versions      = ch_versions.mix(BULK_MARINE.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(BULK_MARINE.out.multiqc)
        }

        // ── SAILOR: all BAMs → Snakemake → ranked BEDs → metaPlotR (3 tiers)
        if (params.run_sailor) {

            // Which BAMs to feed SAILOR depends on whether MARINE ran first
            def ch_sailor_bams
            if (params.run_marine) {
                // BULK_MARINE already aligned / validated all BAMs
                ch_sailor_bams = BULK_MARINE.out.bams
            } else {
                // SAILOR-only: samples must be BAM-start; FASTQ-start requires MARINE first
                ch_by_mode.bulk
                    .filter { meta, _d -> meta.start == 'fastq' }
                    .map { meta, _d ->
                        error("Sample '${meta.id}' is FASTQ-start but --run_marine is false. " +
                              "Cannot run SAILOR without alignment. Enable --run_marine or provide pre-aligned BAMs.")
                    }
                ch_sailor_bams = ch_by_mode.bulk
                    .filter { meta, _d -> meta.start == 'bam' }
                    .map { meta, bam ->
                        def bai = file("${bam}.bai")
                        if (!bai.exists()) bai = file("${bam.toString().replace('.bam', '.bai')}")
                        if (!bai.exists()) {
                            error("Sample '${meta.id}': BAI index not found. " +
                                  "Expected '${bam}.bai' or '${bam.toString().replace('.bam', '.bai')}'.")
                        }
                        [ meta, bam, bai ]
                    }
            }

            // Resolve strandedness: inferred by MARINE or provided by --strandedness
            def ch_sailor_strand
            if (params.run_marine) {
                // Derive consensus from per-sample inferred strand codes
                ch_sailor_strand = BULK_MARINE.out.strand_codes
                    .map { _id, strand -> strand }
                    .collect()
                    .map { strands ->
                        def orig = strands.collect()
                        def unique_strands = orig.unique()
                        if (unique_strands.size() > 1) {
                            def majority = orig.countBy { it }.max { it.value }.key
                            log.warn("SAILOR: samples have mixed strandedness ${orig}. Using most common: ${majority}")
                            majority
                        } else {
                            unique_strands[0]
                        }
                    }
            } else {
                if (params.strandedness == null) {
                    error("--run_sailor without --run_marine requires --strandedness (0, 1, or 2).")
                }
                ch_sailor_strand = channel.value(params.strandedness.toInteger())
            }

            BULK_SAILOR(
                ch_sailor_bams,
                ch_sailor_strand,
                ch_fasta,
                ch_fasta_fai,
                ch_dbsnp_bed,
                ch_snakefile,
                ch_genepred
            )
            ch_versions = ch_versions.mix(BULK_SAILOR.out.versions)

            // ── FLARE: SAILOR outputs → edit-cluster identification (RBP-STAMP)
            // Requires SAILOR (validated in validateInputParameters), so it is
            // nested here where BULK_SAILOR.out is guaranteed to exist.
            if (params.run_flare) {

                // Regions folder: reuse a pre-built one, or generate it from the GTF.
                def ch_flare_regions
                if (params.flare_regions) {
                    ch_flare_regions = channel.value(file(params.flare_regions, checkIfExists: true))
                } else {
                    def ch_flare_scripts = channel.value(
                        file("${projectDir}/assets/workflow_FLARE/scripts", checkIfExists: true)
                    )
                    FLARE_GENERATE_REGIONS(ch_gtf, params.flare_window_size, ch_flare_scripts)
                    ch_versions = ch_versions.mix(FLARE_GENERATE_REGIONS.out.versions)
                    // .first() → value channel so the regions folder broadcasts to every sample
                    ch_flare_regions = FLARE_GENERATE_REGIONS.out.regions.first()
                }

                BULK_FLARE(
                    BULK_SAILOR.out.ranked_beds,
                    BULK_SAILOR.out.output_dir,
                    ch_flare_regions,
                    ch_fasta,
                    ch_fasta_fai,
                    ch_flare_snakefile
                )
                ch_versions = ch_versions.mix(BULK_FLARE.out.versions)
            }
        }
    }

    // ════════════════════════════════════════════════════════════════════════════
    //  SC MODE
    // ════════════════════════════════════════════════════════════════════════════

    if (params.mode == 'sc') {

        SC_MARINE(
            ch_by_mode.sc,
            ch_cellranger_ref,
            ch_fasta,
            ch_gene_bed,
            ch_dbsnp_bed
        )
        ch_versions      = ch_versions.mix(SC_MARINE.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(SC_MARINE.out.multiqc)
    }

    // ── Software versions ─────────────────────────────────────────────────────
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file:  entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'stamp_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    // ── MultiQC ───────────────────────────────────────────────────────────────
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )

    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true)
    )

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'stamp'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList()
    versions       = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
