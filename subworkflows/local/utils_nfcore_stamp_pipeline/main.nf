//
// Subworkflow with functionality specific to the nf-core/stamp pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
// include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/stamp ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/stamp/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --mode <bulk|sc> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input samplesheet
    // Columns (in schema order): sample(meta), fastq_1, fastq_2, bam, library_type, fastq_dir, matrix_dir
    // Mode detection uses ALL columns — any unexpected combination raises an explicit error.
    //
    channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map { meta, fastq_1, fastq_2, bam, library_type, fastq_dir, matrix_dir ->
            classifyAndValidateRow(meta, fastq_1, fastq_2, bam, library_type, fastq_dir, matrix_dir)
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    def mode = params.mode

    if (!mode || !['bulk', 'sc'].contains(mode)) {
        error("--mode must be 'bulk' or 'sc', got: '${mode}'")
    }

    if (mode == 'bulk') {
        if (!params.run_marine && !params.run_sailor) {
            error("--mode bulk: at least one of --run_marine or --run_sailor must be true.")
        }
        // No --strandedness requirement for SAILOR-only runs: BULK_PREPROCESS runs
        // INFER_STRANDEDNESS for every bulk sample regardless of which callers are
        // enabled, and the workflow falls back to that consensus. --strandedness
        // remains available to override it.
        if (params.run_sailor) {
            if (!params.sailor_snakefile) {
                error("--run_sailor true requires --sailor_snakefile (path to SAILOR Snakefile).")
            }
        }
        // MARINE needs the FASTA for its calmd step and SAILOR for its Snakemake
        // rules, so it is required for any bulk run. Both modules declare it as a
        // mandatory input, so a missing value leaves an empty channel and the
        // process is skipped silently rather than failing.
        if (!params.fasta) {
            error("--mode bulk requires --fasta (MARINE calmd step / SAILOR Snakemake rules).")
        }
        if (params.run_flare) {
            // FLARE (RBP-STAMP) operates on SAILOR outputs, so SAILOR must run.
            if (!params.run_sailor) {
                error("--run_flare true requires --run_sailor true (FLARE identifies clusters from SAILOR edit calls).")
            }
            if (!params.flare_snakefile) {
                error("--run_flare true requires --flare_snakefile (path to FLARE Snakefile).")
            }
            // Regions are either supplied pre-built or generated from the GTF.
            if (!params.flare_regions && !params.gtf) {
                error("--run_flare true requires either --flare_regions (pre-built regions folder) or --gtf (to generate one).")
            }
        }
        if (!params.gene_bed) {
            error("--mode bulk requires --gene_bed (BED6 gene model file).")
        }
        if (!params.dbsnp_bed) {
            error("--mode bulk requires --dbsnp_bed (dbSNP BED file).")
        }
        // These were previously gated on --run_marine, from when alignment and
        // quantification lived inside BULK_MARINE. BULK_PREPROCESS now runs them for
        // every bulk run, and both callers end in a metaPlotR step, so gating them
        // let a SAILOR-only run start with an empty channel and skip featureCounts
        // or metaPlotR without reporting anything.
        if (!params.gtf) {
            error("--mode bulk requires --gtf (STAR and featureCounts run for every bulk run).")
        }
        if (!params.genome && !params.genepred) {
            error("--mode bulk requires either --genome (to select a built-in genePred) or --genepred (custom genePred file) for metaPlotR.")
        }
    }

    if (mode == 'sc') {
        if (params.run_flare) {
            error("--run_flare is only supported in --mode bulk (RBP-STAMP is bulk RNA-seq).")
        }
        if (!params.fasta) {
            error("--mode sc requires --fasta (genome FASTA file).")
        }
        if (!params.gene_bed) {
            error("--mode sc requires --gene_bed (BED6 gene model file for MARINE annotation).")
        }
        if (!params.dbsnp_bed) {
            error("--mode sc requires --dbsnp_bed (dbSNP BED file).")
        }
    }
}

//
// Classify a single samplesheet row and validate ALL column combinations explicitly.
// Returns a tuple whose structure depends on the detected mode:
//   bulk + FASTQ  : [ meta, reads ]         reads = [ fastq_1 ] or [ fastq_1, fastq_2 ]
//   bulk + BAM    : [ meta, bam ]
//   sc   + FASTQ  : [ meta, fastq_dir ]
//   sc   + BAM    : [ meta, bam, matrix_dir ]
//
// meta always contains:
//   id         : sample name
//   mode       : 'bulk' | 'sc'
//   start      : 'fastq' | 'bam'
//   single_end : true | false  (bulk only; null for sc)
//
// Resolve a path string relative to the samplesheet file's parent directory.
// Absolute paths and URI schemes (s3://, http://, gs://) pass through unchanged.
// This allows relative paths in samplesheets to work regardless of launchDir.
def resolveFromSamplesheet(String p) {
    if (!p || p == 'null') return null
    if (p.startsWith('/') || p =~ /^[a-zA-Z][a-zA-Z0-9+\-.]*:\/\//) return file(p)
    def base = file(params.input).toAbsolutePath().parent
    return file(base.resolve(p))
}

def classifyAndValidateRow(meta, fastq_1, fastq_2, bam, library_type, fastq_dir, matrix_dir) {
    def sample   = meta.id
    def mode     = params.mode

    // Treat empty strings as absent (nf-schema returns "" for missing optional columns)
    def has_fastq_1    = fastq_1     ? true : false
    def has_fastq_2    = fastq_2     ? true : false
    def has_bam        = bam         ? true : false
    def has_lib        = library_type ? true : false
    def has_fastq_dir  = fastq_dir   ? true : false
    def has_matrix_dir = matrix_dir  ? true : false

    if (mode == 'bulk') {

        // ── bulk FASTQ mode ────────────────────────────────────────────────────
        if (has_fastq_1) {
            // Check no SC or BAM columns are present
            if (has_bam) {
                error("Sample '${sample}' [bulk FASTQ]: 'bam' column must be absent when 'fastq_1' is provided.")
            }
            if (has_fastq_dir) {
                error("Sample '${sample}' [bulk FASTQ]: 'fastq_dir' is a single-cell column and must be absent.")
            }
            if (has_matrix_dir) {
                error("Sample '${sample}' [bulk FASTQ]: 'matrix_dir' is a single-cell column and must be absent.")
            }
            // library_type required
            if (!has_lib) {
                error("Sample '${sample}' [bulk FASTQ]: 'library_type' (SE or PE) is required.")
            }
            // --star_index is only needed for FASTQ-start samples, so it cannot be
            // checked in validateInputParameters() — that runs before the samplesheet
            // is read. Without it STAR_ALIGN gets an empty channel and is skipped, so
            // the sample would silently produce no BAM and vanish from the results.
            if (!params.star_index) {
                error("Sample '${sample}' [bulk FASTQ]: --star_index is required to align FASTQ input. " +
                      "Supply --star_index, or provide a pre-aligned 'bam' instead.")
            }
            // SE/PE consistency
            if (library_type == 'SE' && has_fastq_2) {
                error("Sample '${sample}' [bulk FASTQ]: library_type is 'SE' but 'fastq_2' is also provided.")
            }
            if (library_type == 'PE' && !has_fastq_2) {
                error("Sample '${sample}' [bulk FASTQ]: library_type is 'PE' but 'fastq_2' is missing.")
            }

            def single_end = (library_type == 'SE')
            def reads      = single_end ? [ resolveFromSamplesheet(fastq_1) ] : [ resolveFromSamplesheet(fastq_1), resolveFromSamplesheet(fastq_2) ]
            def new_meta   = meta + [ mode: 'bulk', start: 'fastq', single_end: single_end ]
            return [ new_meta, reads ]
        }

        // ── bulk BAM mode ──────────────────────────────────────────────────────
        else if (has_bam) {
            if (has_fastq_1 || has_fastq_2) {
                error("Sample '${sample}' [bulk BAM]: 'fastq_1'/'fastq_2' must be absent when 'bam' is provided.")
            }
            if (has_fastq_dir) {
                error("Sample '${sample}' [bulk BAM]: 'fastq_dir' is a single-cell column and must be absent.")
            }
            if (has_matrix_dir) {
                error("Sample '${sample}' [bulk BAM]: 'matrix_dir' is a single-cell column and must be absent.")
            }
            if (!has_lib) {
                error("Sample '${sample}' [bulk BAM]: 'library_type' (SE or PE) is required.")
            }

            def single_end = (library_type == 'SE')
            def new_meta   = meta + [ mode: 'bulk', start: 'bam', single_end: single_end ]
            return [ new_meta, resolveFromSamplesheet(bam) ]
        }

        else {
            error(
                "Sample '${sample}' [--mode bulk]: cannot determine start point.\n" +
                "  Provide 'fastq_1' (and optionally 'fastq_2') for FASTQ input, or 'bam' for BAM input.\n" +
                "  Both 'library_type' and either 'fastq_1'/'bam' are required."
            )
        }

    } else if (mode == 'sc') {

        // ── sc FASTQ mode ──────────────────────────────────────────────────────
        if (has_fastq_dir) {
            if (has_bam) {
                error("Sample '${sample}' [sc FASTQ]: 'bam' must be absent when 'fastq_dir' is provided.")
            }
            if (has_fastq_1 || has_fastq_2) {
                error("Sample '${sample}' [sc FASTQ]: 'fastq_1'/'fastq_2' are bulk columns and must be absent. Provide 'fastq_dir' for single-cell FASTQ.")
            }
            if (has_lib) {
                error("Sample '${sample}' [sc FASTQ]: 'library_type' is a bulk column and must be absent.")
            }
            if (has_matrix_dir) {
                error("Sample '${sample}' [sc FASTQ]: 'matrix_dir' must be absent when starting from FASTQ (CellRanger will generate it).")
            }

            def new_meta = meta + [ mode: 'sc', start: 'fastq', single_end: null ]
            return [ new_meta, resolveFromSamplesheet(fastq_dir) ]
        }

        // ── sc BAM mode ────────────────────────────────────────────────────────
        else if (has_bam) {
            if (has_fastq_dir) {
                error("Sample '${sample}' [sc BAM]: 'fastq_dir' must be absent when 'bam' is provided.")
            }
            if (has_fastq_1 || has_fastq_2) {
                error("Sample '${sample}' [sc BAM]: 'fastq_1'/'fastq_2' are bulk columns and must be absent.")
            }
            if (has_lib) {
                error("Sample '${sample}' [sc BAM]: 'library_type' is a bulk column and must be absent.")
            }
            if (!has_matrix_dir) {
                error("Sample '${sample}' [sc BAM]: 'matrix_dir' (CellRanger filtered_feature_bc_matrix directory) is required.")
            }

            def new_meta = meta + [ mode: 'sc', start: 'bam', single_end: null ]
            return [ new_meta, resolveFromSamplesheet(bam), resolveFromSamplesheet(matrix_dir) ]
        }

        else {
            error(
                "Sample '${sample}' [--mode sc]: cannot determine start point.\n" +
                "  Provide 'fastq_dir' for FASTQ input (CellRanger will run), or 'bam' + 'matrix_dir' for BAM input."
            )
        }

    } else {
        error("Unknown --mode '${mode}'. Must be 'bulk' or 'sc'.")
    }
}

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            params.mode == 'bulk' ? "fastp (Chen et al. 2018)," : "",
            params.mode == 'bulk' ? "STAR (Dobin et al. 2013)," : "",
            params.mode == 'sc'   ? "CellRanger (10x Genomics)," : "",
            "MARINE (Brannan et al.),"  ,
            params.mode == 'bulk' && params.run_sailor ? "SAILOR (Vogel et al.)," : "",
            "featureCounts (Liao et al. 2014),"  ,
            "MultiQC (Ewels et al. 2016)",
            "."
        ].findAll { it }.join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048.</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    if (meta.manifest_map.doi) {
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
