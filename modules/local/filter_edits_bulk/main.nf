// Filters MARINE bulk edit calls against dbSNP and applies max editing-fraction
// threshold. For unstranded libraries (strandedness=0) also re-annotates sites
// using the BED6 gene model (MARINE doesn't annotate unstranded runs reliably).
process FILTER_EDITS_BULK {
    tag "${meta.id}"
    label 'process_low'
    publishDir { "${params.outdir}/06_filter_normalize/${meta.id}/filtered" }, mode: params.publish_dir_mode

    conda 'conda-forge::python=3.8 conda-forge::pandas=2.0 bioconda::pybedtools=0.9 conda-forge::matplotlib-base=3.7'
    container "docker.io/aravindsundaravadivelu/ribostamp_utils:1.0.0"

    input:
    tuple val(meta), val(strandedness), path(marine_dir)
    path dbsnp_bed
    path gene_bed

    output:
    tuple val(meta), val(strandedness), path("filtered_edits.tsv"), emit: filtered
    path "versions.yml",                                            emit: versions

    script:
    def prefix = meta.id
    // For strandedness 1/2 prefer the MARINE-annotated file; fall back to unannotated.
    // For strandedness 0 always use unannotated (filter script re-annotates it).
    def marine_input = strandedness == 0
        ? "${marine_dir}/final_filtered_site_info.tsv"
        : "${marine_dir}/final_filtered_site_info_annotated.tsv"

    def annotation_arg = strandedness == 0 ? "--annotation-bed ${gene_bed}" : ""
    def skip_flags = [
        params.filter_bulk_multiallelic ? "" : "--no-filter-multiallelic",
        params.filter_bulk_dbsnp        ? "" : "--no-filter-dbsnp",
        params.filter_bulk_max_frac     ? "" : "--no-filter-max-frac",
        params.filter_bulk_unannotated  ? "" : "--no-filter-unannotated",
    ].findAll { it }.join(" ")
    """
    # Prefer annotated output; fall back to unannotated if it doesn't exist
    MARINE_INPUT="${marine_input}"
    if [ ! -f "\${MARINE_INPUT}" ]; then
        MARINE_INPUT="${marine_dir}/final_filtered_site_info.tsv"
    fi

    helper_filter_edits_bulk.py \\
        --marine-results "\${MARINE_INPUT}" \\
        --strandedness ${strandedness} \\
        --dbsnp-bed ${dbsnp_bed} \\
        --max-frac ${params.max_frac} \\
        ${annotation_arg} \\
        ${skip_flags} \\
        --output-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pybedtools: \$(python -c "import pybedtools; print(pybedtools.__version__)")
    END_VERSIONS
    """
}
