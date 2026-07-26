// Filters MARINE SC edit calls against dbSNP and applies minimum coverage /
// maximum editing-fraction thresholds. Also removes multi-conversion sites
// and unannotated entries (feature_type == -1).
// SC mode always uses the annotated MARINE output (strandedness is fixed at 2).
process FILTER_EDITS_SC {
    tag "${meta.id}"
    label 'process_low'
    publishDir { "${params.outdir}/03_filter_sc/${meta.id}" }, mode: params.publish_dir_mode

    conda 'conda-forge::python>=3.8 conda-forge::pandas>=2.0 bioconda::pybedtools>=0.9 conda-forge::matplotlib-base>=3.7'
    container { params.ribostamp_utils_sif as String ?: 'docker.io/aravindsundaravadivelu/ribostamp_utils:1.0.0' }

    input:
    tuple val(meta), path(marine_dir)
    path dbsnp_bed

    output:
    tuple val(meta), path("filtered_edits.tsv"), emit: filtered
    path "versions.yml",                          emit: versions

    script:
    def prefix = meta.id
    def skip_flags = [
        params.filter_sc_multi_conversion ? "" : "--no-filter-multi-conversion",
        params.filter_sc_dbsnp            ? "" : "--no-filter-dbsnp",
        params.filter_sc_min_count        ? "" : "--no-filter-min-count",
        params.filter_sc_max_frac         ? "" : "--no-filter-max-frac",
        params.filter_sc_unannotated      ? "" : "--no-filter-unannotated",
    ].findAll { it }.join(" ")
    """
    helper_filter_edits_sc.py \\
        --marine-results ${marine_dir}/final_filtered_site_info_annotated.tsv \\
        --dbsnp-bed      ${dbsnp_bed} \\
        --min-count      ${params.min_count} \\
        --max-frac       ${params.max_frac} \\
        ${skip_flags} \\
        --output-dir     .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pybedtools: \$(python -c "import pybedtools; print(pybedtools.__version__)")
    END_VERSIONS
    """
}
