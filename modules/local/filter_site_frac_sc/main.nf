// SC Filter 5: drops sites where more than params.site_max_frac of the reads are
// edited, with depth from SITE_DEPTH_SC and edits restricted to the same cells.
process FILTER_SITE_FRAC_SC {
    tag "${meta.id}"
    label 'process_low'
    publishDir { "${params.outdir}/03_filter_sc/${meta.id}" }, mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    conda 'conda-forge::python>=3.8 conda-forge::pandas>=2.0'
    container { params.ribostamp_utils_sif as String ?: 'docker.io/aravindsundaravadivelu/ribostamp_utils:1.0.0' }

    input:
    tuple val(meta), path(filtered_edits), path(site_depth), path(cells)

    output:
    tuple val(meta), path("filtered_edits_site_frac.tsv"), emit: filtered
    tuple val(meta), path("site_frac.tsv"),                emit: sites
    tuple val(meta), path("site_frac_summary.tsv"),        emit: summary
    path "versions.yml",                                   emit: versions

    script:
    """
    helper_filter_site_frac_sc.py \\
        --filtered-edits ${filtered_edits} \\
        --site-depth     ${site_depth} \\
        --barcodes       ${cells} \\
        --site-max-frac  ${params.site_max_frac} \\
        --min-count      ${params.min_count} \\
        ${params.filter_sc_recount_min_count ? '' : '--no-recount-min-count'} \\
        --output-dir     .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
