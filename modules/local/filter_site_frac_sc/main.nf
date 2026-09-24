// SC Filter 5: drops sites where more than params.site_max_frac of the reads are
// edited, with depth from SITE_DEPTH_SC and edits restricted to the same cells.
// Also draws the plots for F1–F5, so it runs whether or not F5 is on: with
// site_depth/cells empty ([]) it only plots the F1–F4 steps.
process FILTER_SITE_FRAC_SC {
    tag "${meta.id}"
    label 'process_low'
    publishDir { "${params.outdir}/03_filter_sc/${meta.id}" }, mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    conda 'conda-forge::python>=3.8 conda-forge::pandas>=2.0 conda-forge::matplotlib-base>=3.7'
    container { params.ribostamp_utils_sif as String ?: 'docker.io/aravindsundaravadivelu/ribostamp_utils:1.0.0' }

    input:
    tuple val(meta), path(filtered_edits), path(step_stats), path(site_depth), path(cells)

    output:
    tuple val(meta), path("filtered_edits_site_frac.tsv"), emit: filtered, optional: true
    tuple val(meta), path("site_frac.tsv"),                emit: sites,    optional: true
    tuple val(meta), path("site_frac_summary.tsv"),        emit: summary,  optional: true
    tuple val(meta), path("filter_step_stats.tsv"),        emit: stats,    optional: true
    tuple val(meta), path("*.png"),                        emit: plots
    path "versions.yml",                                   emit: versions

    script:
    def f5_args = site_depth ? [
        "--filtered-edits ${filtered_edits}",
        "--site-depth ${site_depth}",
        "--barcodes ${cells}",
        "--site-max-frac ${params.site_max_frac}",
        "--min-count ${params.min_count}",
        params.filter_sc_recount_min_count ? '' : '--no-recount-min-count',
    ].findAll { arg -> arg }.join(' ') : ''
    """
    helper_filter_site_frac_sc.py \\
        --step-stats     ${step_stats} \\
        ${f5_args} \\
        --sample         ${meta.id} \\
        --output-dir     .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
