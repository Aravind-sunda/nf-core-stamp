// Normalizes filtered MARINE bulk edit calls to EPR/EPKM/EPKMR/EPM/EPMR
// metrics using the per-sample featureCounts expression matrix.
// Also outputs a BED6 bedgraph of editing fractions (used by metaPlotR).
process NORMALIZE_EDITS_BULK {
    tag "${meta.id}"
    label 'process_low'
    publishDir { "${params.outdir}/06_filter_normalize/${meta.id}/normalized" }, mode: params.publish_dir_mode

    conda 'conda-forge::python=3.8 conda-forge::pandas=2.0 conda-forge::matplotlib-base=3.7'
    container "docker.io/aravindsundaravadivelu/ribostamp_utils:1.0.0"

    input:
    // combined by .combine() in subworkflow: matrix is broadcast 1-to-many across samples
    tuple val(meta), val(strandedness), path(filtered_edits), path(counts_matrix)

    output:
    tuple val(meta), path("bedgraphs/"), emit: bedgraph_dir
    path "versions.yml",                 emit: versions

    script:
    def prefix = meta.id
    """
    helper_normalize_edits_bulk.py \\
        --input_file ${filtered_edits} \\
        --feature_counts_matrix ${counts_matrix} \\
        --sample_name ${prefix} \\
        --strandedness ${strandedness} \\
        --edit_type "${params.edit_type}" \\
        --output_directory .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
