// Normalizes filtered MARINE SC edits by per-cell UMI counts from the CellRanger
// filtered_feature_bc_matrix. Computes EPR, EPKM, EPKMR, EPM, EPMR metrics.
// matrix_dir and filtered_edits are combined in the subworkflow with .join() so
// they arrive together as a single tuple.
process NORMALIZE_EDITS_SC {
    tag "${meta.id}"
    label 'process_low'
    publishDir { "${params.outdir}/04_normalize_sc/${meta.id}" }, mode: params.publish_dir_mode

    conda 'conda-forge::python>=3.8 conda-forge::pandas>=2.0 conda-forge::scanpy>=1.9'
    container { params.ribostamp_utils_sif as String ?: 'docker.io/aravindsundaravadivelu/ribostamp_utils:1.0.0' }

    input:
    // filtered_edits and matrix_dir joined per-sample in sc_marine subworkflow
    tuple val(meta), path(filtered_edits), path(matrix_dir)
    path gene_bed

    output:
    tuple val(meta), path("normalized_edits.tsv"), emit: normalized
    path "versions.yml",                            emit: versions

    script:
    def prefix = meta.id
    """
    helper_normalize_edits_umi_sc.py \\
        --filtered_edits  ${filtered_edits} \\
        --counts_matrix   ${matrix_dir} \\
        --bed             ${gene_bed} \\
        --output_dir      . \\
        --edit-type       "${params.edit_type}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
