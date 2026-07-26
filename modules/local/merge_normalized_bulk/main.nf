// Merges per-sample normalized edit tables into gene x sample matrices, one per
// metric (EPR/EPKM/EPKMR/EPM/EPMR plus their inputs). Runs once per pipeline
// after all NORMALIZE_EDITS_BULK jobs finish.
process MERGE_NORMALIZED_BULK {
    label 'process_single'
    publishDir "${params.outdir}/06_filter_normalize/combined", mode: params.publish_dir_mode

    conda 'conda-forge::python=3.8 conda-forge::pandas=2.0'
    container "docker.io/aravindsundaravadivelu/ribostamp_utils:1.0.0"

    input:
    path normalized_files  // list of all *.EPR_EPKM_normalized.tsv files

    output:
    path "normalized_matrix_*.tsv", emit: matrices
    path "versions.yml",            emit: versions

    script:
    """
    helper_merge_normalized_bulk.py \\
        --indir . \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
