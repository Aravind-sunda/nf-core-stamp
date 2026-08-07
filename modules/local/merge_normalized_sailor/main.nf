// Merges per-sample normalized SAILOR edit tables into gene x sample matrices, one
// per metric (EPR/EPKM/EPKMR/EPM/EPMR plus their inputs). Runs once per pipeline
// after all NORMALIZE_EDITS_SAILOR jobs finish.
//
// Shares helper_merge_normalized_bulk.py with MERGE_NORMALIZED_BULK — the two
// callers' normalized tables have identical schemas. Kept as a separate process so
// the SAILOR route's arguments and publish path can diverge from MARINE's without
// touching the MARINE side.
process MERGE_NORMALIZED_SAILOR {
    label 'process_single'
    publishDir "${params.outdir}/04_sailor/normalized/combined", mode: params.publish_dir_mode

    conda 'conda-forge::python=3.8 conda-forge::pandas=2.0'
    container { params.ribostamp_utils_sif as String ?: 'docker.io/aravindsundaravadivelu/ribostamp_utils:1.0.0' }

    input:
    path normalized_files  // list of all *.sailor.EPR_EPKM_normalized.tsv files

    output:
    path "normalized_matrix_*.tsv", emit: matrices
    path "versions.yml",            emit: versions

    script:
    // To also emit matrices for the SAILOR-only per-gene columns, add to the call below:
    //     --extra-metrics n_sites mean_confidence \\
    """
    helper_merge_normalized_bulk.py \\
        --indir . \\
        --outdir . \\
        --suffix .sailor.EPR_EPKM_normalized.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
