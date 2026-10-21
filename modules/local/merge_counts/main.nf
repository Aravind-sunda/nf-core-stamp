// Merges per-sample featureCounts outputs into one combined count matrix.
// Runs once per pipeline (not per sample) after all featureCounts jobs finish.
process MERGE_COUNTS {
    label 'process_single'
    publishDir "${params.outdir}/05_featurecounts", mode: params.publish_dir_mode

    conda 'conda-forge::python=3.8 conda-forge::pandas=2.0'
    container { params.ribostamp_utils_sif as String ?: 'docker.io/aravindsundaravadivelu/ribostamp_utils:1.0.0' }

    input:
    path counts_files  // list of all *.featurecounts.txt files

    output:
    path "counts_matrix_combined.tsv", emit: matrix
    path "versions.yml",               emit: versions

    script:
    """
    helper_merge_counts.py \\
        --indir . \\
        --outfile counts_matrix_combined.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
