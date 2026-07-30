// Annotates SAILOR ranked edit sites to genes and normalises them to expression
// using the same EPR/EPKM/EPKMR/EPM/EPMR metrics as NORMALIZE_EDITS_BULK, so the
// MARINE and SAILOR routes produce directly comparable tables.
process NORMALIZE_EDITS_SAILOR {
    tag "${meta.id}"
    label 'process_low'
    publishDir { "${params.outdir}/04_sailor/normalized/${meta.id}" }, mode: params.publish_dir_mode

    conda 'conda-forge::python=3.8 conda-forge::pandas=2.0 bioconda::pybedtools=0.9'
    container { params.ribostamp_utils_sif as String ?: 'docker.io/aravindsundaravadivelu/ribostamp_utils:1.0.0' }

    input:
    // combined by .combine() in the subworkflow: matrix is broadcast 1-to-many across samples
    tuple val(meta), path(ranked_bed), path(counts_matrix)
    val  strandedness  // consensus strandedness (0/1/2) — same value SAILOR itself ran with
    path gene_bed

    output:
    tuple val(meta), path("*.sailor.EPR_EPKM_normalized.tsv"), emit: normalized
    path "versions.yml",                                       emit: versions

    script:
    def prefix = meta.id
    """
    helper_normalize_edits_sailor.py \\
        --input_bed ${ranked_bed} \\
        --feature_counts_matrix ${counts_matrix} \\
        --annotation_bed ${gene_bed} \\
        --sample_name ${prefix} \\
        --strandedness ${strandedness} \\
        --min-confidence ${params.sailor_conf} \\
        --output_directory .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pybedtools: \$(python -c "import pybedtools; print(pybedtools.__version__)")
    END_VERSIONS
    """
}
