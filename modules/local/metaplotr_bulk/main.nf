// Computes metagene distances from the normalized edit-fraction BED6 file.
// The genePred file is selected based on --genome (built-in assets) or
// overridden with --genepred.
process METAPLOTR_BULK {
    tag "${meta.id}"
    label 'process_single'
    publishDir { "${params.outdir}/07_metaplotr/${meta.id}" }, mode: params.publish_dir_mode

    conda 'conda-forge::python=3.8 conda-forge::pandas=2.0'
    container "docker.io/aravindsundaravadivelu/ribostamp_utils:1.0.0"

    input:
    tuple val(meta), path(bedgraph_dir)
    path genepred

    output:
    tuple val(meta), path("${meta.id}.*.dist.measures.txt"), emit: distances
    path "versions.yml",                                     emit: versions

    script:
    def prefix   = meta.id
    def edit_tag = params.edit_type.replace(">", "_")   // "C>T" → "C_T"
    def bed_file = "${bedgraph_dir}/${prefix}.${edit_tag}.edit_fraction.bed"
    """
    # Sort BED6 by chromosome then position (required by metaPlotR)
    sort -k1,1 -k2,2n ${bed_file} > ${prefix}.${edit_tag}.sorted.bed

    helper_calc_metaplot_dist.py \\
        --genePred ${genepred} \\
        --bed ${prefix}.${edit_tag}.sorted.bed \\
        --out ${prefix}.${edit_tag}.dist.measures.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
