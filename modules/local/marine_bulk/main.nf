// Runs MARINE in bulk RNA-seq mode using the provided Singularity SIF.
// Strandedness is taken from the INFER_STRANDEDNESS output (0/1/2).
// Edit type (e.g. "C>T") is converted to MARINE's CT/AG notation for
// the --sailor and --bedgraphs flags.
process MARINE_BULK {
    tag "${meta.id}"
    label 'process_marine_bulk'
    publishDir { "${params.outdir}/04_marine/${meta.id}" }, mode: params.publish_dir_mode

    container { params.marine_sif as String ?: 'docker.io/aravindsundaravadivelu/marine:1.0.2' }

    input:
    tuple val(meta), path(bam), path(bai), val(strandedness)
    path gene_bed

    output:
    tuple val(meta), val(strandedness), path("${meta.id}/"), emit: results
    path "versions.yml",                                      emit: versions

    script:
    def prefix     = meta.id
    def pe_flag    = meta.single_end ? "" : "--paired_end"
    // Convert "C>T" → "CT", "A>G" → "AG" for MARINE's --sailor / --bedgraphs flags
    def edit_code  = params.edit_type.replace(">", "")
    """
    python /opt/MARINE/marine.py \\
        --bam ${bam} \\
        --output_folder ${prefix} \\
        --annotation_bedfile_path ${gene_bed} \\
        --strandedness ${strandedness} \\
        --cores ${task.cpus} \\
        --min_base_quality ${params.min_base_quality} \\
        --min_read_quality ${params.min_read_quality} \\
        --sailor ${edit_code} \\
        --bedgraphs ${edit_code} \\
        ${pe_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        marine: 1.0.2
    END_VERSIONS
    """
}
