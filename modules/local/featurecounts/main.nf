process FEATURECOUNTS {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/05_featurecounts", mode: params.publish_dir_mode

    conda 'bioconda::subread=2.1.1'
    container "community.wave.seqera.io/library/subread:2.1.1--cbb5cb85f59ac813"

    input:
    tuple val(meta), path(bam), val(strandedness)
    path gtf

    output:
    tuple val(meta), path("${meta.id}.featurecounts.txt"),           emit: counts
    tuple val(meta), path("${meta.id}.featurecounts.txt.summary"),   emit: summary
    path "versions.yml",                                             emit: versions

    script:
    def prefix   = meta.id
    def pe_flag  = meta.single_end ? "" : "-p --countReadPairs"
    """
    featureCounts \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -t exon \\
        -g gene_name \\
        -s ${strandedness} \\
        ${pe_flag} \\
        -o ${prefix}.featurecounts.txt \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(featureCounts -v 2>&1 | grep -oP '(?<=featureCounts v)\\S+')
    END_VERSIONS
    """
}
