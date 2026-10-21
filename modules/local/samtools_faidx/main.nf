process SAMTOOLS_FAIDX {
    tag "$fasta"
    label 'process_low'

    conda 'bioconda::samtools'
    container "community.wave.seqera.io/library/samtools_star:3d56ec4ef8fcee61"

    input:
    path fasta

    output:
    path "${fasta}.fai", emit: fai
    path "versions.yml", emit: versions

    script:
    """
    samtools faidx ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
