process FASTP {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/02_fastp", mode: params.publish_dir_mode

    conda 'bioconda::fastp=1.3.4'
    container "community.wave.seqera.io/library/fastp:1.3.4--b75b637bf1c0f4b1"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}*.trimmed*.fastq.gz"), emit: reads
    path "${meta.id}.fastp.html",                           emit: html
    path "${meta.id}.fastp.json",                           emit: json
    path "versions.yml",                                    emit: versions

    script:
    def prefix = meta.id
    if (meta.single_end) {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --out1 ${prefix}.trimmed.fastq.gz \\
            --thread ${task.cpus} \\
            --qualified_quality_phred 6 \\
            --length_required 20 \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
        END_VERSIONS
        """
    } else {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${prefix}_R1.trimmed.fastq.gz \\
            --out2 ${prefix}_R2.trimmed.fastq.gz \\
            --detect_adapter_for_pe \\
            --thread ${task.cpus} \\
            --qualified_quality_phred 6 \\
            --length_required 20 \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
        END_VERSIONS
        """
    }
}
