process STAR_ALIGN {
    tag "${meta.id}"
    label 'process_star'
    publishDir "${params.outdir}/03_star", mode: params.publish_dir_mode

    conda 'bioconda::star bioconda::samtools'
    container "community.wave.seqera.io/library/samtools_star:3d56ec4ef8fcee61"

    input:
    tuple val(meta), path(reads)
    path  star_index
    path  gtf

    output:
    tuple val(meta), path("${meta.id}.Aligned.sortedByCoord.out.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.Aligned.sortedByCoord.out.bam.bai"), emit: bai
    path "${meta.id}.Log.final.out",                                        emit: log_final
    path "${meta.id}.Log.out",                                              emit: log_out
    path "${meta.id}.Log.progress.out",                                     emit: log_progress
    path "versions.yml",                                                    emit: versions

    script:
    def prefix   = meta.id
    def read_cmd = meta.single_end ? "${reads[0]}" : "${reads[0]} ${reads[1]}"
    """
    STAR \\
        --alignEndsType EndToEnd \\
        --genomeDir ${star_index} \\
        --genomeLoad NoSharedMemory \\
        --outBAMcompression 10 \\
        --outFileNamePrefix ${prefix}. \\
        --outFilterMultimapNmax 10 \\
        --outFilterMultimapScoreRange 1 \\
        --outFilterScoreMin 10 \\
        --outReadsUnmapped Fastx \\
        --outSAMattributes All \\
        --outSAMmode Full \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --readFilesCommand zcat \\
        --readFilesIn ${read_cmd} \\
        --runMode alignReads \\
        --runThreadN ${task.cpus} \\
        --sjdbGTFfile ${gtf}

    samtools index -@ ${task.cpus} ${prefix}.Aligned.sortedByCoord.out.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
