// Runs MARINE in bulk RNA-seq mode using the provided Singularity SIF.
// Strandedness is taken from the INFER_STRANDEDNESS output (0/1/2).
// Edit type (e.g. "C>T") is converted to MARINE's CT/AG notation for
// the --sailor and --bedgraphs flags.
// MD tags are required by MARINE to detect mismatches; STAR (--outSAMattributes All)
// emits them, but user-supplied BAMs often lack them — samtools calmd is run
// automatically before MARINE when they are absent (mirrors MARINE_SC).
process MARINE_BULK {
    tag "${meta.id}"
    label 'process_marine_bulk'
    // MARINE names its own output folder after the sample (--output_folder below), so
    // publishing into a further ${meta.id} directory produced 04_marine/<sample>/<sample>/.
    // Publish into 04_marine and let that folder supply the sample level. The pattern is
    // required: without it versions.yml from every sample would collide at 04_marine/.
    // versions.yml is still emitted below and aggregated into pipeline_info/.
    publishDir { "${params.outdir}/04_marine" }, mode: params.publish_dir_mode, pattern: "${meta.id}"

    // No conda directive: MARINE has no Bioconda package.
    // Use -profile singularity or -profile docker; -profile conda is not supported for this process.
    container { params.marine_sif as String ?: 'docker.io/aravindsundaravadivelu/marine:1.0.2' }

    input:
    tuple val(meta), path(bam), path(bai), val(strandedness)
    path fasta
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
    # Check if MD tags are present in the first 100 mapped reads; add them with
    # samtools calmd if absent (MARINE finds no edits on a BAM lacking MD tags).
    MD_COUNT=\$(samtools view -F 4 ${bam} 2>/dev/null | head -100 | grep -c 'MD:Z:' || echo 0)

    if [[ "\${MD_COUNT}" -gt 0 ]]; then
        FINAL_BAM="${bam}"
    else
        samtools calmd -b ${bam} ${fasta} > ${prefix}.md.bam
        FINAL_BAM="${prefix}.md.bam"
    fi

    # Index if neither BAI convention is present
    if [[ ! -f "\${FINAL_BAM}.bai" ]] && [[ ! -f "\${FINAL_BAM%.bam}.bai" ]]; then
        samtools index -@ ${task.cpus} "\${FINAL_BAM}"
    fi

    python /opt/MARINE/marine.py \\
        --bam "\${FINAL_BAM}" \\
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
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
