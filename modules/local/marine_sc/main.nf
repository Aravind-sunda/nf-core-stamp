// Runs MARINE in single-cell mode.
// Strandedness is hardcoded to 2 for 10x STAMP data (riboSTAMP uses reverse strand).
// MD tags are required by MARINE; if the input BAM lacks them (common for CellRanger
// output) samtools calmd is run automatically before MARINE.
process MARINE_SC {
    tag "${meta.id}"
    label 'process_marine_sc'
    publishDir { "${params.outdir}/02_marine_sc/${meta.id}" }, mode: params.publish_dir_mode

    container { params.marine_sif as String ?: 'docker.io/aravindsundaravadivelu/marine:1.0.2' }

    input:
    tuple val(meta), path(bam), path(bai), path(matrix_dir)
    path fasta
    path gene_bed

    output:
    tuple val(meta), path("marine_output/"), emit: results
    path "versions.yml",                     emit: versions

    script:
    def prefix = meta.id
    """
    # Resolve barcodes whitelist from the matrix directory
    if [[ -f "${matrix_dir}/barcodes.tsv.gz" ]]; then
        BARCODES="${matrix_dir}/barcodes.tsv.gz"
    elif [[ -f "${matrix_dir}/barcodes.tsv" ]]; then
        BARCODES="${matrix_dir}/barcodes.tsv"
    else
        echo "ERROR: barcodes file not found in ${matrix_dir}" >&2
        exit 1
    fi

    # Check if MD tags are present in the first 100 mapped reads
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

    # Run MARINE in single-cell mode (strandedness=2 is always correct for 10x STAMP)
    python /opt/MARINE/marine.py \\
        --bam_filepath            "\${FINAL_BAM}" \\
        --output_folder           marine_output \\
        --barcode_whitelist_file  "\${BARCODES}" \\
        --annotation_bedfile_path ${gene_bed} \\
        --barcode_tag             ${params.barcode_tag} \\
        --strandedness            2 \\
        --cores                   ${task.cpus} \\
        --min_read_quality        ${params.min_read_quality} \\
        --min_base_quality        ${params.min_base_quality}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        marine: 1.0.2
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
// do not need the following intermediate files from the flag
// --keep_intermediate_files \\