// Runs CellRanger count for one 10x Genomics single-cell sample.
// IMPORTANT: CellRanger requires a licensed installation. The community container
// at 'nf-core/cellranger' embeds the binary; ensure your institution accepts the
// 10x Genomics license before use. See nf-core/scrnaseq for the container approach.
//
// FASTQ naming requirement: files inside fastq_dir must follow the 10x naming scheme:
//   {sample}_S{n}_L{lane}_R1_001.fastq.gz  (barcodes + UMI)
//   {sample}_S{n}_L{lane}_R2_001.fastq.gz  (cDNA read)
// Rename downloaded SRA/GEO FASTQs to match this convention before running.
process CELLRANGER {
    tag "${meta.id}"
    label 'process_cellranger'
    publishDir { "${params.outdir}/01_cellranger/${meta.id}" }, mode: params.publish_dir_mode,
               saveAs: { fname -> fname.contains("/") ? fname.substring(fname.indexOf("/") + 1) : fname }

    // Community CellRanger container — requires 10x Genomics license acceptance.
    // Override via conf/modules.config if using a locally built or licensed image.
    container "nf-core/cellranger:9.0.1"

    input:
    tuple val(meta), path(fastq_dir)
    path  cellranger_ref

    output:
    tuple val(meta), path("${meta.id}_count/outs/possorted_genome_bam.bam"),
                     path("${meta.id}_count/outs/possorted_genome_bam.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}_count/outs/filtered_feature_bc_matrix/"), emit: matrix_dir
    tuple val(meta), path("${meta.id}_count/outs/"),                            emit: outs_dir
    path "versions.yml",                                                         emit: versions

    script:
    def prefix = meta.id
    def mem_gb = (task.memory.toGiga() as int)
    // --id uses "${prefix}_count" (not "${prefix}") to avoid a name clash with the
    // staged fastq_dir, which Nextflow stages under its basename. If --id equals the
    // fastq_dir name, CellRanger finds the existing directory and errors with
    // "not a pipestance directory". publishDir still uses meta.id for the final path.
    """
    cellranger count \\
        --id="${prefix}_count" \\
        --transcriptome="${cellranger_ref}" \\
        --fastqs="${fastq_dir}" \\
        --sample="${prefix}" \\
        --localcores=${task.cpus} \\
        --localmem=${mem_gb} \\
        --create-bam=true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(cellranger --version 2>&1 | grep -oP 'cellranger \\K[0-9.]+' || cellranger --version 2>&1)
    END_VERSIONS
    """
}
