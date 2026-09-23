// Total read depth at each F1–F4 surviving site, across every cell in the barcode
// list -- including cells with no edit there. This is the denominator for SC
// Filter 5; MARINE's `coverage` column only covers cells that carry an edit.
//
// Read selection mirrors MARINE_SC so numerator and denominator count the same reads:
//   -D CB:<cells>   same cells as the edit table
//   -d xf:25        CellRanger "confidently mapped, UMI-deduplicated", what MARINE keeps
//   -q MAPQ         = MARINE --min_read_quality
//   depth -q BQ     = MARINE --min_base_quality
//   depth -a        listed positions with zero depth are reported, not dropped
// Two view passes because samtools keeps one tag-filter slot: -D CB and -d xf in
// one call fail with 'Different tag "CB" was specified before'.
process SITE_DEPTH_SC {
    tag "${meta.id}"
    label 'process_medium'
    publishDir { "${params.outdir}/03_filter_sc/${meta.id}" }, mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    container { params.marine_sif as String ?: 'docker.io/aravindsundaravadivelu/marine:1.0.2' }

    input:
    tuple val(meta), path(bam), path(bai), path(barcodes, stageAs: 'input_barcodes/*'), path(sites_bed)

    output:
    tuple val(meta), path("site_depth.tsv"), path("cells.txt"), emit: depth
    path "versions.yml",                                        emit: versions

    script:
    """
    zcat -f ${barcodes} | cut -f1 | grep -v '^\$' > cells.txt

    printf 'contig\\tposition\\tdepth\\n' > site_depth.tsv
    samtools view -@ ${task.cpus} -u -q ${params.min_read_quality} -D ${params.barcode_tag}:cells.txt ${bam} \\
        | samtools view -@ ${task.cpus} -u -d xf:25 - \\
        | samtools depth -a -b ${sites_bed} -q ${params.min_base_quality} - \\
        >> site_depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
