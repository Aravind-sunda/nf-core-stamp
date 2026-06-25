// Runs RSeQC infer_experiment.py on each BAM, then parses the output into
// a strand code (0 = unstranded, 1 = forward, 2 = reverse) using the same
// 0.75-threshold logic as the original bash pipeline.
process INFER_STRANDEDNESS {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/04_strandedness", mode: params.publish_dir_mode

    conda 'bioconda::rseqc=5.0.4'
    container "community.wave.seqera.io/library/rseqc:5.0.4--6dbd0838c4d673ae"

    input:
    tuple val(meta), path(bam), path(bai)
    path gene_bed

    output:
    tuple val(meta), path("${meta.id}_strandedness.txt"), emit: strandedness_txt
    tuple val(meta), path("${meta.id}_strand_code.txt"),  emit: strand_code_txt
    path "versions.yml",                                  emit: versions

    script:
    def prefix = meta.id
    """
    infer_experiment.py \\
        -i ${bam} \\
        -r ${gene_bed} \\
        > ${prefix}_strandedness.txt

    # Parse the infer_experiment.py output to a single strand code (0/1/2).
    # Anchors on the "This is" line, then reads +2 (forward) and +3 (reverse)
    # fractions. Threshold: >= 0.75 = stranded; below for both = unstranded.
    awk '
        /This is/{ found=NR }
        found && NR==found+2{ fwd=\$NF }
        found && NR==found+3{ rev=\$NF }
        END{
            if (fwd+0 >= 0.75) print 1
            else if (rev+0 >= 0.75) print 2
            else print 0
        }
    ' ${prefix}_strandedness.txt > ${prefix}_strand_code.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(infer_experiment.py --version 2>&1 | grep -oP '(?<=version )\\S+')
    END_VERSIONS
    """
}
