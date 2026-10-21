// Computes metagene distances from SAILOR ranked BED at three confidence tiers:
//   all sites / confidence ≥ 0.5 / confidence ≥ 0.9
// Unstranded libraries (strandedness=0) use --ignore-strand in metaPlotR.
process METAPLOTR_SAILOR {
    tag "${meta.id}"
    label 'process_single'
    publishDir { "${params.outdir}/04_sailor/metaplotr/${meta.id}" }, mode: params.publish_dir_mode

    conda 'conda-forge::python>=3.8 conda-forge::pandas>=2.0'
    container { params.ribostamp_utils_sif as String ?: 'docker.io/aravindsundaravadivelu/ribostamp_utils:1.0.0' }

    input:
    tuple val(meta), path(ranked_bed)
    val  strandedness  // int: 0/1/2 — same consensus used for SAILOR
    path genepred

    output:
    tuple val(meta), path("*.dist.measures.txt"), emit: distances
    path "versions.yml",                           emit: versions

    script:
    def prefix = meta.id
    // SAILOR ranked BED layout: chr start end confidence reads strand (6 columns)
    // Column 4 (1-indexed) is the Bayesian posterior probability / confidence score (0–1)
    // Column 5 is the reads field (e.g. "1,9" = edited,total) — non-numeric, not the confidence
    // Unstranded libraries: metaPlotR cannot determine 5'/3' direction
    def ignore_strand_opt = strandedness == 0 ? "--ignore-strand" : ""
    """
    # Three confidence tiers: all sites / ≥0.5 / ≥0.9
    cp ${ranked_bed} ${prefix}.all.bed
    awk '\$4 >= 0.5' ${ranked_bed} > ${prefix}.conf0.5.bed
    awk '\$4 >= 0.9' ${ranked_bed} > ${prefix}.conf0.9.bed

    for TIER in all conf0.5 conf0.9; do
        BED_FILE="${prefix}.\${TIER}.bed"
        OUT_FILE="${prefix}.\${TIER}.dist.measures.txt"

        if [[ -s "\${BED_FILE}" ]]; then
            helper_calc_metaplot_dist.py \\
                --genePred ${genepred} \\
                --bed "\${BED_FILE}" \\
                --out "\${OUT_FILE}" \\
                ${ignore_strand_opt}
        else
            # Placeholder so publishDir always has 3 output files per sample
            printf "no sites in %s tier for sample %s\\n" "\${TIER}" "${prefix}" > "\${OUT_FILE}"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
