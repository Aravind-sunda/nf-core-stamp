// Guarantees the dbSNP BED is coordinate-sorted so that downstream bedtools calls
// can use the `-sorted` sweep, which streams both files instead of loading the whole
// dbSNP into an in-memory interval tree (a full hg38 dbSNP costs >100 GB that way).
//
// The sort order here (`sort -k1,1 -k2,2n`) must stay identical to the order the
// filter scripts sort their edit sites with, or bedtools `-sorted` will abort.
process PREPARE_DBSNP {
    tag "${dbsnp_bed.name}"
    label 'process_low'

    conda 'conda-forge::coreutils=9.1'
    container { params.ribostamp_utils_sif as String ?: 'docker.io/aravindsundaravadivelu/ribostamp_utils:1.0.0' }

    input:
    path dbsnp_bed

    output:
    path "dbsnp_sorted.bed", emit: bed
    path "versions.yml",     emit: versions

    script:
    // Leave headroom under the task allocation so sort spills to disk rather than being OOM-killed.
    def sort_mem = Math.max(1, task.memory.toGiga().intValue() - 2)
    """
    # LC_ALL=C is required, not cosmetic. The filter scripts order their edit sites
    # with pandas, which compares strings by Unicode code point. A UTF-8 locale makes
    # GNU sort collate case-insensitively and ignore punctuation, so it orders
    # 'chr1' before 'GL000008.2' while pandas does the reverse. Mixed references —
    # chr-prefixed primaries alongside GL/KI scaffolds and transgene contigs — then
    # disagree between the two sides and bedtools '-sorted' aborts.
    export LC_ALL=C

    # Checking first keeps an already-sorted file free: it is symlinked through
    # rather than copied, which matters when the file is tens of GB.
    if sort -c -k1,1 -k2,2n ${dbsnp_bed} 2>/dev/null; then
        echo "dbSNP BED is already coordinate-sorted — using as-is."
        ln -s \$(readlink -f ${dbsnp_bed}) dbsnp_sorted.bed
    else
        echo "dbSNP BED is not coordinate-sorted — sorting once (cached for later runs)."
        sort -k1,1 -k2,2n -S ${sort_mem}G -T . ${dbsnp_bed} > dbsnp_sorted.bed
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version | head -n1 | sed 's/^sort (GNU coreutils) //')
    END_VERSIONS
    """
}
