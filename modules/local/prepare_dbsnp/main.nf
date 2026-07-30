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
    // Written into sorted/ rather than the task root so the result can never collide
    // with the staged input. That collision is reachable in normal use: the output is
    // named after the input, so feeding a previously sorted BED back in would otherwise
    // abort with "ln: failed to create symbolic link ...: File exists".
    path "sorted/*",     emit: bed
    path "versions.yml", emit: versions

    script:
    // Leave headroom under the task allocation so sort spills to disk rather than being OOM-killed.
    def sort_mem = Math.max(1, task.memory.toGiga().intValue() - 2)
    // Name the result after the input with '.sorted' before the extension
    // (hg38_dbsnp.bed3 -> hg38_dbsnp.sorted.bed3). Any existing '.sorted' is dropped
    // first, so re-feeding a sorted BED does not accumulate '.sorted.sorted'.
    def ext    = dbsnp_bed.name.tokenize('.').last()
    def base   = dbsnp_bed.name.substring(0, dbsnp_bed.name.length() - ext.length() - 1)
    def prefix = base.endsWith('.sorted') ? base.substring(0, base.length() - 7) : base
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
    mkdir -p sorted
    if sort -c -k1,1 -k2,2n ${dbsnp_bed} 2>/dev/null; then
        echo "dbSNP BED is already coordinate-sorted — using as-is."
        ln -s \$(readlink -f ${dbsnp_bed}) sorted/${prefix}.sorted.${ext}
    else
        echo "dbSNP BED is not coordinate-sorted — sorting once (cached for later runs)."
        sort -k1,1 -k2,2n -S ${sort_mem}G -T . ${dbsnp_bed} > sorted/${prefix}.sorted.${ext}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version | head -n1 | sed 's/^sort (GNU coreutils) //')
    END_VERSIONS
    """
}
