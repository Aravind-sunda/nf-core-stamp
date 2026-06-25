// Decompress a single .gz file. Used to normalise compressed reference inputs
// (GTF, BED) before passing them to tools that do not accept gzipped files.
process GUNZIP {
    tag "$archive"
    label 'process_single'

    input:
    path archive

    output:
    path "${archive.baseName}", emit: file
    path "versions.yml",        emit: versions

    script:
    """
    gunzip -c ${archive} > ${archive.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
