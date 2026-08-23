// Runs FLARE cluster-identification mode for one RBP-STAMP sample, downstream
// of SAILOR. Like the SAILOR process, this has NO container directive: FLARE's
// bundled Snakefile drives its own Singularity containers (docker://ekofman/editc:v2)
// via --use-singularity, and nesting containers must be avoided. Snakemake must
// therefore be pre-installed on PATH (or set --sailor_snakemake_path) when using
// -profile singularity — the same requirement as SAILOR.
//
// FLARE consumes this sample's SAILOR outputs (ranked BED, fwd/rev bigwigs,
// filtered-merged BAM) plus the genome-level regions folder.
process FLARE {
    tag "$meta.id"
    label 'process_sailor'
    // FLARE names its own output folder after the sample (--output_folder below), so
    // publishing into a further ${meta.id} directory would give 05_flare/<sample>/<sample>/.
    // Publish into 05_flare and let FLARE's folder supply the sample level — the same
    // arrangement as MARINE_BULK. Previously this folder was the constant 'flare_output',
    // so every sample published to one path and the last to finish replaced the rest.
    // saveAs drops versions.yml, which would otherwise collide across samples at
    // 05_flare/; it is still emitted below and aggregated into pipeline_info/.
    // Must be saveAs rather than `pattern`: pattern is evaluated eagerly, where meta
    // is not in scope, so referencing it there fails with "No such variable: meta".
    publishDir { "${params.outdir}/05_flare" }, mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    // Snakemake is installed via conda when using -profile conda; ignored under
    // -profile singularity (conda.enabled=false globally), where it must be on PATH.
    conda 'bioconda::snakemake-minimal=9.13.4 conda-forge::python>=3.8'

    input:
    tuple val(meta), path(stamp_bed), path(fwd_bw), path(rev_bw), path(bam), path(bai)
    path fasta
    path fasta_fai   // staged as genome.fa.fai so pysam finds the index
    path regions     // FLARE regions folder (cluster-identification mode)
    path snakefile   // FLARE Snakefile

    output:
    tuple val(meta), path("${meta.id}/FLARE/${meta.id}_merged_sorted_peaks.*.scored.tsv"), emit: scored_peaks
    path "${meta.id}/",                                                                    emit: output_dir
    path "versions.yml",                                                                   emit: versions

    script:
    def prefix     = meta.id
    // FLARE expects the 2-char edit code (C>T → CT, A>G → AG), same as SAILOR.
    def edit_flare = params.edit_type.replace(">", "").replaceAll(/\s/, "")
    def sing_cache = params.sailor_singularity_cache ?: ''
    def sing_bind  = params.sailor_singularity_bind  ?: ''
    """
    # Verify Snakemake is available before doing any work
    SNAKEMAKE_BIN="${params.sailor_snakemake_path}"
    command -v "\${SNAKEMAKE_BIN}" > /dev/null 2>&1 || {
        echo "ERROR: snakemake not found at '\${SNAKEMAKE_BIN}'."
        echo "       Pre-install it (e.g. conda install snakemake) or set --sailor_snakemake_path."
        exit 1
    }

    mkdir -p ${prefix}

    # Build the per-sample FLARE (cluster-identification) JSON config. Reference the
    # staged workdir copies (via \$(pwd)) rather than realpath so index sidecars stay
    # adjacent: fasta.fai (which may have been generated into the workdir) sits next
    # to the staged fasta, and the BAM's .bai next to the staged BAM. The workdir is
    # bound into the editc container, and each input's real directory is bound too, so
    # symlinks resolve either way. The regions folder needs no sidecar, so realpath is
    # fine there.
    helper_make_flare_json.py \\
        --label ${prefix} \\
        --output_json ./flare_config.json \\
        --output_folder \$(pwd)/${prefix} \\
        --stamp_sites_file \$(pwd)/${stamp_bed} \\
        --forward_bw \$(pwd)/${fwd_bw} \\
        --reverse_bw \$(pwd)/${rev_bw} \\
        --bam \$(pwd)/${bam} \\
        --fasta \$(pwd)/${fasta} \\
        --regions \$(realpath ${regions}) \\
        --edit_type ${edit_flare} \\
        --fdr_threshold ${params.flare_fdr} \\
        --max_merge_dist ${params.flare_max_merge_dist}

    # Build Snakemake options. --snakefile is realpath'd so FLARE's Snakefile can
    # derive its scripts directory (<assets>/workflow_FLARE/scripts) from the path.
    SNAKEMAKE_OPTS=(
        --snakefile \$(realpath ${snakefile})
        --configfile \$(realpath ./flare_config.json)
        --cores ${task.cpus}
        --latency-wait 120
        --rerun-incomplete
    )

    SNAKEMAKE_OPTS+=(--use-singularity)
    [[ -n "${sing_cache}" ]] && SNAKEMAKE_OPTS+=(--singularity-prefix "${sing_cache}")

    # Resolve bind paths so the editc container can reach every input. Use the
    # user-supplied value, or auto-derive from the real (post-symlink) paths of
    # all inputs, plus the FLARE Snakefile dir (which holds the scripts/ folder).
    if [[ -n "${sing_bind}" ]]; then
        SNAKEMAKE_OPTS+=(--singularity-args "--bind ${sing_bind}")
    else
        _bind_dirs=\$(
            { dirname \$(realpath ${fasta});
              dirname \$(realpath ${fasta_fai});
              dirname \$(realpath ${stamp_bed});
              dirname \$(realpath ${fwd_bw});
              dirname \$(realpath ${rev_bw});
              dirname \$(realpath ${bam});
              realpath ${regions};
              dirname \$(realpath ${snakefile});
              pwd;
            } | sort -u | paste -sd,
        )
        [[ -n "\$_bind_dirs" ]] && SNAKEMAKE_OPTS+=(--singularity-args "--bind \$_bind_dirs")
    fi

    "\${SNAKEMAKE_BIN}" "\${SNAKEMAKE_OPTS[@]}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snakemake: \$("\${SNAKEMAKE_BIN}" --version)
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
