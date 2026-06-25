// Runs the SAILOR Bayesian RNA-editing caller for all bulk samples in one batch.
// Snakemake orchestrates per-sample parallelism and manages its own containers.
// This process intentionally has NO container directive so that Snakemake can
// call Singularity (or conda) for SAILOR rules without nesting containers.
process SAILOR {
    tag "sailor"
    label 'process_sailor'
    publishDir "${params.outdir}/04_sailor", mode: params.publish_dir_mode

    // Snakemake is installed via conda when using -profile conda.
    // Under -profile singularity, this conda directive is ignored (conda.enabled=false globally);
    // Snakemake must be pre-installed and available on PATH (or set --sailor_snakemake_path).
    conda 'bioconda::snakemake-minimal=9.13.4 conda-forge::python>=3.8'

    input:
    path bams           // collected list of BAM files staged into process workdir
    path bais           // collected list of BAI index files staged into process workdir
    val  strandedness   // consensus int: 0=unstranded, 1=forward, 2=reverse
    val  library_type   // 'single' | 'paired'
    path fasta
    path fasta_fai  // staged as genome.fa.fai alongside fasta so bam_to_bw.sh finds $fasta.fai
    path dbsnp_bed
    path snakefile

    output:
    path "sailor_output/*.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed", emit: ranked_beds
    path "sailor_output/",                                                                       emit: output_dir
    path "versions.yml",                                                                         emit: versions

    script:
    // Derive SAILOR parameters from strandedness code:
    //   0 (unstranded) → reverse_stranded=false, mm_tolerance=2
    //   1 (forward)    → reverse_stranded=false, mm_tolerance=1
    //   2 (reverse)    → reverse_stranded=true,  mm_tolerance=1
    def reverse_stranded = strandedness == 2 ? "true" : "false"
    def mm_tolerance     = strandedness == 0 ? 2 : 1
    // SAILOR uses the 2-char code: 'C>T' → 'CT', 'A>G' → 'AG'
    def edit_sailor      = params.edit_type.replace(">", "").replaceAll(/\s/, "")
    // Optional Singularity paths (resolved to empty string if param is null)
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

    mkdir -p sailor_output

    # STAR produces BAMs named {id}.Aligned.sortedByCoord.out.bam. Strip the STAR
    # suffix so SAILOR uses the clean sample ID, not the full filename, as the
    # sample key. BAM-start samples are already cleanly named and unaffected.
    for bam in *.Aligned.sortedByCoord.out.bam; do
        [[ -f "\$bam" ]] || continue
        sample="\${bam%.Aligned.sortedByCoord.out.bam}"
        mv "\$bam" "\${sample}.bam"
        [[ -f "\${bam}.bai" ]] && mv "\${bam}.bai" "\${sample}.bam.bai" || true
    done

    # Generate SAILOR JSON config; all staged BAMs are discoverable in the workdir
    helper_make_sailor_json.py \\
        --samples_path . \\
        --output_dir ./sailor_output \\
        --output_json ./sailor_config.json \\
        --edit_type ${edit_sailor} \\
        --reverse_stranded ${reverse_stranded} \\
        --library ${library_type} \\
        --reference_fasta ${fasta} \\
        --known_snps ${dbsnp_bed} \\
        --mm_tolerance ${mm_tolerance}

    # Build Snakemake options array (bash arrays preserve quoting across whitespace)
    SNAKEMAKE_OPTS=(
        --snakefile \$(realpath ${snakefile})
        --configfile ./sailor_config.json
        --cores ${task.cpus}
        --latency-wait 120
        --rerun-incomplete
    )

    SNAKEMAKE_OPTS+=(--use-singularity)
    [[ -n "${sing_cache}" ]] && SNAKEMAKE_OPTS+=(--singularity-prefix "${sing_cache}")

    # Resolve bind paths: use user-supplied value, or auto-derive from the real
    # (post-symlink) paths of reference files so Singularity containers can reach
    # them regardless of where the user's data lives on the filesystem.
    if [[ -n "${sing_bind}" ]]; then
        SNAKEMAKE_OPTS+=(--singularity-args "--bind ${sing_bind}")
    else
        _bind_dirs=\$(
            { dirname \$(realpath ${fasta});
              dirname \$(realpath ${fasta_fai});
              dirname \$(realpath ${dbsnp_bed});
              dirname \$(realpath ${snakefile});
              for f in *.bam; do [[ -f "\$f" ]] && dirname \$(realpath "\$f"); done;
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
