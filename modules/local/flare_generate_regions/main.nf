// Builds the genome-level FLARE regions folder from a GTF, used by FLARE
// cluster-identification mode. This is expensive (~8-10 GB, slow) and genome-,
// not sample-, specific — so it runs once and is cached. Supply a pre-built
// folder with --flare_regions to skip this entirely.
//
// Runs generate_regions.py directly (not via Snakemake), so it uses a real
// container (the editc image, which ships the pandas/pybedtools stack the script
// needs). Override with --editc_sif to point at a local SIF on offline HPC.
process FLARE_GENERATE_REGIONS {
    tag "flare_regions"
    label 'process_medium'
    publishDir "${params.outdir}/05_flare", mode: params.publish_dir_mode

    container { params.editc_sif as String ?: 'docker.io/ekofman/editc:v2' }

    input:
    path gtf
    val  window_size
    path scripts_dir   // assets/workflow_FLARE/scripts (staged so generate_regions.py is present)

    output:
    path "flare_regions/",            emit: regions
    path "flare_regions_region_map.bed", emit: region_map
    path "versions.yml",              emit: versions

    script:
    """
    # generate_regions.py writes chunk files into <output_dir>/ and a
    # <output_dir>_region_map.bed alongside it. We name the folder flare_regions.
    python ${scripts_dir}/generate_regions.py \\
        ${gtf} \\
        flare_regions \\
        --window_size ${window_size}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
        pybedtools: \$(python -c 'import pybedtools; print(pybedtools.__version__)')
    END_VERSIONS
    """
}
