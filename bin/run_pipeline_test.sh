#!/bin/bash
#SBATCH --job-name=stamp_test
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --time=08:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asundaravadivelu@houstonmethodist.org
#SBATCH --output=slurm_%u_%x_%j.log

# ── Profiles to run — comment out any you want to skip ───────────────────────
PROFILES=(
    test_bulk_bam
    test_bulk_fastq
    test_bulk_bam_flare
    test_bulk_fastq_flare
    test_sc_bam
    test_sc_fastq
)

# ── Environment ───────────────────────────────────────────────────────────────
module load mamba
conda init bash
mamba activate
mamba activate nextflow

export NXF_SINGULARITY_CACHEDIR=/condo/brannanlab/tmhaxs421/singularity_cache
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

PIPELINE_DIR="/condo/brannanlab/tmhaxs421/MARINE_NextFlow/nf-core-stamp"
OUTPUT_DIR="/condo/brannanlab/tmhaxs421/MARINE_NextFlow/nfcore-stamp-run"

# ── System-specific paths — edit these for your cluster ──────────────────────
# Path to the snakemake binary used by SAILOR.
# Use 'snakemake' if it is already on PATH, otherwise provide the full path,
# e.g. '/path/to/conda/envs/snakemake/bin/snakemake'
SNAKEMAKE_PATH="/cm/shared/apps/snakemake/9.13.4/bin/snakemake"

# NXF_SINGULARITY_CACHEDIR (set above) is automatically used as the SAILOR
# Singularity cache via nextflow.config. Only set --sailor_singularity_cache
# explicitly if you need SAILOR to cache containers in a different location.

# ── Loop ─────────────────────────────────────────────────────────────────────
for PROFILE in "${PROFILES[@]}"; do
    echo "======================================================"
    echo "  stamp test: ${PROFILE}"
    echo "  started: $(date)"
    echo "======================================================"

    nextflow run "$PIPELINE_DIR/main.nf" \
        -profile "${PROFILE},singularity" \
        --outdir "${OUTPUT_DIR}/results_${PROFILE}" \
        --sailor_snakemake_path "$SNAKEMAKE_PATH" \
        -resume

    echo "======================================================"
    echo "  stamp test: ${PROFILE} — DONE ($(date))"
    echo "======================================================"
    echo ""
done


# Add the following line to specify the path to the marine_1.0.2.sif file or it will downlaod from docker hub
# --marine_sif "/condo/brannanlab/tmhaxs421/singularity_cache/marine_1.0.2.sif" \

# for PROFILE in "${PROFILES[@]}"; do
# nextflow run /home/tmhaxs421/brannanlab/tmhaxs421/MARINE_NextFlow/nf-core-stamp/main.nf \
#     --profile "${PROFILE},singularity" \
#     -preview -with-dag dag_${PROFILE}.mmd \
#     --sailor_snakemake_path "$SNAKEMAKE_PATH" \
# done