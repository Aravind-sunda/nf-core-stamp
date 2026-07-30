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

# ── Environment ───────────────────────────────────────────────────────────────
module load mamba
conda init bash
mamba activate
mamba activate nextflow

export NXF_SINGULARITY_CACHEDIR=/condo/brannanlab/tmhaxs421/singularity_cache
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

PIPELINE_DIR="/condo/brannanlab/tmhaxs421/MARINE_NextFlow/nf-core-stamp"
SNAKEMAKE_PATH="/cm/shared/apps/snakemake/9.13.4/bin/snakemake"

PROFILE="path/to/your/config_file.config"
OUTPUT_DIR="path/to/your/output_directory"
# results will be stores under ${OUTPUT_DIR}/results

nextflow run "$PIPELINE_DIR/main.nf" \
    -c ${PROFILE} \
    -profile "singularity" \
    --outdir "${OUTPUT_DIR}/results" \
    --sailor_snakemake_path "$SNAKEMAKE_PATH" \
    # -resume