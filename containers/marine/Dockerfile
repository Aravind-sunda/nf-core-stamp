# MARINE — Multi-core Algorithm for Rapid Identification of Nucleotide Edits
# Converted from marine_container.def (Singularity) to Docker format.
#
# Build:
#   docker build --platform linux/amd64 -t aravindsundaravadivelu/marine:1.0.2 .
#
# Push to Docker Hub (then use as params.marine_sif via Singularity pull, or directly in Docker mode):
#   docker push docker.io/aravindsundaravadivelu/marine:1.0.2
#
# Build local Singularity SIF from this Dockerfile (no Docker Hub required):
#   singularity build marine.sif docker-daemon://aravindsundaravadivelu/marine:1.0.2
#
# To use the pre-built SIF directly (current approach):
#   Set params.marine_sif = '/path/to/marine_container.sif' in your nextflow.config

# note this refers to the marine-main repo. either copy that repo here or use the following command to refer to the marine-main repo directly:
# docker build --platform linux/amd64 \
#     -f /path/to/nf-core-ribostamp/containers/marine/Dockerfile \
#     -t aravindsundaravadivelu/marine:1.0.2 \
#     /home/tmhaxs421/brannanlab/tmhaxs421/MARINE_NextFlow/trial_pipeline/containers/MARINE-main

FROM condaforge/miniforge3:latest

LABEL tool="MARINE" \
      version="1.0.2" \
      description="Multi-Core Algorithm for Rapid Identification of Nucleotide Edits" \
      python="3.8.18"

# Copy environment spec and source code
COPY marine_environment_mamba_py38.yaml /opt/marine_environment_mamba_py38.yaml
COPY marine.py                          /opt/MARINE/marine.py
COPY src/                               /opt/MARINE/src/
COPY annotations/                       /opt/MARINE/annotations/

# Build conda environment
RUN set -e && \
    conda config --system --set ssl_verify false && \
    conda config --system --set channel_priority flexible && \
    mamba env create -f /opt/marine_environment_mamba_py38.yaml && \
    # Regenerate SONAME symlinks (some conda packages omit them)
    /sbin/ldconfig -n /opt/conda/envs/marine_environment/lib 2>/dev/null || true && \
    /sbin/ldconfig -n /opt/conda/lib 2>/dev/null || true && \
    # Pre-warm Numba JIT cache — scanpy uses @njit(cache=True); bake cache here
    # while filesystem is still writable (container squashfs is read-only at runtime)
    NUMBA_CACHE_DIR=/tmp/numba_marine_cache \
        /opt/conda/envs/marine_environment/bin/python -c "import scanpy" 2>&1 || true && \
    # Reduce image size
    mamba clean -afy && \
    find /opt/conda -type f -name "*.pyc" -delete

ENV PATH=/opt/conda/envs/marine_environment/bin:$PATH \
    LD_LIBRARY_PATH=/opt/conda/envs/marine_environment/lib:/opt/conda/lib \
    CONDA_DEFAULT_ENV=marine_environment \
    CONDA_PREFIX=/opt/conda/envs/marine_environment \
    # Numba and Matplotlib caches must be on a writable path at runtime
    NUMBA_CACHE_DIR=/tmp/numba_marine_cache \
    MPLCONFIGDIR=/tmp/matplotlib_cache

# Sanity check: marine.py --help must exit 0
RUN conda env list && \
    NUMBA_CACHE_DIR=/tmp/numba_marine_cache \
    python /opt/MARINE/marine.py --help

# Bundled annotations (pass to --annotation_bedfile_path):
#   /opt/MARINE/annotations/hg38_gencode.v35.annotation.genes.bed
#   /opt/MARINE/annotations/cellranger-GRCh38-3.0.0.annotation.genes.bed
#   /opt/MARINE/annotations/cellranger-GRCh37.annotation.genes.bed
#   /opt/MARINE/annotations/cellranger-mm10-3.0.0.annotation.genes.bed

ENTRYPOINT ["python", "/opt/MARINE/marine.py"]
