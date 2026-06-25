FROM condaforge/miniforge3:latest

# docker build --platform linux/amd64 -t aravindsundaravadivelu/ribostamp_utils:1.0.0 .
# docker push docker.io/aravindsundaravadivelu/ribostamp_utils:1.0.0

LABEL tool="ribostamp_utils" \
      version="1.0" \
      description="Python helper scripts for the ribostamp nf-core pipeline" \
      python="3.10"

# Disable SSL verification for corporate/hospital networks that use SSL inspection
# proxies presenting self-signed certificates (safe for internal build-time use).
RUN conda config --system --set ssl_verify false

# Single mamba install to minimize layers
# bedtools is a runtime dependency of pybedtools (binary must be on PATH)
RUN mamba install -y -c conda-forge -c bioconda \
        python=3.10 \
        pandas>=2.0 \
        numpy \
        matplotlib-base>=3.7 \
        pybedtools>=0.9 \
        bedtools \
        scanpy>=1.9 \
    && mamba clean -afy \
    && find /opt/conda -type f -name "*.pyc" -delete

# Pre-warm matplotlib and numba caches into a writable location that
# can be overridden at runtime.  MPLCONFIGDIR must be writable; /tmp is
# always writable inside containers (and in Nextflow task workdirs).
ENV MPLCONFIGDIR=/tmp/matplotlib_cache \
    NUMBA_CACHE_DIR=/tmp/numba_cache \
    PYTHONNOUSERSITE=1

# Sanity check — these are the imports used by the bin/ scripts
RUN python -c "import pandas, numpy, matplotlib, pybedtools, scanpy; print('ribostamp_utils OK')"
