# Running `nf-core/ribostamp` on an Offline HPC Cluster

A step-by-step guide for running the STAMP / Ribo-STAMP Nextflow pipeline on no-internet HPC cluster using Singularity/Apptainer.

> **Context.** This guide was written while testing `nf-core-ribostamp-dev`
> (v1.0.0dev) on an HPC cluster whose compute nodes have **no outbound internet
> access** and whose egress proxy **blocks the Seqera Wave registry**
> (`community.wave.seqera.io`, `community-cr-prod.seqera.io` → HTTP 403) while
> permitting `quay.io`, `docker.io`, and `depot.galaxyproject.org`. If your
> cluster has open egress, most of the offline workarounds below are unnecessary.

---

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Install Nextflow](#2-install-nextflow-conda)
3. [Get the pipeline and test data](#3-get-the-pipeline-and-test-data)
4. [Configure paths in the run script](#4-configure-paths-in-the-run-script)
5. [Install the nf-schema plugin offline](#5-install-the-nf-schema-plugin-offline)
6. [Install SAILOR's Snakemake dependency](#6-install-sailors-snakemake-dependency)
7. [Pre-pull blocked containers on a connected machine](#7-pre-pull-blocked-containers)
8. [Transfer images to the cluster cache](#8-transfer-images-to-the-cluster-cache)
9. [The offline run script](#9-the-offline-run-script)
10. [Troubleshooting reference](#10-troubleshooting-reference)

---

## 1. Prerequisites

| Requirement | Notes |
|---|---|
| HPC access | SLURM scheduler assumed; login node with `apptainer` module |
| Apptainer/Singularity | `apptainer/1.3.4` or `apptainer/1.4.2` via `module load` |
| Conda/Miniconda | Used for host-side tools (Nextflow, Snakemake) since the cluster ships no native Nextflow/Java |
| A machine with open internet | Needed to pre-pull Wave-blocked containers (a local Mac/Linux box works) |

> **Cluster note.** On this cluster there is no Singularity module, no
> system Java/OpenJDK, and no native Nextflow. Conda (`bioconda` +
> `conda-forge`) is the practical path for host-side dependencies.

---

## 2. Install Nextflow (conda)

First-time users should start here. Create a dedicated conda environment:

```bash
conda create -n nextflow -c bioconda -c conda-forge nextflow
conda activate nextflow

# verify
nextflow -version
java -version
```

---

## 3. Get the pipeline and test data

```bash
# Copy/clone the pipeline repo into your cluster workspace
cd /path/to/your/workspace
git clone <pipeline-repo-url> nf-core-ribostamp-dev

# Download the test data into the assets folder (per pipeline docs)
```

> ⚠️ **Edit the run script before launching.** The shipped
> `bin/run_pipeline_test.sh` contains hard-coded paths from the original
> author's environment (e.g. `/condo/brannanlab/...`) and calls `mamba`, which
> may not be on your `PATH`. First-time users **must** update these paths and
> the `main.nf` location for their own cluster, or you will see errors like:
>
> ```
> mamba: command not found
> mkdir: cannot create directory '/condo': Permission denied
> curl: (22) The requested URL returned error: 403
> Cannot find script file: /condo/brannanlab/.../main.nf
> ```

---

## 4. Configure paths in the run script

Set these environment variables at the top of your run script. They are
**required for offline operation**:

```bash
# Point NXF_HOME at the pipeline's local .nextflow dir (so the offline
# plugin is found — see step 5)
export NXF_HOME=/path/to/nf-core-ribostamp-dev/.nextflow

# Stops Nextflow's OWN network calls (config/plugin fetches).
# NOTE: this does NOT stop Singularity/Apptainer image pulls — those are
# a separate concern handled in steps 7–8.
export NXF_OFFLINE=true

# Singularity image cache — pre-pulled images live here
export NXF_SINGULARITY_CACHEDIR=/path/to/singularity_cache
mkdir -p "$NXF_SINGULARITY_CACHEDIR"
```

---

## 5. Install the nf-schema plugin offline

Because the cluster can't reach the internet, install the plugin on a
connected machine and copy it in:

```bash
# On a machine WITH internet:
nextflow plugin install nf-schema@2.5.1
```

Then transfer the plugin into the pipeline's plugin path on the cluster:

```
/path/to/nf-core-ribostamp-dev/.nextflow/plugins/nf-schema-2.5.1/
```

Make sure `NXF_HOME` (step 4) points at that `.nextflow` directory so
Nextflow discovers the plugin locally instead of trying to fetch it.

---

## 6. Install SAILOR's Snakemake dependency

SAILOR requires **Snakemake 8.x** (it uses `--use-singularity` /
`--singularity-*`, which were removed in Snakemake 9). Create a dedicated env:

```bash
conda update conda
conda install -c conda-forge mamba

conda create --name snakemake8
conda activate snakemake8
mamba install -c conda-forge -c bioconda 'snakemake>=8,<9'
```

Then point the run script at the **absolute** Snakemake path (so SAILOR execs
the right binary regardless of the active conda env):

```bash
SNAKEMAKE_PATH="/path/to/miniconda3/envs/snakemake8/bin/snakemake"
```

---

## 7. Pre-pull blocked containers

The HPC egress proxy returns **HTTP 403** for the Seqera Wave registry, so
these images can never be pulled from the compute nodes. Pull them **once** on
a machine with open egress, then copy them to the cluster cache.

### 7a. Set up Apptainer locally (macOS via Lima)

Apptainer has no native macOS build, so run it inside a Linux VM using
[Lima](https://apptainer.org/docs/admin/main/installation.html#mac):

```bash
# Install Homebrew (if needed), then Lima
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew install lima

# Start and enter the Apptainer VM
limactl start template://apptainer
limactl shell apptainer
```

### 7b. Prepare working + temp folders

The VM can run out of space during pulls, so set explicit temp dirs:

```bash
mkdir -p ~/ribostamp_images
cd ~/ribostamp_images

mkdir -p ~/ribostamp_images/tmp
export APPTAINER_TMPDIR=~/ribostamp_images/tmp
export TMPDIR=~/ribostamp_images/tmp
```

### 7c. Container inventory

| Image | Used by (module) | Status |
|---|---|---|
| `fastp:1.3.4--b75b637bf1c0f4b1` | `local/fastp` | Wave — blocked |
| `samtools_star:3d56ec4ef8fcee61` | `local/samtools_faidx` + `local/star_align` (same image, twice) | Wave — blocked |
| `fastqc:0.12.1--5cb1a2fa2f18c7c2` | MULTIQC/FASTQC config | Wave — blocked |
| `multiqc:1.34--db7c73dae76bc9e6` | MULTIQC config | Wave — blocked |
| `subread:2.1.1--cbb5cb85f59ac813` | `local/featurecounts` | Wave — blocked |
| `rseqc:5.0.4--6dbd0838c4d673ae` | `local/infer_strandedness` | Wave — blocked |

### 7d. Pull the images (amd64)

```bash
apptainer pull --arch amd64 multiqc-1.34.img \
  docker://community.wave.seqera.io/library/multiqc:1.34--db7c73dae76bc9e6

apptainer pull --arch amd64 fastqc-0.12.1.img \
  docker://community.wave.seqera.io/library/fastqc:0.12.1--5cb1a2fa2f18c7c2

apptainer pull --arch amd64 subread-2.1.1.img \
  docker://community.wave.seqera.io/library/subread:2.1.1--cbb5cb85f59ac813

apptainer pull --arch amd64 rseqc-5.0.4.img \
  docker://community.wave.seqera.io/library/rseqc:5.0.4--6dbd0838c4d673ae

# fastp and samtools_star (pull with their exact Wave names preserved)
apptainer pull --arch amd64 \
  community.wave.seqera.io-library-fastp-1.3.4--b75b637bf1c0f4b1.img \
  docker://community.wave.seqera.io/library/fastp:1.3.4--b75b637bf1c0f4b1

apptainer pull --arch amd64 \
  community.wave.seqera.io-library-samtools_star-3d56ec4ef8fcee61.img \
  docker://community.wave.seqera.io/library/samtools_star:3d56ec4ef8fcee61
```

> 💡 **Preserve exact filenames.** Nextflow resolves a cached image by a
> filename derived from the container URI. Keeping the exact
> `community.wave.seqera.io-library-<ref>.img` names lets Nextflow find them
> without a container override. For the images you rename (multiqc, fastqc,
> subread, rseqc), you must map them explicitly in `offline.config` (step 9).

---

## 8. Transfer images to the cluster cache

### 8a. Copy from the Lima VM to your local disk

In a **new** terminal on the host:

```bash
mkdir -p /path/to/local/nf-stamp/ribostamp_images
cd /path/to/local/nf-stamp/ribostamp_images

limactl copy apptainer:/home/lima.guest/ribostamp_images/multiqc-1.34.img .
limactl copy apptainer:/home/lima.guest/ribostamp_images/fastqc-0.12.1.img .
limactl copy apptainer:/home/lima.guest/ribostamp_images/subread-2.1.1.img .
limactl copy apptainer:/home/lima.guest/ribostamp_images/rseqc-5.0.4.img .
limactl copy "apptainer:/home/lima.guest/ribostamp_images/community.wave.seqera.io-library-fastp-1.3.4--b75b637bf1c0f4b1.img" .
limactl copy "apptainer:/home/lima.guest/ribostamp_images/community.wave.seqera.io-library-samtools_star-3d56ec4ef8fcee61.img" .
```

### 8b. Transfer to the cluster cache

```bash
rsync -avP /path/to/local/nf-stamp/ribostamp_images/*.img \
  user@cluster:/path/to/singularity_cache/
```

> 💡 `rsync` is preferred over `scp` for large image transfers (resumable,
> shows progress). On older macOS `rsync` (2.6.9-era), `--mkpath` is
> unavailable — create the destination directory on the cluster manually first.

### 8c. (Optional) Verify tool locations inside each image

Wave images vary by build toolchain and place their tool binary in different
locations. This inventory loop helps confirm where each tool lives (informs the
`PATH` fix in step 9):

```bash
cd /path/to/singularity_cache
for img in *.img *.simg; do
  [[ -f "$img" ]] || continue
  echo "=============================================================="
  echo "IMAGE: $img"
  apptainer exec "$img" sh -c '
    for d in /opt/wave/.pixi/envs/default/bin /opt/conda/bin /usr/local/bin /usr/bin; do
      [ -d "$d" ] || continue
      echo "  [$d]:"
      ls "$d" 2>/dev/null | grep -iE "fastqc|multiqc|fastp|samtools|star|featurecounts|subread|rseqc|infer_experiment|python" | sed "s/^/    /"
    done
  ' 2>/dev/null

  flavor="unknown/oras-native"
  apptainer exec "$img" sh -c '[ -d /opt/wave/.pixi ]' 2>/dev/null && flavor="pixi"
  apptainer exec "$img" sh -c '[ -d /opt/conda ]' 2>/dev/null && flavor="conda"
  echo "  FLAVOR: $flavor"
  echo ""
done
```

---

## 9. The offline run script

Save as `bin/run_test_offline.sh`. It regenerates `offline.config` on every run
so the script is the single source of truth, guards against missing images, and
loops over all four test profiles.

```bash
#!/usr/bin/env bash
set -euo pipefail

# ── Profiles to run — comment out any you want to skip ───────────────
PROFILES=(
  test_bulk_bam
  test_bulk_fastq
  test_sc_bam
  test_sc_fastq
)

# ── Environment ─────────────────────────────────────────────────────
export NXF_HOME=/path/to/nf-core-ribostamp-dev/.nextflow
# NXF_OFFLINE stops Nextflow's own network calls only — it does NOT stop
# Singularity image pulls. Images must already exist in the cache (step 8).
export NXF_OFFLINE=true

export NXF_SINGULARITY_CACHEDIR=/path/to/singularity_cache
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

PIPELINE_DIR="/path/to/nf-core-ribostamp-dev"
export OUTPUT_DIR=/path/to/nfcore-ribostamp-run
mkdir -p "$OUTPUT_DIR"

# ── System-specific paths — edit for your cluster ───────────────────
# SAILOR needs Snakemake 8.x. Absolute path so it runs regardless of the
# active conda env.
SNAKEMAKE_PATH="/path/to/miniconda3/envs/snakemake8/bin/snakemake"

# Wave "pixi" and "conda" container builds put the tool in different dirs.
# Nextflow's `apptainer exec` bypasses the container entrypoint that would
# normally add these to PATH, causing exit-127 "command not found". Prepending
# BOTH dirs (plus standard system dirs) makes tool resolution deterministic for
# every image flavor. Absent dirs are simply skipped.
WAVE_PIXI_BIN="/opt/wave/.pixi/envs/default/bin"
WAVE_CONDA_BIN="/opt/conda/bin"
CONTAINER_PATH="${WAVE_PIXI_BIN}:${WAVE_CONDA_BIN}:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

# ── Offline config (regenerated every run) ──────────────────────────
OFFLINE_CONFIG="$(dirname "$0")/offline.config"
cat > "$OFFLINE_CONFIG" <<EOF
// AUTO-GENERATED by $(basename "$0") — do not edit by hand.

// Prevent FileAlreadyExistsException on re-runs into the same --outdir
report.overwrite   = true
timeline.overwrite = true
trace.overwrite    = true
dag.overwrite      = true

// Silence "Access to undefined parameter" warnings
params.marine_sif          = params.containsKey('marine_sif')          ? params.marine_sif          : null
params.ribostamp_utils_sif = params.containsKey('ribostamp_utils_sif') ? params.ribostamp_utils_sif : null

// Container overrides: bind every Wave/Seqera-blocked process to a local .img.
// When a container directive is an existing file path, Apptainer uses it
// directly and never contacts the blocked registry.
params.sing_cache     = '${NXF_SINGULARITY_CACHEDIR}'
params.container_path = '${CONTAINER_PATH}'

process {
  withName: 'FASTP' {
    container        = "\${params.sing_cache}/community.wave.seqera.io-library-fastp-1.3.4--b75b637bf1c0f4b1.img"
    containerOptions = "--env PATH=\${params.container_path}"
  }
  withName: 'SAMTOOLS_FAIDX' {
    container        = "\${params.sing_cache}/community.wave.seqera.io-library-samtools_star-3d56ec4ef8fcee61.img"
    containerOptions = "--env PATH=\${params.container_path}"
  }
  withName: 'STAR_ALIGN' {
    container        = "\${params.sing_cache}/community.wave.seqera.io-library-samtools_star-3d56ec4ef8fcee61.img"
    containerOptions = "--env PATH=\${params.container_path}"
  }
  withName: 'FEATURECOUNTS' {
    container        = "\${params.sing_cache}/subread-2.1.1.img"
    containerOptions = "--env PATH=\${params.container_path}"
  }
  withName: 'INFER_STRANDEDNESS' {
    container        = "\${params.sing_cache}/rseqc-5.0.4.img"
    containerOptions = "--env PATH=\${params.container_path}"
  }
  withName: 'FASTQC' {
    container        = "\${params.sing_cache}/fastqc-0.12.1.img"
    containerOptions = "--env PATH=\${params.container_path}"
  }
  withName: 'MULTIQC' {
    container        = "\${params.sing_cache}/multiqc-1.34.img"
    containerOptions = "--env PATH=\${params.container_path}"
  }
}
EOF
echo "Wrote $OFFLINE_CONFIG"

# ── Guard: refuse to launch if any required image is missing ────────
REQUIRED_IMAGES=(
  "community.wave.seqera.io-library-fastp-1.3.4--b75b637bf1c0f4b1.img"
  "community.wave.seqera.io-library-samtools_star-3d56ec4ef8fcee61.img"
  "fastqc-0.12.1.img"
  "multiqc-1.34.img"
  "subread-2.1.1.img"
  "rseqc-5.0.4.img"
)
missing=0
for img in "${REQUIRED_IMAGES[@]}"; do
  if [[ ! -e "$NXF_SINGULARITY_CACHEDIR/$img" ]]; then
    echo "MISSING cached image: $img"
    missing=1
  fi
done
if [[ "$missing" -eq 1 ]]; then
  echo ""
  echo "One or more blocked images are not in the cache. Pre-pull them from a"
  echo "connected host (step 7), copy into \$NXF_SINGULARITY_CACHEDIR (step 8),"
  echo "then re-run. To bypass (e.g. you overrode containers to quay.io), set"
  echo "SKIP_IMAGE_CHECK=1 before running."
  [[ "${SKIP_IMAGE_CHECK:-0}" != "1" ]] && exit 1
fi
echo "All required container images present."

# ── Loop ────────────────────────────────────────────────────────────
# Run this script FROM the same directory every time (work/ lives in
# scripts/), otherwise -resume won't find the cache.
for PROFILE in "${PROFILES[@]}"; do
  echo "======================================================"
  echo " ribostamp test: ${PROFILE}   started: $(date)"
  echo "======================================================"

  # -c "$OFFLINE_CONFIG" passed LAST so container overrides win.
  nextflow run "$PIPELINE_DIR/main.nf" \
    -profile "${PROFILE},singularity" \
    -c "$OFFLINE_CONFIG" \
    --outdir "${OUTPUT_DIR}/results_off_${PROFILE}" \
    --sailor_snakemake_path "$SNAKEMAKE_PATH" \
    -resume

  # Optional: supply local .sif files to skip Docker Hub pulls:
  #   --marine_sif          "${NXF_SINGULARITY_CACHEDIR}/docker.io-aravindsundaravadivelu-marine-1.0.2.img"
  #   --ribostamp_utils_sif "${NXF_SINGULARITY_CACHEDIR}/docker.io-aravindsundaravadivelu-ribostamp_utils-1.0.0.img"

  echo "======================================================"
  echo " ribostamp test: ${PROFILE} — DONE ($(date))"
  echo "======================================================"
  echo ""
done
```

---

## 10. Troubleshooting reference

| Symptom | Cause | Fix |
|---|---|---|
| `curl: (22) ... error: 403` on `*.wave.seqera.io` | Egress proxy blocks the Seqera Wave registry | Pre-pull images on a connected host and cache them (steps 7–8) |
| `mamba: command not found` / hard-coded `/condo/...` paths | Shipped test script has the original author's paths | Edit all paths + `main.nf` location before launching (step 3–4) |
| Plugin/config fetch fails on compute node | No internet; Nextflow tries a remote fetch | `export NXF_OFFLINE=true`; install nf-schema offline (steps 4–5) |
| Image still pulled despite `NXF_OFFLINE=true` | `NXF_OFFLINE` only stops **Nextflow's** calls, not Apptainer pulls | Ensure all images are pre-cached with exact filenames (steps 7–8) |
| Exit code 127 / "command not found" inside a container | Nextflow's `apptainer exec` bypasses the entrypoint that sets `PATH` | Inject `PATH` via `containerOptions` covering pixi + conda bin dirs (step 9) |
| `FileAlreadyExistsException` on report HTML during re-runs | Report/timeline/trace/dag files already exist | `*.overwrite = true` in `offline.config` (step 9) |
| Snakemake errors from SAILOR | Snakemake 9 removed `--use-singularity` | Install Snakemake 8.x and point `--sailor_snakemake_path` at it (step 6) |
| `-resume` doesn't reuse cache | Script launched from a different directory | Always run from the same directory (Nextflow logs/`work/` live in `scripts/` on this cluster) |

---

## Recommended durable fix

The offline workarounds above are robust but maintenance-heavy — **every**
future pipeline that uses Wave-hosted containers will hit the same 403. The
durable solution is to **request that HPC admins allowlist**
`community.wave.seqera.io` and `community-cr-prod.seqera.io` on the egress
proxy. Once allowlisted, the pre-pull/transfer steps (7–8) and the container
overrides in `offline.config` become unnecessary.
