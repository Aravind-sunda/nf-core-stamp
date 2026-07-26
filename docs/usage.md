# nf-core/stamp: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/stamp/usage](https://nf-co.re/stamp/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

nf-core/stamp analyses STAMP experiments, in which an RNA base editor deposits C-to-U (or A-to-I) edits on the transcripts it contacts. Which analysis you get is controlled by two things: `--mode`, which selects bulk or single-cell processing, and a set of `--run_*` toggles that switch individual edit callers on and off.

| Assay              | Parameters                                    | What you get                                             |
| ------------------ | --------------------------------------------- | -------------------------------------------------------- |
| Ribo-STAMP (bulk)  | `--mode bulk` (defaults)                      | Edits per gene, normalised to expression                  |
| RBP-STAMP (bulk)   | `--mode bulk --run_flare`                     | The above, plus FDR-scored edit clusters (binding sites)  |
| Single-cell STAMP  | `--mode sc`                                   | Per-cell-barcode edit calls, normalised to UMI counts     |

The bulk toggles are `--run_marine` (default `true`) and `--run_sailor` (default `true`). MARINE produces per-site edit calls and gene-level quantification; SAILOR produces confidence-ranked BED files. `--run_flare` builds on SAILOR's output and therefore requires `--run_sailor true`; it is bulk-only. Leave `--run_flare false` (the default) for Ribo-STAMP.

`--edit_type` must match your base editor: `C>T` for APOBEC1 (the default), `A>G` for ADAR.

## Samplesheet input

The pipeline takes a single comma-separated samplesheet via `--input`. There is one samplesheet format, and the **columns you fill in determine the start point**. Unused columns may be omitted entirely or left empty.

```bash
--input '[path to samplesheet file]'
```

Relative paths inside the samplesheet are resolved against the samplesheet's own directory, so a samplesheet and its data can be moved together.

### Accepted column combinations

Each row must match exactly one of the four combinations below. Any other combination is rejected with an explicit error naming the offending sample, so a typo fails immediately rather than silently changing what runs.

| Mode   | Start point | Required columns                        | Must be absent                        |
| ------ | ----------- | --------------------------------------- | ------------------------------------- |
| `bulk` | FASTQ       | `sample`, `fastq_1`, `library_type`      | `bam`, `fastq_dir`, `matrix_dir`      |
| `bulk` | BAM         | `sample`, `bam`, `library_type`           | `fastq_1`, `fastq_2`, `fastq_dir`, `matrix_dir` |
| `sc`   | FASTQ       | `sample`, `fastq_dir`                    | `bam`, `fastq_1`, `fastq_2`, `library_type`, `matrix_dir` |
| `sc`   | BAM         | `sample`, `bam`, `matrix_dir`             | `fastq_1`, `fastq_2`, `library_type`, `fastq_dir` |

For bulk FASTQ, `fastq_2` is required when `library_type` is `PE` and must be absent when it is `SE`.

| Column         | Description                                                                                                                       |
| -------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| `sample`       | Sample name. Must be unique. Spaces are converted to underscores (`_`).                                                            |
| `fastq_1`      | Path to R1 FASTQ. Must be gzipped, ending `.fastq.gz` or `.fq.gz`.                                                                 |
| `fastq_2`      | Path to R2 FASTQ, for paired-end bulk libraries only.                                                                              |
| `bam`          | Path to a coordinate-sorted, indexed BAM. The `.bai` must sit next to it.                                                          |
| `library_type` | `SE` or `PE`. Bulk only.                                                                                                          |
| `fastq_dir`    | Directory of 10x FASTQs following the `SAMPLE_S1_L001_R1_001.fastq.gz` convention. Single-cell FASTQ start only.                    |
| `matrix_dir`   | Path to a Cell Ranger `filtered_feature_bc_matrix` directory. Single-cell BAM start only; Cell Ranger generates it on FASTQ start. |

### Examples

Bulk, FASTQ start, mixed single- and paired-end:

```csv title="samplesheet_bulk_fastq.csv"
sample,fastq_1,fastq_2,library_type
CONTROL_REP1,ctrl_rep1_R1.fastq.gz,ctrl_rep1_R2.fastq.gz,PE
CONTROL_REP2,ctrl_rep2_R1.fastq.gz,,SE
DOX_REP1,dox_rep1_R1.fastq.gz,dox_rep1_R2.fastq.gz,PE
```

Bulk, BAM start:

```csv title="samplesheet_bulk_bam.csv"
sample,bam,library_type
CONTROL_REP1,ctrl_rep1.sorted.bam,SE
DOX_REP1,dox_rep1.sorted.bam,SE
```

Single-cell, FASTQ start (Cell Ranger runs):

```csv title="samplesheet_sc_fastq.csv"
sample,fastq_dir
STAMP_10X_1,/data/fastqs/STAMP_10X_1/
```

Single-cell, BAM start (Cell Ranger already run):

```csv title="samplesheet_sc_bam.csv"
sample,bam,matrix_dir
STAMP_10X_1,/data/cr/STAMP_10X_1/outs/possorted_genome_bam.bam,/data/cr/STAMP_10X_1/outs/filtered_feature_bc_matrix
```

Example samplesheets are bundled under [`assets/samplesheets/`](../assets/samplesheets/).

> [!NOTE]
> The pipeline does not merge multiple rows sharing a `sample` name. If you sequenced a library across several lanes, concatenate the FASTQs before running.

## Reference files

| Parameter      | Required for                          | Notes                                                                        |
| -------------- | ------------------------------------- | ---------------------------------------------------------------------------- |
| `--fasta`      | All modes                             | Genome FASTA. Used for `samtools calmd`, SAILOR and FLARE.                    |
| `--gtf`        | `--run_marine`, and FLARE region generation | Gene annotation, plain or gzipped.                                        |
| `--gene_bed`   | All modes                             | BED6 gene models, for RSeQC and edit annotation.                              |
| `--dbsnp_bed`  | All modes                             | dbSNP BED, for filtering known germline variants.                             |
| `--star_index` | Bulk FASTQ start                      | Pre-built STAR index directory.                                              |
| `--cellranger_ref` | Single-cell FASTQ start           | Pre-built Cell Ranger reference.                                             |
| `--genome`     | `--run_marine`                        | Selects a bundled metaPlotR genePred (`hg19`, `hg38`, `hg38_V44`, `mm10`, `mm39`). Override with `--genepred`. |

BAMs must carry `MD` tags for MARINE to detect edits. If they do not, the pipeline runs `samtools calmd` automatically, which needs `--fasta` to match the reference the BAM was aligned to.

## Running SAILOR and FLARE

SAILOR and FLARE are Snakemake workflows bundled under `assets/`. Unlike every other step, they are **not** wrapped in a Nextflow container: Snakemake launches its own Singularity containers, and nesting containers does not work. This means the host running these processes needs `snakemake` and `singularity` available. The pipeline is tested against Snakemake 9.13.4.

| Parameter                    | Purpose                                                                                       |
| ---------------------------- | --------------------------------------------------------------------------------------------- |
| `--sailor_snakemake_path`    | Path to the `snakemake` binary if it is not on `PATH`. Also used by FLARE.                     |
| `--sailor_singularity_cache` | Directory for Snakemake's Singularity images. Defaults to `$NXF_SINGULARITY_CACHEDIR`.         |
| `--sailor_singularity_bind`  | Comma-separated bind paths. Must cover the locations of `--fasta` and `--dbsnp_bed`.           |

### FLARE regions

FLARE scores edit clusters against a genome-wide set of sliding windows. Build these once with `--flare_regions` pointing at a folder produced by FLARE's `generate_regions.py`; if you omit it, the pipeline generates the folder from `--gtf` on the fly. The folder is large (roughly 8–10 GB for a human genome), so generating it once and reusing it across runs is worth the setup.

Cluster scoring is tuned with `--flare_fdr` (default `0.1`) and `--flare_max_merge_dist` (default `15` bp).

> [!NOTE]
> Only FLARE's cluster-identification mode is wired into this pipeline. FLARE's edit-fraction mode depends heavily on how you define your regions of interest, so run it yourself against the [upstream FLARE workflow](https://github.com/YeoLab/FLARE) if you need it.

## Running on an HPC without internet access

If your compute nodes cannot reach the internet, pre-download all containers into a local cache and point the pipeline at it. Follow [`RUNNING_ON_HPC.md`](../RUNNING_ON_HPC.md), then use `--marine_sif`, `--ribostamp_utils_sif` and `--editc_sif` to override the Docker Hub images with local `.sif` files.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/stamp \
    -profile docker \
    --mode bulk \
    --input ./samplesheet.csv \
    --outdir ./results \
    --fasta hg38.fa \
    --gtf gencode.v44.annotation.gtf.gz \
    --gene_bed hg38_genes.bed \
    --dbsnp_bed hg38_dbsnp.bed \
    --star_index star_index_hg38/ \
    --genome hg38 \
    --edit_type 'C>T'
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/stamp -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
mode: 'bulk'
input: './samplesheet.csv'
outdir: './results/'
genome: 'hg38'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/stamp
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/stamp releases page](https://github.com/nf-core/stamp/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
