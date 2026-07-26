# nf-core/stamp: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory. Directories are numbered in the order the steps run, and which ones appear depends on `--mode` and the `--run_*` toggles.

## Pipeline overview

**Bulk mode** (`--mode bulk`):

- [FastQC](#fastqc) - Raw read QC
- [fastp](#fastp) - Adapter and quality trimming
- [STAR](#star) - Alignment
- [Strandedness](#strandedness) - Library strandedness inference
- [MARINE](#marine-bulk) - Edit site calling
- [SAILOR](#sailor) - Confidence-ranked edit calls
- [FLARE](#flare) - Edit-cluster identification for RBP-STAMP
- [featureCounts](#featurecounts) - Gene-level quantification
- [Filtered and normalised edits](#filtered-and-normalised-edits) - Final edit calls
- [metaPlotR](#metaplotr) - Metagene distribution of edit sites

**Single-cell mode** (`--mode sc`):

- [Cell Ranger](#cell-ranger) - Alignment and cell calling
- [MARINE](#marine-single-cell) - Per-barcode edit site calling
- [Filtered and normalised edits (single-cell)](#filtered-and-normalised-edits-single-cell)

**All modes**:

- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences.

### fastp

<details markdown="1">
<summary>Output files</summary>

- `02_fastp/`
  - `*.fastp.html`, `*.fastp.json`: Trimming reports.
  - `*.trimmed.fastq.gz`: Trimmed reads passed to STAR.

</details>

[fastp](https://github.com/OpenGene/fastp) removes adapters and low-quality bases. Only produced for FASTQ-start samples.

### STAR

<details markdown="1">
<summary>Output files</summary>

- `03_star/`
  - `*.Aligned.sortedByCoord.out.bam(.bai)`: Coordinate-sorted, indexed alignments.
  - `*.Log.final.out`, `*.Log.out`, `*.Log.progress.out`: Alignment logs; the final log is parsed by MultiQC.

</details>

[STAR](https://github.com/alexdobin/STAR) aligns reads to the genome. BAMs are written with `MD` tags, which MARINE requires to detect edits.

### Strandedness

<details markdown="1">
<summary>Output files</summary>

- `04_strandedness/`
  - `*_strandedness.txt`: Raw `infer_experiment.py` output.
  - `*_strand_code.txt`: The inferred code (`0` unstranded, `1` forward, `2` reverse).

</details>

[RSeQC](https://rseqc.sourceforge.net/) infers library strandedness, which determines how edits are assigned to strands downstream. Supply `--strandedness` to skip inference.

### MARINE (bulk)

<details markdown="1">
<summary>Output files</summary>

- `04_marine/<sample>/`
  - `final_filtered_site_info.tsv`: All called edit sites.
  - `final_filtered_site_info_annotated.tsv`: The same sites annotated against `--gene_bed` (stranded libraries only).

</details>

[MARINE](https://github.com/YeoLab/MARINE) reads `MD` tags to call nucleotide edits of the type given by `--edit_type`. These are raw calls; filtering happens later. BAMs lacking `MD` tags are passed through `samtools calmd` first.

### SAILOR

<details markdown="1">
<summary>Output files</summary>

- `04_sailor/sailor_output/`
  - `*.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed`: Confidence-ranked edit sites, one BED per sample.
  - `8_bw_and_bam/`: Per-sample coverage bigWigs and filtered BAMs, consumed by FLARE.
- `04_sailor/metaplotr/<sample>/`
  - `*.dist.measures.txt`: Metagene coordinates of the ranked edit sites.

</details>

[SAILOR](https://github.com/YeoLab/FLARE) scores each candidate edit site with a beta-binomial confidence value. The ranked BED is the input to FLARE. Produced when `--run_sailor true` (the default).

### FLARE

<details markdown="1">
<summary>Output files</summary>

- `05_flare/flare_output/FLARE/`
  - `<sample>_merged_sorted_peaks.*.scored.tsv`: Edit clusters with FDR scores. **This is the primary RBP-STAMP result** — each row is a candidate RBP binding site.
- `05_flare/flare_regions/`
  - Sliding-window region files, generated from `--gtf` when `--flare_regions` is not supplied.

</details>

[FLARE](https://github.com/YeoLab/FLARE) groups SAILOR's edit sites into clusters and tests each against a background model, reporting those below `--flare_fdr`. Adjacent significant windows within `--flare_max_merge_dist` bp are merged. Produced only when `--run_flare true`.

Only cluster-identification mode runs here. FLARE's edit-fraction mode is left to the user; see the [FLARE documentation](https://github.com/YeoLab/FLARE).

### featureCounts

<details markdown="1">
<summary>Output files</summary>

- `05_featurecounts/`
  - `*.featurecounts.txt`: Per-sample gene counts.
  - `*.featurecounts.txt.summary`: Assignment summary, parsed by MultiQC.
  - `counts_matrix_combined.tsv`: All samples merged into one matrix.

</details>

[featureCounts](https://subread.sourceforge.net/) quantifies reads per gene. The merged matrix is used to normalise edit counts to expression, so that a heavily edited gene is not simply a heavily expressed one.

### Filtered and normalised edits

<details markdown="1">
<summary>Output files</summary>

- `06_filter_normalize/<sample>/filtered/`
  - `filtered_edits.tsv`: Edit sites surviving all enabled filters.
- `06_filter_normalize/<sample>/normalized/bedgraphs/`
  - `<sample>.<edit_type>.edit_fraction.bed`: BED6 of edit fractions, used by metaPlotR.
  - Per-gene edit rates normalised to expression.

</details>

Edits pass through four filters, each individually switchable: multiallelic sites (`--filter_bulk_multiallelic`), dbSNP overlap (`--filter_bulk_dbsnp`), editing fraction above `--max_frac` (`--filter_bulk_max_frac`), and sites without gene annotation (`--filter_bulk_unannotated`). Surviving sites are then normalised against the featureCounts matrix.

**For Ribo-STAMP these are the primary results**: the normalised edit rate per gene is the readout of ribosome association.

### metaPlotR

<details markdown="1">
<summary>Output files</summary>

- `07_metaplotr/<sample>/`
  - `*.dist.measures.txt`: Position of each edit site in normalised transcript coordinates (5' UTR / CDS / 3' UTR).

</details>

[metaPlotR](https://github.com/olarerin/metaPlotR) maps edit sites onto a metagene so their transcript-level distribution can be plotted. The genePred annotation is chosen from `--genome` or supplied with `--genepred`.

### Cell Ranger

<details markdown="1">
<summary>Output files</summary>

- `01_cellranger/<sample>/`
  - `outs/possorted_genome_bam.bam(.bai)`: Aligned reads tagged with cell barcodes.
  - `outs/filtered_feature_bc_matrix/`: Cell-by-gene count matrix.
  - `outs/web_summary.html`: Cell Ranger QC report.

</details>

[Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger) aligns 10x reads and calls cells. Only run for single-cell FASTQ-start samples; BAM-start samples supply these files directly.

> [!NOTE]
> Cell Ranger requires a 10x Genomics licence, which you must accept before use.

### MARINE (single-cell)

<details markdown="1">
<summary>Output files</summary>

- `02_marine_sc/<sample>/marine_output/`
  - Per-barcode edit site calls, keyed by the BAM tag given in `--barcode_tag`.

</details>

### Filtered and normalised edits (single-cell)

<details markdown="1">
<summary>Output files</summary>

- `03_filter_sc/<sample>/`
  - `filtered_edits.tsv`: Edit sites surviving all enabled filters.
- `04_normalize_sc/<sample>/`
  - `normalized_edits.tsv`: Edit counts normalised to per-cell UMI totals. **This is the primary single-cell result.**

</details>

Five filters apply here, each switchable: multiple conversion types per barcode (`--filter_sc_multi_conversion`), dbSNP overlap (`--filter_sc_dbsnp`), fewer than `--min_count` edited reads (`--filter_sc_min_count`), editing fraction above `--max_frac` (`--filter_sc_max_frac`), and unannotated sites (`--filter_sc_unannotated`).

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) collates QC into a single report, alongside the software versions used in the run. What it contains depends on the route taken:

| Route             | Sections in the report                                                  |
| ----------------- | ----------------------------------------------------------------------- |
| Bulk, FASTQ start | FastQC, fastp, STAR, RSeQC (strandedness), featureCounts                 |
| Bulk, BAM start   | RSeQC, featureCounts                                                     |
| `sc`, FASTQ start | Cell Ranger                                                              |
| `sc`, BAM start   | Software versions only                                                   |

The edit-calling steps — MARINE, SAILOR, FLARE — and the filter, normalise and metaPlotR steps produce no MultiQC-parsable output, so they do not appear in the report.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
