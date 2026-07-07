#!/bin/bash

#SBATCH --job-name=test_dataset_flare
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your@email.com  # ← change to your email
#SBATCH --output=slurm_%u_%x_%j.log


# make_test_dataset_flare.sh
# Builds a small bulk RBP-STAMP (RBFOX2-BAMR) test dataset for exercising the
# ribostamp SAILOR → FLARE path start to finish, in both BAM-start and
# FASTQ-start modes. Subsamples the source BAMs to chr1, then derives matching
# FASTQs. All chr1 reference/annotation/dbSNP/STAR-index files are REUSED from
# assets/test_data/refs/ (built by bin/make_test_dataset.sh) — this script only
# produces the RBFOX2 BAMs, FASTQs, and samplesheets.
#
# Companion configs:
#   conf/test_bulk_bam_flare.config
#   conf/test_bulk_fastq_flare.config
#
# Usage:  sbatch bin/make_test_dataset_flare.sh
#         # or, for a quick local build:  THREADS=8 bash bin/make_test_dataset_flare.sh
#
# Requirements: samtools (>=1.15) on PATH.
module load samtools 2>/dev/null || true

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION — edit these paths before running
# ─────────────────────────────────────────────────────────────────────────────

# -- Source bulk BAMs (sorted; indexed automatically if .bai missing) ---------
#    control = -DOX (APOBEC fusion not induced), exp = +DOX (editing induced)
BULK_CTRL_BAM="/home/tmhaxs421/brannanlab/Vrutant/bulk-RNAseq/Project1/bam_index/293XTBAMRRBFOX2-DOX1_S1_merge_trimAligned.sorted.bam"
BULK_EXP_BAM="/home/tmhaxs421/brannanlab/Vrutant/bulk-RNAseq/Project1/bam_index/293XTBAMRRBFOX2plusDOX1_S7_merge_trimAligned.sorted.bam"

# -- Sample names (used in samplesheet and output file names) -----------------
SAMPLE_CTRL="minusdox"
SAMPLE_EXP="plusdox"

# -- Subsampling settings -----------------------------------------------------
OUTDIR="../assets/test_data/flare"   # New RBFOX2 data (refs are shared, see below)
REFDIR="../assets/test_data/refs"    # Existing chr1 refs (reused, not recreated)
CHROM="chr1"
# Target mapped reads per sample. Higher than the MARINE-only test (150k) so that
# editing depth survives for SAILOR to call C>T edits and FLARE to cluster them.
BULK_TARGET=500000
SEED=42                # Random seed — keep constant so datasets are reproducible
THREADS="${THREADS:-16}"   # samtools thread count (override: THREADS=8 bash ...)

# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

chr_mapped_count() { samtools idxstats "$1" | awk -v c="$CHROM" '$1==c { print $3 }'; }

detect_lib() {
    local n
    n=$(samtools view -f 1 -c --threads "$THREADS" "$1" 2>/dev/null || echo 0)
    [[ "$n" -gt 0 ]] && echo "PE" || echo "SE"
}

ensure_index() {
    local bam="$1"
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
        log "  Indexing $(basename "$bam")..."
        samtools index --threads "$THREADS" "$bam"
    fi
}

# Extract CHROM reads from INPUT, subsample to TARGET, write to OUTPUT
subset_bam() {
    local input="$1" output="$2" target="$3" label="$4"
    log "--- ${label} ---"
    log "  Input:  $input"
    ensure_index "$input"

    local total
    total=$(chr_mapped_count "$input")
    if [[ -z "$total" || "$total" -eq 0 ]]; then
        die "No reads found on '$CHROM' in $(basename "$input").
  Check the chromosome name — run: samtools idxstats $input | head"
    fi
    log "  Reads on ${CHROM}: ${total}"

    if [[ "$total" -le "$target" ]]; then
        log "  Total (${total}) <= target (${target}) — keeping all ${CHROM} reads"
        samtools view -b --threads "$THREADS" "$input" "$CHROM" > "$output"
    else
        local frac frac_digits s_arg
        frac=$(awk -v t="$target" -v n="$total" 'BEGIN { printf "%.4f", t/n }')
        frac_digits="${frac#0.}"
        s_arg="${SEED}.${frac_digits}"
        log "  Subsampling ${frac} of reads (samtools -s ${s_arg})"
        samtools view -b --threads "$THREADS" "$input" "$CHROM" \
            | samtools view -b -s "$s_arg" --threads "$THREADS" \
            > "$output"
    fi
    samtools index "$output"
    local final
    final=$(samtools flagstat "$output" | awk '/mapped \(/{print $1; exit}')
    log "  Output: $output  (final mapped reads: ${final})"
    echo ""
}

bam_to_fastq() {
    local bam="$1" out_prefix="$2" lib="$3" label="$4"
    log "--- ${label} ---"
    if [[ "$lib" == "PE" ]]; then
        samtools fastq --threads "$THREADS" \
            -1 "${out_prefix}_R1.fastq.gz" -2 "${out_prefix}_R2.fastq.gz" \
            -0 /dev/null -s /dev/null "$bam"
        log "  R1: ${out_prefix}_R1.fastq.gz   R2: ${out_prefix}_R2.fastq.gz"
    else
        samtools fastq --threads "$THREADS" "$bam" | gzip > "${out_prefix}_R1.fastq.gz"
        log "  R1: ${out_prefix}_R1.fastq.gz"
    fi
    echo ""
}

validate_config() {
    local missing=0
    check() { [[ -e "$1" ]] || { echo "  Not found: $2 = $1" >&2; missing=1; }; }
    check "$BULK_CTRL_BAM"  "BULK_CTRL_BAM"
    check "$BULK_EXP_BAM"   "BULK_EXP_BAM"
    # Shared chr1 refs (created by make_test_dataset.sh) — must already exist
    check "${REFDIR}/hg38_${CHROM}.fa"                "refs/hg38_${CHROM}.fa"
    check "${REFDIR}/${CHROM}_hg38_genes.gtf.gz"      "refs/${CHROM}_hg38_genes.gtf.gz"
    check "${REFDIR}/${CHROM}_hg38_genes.bed.gz"      "refs/${CHROM}_hg38_genes.bed.gz"
    check "${REFDIR}/${CHROM}_hg38_dbsnp.bed"         "refs/${CHROM}_hg38_dbsnp.bed"
    check "${REFDIR}/star_index_hg38_${CHROM}/SA"     "refs/star_index_hg38_${CHROM}"
    if [[ "$missing" -eq 1 ]]; then
        die "Missing inputs. Run bin/make_test_dataset.sh first to build the chr1 refs."
    fi
    command -v samtools >/dev/null || die "samtools not on PATH"
}

# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

log "ribostamp FLARE (RBP-STAMP) test dataset builder"
log "Chromosome: ${CHROM} | Bulk target: ${BULK_TARGET} reads/sample | Threads: ${THREADS}"
log "Reusing shared chr1 refs from: ${REFDIR}"
echo ""

validate_config

mkdir -p "${OUTDIR}/bulk_bams" "${OUTDIR}/bulk_fastqs" "${OUTDIR}/samplesheets"

# ── Subset + subsample BAMs ──────────────────────────────────────────────────
log "====== BULK BAMs → ${CHROM} ======"
subset_bam "$BULK_CTRL_BAM" "${OUTDIR}/bulk_bams/${SAMPLE_CTRL}_${CHROM}.bam" "$BULK_TARGET" "Control (-DOX)"
subset_bam "$BULK_EXP_BAM"  "${OUTDIR}/bulk_bams/${SAMPLE_EXP}_${CHROM}.bam"  "$BULK_TARGET" "Experimental (+DOX)"

BULK_LIB=$(detect_lib "${OUTDIR}/bulk_bams/${SAMPLE_CTRL}_${CHROM}.bam")
log "Detected bulk library type: ${BULK_LIB}"
echo ""

# ── BAM-start samplesheet ────────────────────────────────────────────────────
BAM_SHEET="${OUTDIR}/samplesheets/samplesheet_bulk_bam_flare.csv"
{
    printf "sample,bam,library_type\n"
    printf "%s,%s,%s\n" "$SAMPLE_CTRL" "../bulk_bams/${SAMPLE_CTRL}_${CHROM}.bam" "$BULK_LIB"
    printf "%s,%s,%s\n" "$SAMPLE_EXP"  "../bulk_bams/${SAMPLE_EXP}_${CHROM}.bam"  "$BULK_LIB"
} > "$BAM_SHEET"
log "BAM-start samplesheet: ${BAM_SHEET}"
echo ""

# ── FASTQ-start (derive FASTQ from the subset BAMs → identical reads) ─────────
log "====== BULK ${CHROM} BAM → FASTQ ======"
bam_to_fastq "${OUTDIR}/bulk_bams/${SAMPLE_CTRL}_${CHROM}.bam" "${OUTDIR}/bulk_fastqs/${SAMPLE_CTRL}_${CHROM}" "$BULK_LIB" "Control → FASTQ"
bam_to_fastq "${OUTDIR}/bulk_bams/${SAMPLE_EXP}_${CHROM}.bam"  "${OUTDIR}/bulk_fastqs/${SAMPLE_EXP}_${CHROM}"  "$BULK_LIB" "Experimental → FASTQ"

FASTQ_SHEET="${OUTDIR}/samplesheets/samplesheet_bulk_fastq_flare.csv"
{
    printf "sample,fastq_1,fastq_2,library_type\n"
    if [[ "$BULK_LIB" == "PE" ]]; then
        printf "%s,%s,%s,PE\n" "$SAMPLE_CTRL" "../bulk_fastqs/${SAMPLE_CTRL}_${CHROM}_R1.fastq.gz" "../bulk_fastqs/${SAMPLE_CTRL}_${CHROM}_R2.fastq.gz"
        printf "%s,%s,%s,PE\n" "$SAMPLE_EXP"  "../bulk_fastqs/${SAMPLE_EXP}_${CHROM}_R1.fastq.gz"  "../bulk_fastqs/${SAMPLE_EXP}_${CHROM}_R2.fastq.gz"
    else
        printf "%s,%s,,SE\n" "$SAMPLE_CTRL" "../bulk_fastqs/${SAMPLE_CTRL}_${CHROM}_R1.fastq.gz"
        printf "%s,%s,,SE\n" "$SAMPLE_EXP"  "../bulk_fastqs/${SAMPLE_EXP}_${CHROM}_R1.fastq.gz"
    fi
} > "$FASTQ_SHEET"
log "FASTQ-start samplesheet: ${FASTQ_SHEET}"
echo ""

# ── Summary ──────────────────────────────────────────────────────────────────
log "====== SUMMARY ======"
echo ""
find "${OUTDIR}" -type f | sort | while read -r f; do
    printf "  %-64s %s\n" "$f" "$(du -sh "$f" 2>/dev/null | cut -f1)"
done
echo ""
log "New FLARE test-data size: $(du -sh "$OUTDIR" | cut -f1)  (refs are shared from ${REFDIR})"
echo ""
log "Detected library type: ${BULK_LIB}  (edit_type is C>T — RBFOX2-BAMR / APOBEC1)"
echo ""
log "Run the pipeline tests with:"
echo "    nextflow run . -profile singularity,test_bulk_bam_flare   --outdir results_bam_flare"
echo "    nextflow run . -profile singularity,test_bulk_fastq_flare --outdir results_fastq_flare"
echo ""
log "Done."
