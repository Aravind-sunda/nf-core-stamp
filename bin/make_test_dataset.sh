#!/bin/bash

#SBATCH --job-name=test_dataset
#SBATCH --nodes=1
#SBATCH --partition=defq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=0G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your@email.com  # ← change to your email
#SBATCH --output=slurm_%u_%x_%j.log


# make_test_dataset.sh
# Subsamples bulk and single-cell BAMs to a single chromosome for ribostamp
# test dataset creation. Produces samplesheets ready for conf/test.config.
#
# Usage: bash make_test_dataset.sh
#
# Requirements: samtools (>=1.15) and STAR on PATH.
# ── Adjust these module-load lines for your cluster (or remove if tools are
#    already on PATH via conda/module system):
module load samtools
module load star/2.7.10b
module load cellranger/9.0.1

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION — edit these paths before running
# ─────────────────────────────────────────────────────────────────────────────

# -- Bulk BAMs (sorted + indexed) ---------------------------------------------
BULK_CTRL_BAM="/home/tmhaxs421/brannanlab/tmhaxs421/TERT/test/03_star/293xt_BAMR_RPS2_DOX_rep1.Aligned.sortedByCoord.out.bam"   # e.g. /condo/brannanlab/data/bulk/ctrl_rep1.bam
BULK_EXP_BAM="/home/tmhaxs421/brannanlab/tmhaxs421/TERT/test/03_star/293xt_BAMR_RPS2_plusDOX_rep1.Aligned.sortedByCoord.out.bam"    # e.g. /condo/brannanlab/data/bulk/exp_rep1.bam

# -- Single-cell BAMs (CellRanger possorted_genome_bam.bam) -------------------
SC_CTRL_BAM="/home/tmhaxs421/brannanlab/tmhaxs421/TCA_timepoint/data/t1_mESC_wt_a_S1/possorted_genome_bam.bam"     # e.g. /condo/brannanlab/data/sc/ctrl/outs/possorted_genome_bam.bam
SC_EXP_BAM="/home/tmhaxs421/brannanlab/tmhaxs421/TCA_timepoint/data/t1_shCTRL/outs/possorted_genome_bam.bam"      # e.g. /condo/brannanlab/data/sc/exp/outs/possorted_genome_bam.bam
# Note: filtered_feature_bc_matrix is generated automatically from the subsampled
# sc BAMs using CB/UB/GX tags — no need to supply an existing matrix directory.

# -- Sample names (used in samplesheet and output file names) -----------------
SAMPLE_CTRL="ctrl"
SAMPLE_EXP="exp"

# -- Subsampling settings -----------------------------------------------------
OUTDIR="../assets/test_data"  # Root output directory (created if absent)
CHROM="chr1"        # Chromosome to subset — must match BAM header naming

# Target mapped read counts per sample.
# Bulk 150k:  ~15-30 MB BAM, MARINE finishes in ~10-20 min on 4 CPUs.
# SC   500k:  spread across thousands of cells; gives ~50-100 reads/cell,
#             enough for MARINE SC to call editing in the most-covered cells.
BULK_TARGET=150000
SC_TARGET=500000
SEED=42       # Random seed — keep constant so datasets are reproducible
THREADS=36     # samtools / STAR thread count

# -- hg38 references for FASTQ-start bulk test (chr1 STAR index) --------------
HG38_FASTA="/home/tmhaxs421/brannanlab/tmhaxs421/reference/cellranger_reference/custom_cellranger/GRCh38-2024-A_mruby_gfp/genome.fa"
HG38_GTF="/home/tmhaxs421/brannanlab/tmhaxs421/reference/cellranger_reference/custom_cellranger/GRCh38-2024-A_mruby_gfp/genes.gtf"
HG38_GENE_BED="/condo/brannanlab/tmhaxs421/reference/cellranger_reference/custom_cellranger/GRCh38-2024-A_mruby_gfp/genes_mruby_gfp.bed"
STAR_INDEX_DIR="${OUTDIR}/refs/star_index_hg38_chr1"   # built once, reused

# -- hg38 dbSNP for bulk BAM/FASTQ-start tests --------------------------------
HG38_DBSNP_BED3="/home/tmhaxs421/brannanlab/tmhaxs421/reference/dbSNP/dbsnp-hg38/hg38_bed/hg38_dbsnp_combined.bed3"

# -- mm10 dbSNP for SC BAM-start test -----------------------------------------
MM10_DBSNP_BED3="/home/tmhaxs421/brannanlab/tmhaxs421/reference/dbSNP/dbsnp142-mm10/mm10_dbsnp_combined.bed3"

# -- mm10 references for SC annotation and chr1-only CellRanger reference ------
MM10_FASTA="/home/tmhaxs421/brannanlab/tmhaxs421/reference/cellranger_reference/custom_cellranger/mRuby3_mm10/mRUBY3_mm10/fasta/genome.fa"
MM10_GTF="/home/tmhaxs421/brannanlab/tmhaxs421/reference/cellranger_reference/custom_cellranger/mRuby3_mm10/mRUBY3_mm10/genes/genes.gtf.gz"
MM10_GENE_BED="/home/tmhaxs421/brannanlab/tmhaxs421/reference/cellranger_reference/custom_cellranger/mRuby3_mm10/mRUBY3_mm10/genes/genes.bed"
CELLRANGER_PATH="cellranger"    # set to full path if cellranger is not on PATH

# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

log() { echo "[$(date '+%H:%M:%S')] $*"; }

die() { echo "ERROR: $*" >&2; exit 1; }

# Count mapped reads on CHROM in a BAM (requires index)
chr_mapped_count() {
    samtools idxstats "$1" | awk -v c="$CHROM" '$1==c { print $3 }'
}

# Detect library type from BAM flags: PE or SE
detect_lib() {
    local n
    n=$(samtools view -f 1 -c --threads "$THREADS" "$1" 2>/dev/null || echo 0)
    [[ "$n" -gt 0 ]] && echo "PE" || echo "SE"
}

# Ensure BAM is indexed — index in-place if .bai is missing
ensure_index() {
    local bam="$1"
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
        log "  Indexing $(basename "$bam")..."
        samtools index --threads "$THREADS" "$bam"
    fi
}

# Extract CHROM reads from INPUT, subsample to TARGET, write to OUTPUT
subset_bam() {
    local input="$1"
    local output="$2"
    local target="$3"
    local label="$4"

    log "--- ${label} ---"
    log "  Input:  $input"

    ensure_index "$input"

    local total
    total=$(chr_mapped_count "$input")

    if [[ -z "$total" || "$total" -eq 0 ]]; then
        die "No reads found on '$CHROM' in $(basename "$input").
  Check the chromosome name — BAM may use '1' instead of 'chr1'.
  Run: samtools idxstats $input | head"
    fi

    log "  Reads on ${CHROM}: ${total}"

    if [[ "$total" -le "$target" ]]; then
        log "  Total (${total}) <= target (${target}) — keeping all ${CHROM} reads"
        samtools view -b --threads "$THREADS" "$input" "$CHROM" > "$output"
    else
        # samtools -s format: {seed}.{4-digit fraction}
        # e.g. seed=42, fraction=0.03  →  -s 42.0300
        local frac
        frac=$(awk -v t="$target" -v n="$total" 'BEGIN { printf "%.4f", t/n }')
        local frac_digits="${frac#0.}"   # strip leading "0." → "0300"
        local s_arg="${SEED}.${frac_digits}"
        log "  Subsampling ${frac} of reads (samtools -s ${s_arg})"

        samtools view -b --threads "$THREADS" "$input" "$CHROM" \
            | samtools view -b -s "$s_arg" --threads "$THREADS" \
            > "$output"
    fi

    samtools index "$output"

    local final
    final=$(samtools flagstat "$output" | awk '/mapped \(/{print $1; exit}')
    log "  Output: $output"
    log "  Final mapped reads: ${final}"
    echo ""
}

# ─────────────────────────────────────────────────────────────────────────────
# VALIDATION
# ─────────────────────────────────────────────────────────────────────────────

validate_config() {
    local missing=0

    check_var() {
        local val="$1" name="$2"
        if [[ -z "$val" ]]; then
            echo "  Not set: $name" >&2
            missing=1
        elif [[ ! -e "$val" ]]; then
            echo "  Not found: $name = $val" >&2
            missing=1
        fi
    }

    check_var "$BULK_CTRL_BAM"   "BULK_CTRL_BAM"
    check_var "$BULK_EXP_BAM"    "BULK_EXP_BAM"
    check_var "$SC_CTRL_BAM"     "SC_CTRL_BAM"
    check_var "$SC_EXP_BAM"      "SC_EXP_BAM"
    check_var "$HG38_FASTA"      "HG38_FASTA"
    check_var "$HG38_GTF"        "HG38_GTF"
    check_var "$HG38_GENE_BED"   "HG38_GENE_BED"
    check_var "$HG38_DBSNP_BED3" "HG38_DBSNP_BED3"
    check_var "$MM10_DBSNP_BED3" "MM10_DBSNP_BED3"
    check_var "$MM10_FASTA"      "MM10_FASTA"
    check_var "$MM10_GTF"        "MM10_GTF"
    check_var "$MM10_GENE_BED"   "MM10_GENE_BED"

    if [[ "$missing" -eq 1 ]]; then die "Fix the paths above and re-run."; fi
}

# Build filtered_feature_bc_matrix from a CellRanger-tagged BAM.
# Extracts CB (corrected barcode), UB (corrected UMI), GX (gene ID), GN (gene name).
# One UMI = one count; duplicate UBs per (barcode, gene) are deduplicated.
# Output: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz in 10x MTX format.
generate_sc_matrix() {
    local bam="$1"
    local out_dir="$2"
    local label="$3"
    log "--- ${label} ---"
    log "  BAM: $bam"
    log "  Output: $out_dir"

    python3 - "$bam" "$out_dir" <<'PYEOF'
import sys, gzip, subprocess, collections, os

bam_path, out_dir = sys.argv[1], sys.argv[2]

umi_sets  = collections.defaultdict(set)   # (cb, gx) -> set of UBs
gene_names = {}                            # gx -> gn

proc = subprocess.Popen(
    ['samtools', 'view', '-F', '2048', bam_path],
    stdout=subprocess.PIPE, text=True
)
for line in proc.stdout:
    f = line.rstrip('\n').split('\t')
    if len(f) < 12:
        continue
    cb = ub = gx = gn = None
    for tag in f[11:]:
        if tag[:3] == 'CB:':   cb = tag[5:]
        elif tag[:3] == 'UB:': ub = tag[5:]
        elif tag[:3] == 'GX:': gx = tag[5:]
        elif tag[:3] == 'GN:': gn = tag[5:]
    if cb and ub and gx:
        umi_sets[(cb, gx)].add(ub)
        if gx not in gene_names and gn:
            gene_names[gx] = gn
proc.wait()
if proc.returncode != 0:
    sys.exit(proc.returncode)

barcodes = sorted({cb for cb, _ in umi_sets})
genes    = sorted({gx for _, gx in umi_sets})
bc_idx   = {cb: i + 1 for i, cb in enumerate(barcodes)}
gx_idx   = {gx: i + 1 for i, gx in enumerate(genes)}

entries = sorted(
    (gx_idx[gx], bc_idx[cb], len(ubs))
    for (cb, gx), ubs in umi_sets.items()
)

with gzip.open(f'{out_dir}/barcodes.tsv.gz', 'wt') as fh:
    fh.write('\n'.join(barcodes) + '\n')

with gzip.open(f'{out_dir}/features.tsv.gz', 'wt') as fh:
    for gx in genes:
        fh.write(f'{gx}\t{gene_names.get(gx, gx)}\tGene Expression\n')

with gzip.open(f'{out_dir}/matrix.mtx.gz', 'wt') as fh:
    fh.write('%%MatrixMarket matrix coordinate integer general\n')
    fh.write('%metadata_json: {"format_version": 2}\n')
    fh.write(f'{len(genes)} {len(barcodes)} {len(entries)}\n')
    for gi, bi, count in entries:
        fh.write(f'{gi} {bi} {count}\n')

print(f'  {len(barcodes):,} barcodes, {len(genes):,} genes, {len(entries):,} nonzero entries',
      file=sys.stderr)
PYEOF

    log "  Done."
    echo ""
}

# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

log "ribostamp test dataset builder"
log "Chromosome: ${CHROM} | Bulk target: ${BULK_TARGET} reads | SC target: ${SC_TARGET} reads"
echo ""

validate_config

mkdir -p \
    "${OUTDIR}/bulk_bams" \
    "${OUTDIR}/bulk_fastqs" \
    "${OUTDIR}/sc_bams" \
    "${OUTDIR}/sc_fastqs/${SAMPLE_CTRL}/fastqs" \
    "${OUTDIR}/sc_fastqs/${SAMPLE_EXP}/fastqs" \
    "${OUTDIR}/sc_matrix/${SAMPLE_CTRL}/filtered_feature_bc_matrix" \
    "${OUTDIR}/sc_matrix/${SAMPLE_EXP}/filtered_feature_bc_matrix" \
    "${OUTDIR}/samplesheets" \
    "${OUTDIR}/refs"

# ── Bulk ──────────────────────────────────────────────────────────────────────
log "====== BULK ======"

subset_bam \
    "$BULK_CTRL_BAM" \
    "${OUTDIR}/bulk_bams/${SAMPLE_CTRL}_${CHROM}.bam" \
    "$BULK_TARGET" \
    "Bulk control"

subset_bam \
    "$BULK_EXP_BAM" \
    "${OUTDIR}/bulk_bams/${SAMPLE_EXP}_${CHROM}.bam" \
    "$BULK_TARGET" \
    "Bulk experimental"

BULK_LIB=$(detect_lib "${OUTDIR}/bulk_bams/${SAMPLE_CTRL}_${CHROM}.bam")
log "Detected bulk library type: ${BULK_LIB}"

BULK_SHEET="${OUTDIR}/samplesheets/samplesheet_bulk_bam.csv"
{
    printf "sample,bam,library_type\n"
    printf "%s,%s,%s\n" "$SAMPLE_CTRL" "../bulk_bams/${SAMPLE_CTRL}_${CHROM}.bam" "$BULK_LIB"
    printf "%s,%s,%s\n" "$SAMPLE_EXP"  "../bulk_bams/${SAMPLE_EXP}_${CHROM}.bam"  "$BULK_LIB"
} > "$BULK_SHEET"
log "Bulk samplesheet written: ${BULK_SHEET}"
echo ""

# ── Single-cell ───────────────────────────────────────────────────────────────
log "====== SINGLE-CELL ======"

subset_bam \
    "$SC_CTRL_BAM" \
    "${OUTDIR}/sc_bams/${SAMPLE_CTRL}_${CHROM}.bam" \
    "$SC_TARGET" \
    "SC control"

subset_bam \
    "$SC_EXP_BAM" \
    "${OUTDIR}/sc_bams/${SAMPLE_EXP}_${CHROM}.bam" \
    "$SC_TARGET" \
    "SC experimental"

SC_SHEET="${OUTDIR}/samplesheets/samplesheet_sc_bam.csv"
{
    printf "sample,bam,matrix_dir\n"
    printf "%s,%s,%s\n" \
        "$SAMPLE_CTRL" \
        "../sc_bams/${SAMPLE_CTRL}_${CHROM}.bam" \
        "../sc_matrix/${SAMPLE_CTRL}/filtered_feature_bc_matrix"
    printf "%s,%s,%s\n" \
        "$SAMPLE_EXP" \
        "../sc_bams/${SAMPLE_EXP}_${CHROM}.bam" \
        "../sc_matrix/${SAMPLE_EXP}/filtered_feature_bc_matrix"
} > "$SC_SHEET"
log "SC samplesheet written: ${SC_SHEET}"
echo ""

# ── Generate filtered_feature_bc_matrix from sc BAMs ─────────────────────────
log "====== SC FILTERED FEATURE BC MATRIX ======"

generate_sc_matrix \
    "${OUTDIR}/sc_bams/${SAMPLE_CTRL}_${CHROM}.bam" \
    "${OUTDIR}/sc_matrix/${SAMPLE_CTRL}/filtered_feature_bc_matrix" \
    "SC control → matrix"

generate_sc_matrix \
    "${OUTDIR}/sc_bams/${SAMPLE_EXP}_${CHROM}.bam" \
    "${OUTDIR}/sc_matrix/${SAMPLE_EXP}/filtered_feature_bc_matrix" \
    "SC experimental → matrix"

# ── hg38 chr1 STAR index (for FASTQ-start bulk test) ─────────────────────────
log "====== HG38 CHR1 STAR INDEX ======"

CHR1_FASTA="${OUTDIR}/refs/hg38_chr1.fa"
if [[ ! -f "$CHR1_FASTA" ]]; then
    log "Extracting hg38 chr1 FASTA..."
    samtools faidx "$HG38_FASTA" chr1 > "$CHR1_FASTA"
    log "  Written: $CHR1_FASTA"
else
    log "  Already exists: $CHR1_FASTA"
fi

if [[ ! -f "${STAR_INDEX_DIR}/SA" ]]; then
    log "Building chr1 STAR index (~5 min)..."
    mkdir -p "$STAR_INDEX_DIR"
    STAR \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX_DIR" \
        --genomeFastaFiles "$CHR1_FASTA" \
        --sjdbGTFfile "$HG38_GTF" \
        --genomeSAindexNbases 12 \
        --runThreadN "$THREADS"
    log "  STAR index built: $STAR_INDEX_DIR"
else
    log "  STAR index already exists: $STAR_INDEX_DIR"
fi
echo ""

# ── hg38 chr1 GTF and gene BED (bundled test references) ─────────────────────
log "====== HG38 CHR1 ANNOTATION FILES ======"

CHR1_GTF="${OUTDIR}/refs/chr1_hg38_genes.gtf.gz"
if [[ ! -f "$CHR1_GTF" ]]; then
    log "Extracting chr1-only GTF..."
    awk 'substr($0,1,1)=="#" || $1=="chr1"' "$HG38_GTF" | gzip > "$CHR1_GTF"
    log "  Written: $CHR1_GTF ($(zcat "$CHR1_GTF" | wc -l) lines)"
else
    log "  Already exists: $CHR1_GTF"
fi

CHR1_GENE_BED="${OUTDIR}/refs/chr1_hg38_genes.bed.gz"
if [[ ! -f "$CHR1_GENE_BED" ]]; then
    log "Extracting chr1-only gene BED..."
    awk '$1=="chr1"' "$HG38_GENE_BED" | gzip > "$CHR1_GENE_BED"
    log "  Written: $CHR1_GENE_BED ($(zcat "$CHR1_GENE_BED" | wc -l) lines)"
else
    log "  Already exists: $CHR1_GENE_BED"
fi

# ── hg38 chr1 dbSNP (for bulk BAM/FASTQ-start tests) ─────────────────────────
log "====== HG38 CHR1 DBSNP ======"

CHR1_HG38_DBSNP="${OUTDIR}/refs/chr1_hg38_dbsnp.bed"
if [[ ! -f "$CHR1_HG38_DBSNP" ]]; then
    log "Subsetting hg38 dbSNP to chr1 (~2 min)..."
    grep $'^chr1\t' "$HG38_DBSNP_BED3" > "$CHR1_HG38_DBSNP"
    log "  Written: $CHR1_HG38_DBSNP ($(wc -l < "$CHR1_HG38_DBSNP") lines)"
else
    log "  Already exists: $CHR1_HG38_DBSNP ($(wc -l < "$CHR1_HG38_DBSNP") lines)"
fi
echo ""
echo ""

# ── mm10 chr1 dbSNP (for SC BAM-start test) ───────────────────────────────────
log "====== MM10 CHR1 DBSNP ======"

CHR1_MM10_DBSNP="${OUTDIR}/refs/chr1_mm10_dbsnp.bed"
if [[ ! -f "$CHR1_MM10_DBSNP" ]]; then
    log "Subsetting mm10 dbSNP to chr1 (~2 min)..."
    grep $'^chr1\t' "$MM10_DBSNP_BED3" > "$CHR1_MM10_DBSNP"
    log "  Written: $CHR1_MM10_DBSNP ($(wc -l < "$CHR1_MM10_DBSNP") lines)"
else
    log "  Already exists: $CHR1_MM10_DBSNP ($(wc -l < "$CHR1_MM10_DBSNP") lines)"
fi
echo ""

# ── mm10 chr1 FASTA (for SC BAM-start test) ───────────────────────────────────
log "====== MM10 CHR1 FASTA ======"

MM10_CHR1_FASTA="${OUTDIR}/refs/mm10_chr1.fa"
if [[ ! -f "$MM10_CHR1_FASTA" ]]; then
    log "Extracting mm10 chr1 FASTA..."
    samtools faidx "$MM10_FASTA" chr1 > "$MM10_CHR1_FASTA"
    samtools faidx "$MM10_CHR1_FASTA"
    log "  Written: $MM10_CHR1_FASTA"
else
    log "  Already exists: $MM10_CHR1_FASTA"
fi
echo ""

# ── mm10 chr1 annotation files (GTF + gene BED) ──────────────────────────────
log "====== MM10 CHR1 ANNOTATION FILES ======"

CHR1_MM10_GTF="${OUTDIR}/refs/chr1_mm10_genes.gtf.gz"
if [[ ! -f "$CHR1_MM10_GTF" ]]; then
    log "Extracting chr1-only mm10 GTF..."
    zcat "$MM10_GTF" | awk 'substr($0,1,1)=="#" || $1=="chr1"' | gzip > "$CHR1_MM10_GTF"
    log "  Written: $CHR1_MM10_GTF ($(zcat "$CHR1_MM10_GTF" | wc -l) lines)"
else
    log "  Already exists: $CHR1_MM10_GTF"
fi

CHR1_MM10_GENE_BED="${OUTDIR}/refs/chr1_mm10_genes.bed.gz"
if [[ ! -f "$CHR1_MM10_GENE_BED" ]]; then
    log "Extracting chr1-only mm10 gene BED..."
    awk '$1=="chr1"' "$MM10_GENE_BED" | gzip > "$CHR1_MM10_GENE_BED"
    log "  Written: $CHR1_MM10_GENE_BED ($(zcat "$CHR1_MM10_GENE_BED" | wc -l) lines)"
else
    log "  Already exists: $CHR1_MM10_GENE_BED"
fi
echo ""

# ── chr1-only CellRanger reference (for SC FASTQ-start test) ─────────────────
log "====== MM10 CHR1 CELLRANGER REFERENCE ======"

CELLRANGER_REF_DIR="${OUTDIR}/refs/cellranger_ref_mm10"
if [[ ! -f "${CELLRANGER_REF_DIR}/reference.json" ]]; then
    log "Building chr1-only CellRanger reference (~15-20 min)..."
    # Resolve absolute paths before pushd — relative paths break after cd
    _fasta_abs=$(realpath "$MM10_CHR1_FASTA")
    _gtf_abs=$(realpath "$CHR1_MM10_GTF")
    # cellranger mkref writes output to ./<genome>/ in the cwd
    pushd "${OUTDIR}/refs" > /dev/null
    "${CELLRANGER_PATH}" mkref \
        --genome=cellranger_ref_mm10 \
        --fasta="$_fasta_abs" \
        --genes="$_gtf_abs" \
        --nthreads="$THREADS"
    popd > /dev/null
    log "  Written: $CELLRANGER_REF_DIR"
else
    log "  Already exists: $CELLRANGER_REF_DIR"
fi
echo ""

# ── FASTQ-start (convert hg38 bulk chr1 BAMs → FASTQ) ────────────────────────
log "====== BULK FASTQ-START ======"

bam_to_fastq() {
    local bam="$1" out_prefix="$2" lib="$3" label="$4"
    log "--- ${label} ---"
    if [[ "$lib" == "PE" ]]; then
        samtools fastq --threads "$THREADS" \
            -1 "${out_prefix}_R1.fastq.gz" \
            -2 "${out_prefix}_R2.fastq.gz" \
            -0 /dev/null -s /dev/null "$bam"
        log "  R1: ${out_prefix}_R1.fastq.gz"
        log "  R2: ${out_prefix}_R2.fastq.gz"
    else
        samtools fastq --threads "$THREADS" "$bam" | gzip > "${out_prefix}_R1.fastq.gz"
        log "  R1: ${out_prefix}_R1.fastq.gz"
    fi
    echo ""
}

bam_to_fastq \
    "${OUTDIR}/bulk_bams/${SAMPLE_CTRL}_${CHROM}.bam" \
    "${OUTDIR}/bulk_fastqs/${SAMPLE_CTRL}_${CHROM}" \
    "$BULK_LIB" \
    "Bulk control → FASTQ"

bam_to_fastq \
    "${OUTDIR}/bulk_bams/${SAMPLE_EXP}_${CHROM}.bam" \
    "${OUTDIR}/bulk_fastqs/${SAMPLE_EXP}_${CHROM}" \
    "$BULK_LIB" \
    "Bulk experimental → FASTQ"

FASTQ_SHEET="${OUTDIR}/samplesheets/samplesheet_bulk_fastq.csv"
{
    printf "sample,fastq_1,fastq_2,library_type\n"
    if [[ "$BULK_LIB" == "PE" ]]; then
        printf "%s,%s,%s,%s\n" \
            "$SAMPLE_CTRL" \
            "../bulk_fastqs/${SAMPLE_CTRL}_${CHROM}_R1.fastq.gz" \
            "../bulk_fastqs/${SAMPLE_CTRL}_${CHROM}_R2.fastq.gz" \
            "PE"
        printf "%s,%s,%s,%s\n" \
            "$SAMPLE_EXP" \
            "../bulk_fastqs/${SAMPLE_EXP}_${CHROM}_R1.fastq.gz" \
            "../bulk_fastqs/${SAMPLE_EXP}_${CHROM}_R2.fastq.gz" \
            "PE"
    else
        printf "%s,%s,,%s\n" "$SAMPLE_CTRL" "../bulk_fastqs/${SAMPLE_CTRL}_${CHROM}_R1.fastq.gz" "SE"
        printf "%s,%s,,%s\n" "$SAMPLE_EXP"  "../bulk_fastqs/${SAMPLE_EXP}_${CHROM}_R1.fastq.gz"  "SE"
    fi
} > "$FASTQ_SHEET"
log "FASTQ samplesheet written: ${FASTQ_SHEET}"
echo ""

# ── SC FASTQ-start (reconstruct 10x R1/R2 from chr1-subsampled CellRanger BAMs) ──
log "====== SC FASTQ-START ======"
#
# CellRanger BAMs encode the original 10x reads as SAM tags:
#   CR / CY  — raw cell barcode sequence / quality (R1 first 16bp)
#   UR / UY  — raw UMI sequence / quality         (R1 last 12bp)
#   SEQ/QUAL — cDNA read                          (R2)
# We reconstruct R1 and R2 from these tags so CellRanger count can re-process
# the chr1-subsampled data. Output follows CellRanger naming convention:
#   {sample}_S1_L001_R{1,2}_001.fastq.gz

sc_bam_to_10x_fastq() {
    local bam="$1" outdir="$2" sample="$3" label="$4"
    log "--- ${label} ---"
    mkdir -p "$outdir"

    local r1="${outdir}/${sample}_S1_L001_R1_001.fastq.gz"
    local r2="${outdir}/${sample}_S1_L001_R2_001.fastq.gz"

    # Python inline: reads BAM via samtools view, writes gzipped paired FASTQs.
    # -F 2048 excludes supplementary alignments (same QNAME, would duplicate entries).
    python3 - "$bam" "$r1" "$r2" <<'PYEOF'
import sys, gzip, subprocess

bam_path, r1_path, r2_path = sys.argv[1:]

def parse_tags(fields):
    tags = {}
    for f in fields:
        parts = f.split(':', 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags

with gzip.open(r1_path, 'wt', compresslevel=1) as r1_fh, \
     gzip.open(r2_path, 'wt', compresslevel=1) as r2_fh:
    proc = subprocess.Popen(
        ['samtools', 'view', '-F', '2048', bam_path],
        stdout=subprocess.PIPE, text=True
    )
    n = 0
    for line in proc.stdout:
        f = line.rstrip('\n').split('\t')
        if len(f) < 11 or f[9] == '*':
            continue
        qname = f[0]
        seq   = f[9]
        qual  = f[10] if f[10] != '*' else 'I' * len(f[9])
        tags  = parse_tags(f[11:])
        # Fallback to Ns if barcode/UMI tags are absent
        cr = tags.get('CR', 'N' * 16)
        ur = tags.get('UR', 'N' * 12)
        cy = tags.get('CY', 'I' * len(cr))
        uy = tags.get('UY', 'I' * len(ur))
        r1_fh.write(f'@{qname}\n{cr}{ur}\n+\n{cy}{uy}\n')
        r2_fh.write(f'@{qname}\n{seq}\n+\n{qual}\n')
        n += 1
    proc.wait()
    print(f'  Wrote {n:,} read pairs', file=sys.stderr)
    if proc.returncode != 0:
        sys.exit(proc.returncode)
PYEOF

    log "  R1: ${r1}"
    log "  R2: ${r2}"
    echo ""
}

sc_bam_to_10x_fastq \
    "${OUTDIR}/sc_bams/${SAMPLE_CTRL}_${CHROM}.bam" \
    "${OUTDIR}/sc_fastqs/${SAMPLE_CTRL}" \
    "$SAMPLE_CTRL" \
    "SC control → 10x FASTQ"

sc_bam_to_10x_fastq \
    "${OUTDIR}/sc_bams/${SAMPLE_EXP}_${CHROM}.bam" \
    "${OUTDIR}/sc_fastqs/${SAMPLE_EXP}" \
    "$SAMPLE_EXP" \
    "SC experimental → 10x FASTQ"

SC_FASTQ_SHEET="${OUTDIR}/samplesheets/samplesheet_sc_fastq.csv"
{
    printf "sample,fastq_dir\n"
    printf "%s,%s\n" "$SAMPLE_CTRL" "../sc_fastqs/${SAMPLE_CTRL}"
    printf "%s,%s\n" "$SAMPLE_EXP"  "../sc_fastqs/${SAMPLE_EXP}"
} > "$SC_FASTQ_SHEET"
log "SC FASTQ samplesheet written: ${SC_FASTQ_SHEET}"
echo ""

# ── Summary ───────────────────────────────────────────────────────────────────
log "====== SUMMARY ======"
echo ""
log "All output files:"
find "${OUTDIR}" -type f | sort | while read -r f; do
    printf "  %-60s %s\n" "$f" "$(du -sh "$f" 2>/dev/null | cut -f1)"
done
echo ""
log "Config hints:"
echo ""
echo "  conf/test.config (bulk BAM-start, hg38) — already configured"
echo "    input = '$(realpath "$BULK_SHEET")'"
echo ""
echo "  conf/test_fastq.config (bulk FASTQ-start, hg38)"
echo "    input      = '$(realpath "$FASTQ_SHEET")'"
echo "    star_index = '$(realpath "$STAR_INDEX_DIR")'"
echo ""
echo "  conf/test_sc_bam.config (SC BAM-start, mm10)"
echo "    input     = '$(realpath "$SC_SHEET")'"
echo "    dbsnp_bed = '$(realpath "$CHR1_MM10_DBSNP")'"
echo ""
echo "  conf/test_sc_fastq.config (SC FASTQ-start, mm10) — update input path:"
echo "    input = '$(realpath "$SC_FASTQ_SHEET")'"
echo ""
log "Done."
