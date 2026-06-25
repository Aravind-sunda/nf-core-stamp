#!/usr/bin/env python3
"""
metaplot_dist.py — a single-file Python reimplementation of the metaPlotR
(Olarerin-George & Jaffrey, Bioinformatics 2017) distance-measures pipeline.

It takes a gene-prediction (genePred) annotation and a BED6 file of sites
(e.g. RNA-editing positions, m6A sites, RBP crosslink sites) and produces the
"distance measures" table that drives metagene plots:

    chr  coord  gene_name  refseqID  rel_location
    utr5_st  utr5_end  cds_st  cds_end  utr3_st  utr3_end
    utr5_size  cds_size  utr3_size

It reproduces the exact arithmetic of the original Perl scripts
(make_annot_bed.pl + size_of_cds_utrs.pl + annotate_bed_file.pl +
rel_and_abs_dist_calc.pl) but skips the genome-sized intermediate "master
annotation BED": all mRNA-coordinate math is done analytically per transcript,
so no genome FASTA and no bedtools are required.

Faithfulness notes (matching the original behaviour by default):
  * Sites are assigned to a transcript only if strands match (like
    `intersectBed -s`). Use --ignore-strand to disable.
  * By default only the longest transcript per gene is used (one isoform per
    gene). Use --all-isoforms to instead emit one row per (site, overlapping
    isoform) — the redundant table produced by the original metaPlotR.
  * Non-coding transcripts (cdsStart == cdsEnd) are excluded.
  * Transcripts lacking ANY of the three regions (no 5'UTR, no CDS, or no
    3'UTR in mRNA space) are excluded — this mirrors the original, where a
    missing region writes "NA" and the site is dropped. Use
    --keep-incomplete to relax this (regions present are still reported;
    rel_location is only defined for sites in an annotated region).
  * `*_size` columns are max-min of mRNA endpoints (NOT +1), matching the
    original size_features().

Author: generated for Wall-e; algorithm faithful to olarerin/metaPlotR.
"""

import argparse
import gzip
import sys
from collections import defaultdict

BIN = 16384  # genomic binning window for fast single-nucleotide overlap lookup


def open_file(path):
    """Open a plain-text or gzip-compressed file transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def log(msg):
    print(msg, file=sys.stderr, flush=True)


# --------------------------------------------------------------------------- #
# genePred parsing
# --------------------------------------------------------------------------- #
def detect_layout(fields):
    """Return (has_bin, name2_idx_or_None) by inspecting one data row.

    Standard genePred:      name chrom strand txStart ... (strand at idx 2)
    genePred with UCSC bin: bin name chrom strand ...     (strand at idx 3)
    name2 (gene symbol) lives at idx 11 (no bin) or 12 (bin) in genePredExt.
    """
    if len(fields) > 3 and fields[3] in ("+", "-"):
        has_bin = True
    elif len(fields) > 2 and fields[2] in ("+", "-"):
        has_bin = False
    else:
        raise ValueError(
            "Could not locate the strand column; is this a genePred file? "
            "First row: " + "\t".join(fields[:6])
        )
    # name2 index for the gene symbol, if the extended columns are present
    name2_idx = None
    base = 1 if has_bin else 0
    # genePredExt order after strand: txStart txEnd cdsStart cdsEnd exonCount
    # exonStarts exonEnds score name2 ...
    candidate = base + 11  # name + chrom + strand + 8 fields + score = name2
    if len(fields) > candidate:
        name2_idx = candidate
    return has_bin, name2_idx


def parse_genepred(path, force_bin=None):
    """Yield dicts: name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd,
    exonStarts (list, 0-based), exonEnds (list, 0-based half-open), gene_name."""
    with open_file(path) as fh:
        has_bin = force_bin
        name2_idx = None
        layout_done = has_bin is not None
        for raw in fh:
            if not raw.strip():
                continue
            f = raw.rstrip("\n").split("\t")
            # Skip an obvious header line
            if f[0].startswith("#") or f[0].lower() in ("bin", "name"):
                continue
            if not layout_done:
                has_bin, name2_idx = detect_layout(f)
                layout_done = True
            elif name2_idx is None and force_bin is not None:
                # layout forced; still try to find name2
                _, name2_idx = detect_layout(f)
            o = 1 if has_bin else 0
            try:
                name = f[o + 0]
                chrom = f[o + 1]
                strand = f[o + 2]
                txStart = int(f[o + 3])
                txEnd = int(f[o + 4])
                cdsStart = int(f[o + 5])
                cdsEnd = int(f[o + 6])
                exonStarts = [int(x) for x in f[o + 8].rstrip(",").split(",") if x != ""]
                exonEnds = [int(x) for x in f[o + 9].rstrip(",").split(",") if x != ""]
            except (IndexError, ValueError) as e:
                log(f"WARN: skipping malformed genePred row: {e}")
                continue
            gene = name
            if name2_idx is not None and len(f) > name2_idx and f[name2_idx]:
                gene = f[name2_idx]
            yield {
                "name": name, "chrom": chrom, "strand": strand,
                "txStart": txStart, "txEnd": txEnd,
                "cdsStart": cdsStart, "cdsEnd": cdsEnd,
                "exonStarts": exonStarts, "exonEnds": exonEnds,
                "gene": gene,
            }


# --------------------------------------------------------------------------- #
# Per-transcript geometry
# --------------------------------------------------------------------------- #
def exonic_len_in(exons, lo, hi):
    """Total exonic length intersecting genomic half-open interval [lo, hi)."""
    total = 0
    for s, e in exons:
        a = s if s > lo else lo
        b = e if e < hi else hi
        if b > a:
            total += (b - a)
    return total


def build_transcript(tx, keep_incomplete=False):
    """Compute mRNA geometry for one transcript.

    Returns a dict with sorted exons, cumulative offsets, total length L,
    strand, gene, chrom, and the six mRNA endpoints
    [u5_st,u5_end,cds_st,cds_end,u3_st,u3_end] (1-based mRNA positions),
    or None if the transcript should be excluded.
    """
    if tx["cdsStart"] >= tx["cdsEnd"]:
        return None  # non-coding (matches cdsStart==cdsEnd ncRNA exclusion)

    exons = sorted(zip(tx["exonStarts"], tx["exonEnds"]))
    L = sum(e - s for s, e in exons)
    if L == 0:
        return None

    # Exonic lengths in the three genomic partitions (standard genePred split)
    L_left = exonic_len_in(exons, tx["txStart"], tx["cdsStart"])   # genomic 5' side
    L_cds = exonic_len_in(exons, tx["cdsStart"], tx["cdsEnd"])
    L_right = exonic_len_in(exons, tx["cdsEnd"], tx["txEnd"])      # genomic 3' side

    if tx["strand"] == "+":
        a, b, c = L_left, L_cds, L_right       # 5'UTR, CDS, 3'UTR exonic lengths
    elif tx["strand"] == "-":
        a, b, c = L_right, L_cds, L_left
    else:
        return None  # require explicit strand

    if not keep_incomplete and (a == 0 or b == 0 or c == 0):
        return None  # mirrors original: any missing region -> excluded

    # mRNA endpoints (1-based). Regions are contiguous in mRNA space.
    endpts = [1, a, a + 1, a + b, a + b + 1, a + b + c]

    # Cumulative exonic length preceding each exon (ascending genomic order)
    cum = []
    acc = 0
    for s, e in exons:
        cum.append(acc)
        acc += (e - s)

    return {
        "name": tx["name"], "gene": tx["gene"], "chrom": tx["chrom"],
        "strand": tx["strand"], "exons": exons, "cum": cum, "L": L,
        "endpts": endpts,
        "trx_len": (a + b + c),  # for longest-isoform selection
    }


def mrna_pos_of(t, g):
    """1-based mRNA (5'->3') position of genomic 0-based coordinate g in
    transcript t, or None if g is not exonic in t."""
    exons = t["exons"]
    cum = t["cum"]
    for i, (s, e) in enumerate(exons):
        if s <= g < e:
            genomic_rank = cum[i] + (g - s) + 1   # rank from low genomic coord
            if t["strand"] == "+":
                return genomic_rank
            else:
                return t["L"] - genomic_rank + 1
    return None


def rel_distance(mrna_pos, endpts):
    u5_st, u5_end, cds_st, cds_end, u3_st, u3_end = endpts
    if u5_st <= mrna_pos <= u5_end:
        return (mrna_pos - u5_st + 1) / (u5_end - u5_st + 1)
    if cds_st <= mrna_pos <= cds_end:
        return (mrna_pos - cds_st + 1) / (cds_end - cds_st + 1) + 1
    if u3_st <= mrna_pos <= u3_end:
        return (mrna_pos - u3_st + 1) / (u3_end - u3_st + 1) + 2
    return None


def abs_distances(mrna_pos, endpts):
    return [mrna_pos - p for p in endpts]


def feature_sizes(endpts):
    return (endpts[1] - endpts[0], endpts[3] - endpts[2], endpts[5] - endpts[4])


# --------------------------------------------------------------------------- #
# Query / annotation consistency
# --------------------------------------------------------------------------- #
def check_chrom_consistency(bed_path, genepred_chroms):
    """Scan query-BED chromosome names and compare to the genePred's.

    Warns about BED chromosomes absent from the genePred (those sites are
    silently skipped during overlap), and aborts if there is NO overlap at
    all -- almost always a 'chr1' vs '1' naming mismatch. Returns the set of
    BED chromosomes seen."""
    bed_chroms = set()
    with open(bed_path) as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            bed_chroms.add(raw.split("\t", 1)[0])

    overlap = bed_chroms & genepred_chroms
    missing = bed_chroms - genepred_chroms

    if missing:
        log(f"WARN: {len(missing)} BED chromosome name(s) absent from the "
            f"genePred (their sites will be skipped): {sorted(missing)[:5]}"
            + (" ..." if len(missing) > 5 else ""))

    if not overlap:
        bed_has_chr = any(c.startswith("chr") for c in bed_chroms)
        gp_has_chr = any(c.startswith("chr") for c in genepred_chroms)
        hint = ""
        if gp_has_chr and not bed_has_chr:
            hint = "  (genePred uses a 'chr' prefix; the BED does not)"
        elif bed_has_chr and not gp_has_chr:
            hint = "  (the BED uses a 'chr' prefix; the genePred does not)"
        log("ERROR: no chromosome names are shared between the BED and the "
            "genePred." + hint)
        log(f"       BED examples:      {sorted(bed_chroms)[:4]}")
        log(f"       genePred examples: {sorted(genepred_chroms)[:4]}")
        sys.exit(1)

    log(f"  chromosome check: {len(overlap)} shared, "
        f"{len(missing)} BED-only (skipped)")
    return bed_chroms


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description="metaPlotR distance-measures, single-script Python port.")
    ap.add_argument("--genePred", required=True,
                    help="UCSC genePred / genePredExt annotation (optionally "
                         "with leading bin column; auto-detected). "
                         "Accepts plain text or gzip-compressed (.gz) files.")
    ap.add_argument("--bed", required=True,
                    help="BED6 of sites: 0-based, single-nucleotide, stranded "
                         "(column 6). Columns 4/5 are ignored.")
    ap.add_argument("-o", "--out", default="-",
                    help="Output distance-measures file (default: stdout).")
    ap.add_argument("--ignore-strand", action="store_true",
                    help="Assign sites to transcripts regardless of strand "
                         "(default: require matching strand, like bedtools -s).")
    ap.add_argument("--keep-incomplete", action="store_true",
                    help="Keep transcripts missing a 5'UTR/CDS/3'UTR region "
                         "(default: exclude, matching the original).")
    ap.add_argument("--all-isoforms", action="store_true",
                    help="Keep every overlapping isoform (one row per site x "
                         "transcript), like the original metaPlotR. Default: "
                         "use only the longest transcript per gene.")
    ap.add_argument("--force-bin", choices=["0", "1"], default=None,
                    help="Override genePred layout detection: 1 = has leading "
                         "bin column, 0 = no bin column.")
    args = ap.parse_args()

    force_bin = None if args.force_bin is None else (args.force_bin == "1")

    # --- Load + build transcripts -----------------------------------------
    log("Parsing genePred and building transcript geometry...")
    transcripts = []
    n_raw = n_kept = 0
    for tx in parse_genepred(args.genePred, force_bin=force_bin):
        n_raw += 1
        built = build_transcript(tx, keep_incomplete=args.keep_incomplete)
        if built is not None:
            transcripts.append(built)
            n_kept += 1
    log(f"  transcripts read: {n_raw}; usable: {n_kept}")
    genepred_chroms = {t["chrom"] for t in transcripts}

    # --- Check BED vs genePred chromosome naming up front -----------------
    check_chrom_consistency(args.bed, genepred_chroms)

    # --- Isoform selection: longest per gene by default -------------------
    if not args.all_isoforms:
        best = {}
        for t in transcripts:
            cur = best.get(t["gene"])
            if cur is None or t["trx_len"] > cur["trx_len"]:
                best[t["gene"]] = t
        transcripts = list(best.values())
        log(f"  longest-isoform selection: {len(transcripts)} transcripts "
            f"(one per gene); use --all-isoforms to keep all")
    else:
        log(f"  keeping all isoforms: {len(transcripts)} transcripts")

    # --- Index exons by genomic bin for fast single-nt overlap ------------
    index = defaultdict(list)  # (chrom, bin) -> list of transcript objects
    for t in transcripts:
        seen_bins = set()
        for s, e in t["exons"]:
            for b in range(s // BIN, (e - 1) // BIN + 1):
                key = (t["chrom"], b)
                if key not in seen_bins:
                    index[key].append(t)
                    seen_bins.add(key)

    # --- Stream the query BED ---------------------------------------------
    out = sys.stdout if args.out == "-" else open(args.out, "w")
    out.write("chr\tcoord\tgene_name\trefseqID\trel_location\t"
              "utr5_st\tutr5_end\tcds_st\tcds_end\tutr3_st\tutr3_end\t"
              "utr5_size\tcds_size\tutr3_size\n")

    total = matched = emitted = 0
    with open(args.bed) as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            f = raw.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            chrom = f[0]
            start = int(f[1])
            end = int(f[2])
            strand = f[5] if len(f) >= 6 else "."
            total += 1
            hit_any = False
            # iterate each genomic position covered by the site (1 for editing)
            for g in range(start, end):
                key = (chrom, g // BIN)
                cands = index.get(key)
                if not cands:
                    continue
                for t in cands:
                    if not args.ignore_strand and strand in ("+", "-") \
                            and t["strand"] != strand:
                        continue
                    mp = mrna_pos_of(t, g)
                    if mp is None:
                        continue
                    rel = rel_distance(mp, t["endpts"])
                    if rel is None:
                        continue  # site exonic but outside annotated regions
                    ad = abs_distances(mp, t["endpts"])
                    s5, sc, s3 = feature_sizes(t["endpts"])
                    out.write(
                        f"{chrom}\t{g + 1}\t{t['gene']}\t{t['name']}\t"
                        f"{rel:.6f}\t"
                        f"{ad[0]}\t{ad[1]}\t{ad[2]}\t{ad[3]}\t{ad[4]}\t{ad[5]}\t"
                        f"{s5}\t{sc}\t{s3}\n")
                    emitted += 1
                    hit_any = True
            if hit_any:
                matched += 1

    if out is not sys.stdout:
        out.close()
    log(f"** query sites: {total}; sites with >=1 transcript hit: {matched}; "
        f"rows emitted: {emitted}")


if __name__ == "__main__":
    main()
