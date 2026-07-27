# Local modifications to vendored SAILOR

This is a patched copy of upstream [YeoLab/FLARE](https://github.com/YeoLab/FLARE) (`workflow_sailor`).
Re-apply these when re-vendoring a newer version. Find edits in code with:
`grep -rn "LOCAL PATCH (stamp)" assets/workflow_sailor/`

| Date | File:line | Change | Why |
|------|-----------|--------|-----|
| 2026-07-26 | `Snakefile:284` | Added `threads: max(1, workflow.cores // 2)` to `rule filter_known_snp` | `filter_known_snp.py` holds the entire dbSNP BED in memory (50-100 GB for hg38). No rule in this Snakefile declares `threads:`, so Snakemake scheduled every job as 1 thread and ran `--cores` of them at once — 17 OOM kills on a 16-sample hg38 run at `--cores 36`. Capping this rule at 2 concurrent jobs leaves all other rules fully parallel. Behaviour unchanged. |
| 2026-07-26 | `scripts/filter_known_snp.py:88` | Comment only — records the inverted-join rewrite, not enabled | Documents the real fix for the memory blow-up above so it isn't lost. No behavioural change. |
| 2026-07-27 | `Snakefile:158` | Added `threads: workflow.cores` to `rule make_bigwigs` | `scripts/bam_to_bw.sh` writes a shared `8_bw_and_bam/chrom.sizes` and reads it back immediately, so concurrent jobs let one truncate the file mid-read. `bedGraphToBigWig` then fails with "chr1 is not found in chromosome sizes file" but, with no `set -e` in the script, still exits 0 — Snakemake reports success with missing output. Serialising the rule removes the overlap. Behaviour unchanged; see the note below on the cost. |

## Deferred: inverting the join in `filter_known_snp.py`

The `threads:` cap above is a mitigation, not a fix — it bounds concurrency but each
job still materialises ~600M dbSNP rows (~15 GB on disk, 50-100 GB resident) purely to
left-join them against a ~150k-row VCF.

The fix is to build the lookup from the small side and stream the large one, which
bounds peak memory to one chunk plus the VCF. The proposed implementation is kept as a
comment block at `scripts/filter_known_snp.py:88`. It is semantically equivalent to the
current merge: a VCF row matching *k* duplicate dbSNP entries expands to *k* rows that
all carry `KNOWN=1` and are dropped regardless, so "drop if present in dbSNP" is exactly
set membership.

Deliberately not enabled yet — deferred so the vendored script stays byte-compatible
with upstream behaviour until there is time to validate it against a full run.

Note it would still stream the full dbSNP once per job. Removing that repetition across
samples means restructuring the rule itself (filter once for all samples rather than per
sample-strand), which is a larger change to vendored code.

## Deferred: per-sample `chrom.sizes` in `bam_to_bw.sh`

Serialising `make_bigwigs` above removes the race but not its cause: every job still
writes and reads the same `8_bw_and_bam/chrom.sizes`. The root fix is to make the file
per-sample — three references in `scripts/bam_to_bw.sh` (one write at line 60, two reads
at lines 63 and 65) change from `$output_dir/chrom.sizes` to
`$output_dir/$samplename.chrom.sizes`.

Deliberately not done yet: serialisation was preferred as the smaller change to vendored
code. Worth revisiting if bigwig generation becomes a bottleneck, because the current
mitigation costs wall-clock time proportional to sample count — the jobs took roughly
2-3 minutes each in testing, so a 16-sample run serialises to ~40 minutes of what was
previously parallel work. It would also make the script correct for anyone running
SAILOR outside Nextflow, where the `threads` directive does not apply.

## Related, outside this directory

`filter_known_snp.py` matches dbSNP by **exact position equality**, whereas the MARINE-side
filter (`bin/helper_filter_edits_bulk.py`) uses **interval overlap** via bedtools. For
multi-base dbSNP entries — hg38 contains indels up to ~900 bp — a site inside such an entry
is dropped by the MARINE filter but kept by SAILOR's. This is upstream SAILOR behaviour and
is preserved deliberately; changing it would silently alter SAILOR results and break
comparability with previously published SAILOR runs.
