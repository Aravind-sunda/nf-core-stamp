# Local modifications to vendored FLARE

This is a patched copy of upstream [YeoLab/FLARE](https://github.com/YeoLab/FLARE).
Re-apply these when re-vendoring a newer version. Find edits in code with:
`grep -rn "LOCAL PATCH (ribostamp)" assets/workflow_FLARE/`

| Date | File:line | Change | Why |
|------|-----------|--------|-----|
| 2026-07-07 | `Snakefile:75` | `regions_file = ''` → `regions_file = []` | Snakemake ≥9 rejects empty-string rule inputs (`rule split_regions`). No behavioural change. |
