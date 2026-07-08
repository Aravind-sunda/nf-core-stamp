#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
#  Render the nf-core/stamp metro map with nf-metro 1.1.0
#  Docs / directive reference: https://github.com/seqeralabs/nf-metro
#
#  nf-metro 1.1.0 CLI changes vs 0.7.2 used previously:
#    • --no-straight-diamonds        → replaced by  --diamond-style [straight|symmetric]
#      (straight = keep top branch on the main track; symmetric = fan branches evenly)
#    • --center-ports                → unchanged (still valid)
#    • --theme / --logo / --format   → unchanged
#    • --validate (on render)        → new: runs the geometry oracle on the drawn SVG
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

module load mamba
mamba activate
mamba activate nextflow

cd "$(dirname "${BASH_SOURCE[0]}")"

nf-metro validate metro_map.mmd

nf-metro render metro_map.mmd \
    --theme nfcore \
    --diamond-style symmetric \
    --center-ports \
    --validate \
    -o images/stamp_metro.svg \
    --logo images/nf-core-stamp_logo_dark.png

# Render as interactive HTML (pan/zoom, hover stations, click legend to isolate a line).
nf-metro render metro_map.mmd \
    --theme nfcore \
    --diamond-style straight \
    --center-ports \
    --format html \
    -o images/stamp_metro.html \
    --logo images/nf-core-stamp_logo_dark.png
    