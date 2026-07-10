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
#
#  Layout flags shared by every render below:
#    • --line-spread centered  balance bundles about the midline; cuts the canvas
#      ~20% versus the default 'bundle' with no loss of clarity.
#    • --directional           chevrons along every route, so the FASTQ/BAM inputs
#      and every downstream hop show flow direction.
#
#  Deliberately NOT used:
#    • --compact-offsets  saves no canvas here (measured: identical 1568x689 with
#      and without) and opens a one-track gap between the SAILOR and FLARE lines
#      in the pre-processing bundle, because each station is sized only for the
#      lines crossing it.
#    • --line-spread rails  rejected by nf-metro's routing invariants on this map
#      (bundle-order flip at bam_bulk_in -> infer_strand).
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

module load mamba
mamba activate
mamba activate nextflow

cd "$(dirname "${BASH_SOURCE[0]}")"

nf-metro validate metro_map.mmd

nf-metro render metro_map.mmd \
    --theme nfcore \
    --diamond-style straight \
    --center-ports \
    --line-spread centered \
    --directional \
    --validate \
    -o images/stamp_metro.svg \
    --logo images/nf-core-stamp_logo_dark.png

# Animated variant embedded in the README. --embed-font inlines Inter so the SVG
# renders identically on GitHub, which serves it without access to local fonts.
nf-metro render metro_map.mmd \
    --theme nfcore \
    --diamond-style straight \
    --center-ports \
    --line-spread centered \
    --directional \
    --animate \
    --embed-font \
    -o images/stamp_metro_animated.svg \
    --logo images/nf-core-stamp_logo_dark.png

# Render as interactive HTML (pan/zoom, hover stations, click legend to isolate a line).
nf-metro render metro_map.mmd \
    --theme nfcore \
    --diamond-style straight \
    --center-ports \
    --line-spread centered \
    --directional \
    --format html \
    -o images/stamp_metro.html \
    --logo images/nf-core-stamp_logo_dark.png
    