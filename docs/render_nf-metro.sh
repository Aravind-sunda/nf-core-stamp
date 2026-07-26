#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
#  Render the nf-core/stamp metro map with nf-metro 1.1.0
#  Docs / directive reference: https://github.com/seqeralabs/nf-metro
#
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
#    • --mode                  bakes that mode's palette in rather than leaving it
#      to the viewer's OS preference. Passed explicitly on every render so neither
#      variant depends on the '%%metro style: dark' default in the .mmd.
#    • --logo                  must match the mode: the _dark logo is the
#      white-text lockup for dark backgrounds, _light the black-text one for
#      light backgrounds.
#
#  Deliberately NOT used:
#    • --compact-offsets  saves no canvas here (measured: identical 1568x689 with
#      and without) and opens a one-track gap between the SAILOR and FLARE lines
#      in the pre-processing bundle, because each station is sized only for the
#      lines crossing it.
#    • --line-spread rails  rejected by nf-metro's routing invariants on this map
#      (bundle-order flip at bam_bulk_in -> infer_strand).
# ─────────────────────────────────────────────────────────────────────────────
# Activate before 'set -u': the cluster's MKL activation hook reads unset
# variables, which would abort the script.
module load mamba
mamba activate
mamba activate nextflow

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

nf-metro validate metro_map.mmd

# Static SVG, dark — non-animated fallback for dark-theme readers.
nf-metro render metro_map.mmd \
    --theme nfcore \
    --diamond-style straight \
    --center-ports \
    --line-spread centered \
    --directional \
    --mode dark \
    --validate \
    -o images/stamp_metro_dark.svg \
    --logo images/nf-core-stamp_logo_dark.png

# Static SVG, light — same map and layout, light palette.
nf-metro render metro_map.mmd \
    --theme nfcore \
    --diamond-style straight \
    --center-ports \
    --line-spread centered \
    --directional \
    --mode light \
    --validate \
    -o images/stamp_metro_light.svg \
    --logo images/nf-core-stamp_logo_light.png

# Animated SVG, dark — embedded in the README for dark-theme readers.
# --embed-font inlines Inter so the SVG renders identically on GitHub, which
# serves it without access to local fonts.
nf-metro render metro_map.mmd \
    --theme nfcore \
    --diamond-style straight \
    --center-ports \
    --line-spread centered \
    --directional \
    --mode dark \
    --animate \
    --embed-font \
    -o images/stamp_metro_animated_dark.svg \
    --logo images/nf-core-stamp_logo_dark.png

# Animated SVG, light — embedded in the README for light-theme readers.
nf-metro render metro_map.mmd \
    --theme nfcore \
    --diamond-style straight \
    --center-ports \
    --line-spread centered \
    --directional \
    --mode light \
    --animate \
    --embed-font \
    -o images/stamp_metro_animated_light.svg \
    --logo images/nf-core-stamp_logo_light.png

# Interactive HTML, dark — pan/zoom, hover stations, click legend to isolate a
# line. --validate is SVG-only, so the HTML renders cannot carry it.
nf-metro render metro_map.mmd \
    --theme nfcore \
    --diamond-style straight \
    --center-ports \
    --line-spread centered \
    --directional \
    --mode dark \
    --format html \
    -o images/stamp_metro_dark.html \
    --logo images/nf-core-stamp_logo_dark.png

# Interactive HTML, light — same interactions, light palette.
nf-metro render metro_map.mmd \
    --theme nfcore \
    --diamond-style straight \
    --center-ports \
    --line-spread centered \
    --directional \
    --mode light \
    --format html \
    -o images/stamp_metro_light.html \
    --logo images/nf-core-stamp_logo_light.png
