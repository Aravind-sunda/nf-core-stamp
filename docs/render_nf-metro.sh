module load mamba 
mamba activate
mamba activate nextflow

cd /home/tmhaxs421/brannanlab/tmhaxs421/MARINE_NextFlow/nf-core-ribostamp/docs

# Validate first (catches syntax errors without rendering)
nf-metro validate metro_map.mmd

# Render as SVG (static, for README/docs)
nf-metro render metro_map.mmd \
    --theme nfcore \
    -o images/ribostamp_metro.svg

# Render as interactive HTML (pan/zoom, click legend to isolate lines)
nf-metro render metro_map.mmd \
    --theme nfcore \
    --format html \
    -o images/ribostamp_metro.html


nf-metro validate metro_map_v2.mmd

nf-metro render metro_map_v2.mmd --theme nfcore -o images/ribostamp_metro_v2.svg
