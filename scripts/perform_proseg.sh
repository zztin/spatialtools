#!/bin/bash

# Check if output directory is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <output_directory>"
    exit 1
fi

# Define output directory and create it if it doesn't exist
OUTPUT_DIR="$1"
mkdir -p "$OUTPUT_DIR"

# Set input directory (modify if needed)
INPUT_DIR="/hpc/pmc_holstege/lchen/projects/spatial/data/02_proseg/split_cores/transcripts"
cd "$INPUT_DIR" || exit 1


# Perform on only a subset of control donors
#for file in transcripts_D3*.csv.gz transcripts_D5*.csv.gz transcripts_D6*.csv.gz transcripts_D11*.csv.gz; do
    # Skip if no match (e.g., no D11 files)
#    [ -e "$file" ] || continue
#    SAMPLE=$(basename "$file" .csv.gz | sed 's/^transcripts_//')


# Loop over each transcript file
for file in transcripts_*.csv.gz; do
    # Extract SAMPLE name by removing prefix and suffix
    SAMPLE=$(basename "$file" .csv.gz | sed 's/^transcripts_//')

    LOGFILE="${OUTPUT_DIR}/${SAMPLE}.log"

    echo "[$(date)] Processing $SAMPLE..." | tee "$LOGFILE"

    proseg --xenium "$file" \
        --output-expected-counts "${OUTPUT_DIR}/${SAMPLE}_expected-counts.csv.gz" \
        --output-cell-metadata "${OUTPUT_DIR}/${SAMPLE}_cell-metadata.csv.gz" \
        --output-transcript-metadata "${OUTPUT_DIR}/${SAMPLE}_transcript-metadata.csv.gz" \
        --output-gene-metadata "${OUTPUT_DIR}/${SAMPLE}_gene-metadata.csv.gz" \
        --output-rates "${OUTPUT_DIR}/${SAMPLE}_rates.csv.gz" \
        --output-cell-polygons "${OUTPUT_DIR}/${SAMPLE}_cell-polygons.geojson.gz" \
        --output-cell-polygon-layers "${OUTPUT_DIR}/${SAMPLE}_cell-polygon-layers.geojson.gz" \
        --output-union-cell-polygons "${OUTPUT_DIR}/${SAMPLE}_union-cell-polygons.geojson.gz"
        >> "$LOGFILE" 2>&1

    echo "[$(date)] Done with $SAMPLE" | tee -a "$LOGFILE"


done

