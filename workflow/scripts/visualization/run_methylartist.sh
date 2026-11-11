#!/bin/bash
set -euo pipefail

BAM=""
BED_FILE=""
REFERENCE=""
OUTPUT_DIR=""
THREADS="1"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam)
            BAM="$2"
            shift 2
            ;;
        --bed-file)
            BED_FILE="$2"
            shift 2
            ;;
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$BAM" || -z "$BED_FILE" || -z "$REFERENCE" || -z "$OUTPUT_DIR" ]]; then
    echo "Usage: $0 --bam BAM1,BAM2,... --bed-file BED_FILE --reference REFERENCE --output-dir OUTPUT_DIR [--threads N]" >&2
    exit 1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
    echo "--threads must be a positive integer" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Log the BAMs being used
echo "Using BAMs:"
IFS=',' read -r -a BAMS <<< "$BAM"
printf '  %s\n' "${BAMS[@]}"
echo "Total: ${#BAMS[@]} sample(s)"

while IFS= read -r LINE; do
    [[ -z "${LINE:-}" ]] && continue
    [[ "${LINE:0:1}" == "#" ]] && continue

    read -r CHROM START END LABEL _ <<< "$LINE"

    if [[ -z "${CHROM:-}" || -z "${START:-}" || -z "${END:-}" ]]; then
        echo "Skipping malformed BED row (expected >=3 columns): $LINE" >&2
        continue
    fi
    if ! [[ "$START" =~ ^[0-9]+$ && "$END" =~ ^[0-9]+$ ]]; then
        echo "Skipping BED row with non-numeric coordinates: $LINE" >&2
        continue
    fi

    if [[ -z "${LABEL:-}" ]]; then
        LABEL="${CHROM}_${START}_${END}"
    fi

    REGION="${CHROM}:${START}-${END}"
    OUTPUT_PREFIX="${OUTPUT_DIR}/${LABEL}_methylartist_plot"

    echo "Processing ${LABEL} (${REGION})"

    tmp_log=$(mktemp)
    if ! methylartist locus \
        -b "$BAM" \
        -i "$REGION" \
        --samplepalette viridis \
        --labelgenes \
        --ref "$REFERENCE" \
        --motif CG \
        -m m \
        --show_transcripts \
        --exonheight 0.6 \
        --panelratios 2,5,1,3,3 \
        -o "$OUTPUT_PREFIX" >"$tmp_log" 2>&1; then
        if grep -qi "no mods found" "$tmp_log" || grep -qi "missing MM/ML" "$tmp_log"; then
            echo "Warning: no modification tags found for region $REGION; skipping $LABEL" >&2
            rm -f "$tmp_log"
            continue
        else
            cat "$tmp_log" >&2
            rm -f "$tmp_log"
            exit 1
        fi
    fi
    rm -f "$tmp_log"

    echo "Finished ${LABEL}"
done < "$BED_FILE"