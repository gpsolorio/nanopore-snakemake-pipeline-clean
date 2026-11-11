#!/bin/bash
set -euo pipefail

##############################################################################
# run_methylartist_mini.sh
#
# Mini-data variant of run_methylartist.sh. Adaptive/targeted sequencing at
# a single small locus produces full production-scale READ DEPTH even
# though the mini reference is tiny (e.g. ~70x over 186kb -> ~8-15k reads
# per sample). methylartist draws one track per read, so plotting cost
# scales with read COUNT, not genome/window size -- "mini" alone doesn't
# make this step fast.
#
# This variant:
#   1. Subsamples each input BAM down to a target read count before
#      plotting (plotting-only; never touches the real variant-calling
#      BAMs -- subsampled copies are written to a temp dir and discarded).
#   2. Runs independent BED rows (e.g. multiple overlapping transcript
#      isoforms of the same gene) in parallel, bounded by --threads,
#      since methylartist locus itself has no general multithreading flag
#      (only --coverprocs, which is scoped to -c/--plot_coverage).
#
# Use the original run_methylartist.sh for production/full-depth runs.
##############################################################################

BAM=""
BED_FILE=""
REFERENCE=""
OUTPUT_DIR=""
THREADS="1"
MAX_READS="2500"   # target reads per sample after subsampling
SEED="42"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam) BAM="$2"; shift 2 ;;
        --bed-file) BED_FILE="$2"; shift 2 ;;
        --reference) REFERENCE="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --max-reads) MAX_READS="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$BAM" || -z "$BED_FILE" || -z "$REFERENCE" || -z "$OUTPUT_DIR" ]]; then
    echo "Usage: $0 --bam BAM1,BAM2,... --bed-file BED_FILE --reference REFERENCE --output-dir OUTPUT_DIR [--threads N] [--max-reads N] [--seed N]" >&2
    exit 1
fi
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
    echo "--threads must be a positive integer" >&2; exit 1
fi
if ! [[ "$MAX_READS" =~ ^[0-9]+$ ]] || [[ "$MAX_READS" -lt 1 ]]; then
    echo "--max-reads must be a positive integer" >&2; exit 1
fi

mkdir -p "$OUTPUT_DIR"
TMP_BAM_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_BAM_DIR"' EXIT

echo "Using BAMs:"
IFS=',' read -r -a BAMS <<< "$BAM"
printf '  %s\n' "${BAMS[@]}"
echo "Total: ${#BAMS[@]} sample(s)"
echo "Subsampling target: ~${MAX_READS} reads/sample (plotting only; source BAMs untouched)"

SUBSAMPLED_BAMS=()
for src_bam in "${BAMS[@]}"; do
    sample_name="$(basename "$src_bam" .bam)"
    n_reads=$(samtools view -c -F 0x900 "$src_bam")

    if [[ "$n_reads" -le "$MAX_READS" ]]; then
        echo "  [${sample_name}] ${n_reads} reads <= target, using as-is"
        SUBSAMPLED_BAMS+=("$src_bam")
        continue
    fi

    frac=$(awk -v n="$n_reads" -v target="$MAX_READS" 'BEGIN { printf "%.6f", target / n }')
    out_bam="${TMP_BAM_DIR}/${sample_name}_sub.bam"

    echo "  [${sample_name}] ${n_reads} reads -> subsampling to fraction ${frac} (~${MAX_READS} reads)"
    samtools view -b -s "${SEED}.${frac#0.}" "$src_bam" > "$out_bam"
    samtools index "$out_bam"
    SUBSAMPLED_BAMS+=("$out_bam")
done

SUBSAMPLED_BAM_ARG=$(IFS=,; echo "${SUBSAMPLED_BAMS[*]}")

pids=()

run_one_region() {
    local CHROM="$1" START="$2" END="$3" LABEL="$4"
    local REGION="${CHROM}:${START}-${END}"
    local OUTPUT_PREFIX="${OUTPUT_DIR}/${LABEL}_methylartist_plot"

    echo "Processing ${LABEL} (${REGION})"

    local tmp_log
    tmp_log=$(mktemp)
    if ! methylartist locus \
        -b "$SUBSAMPLED_BAM_ARG" \
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
            rm -f "$tmp_log"; return 0
        else
            cat "$tmp_log" >&2; rm -f "$tmp_log"; return 1
        fi
    fi
    rm -f "$tmp_log"
    echo "Finished ${LABEL}"
}

while IFS= read -r LINE; do
    [[ -z "${LINE:-}" ]] && continue
    [[ "${LINE:0:1}" == "#" ]] && continue
    read -r CHROM START END LABEL _ <<< "$LINE"

    if [[ -z "${CHROM:-}" || -z "${START:-}" || -z "${END:-}" ]]; then
        echo "Skipping malformed BED row (expected >=3 columns): $LINE" >&2; continue
    fi
    if ! [[ "$START" =~ ^[0-9]+$ && "$END" =~ ^[0-9]+$ ]]; then
        echo "Skipping BED row with non-numeric coordinates: $LINE" >&2; continue
    fi
    [[ -z "${LABEL:-}" ]] && LABEL="${CHROM}_${START}_${END}"

    run_one_region "$CHROM" "$START" "$END" "$LABEL" &
    pids+=("$!")

    # Throttle to --threads concurrent jobs. Deliberately avoids `wait -n`
    # (bash >=4.3 only) so this works under macOS's stock bash 3.2 as well
    # as Linux. A short poll is fine here: job counts are small (a handful
    # of BED rows) and this loop only spins while at the concurrency cap.
    while (( $(jobs -r -p | wc -l) >= THREADS )); do
        sleep 0.2
    done
done < "$BED_FILE"

fail=0
for pid in "${pids[@]}"; do
    wait "$pid" || fail=1
done
exit "$fail"