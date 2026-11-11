#!/bin/bash
set -euo pipefail

##############################################################################
# extract_ar_reads.sh
#
# Merges only the per-batch BAM chunks in a sample's pass/ directory that
# overlap the region of interest, sorts and indexes the result.
#
# Usage: extract_ar_reads.sh <sample_pass_dir> <region> <label> [batch_size]
#   sample_pass_dir : e.g. .../sample_001_5mC-5hmC/pass
#   region          : ORIGINAL hg38 coordinates, e.g. chrX:67544032-67730619
#   label           : short name for output, e.g. sample_001_AR
#   batch_size      : optional, default 500
#
# Output (relative to this script's location):
#   ../mini-data/bam/<label>_sorted.bam(.bai)
#
# NOTE ON BATCHING: this cluster's hard open-file-descriptor ceiling is
# fixed at 16384 (confirmed via `ulimit -Hn`), and a single sample's pass/
# directory can contain 30k+ chunk BAMs. A single `samtools merge -b filelist`
# call needs close to one fd per listed input, so above ~16k inputs it WILL
# fail with "Too many open files" no matter what -- there is no ulimit value
# available on this cluster that fixes it. Merging in batches keeps peak
# open files at ~batch_size regardless of how many total chunks exist.
##############################################################################

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <sample_pass_dir> <region> <label> [batch_size]" >&2
    exit 1
fi

PASS_DIR="$1"
REGION="$2"
LABEL="$3"
BATCH_SIZE="${4:-500}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"   # workflow/scripts/mini_genomes/ar_test -> repo root
# MINI_OUT_ROOT controls where bam/ (and its batches) get written -- defaults
# to results/mini-test (disposable, safe to rm -rf) for manual/ad hoc runs,
# same as what Snakemake sets explicitly. Never defaults into resources/,
# which holds only the tracked, curated genome/beds.
MINI_OUT_ROOT="${MINI_OUT_ROOT:-${REPO_ROOT}/results/mini-test}"
BAM_DIR="${MINI_OUT_ROOT}/bam"
mkdir -p "$BAM_DIR"

OUT_PREFIX="${BAM_DIR}/${LABEL}"
TMP_LIST=$(mktemp)
BATCH_BAM_LIST=$(mktemp)
BATCH_DIR="${OUT_PREFIX}_batches"
trap 'rm -f "$TMP_LIST" "$BATCH_BAM_LIST"; rm -rf "$BATCH_DIR"' EXIT

echo "Listing chunks in: $PASS_DIR"
find "$PASS_DIR" -maxdepth 1 -name '*.bam' > "$TMP_LIST"
n_total=$(wc -l < "$TMP_LIST")
echo "Chunks found: $n_total"

if [[ "$n_total" -eq 0 ]]; then
    echo "No BAM files found in $PASS_DIR — aborting." >&2
    exit 1
fi

# Merge in batches so peak simultaneously-open files stays ~batch_size,
# regardless of n_total. samtools merge -R still consults each input's
# .bai to determine region overlap -- that part isn't free -- but batching
# avoids ever exceeding the fd ceiling.
mkdir -p "$BATCH_DIR"
split -l "$BATCH_SIZE" -d -a 5 "$TMP_LIST" "${BATCH_DIR}/chunklist_"

n_batches=0
for chunklist in "${BATCH_DIR}"/chunklist_*; do
    n_batches=$((n_batches + 1))
    batch_out="${chunklist}.bam"
    samtools merge -R "$REGION" -b "$chunklist" -f "$batch_out"
    echo "$batch_out" >> "$BATCH_BAM_LIST"
    rm -f "$chunklist"
done
echo "Merged $n_total chunks into $n_batches batch BAM(s) (batch size ${BATCH_SIZE})"

echo "Combining $n_batches batch BAM(s) into final merged BAM..."
samtools merge -f -b "$BATCH_BAM_LIST" "${OUT_PREFIX}_merged.bam"

samtools sort -o "${OUT_PREFIX}_sorted.bam" "${OUT_PREFIX}_merged.bam"
samtools index "${OUT_PREFIX}_sorted.bam"
rm -f "${OUT_PREFIX}_merged.bam"

echo "Done -> ${OUT_PREFIX}_sorted.bam"
samtools flagstat "${OUT_PREFIX}_sorted.bam"