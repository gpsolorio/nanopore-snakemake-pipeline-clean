#!/usr/bin/env bash
#
# extract_mini_bam.sh
#
# Build a "mini" BAM for one sample by pulling only reads that overlap
# the mini-genome regions, using indexed random access on each chunk
# BAM individually -- NOT by merging everything first.
#
# Why this matters: the private cohort sample dirs can be 258GB+ across hundreds of
# small per-runid chunk BAMs; large cohort barcode dirs are ~14GB across a
# dozen large chunks. Every chunk already has a .bai next to it, so
# `samtools view <region args>` on each chunk uses the BAM index (BAI)
# to seek directly to the relevant blocks instead of decompressing the
# whole file. Merging first (as in the original plan) would force you
# to stream the entire dataset through `samtools merge` before ever
# discarding the >99% you don't need -- exactly backwards.
#
# Usage:
#   ./extract_mini_bam.sh <input_bam_glob> <output_bam> [threads]
#
# Example (one sample):
#   ./extract_mini_bam.sh \
#     "/path/to/run1/pass/*.bam" \
#     mini_data/sample_001.bam
#
# Example (one barcode):
#   ./extract_mini_bam.sh \
#     "/path/to/run1/pass/*.bam" \
#     mini_data/barcode01.bam

set -euo pipefail

INPUT_GLOB="${1:?Usage: $0 <input_bam_glob> <output_bam> [threads]}"
OUTPUT_BAM="${2:?Usage: $0 <input_bam_glob> <output_bam> [threads]}"
THREADS="${3:-4}"

# Regions are generated from workflow/scripts/mini_genomes/common/window_spec.py, the single
# source of truth also used to build mini_genome.fa -- this guarantees the
# BAM slice always matches the reference windows, including any locus
# windows (AR_locus, BCOR_locus, etc.) added after the fact.
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
WINDOW_SPEC_DIR="${WINDOW_SPEC_DIR:-$REPO_ROOT/workflow/scripts/mini_genomes/common}"

REGIONS_FILE=$(mktemp)

python3 "$WINDOW_SPEC_DIR/print_hg38_regions.py" > "$REGIONS_FILE"

if [ ! -s "$REGIONS_FILE" ]; then
  echo "ERROR: region generation failed" >&2
  exit 1
fi

mapfile -t MINI_REGIONS < "$REGIONS_FILE"

if [ "${#MINI_REGIONS[@]}" -eq 0 ]; then
  echo "ERROR: no regions generated -- check WINDOW_SPEC_DIR=$WINDOW_SPEC_DIR" >&2
  exit 1
fi
echo "[extract_mini_bam] using ${#MINI_REGIONS[@]} regions from window_spec.py"

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

mkdir -p "$(dirname "$OUTPUT_BAM")"

echo "[extract_mini_bam] scanning: $INPUT_GLOB"
mapfile -t CHUNKS < <(compgen -G "$INPUT_GLOB" || true)

if [ "${#CHUNKS[@]}" -eq 0 ]; then
  echo "ERROR: no BAM files matched $INPUT_GLOB" >&2
  exit 1
fi
echo "[extract_mini_bam] found ${#CHUNKS[@]} chunk BAM(s)"

extract_count=0
extracted_files=()

for chunk in "${CHUNKS[@]}"; do
  # Skip files that have no index -- region-based random access needs one.
  if [ ! -f "${chunk}.bai" ] && [ ! -f "${chunk%.bam}.bai" ]; then
    echo "  [skip, no index] $chunk" >&2
    continue
  fi

  out="$WORKDIR/$(basename "${chunk%.bam}").mini.bam"

  # Indexed random-access extraction: only reads the blocks that
  # overlap these regions, does not decompress the whole chunk.
  if samtools view -@ "$THREADS" -b "$chunk" "${MINI_REGIONS[@]}" -o "$out"; then
    n_reads=$(samtools view -c "$out")
    if [ "$n_reads" -gt 0 ]; then
      extracted_files+=("$out")
      extract_count=$((extract_count + n_reads))
    fi
  fi
done

echo "[extract_mini_bam] extracted $extract_count reads from ${#extracted_files[@]} chunk(s)"

if [ "${#extracted_files[@]}" -eq 0 ]; then
  echo "ERROR: no reads found in any chunk for the mini regions." >&2
  echo "Check that MINI_REGIONS matches make_mini_genome.py's REGIONS spec." >&2
  exit 1
fi

echo "[extract_mini_bam] merging + sorting -> $OUTPUT_BAM (batched, ${#extracted_files[@]} files)"

# Raise the fd limit as far as the hard limit allows -- helps, but batching
# below is what actually makes this scale past ulimit regardless of value.
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true

BATCH_SIZE=500
BATCH_DIR="$WORKDIR/batches"
mkdir -p "$BATCH_DIR"
batch_files=()
batch_num=0
i=0

while [ "$i" -lt "${#extracted_files[@]}" ]; do
  batch=("${extracted_files[@]:i:BATCH_SIZE}")
  batch_out="$BATCH_DIR/batch_${batch_num}.bam"
  samtools cat -o "$batch_out" "${batch[@]}"
  batch_files+=("$batch_out")
  batch_num=$((batch_num + 1))
  i=$((i + BATCH_SIZE))
done

echo "[extract_mini_bam] merged into ${#batch_files[@]} batch file(s), now final cat+sort"
samtools cat "${batch_files[@]}" \
  | samtools sort -@ "$THREADS" -o "$OUTPUT_BAM" -

samtools index "$OUTPUT_BAM"

echo "[extract_mini_bam] done."
samtools flagstat "$OUTPUT_BAM"
