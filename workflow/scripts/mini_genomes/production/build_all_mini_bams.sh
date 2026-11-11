#!/usr/bin/env bash
#
# build_all_mini_bams.sh
#
# Runs the FULL mini-BAM chain for every sample in your cohort:
#
#   1. extract_mini_bam.sh   -- indexed region slice out of hg38 chunks
#                                (output still has hg38 RNAME/POS + header)
#   2. remap_bam_to_mini.py  -- rewrites header + RNAME/POS to match
#                                resources/mini_genome.fa (REQUIRED --
#                                without this, modkit/samtools reject the
#                                header/reference-length mismatch)
#   3. samtools sort + index -- remap's output is coordinate-shifted but
#                                not guaranteed sorted; this produces the
#                                final BAM the pipeline actually reads
#
# Final output lands at workflow/data/bam/{sample}/{sample}.bam, matching
# the layout the pipeline expects.
#
# Usage: ./build_all_mini_bams.sh samples.tsv
# samples.tsv: two tab-separated columns, no header: sample_id <TAB> bam_glob
# (see examples/samples.example.tsv)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_ROOT="workflow/data/bam"
WINDOW_SPEC_DIR="${WINDOW_SPEC_DIR:-workflow/scripts/mini_genomes/common}"
THREADS="${THREADS:-4}"


run_one() {
  local glob="$1"
  local sample="$2"
  local out_dir="$OUT_ROOT/${sample}"
  local final_out="${out_dir}/${sample}.bam"
  local raw_slice="${out_dir}/${sample}.hg38_slice.bam"
  local remapped_unsorted="${out_dir}/${sample}.remapped.unsorted.bam"

  echo "=========================================================="
  echo "Sample: $sample"
  echo "Source: $glob"
  echo "Target: $final_out"
  echo "=========================================================="

  mkdir -p "$out_dir"

  # Step 1: indexed region slice (still hg38 coordinates/header)
  echo "[1/3] extract_mini_bam.sh -> $raw_slice"
  "$SCRIPT_DIR/extract_mini_bam.sh" "$glob" "$raw_slice" "$THREADS"

  # Step 2: rewrite RNAME/POS/header to mini-genome coordinates
  echo "[2/3] remap_bam_to_mini.py -> $remapped_unsorted"
  python3 "$SCRIPT_DIR/remap_bam_to_mini.py" \
    "$WINDOW_SPEC_DIR" "$raw_slice" "$remapped_unsorted"

  # Step 3: sort + index the remapped BAM -- this is the real final output
  echo "[3/3] sort + index -> $final_out"
  samtools sort -@ "$THREADS" -o "$final_out" "$remapped_unsorted"
  samtools index "$final_out"

  echo "[$sample] done. Validating against mini reference header..."
  if command -v diff >/dev/null 2>&1 && [ -f resources/mini_genome.fa.fai ]; then
    diff \
      <(samtools view -H "$final_out" | grep '^@SQ' | sed 's/.*SN:\([^\t]*\).*LN:\([0-9]*\).*/\1\t\2/' | sort) \
      <(cut -f1,2 resources/mini_genome.fa.fai | sort) \
      && echo "[$sample] header matches mini_genome.fa.fai exactly." \
      || echo "[$sample] WARNING: header mismatch vs mini_genome.fa.fai -- investigate before trusting this BAM."
  fi

  samtools flagstat "$final_out"

  # Clean up intermediates -- comment these out if you want to keep them
  # for debugging a failed remap/sort step.
  rm -f "$raw_slice" "${raw_slice}.bai" "$remapped_unsorted"
}

MANIFEST="${1:?Usage: $0 samples.tsv  (columns: sample_id<TAB>bam_glob)}"
[[ -f "$MANIFEST" ]] || { echo "Manifest not found: $MANIFEST" >&2; exit 1; }

while IFS=$'\t' read -r sample glob; do
  [[ -z "$sample" || "$sample" == \#* ]] && continue
  run_one "$glob" "$sample"
done < "$MANIFEST"
