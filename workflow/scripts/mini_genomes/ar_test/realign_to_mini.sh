#!/bin/bash
set -euo pipefail

##############################################################################
# realign_to_mini.sh
#
# Converts a subset BAM back to fastq (preserving MM/ML mod tags), then
# realigns those reads to the mini genome with multiple aligners for
# side-by-side comparison.
#
# Usage: realign_to_mini.sh <label> [mini_genome_fa]
#   label          : must match the label used in extract_ar_reads.sh
#   mini_genome_fa : optional override; defaults to the tracked genome at
#                    resources/mini_ar_test/genome/mini_genome.fa
#
# Output (under $MINI_OUT_ROOT, default results/mini-test/):
#   aligned/<label>/reads.fastq
#   aligned/<label>/minimap2.bam(.bai)
#   aligned/<label>/winnowmap2.bam(.bai)
#   genome_cache/                      (meryl index, regenerated as needed)
##############################################################################

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <label> [mini_genome_fa]" >&2
    exit 1
fi

LABEL="$1"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"

# MINI_OUT_ROOT: where bam/ and aligned/ get written -- defaults to
# results/mini-test (disposable, safe to rm -rf) for manual runs, same as
# what Snakemake sets explicitly.
MINI_OUT_ROOT="${MINI_OUT_ROOT:-${REPO_ROOT}/results/mini-test}"
BAM_DIR="${MINI_OUT_ROOT}/bam"
# GENOME_DIR: the curated, tracked reference -- separate override, defaults
# to the tracked location under resources/, never the disposable output root.
GENOME_DIR="${MINI_GENOME_DIR:-${REPO_ROOT}/resources/mini_ar_test/genome}"
ALIGN_DIR="${MINI_OUT_ROOT}/aligned/${LABEL}"

MINI_GENOME="${2:-${GENOME_DIR}/mini_genome.fa}"
SUBSET_BAM="${BAM_DIR}/${LABEL}_sorted.bam"

if [[ ! -f "$SUBSET_BAM" ]]; then
    echo "Expected subset BAM not found: $SUBSET_BAM" >&2
    echo "(run extract_ar_reads.sh with label '${LABEL}' first)" >&2
    exit 1
fi
if [[ ! -f "$MINI_GENOME" ]]; then
    echo "Mini genome not found: $MINI_GENOME" >&2
    exit 1
fi

mkdir -p "$ALIGN_DIR"

echo "Converting $SUBSET_BAM -> fastq (preserving MM/ML tags)"
samtools fastq -T MM,ML "$SUBSET_BAM" > "${ALIGN_DIR}/reads.fastq"

# --- minimap2 ---
echo "Aligning with minimap2..."
minimap2 -ax map-ont -y --MD "$MINI_GENOME" "${ALIGN_DIR}/reads.fastq" \
    | samtools sort -o "${ALIGN_DIR}/minimap2.bam"
samtools index "${ALIGN_DIR}/minimap2.bam"

# --- winnowmap2 (repetitive k-mer file is cached under the disposable
# output root, NOT resources/mini_ar_test/genome/, which stays tracked-only:
# just mini_genome.fa(.fai). Keeping generated meryl artifacts out of the
# tracked genome dir avoids mixing curated + generated content again.)
MERYL_CACHE_DIR="${MINI_OUT_ROOT}/genome_cache"
mkdir -p "$MERYL_CACHE_DIR"
REPETITIVE_KMERS="${MERYL_CACHE_DIR}/mini_genome.repetitive_k15.txt"
if [[ ! -f "$REPETITIVE_KMERS" ]]; then
    echo "Building repetitive k-mer profile for winnowmap2 (one-time per genome)..."
    meryl count k=15 output "${MERYL_CACHE_DIR}/mini_genome.meryl" "$MINI_GENOME"
    meryl print greater-than distinct=0.9998 "${MERYL_CACHE_DIR}/mini_genome.meryl" \
        > "$REPETITIVE_KMERS"
fi

echo "Aligning with winnowmap2..."
winnowmap -W "$REPETITIVE_KMERS" -ax map-ont -y --MD \
    "$MINI_GENOME" "${ALIGN_DIR}/reads.fastq" \
    | samtools sort -o "${ALIGN_DIR}/winnowmap2.bam"
samtools index "${ALIGN_DIR}/winnowmap2.bam"

echo ""
echo "=== Comparison: ${LABEL} ==="
for bam in "${ALIGN_DIR}"/*.bam; do
    echo "--- $(basename "$bam") ---"
    samtools flagstat "$bam"
done