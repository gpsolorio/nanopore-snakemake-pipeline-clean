#!/bin/bash
set -euo pipefail

##############################################################################
# make_mini_genome.sh
#
# Extracts a region from hg38.fa into its own mini-genome contig, written
# to the tracked resources/mini_ar_test/genome/ location.
#
# Usage: make_mini_genome.sh <chrom> <start> <end> <contig_name> [hg38_fa]
#   chrom/start/end : BED-style (0-based, half-open) coordinates
#   contig_name     : name for the new contig, e.g. AR_locus
#   hg38_fa         : required: path to hg38.fa
#
# NOTE: this OVERWRITES mini_genome.fa rather than appending, to avoid
# duplicate-contig corruption when rebuilding the same locus with corrected
# coordinates. If you need a genuine multi-locus genome, call this once per
# locus against a fresh empty file and manage appends deliberately.
##############################################################################

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 <chrom> <start> <end> <contig_name> [hg38_fa]" >&2
    exit 1
fi

CHROM="$1"
START="$2"
END="$3"
NAME="$4"
HG38_FA="${5:?Pass the path to hg38.fa as argument 5}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"   # workflow/scripts/mini_genomes/ar_test -> repo root
GENOME_DIR="${REPO_ROOT}/resources/mini_ar_test/genome"
OUT_FA="${GENOME_DIR}/mini_genome.fa"

mkdir -p "$GENOME_DIR"

echo "Extracting ${CHROM}:${START}-${END} -> contig '${NAME}'"
samtools faidx "$HG38_FA" "${CHROM}:$((START+1))-${END}" \
    | sed "1s/.*/>${NAME}/" > "$OUT_FA"

rm -f "${OUT_FA}.fai"
samtools faidx "$OUT_FA"

echo "Wrote: ${OUT_FA}"
cat "${OUT_FA}.fai"