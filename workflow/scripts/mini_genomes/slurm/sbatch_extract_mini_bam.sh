#!/usr/bin/env bash
#SBATCH --job-name=extract_mini_bam
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=./logs/extract_mini_bam_%j.out
#SBATCH --error=./logs/extract_mini_bam_%j.err

# sbatch_extract_mini_bam.sh
#
# SLURM wrapper for the FULL mini-BAM chain for one sample:
#
#   0. pre-flight checks     -- confirm the sample directory exists and
#                                actually contains *.bam files before
#                                doing any work
#   1. extract_mini_bam.sh   -- indexed region slice out of hg38 chunks
#                                (output still has hg38 RNAME/POS + header)
#   2. remap_bam_to_mini.py  -- REQUIRED coordinate/header rewrite to match
#                                resources/mini_genome.fa -- extraction alone
#                                does NOT do this (samtools view <region> is
#                                a read filter, not a coordinate transform)
#   3. samtools sort + index -- remap's output is coordinate-shifted but
#                                not guaranteed sorted; this produces the
#                                final BAM the pipeline actually reads
#   4. header validation      -- diff against resources/mini_genome.fa.fai,
#                                fail loud on any mismatch rather than
#                                silently producing a bad BAM
#
# Known directory shapes (confirmed via ls in-session, NOT assumed):
#   layout A:   {sample}_5mC-5hmC/pass/*.bam
#   layout B:   barcode0N_{sample}_5mC-5hmC/*.bam   (no pass/)
#   layout C:    {sample}_5mC-5hmC/*.bam             (no pass/, no barcode prefix)
#
# This script does NOT guess which shape applies -- it is handed the
# already-resolved directory that DIRECTLY contains *.bam files (the
# submission wrapper resolves pass/ vs. no-pass/ before calling this),
# and re-verifies BAMs exist there before doing any work.
#
# Usage:
#   sbatch --job-name=<sample> workflow/scripts/mini_genomes/slurm/sbatch_extract_mini_bam.sh \
#       "<bam_source_dir>" <sample_id>

set -euo pipefail

BAM_SOURCE_DIR="${1:?Usage: sbatch ... $0 <bam_source_dir> <sample_id>}"
SAMPLE="${2:?Usage: sbatch ... $0 <bam_source_dir> <sample_id>}"

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$REPO_ROOT"

# Explicitly activate conda env -- SLURM batch jobs do NOT reliably
# inherit the interactive shell's active environment.
module load conda 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ngs-utils-env

if ! command -v samtools &> /dev/null; then
    echo "FAIL: samtools not found on PATH after activating ngs-utils-env." >&2
    exit 1
fi

SCRIPT_DIR="workflow/scripts/mini_genomes/production"
WINDOW_SPEC_DIR="${WINDOW_SPEC_DIR:-workflow/scripts/mini_genomes/common}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

OUT_DIR="workflow/data/bam/${SAMPLE}"
FINAL_OUT="${OUT_DIR}/${SAMPLE}.bam"
RAW_SLICE="${OUT_DIR}/${SAMPLE}.hg38_slice.bam"
REMAPPED_UNSORTED="${OUT_DIR}/${SAMPLE}.remapped.unsorted.bam"

LOGDIR="logs"
mkdir -p "$LOGDIR" "$OUT_DIR"

echo "=========================================================="
echo "Sample:     $SAMPLE"
echo "Source dir: $BAM_SOURCE_DIR"
echo "Target:     $FINAL_OUT"
echo "Threads:    $THREADS"
echo "Node:       $(hostname)"
echo "=========================================================="

# ---- Pre-flight checks -----------------------------------------------
echo "[0/4] pre-flight checks"

if [ ! -d "$BAM_SOURCE_DIR" ]; then
    echo "FAIL: $BAM_SOURCE_DIR does not exist -- check the path passed in." >&2
    exit 1
fi

n_bams=$(find "$BAM_SOURCE_DIR" -maxdepth 1 -name "*.bam" | wc -l | tr -d ' ')
if [ "$n_bams" -eq 0 ]; then
    echo "FAIL: no *.bam files directly in $BAM_SOURCE_DIR" >&2
    echo "      If this sample uses a pass/ subfolder, the submission" >&2
    echo "      wrapper should have resolved that already -- check" >&2
    echo "      submit_mini_bam_batch.sh's directory detection." >&2
    exit 1
fi
echo "  found $n_bams BAM file(s) in $BAM_SOURCE_DIR"

if [ ! -f resources/mini_genome.fa.fai ]; then
    echo "FAIL: resources/mini_genome.fa.fai not found -- mini reference" >&2
    echo "      must exist before extraction/remapping can be validated." >&2
    exit 1
fi

if [ ! -f "$SCRIPT_DIR/extract_mini_bam.sh" ] || [ ! -f "$SCRIPT_DIR/remap_bam_to_mini.py" ]; then
    echo "FAIL: expected scripts not found under $SCRIPT_DIR" >&2
    exit 1
fi

# ---- Step 1: extract -----------------------------------------------------
echo "[1/4] extract_mini_bam.sh -> $RAW_SLICE"
"$SCRIPT_DIR/extract_mini_bam.sh" "${BAM_SOURCE_DIR}/*.bam" "$RAW_SLICE" "$THREADS"

# ---- Step 2: remap --------------------------------------------------------
echo "[2/4] remap_bam_to_mini.py -> $REMAPPED_UNSORTED"
python3 "$SCRIPT_DIR/remap_bam_to_mini.py" \
    "$WINDOW_SPEC_DIR" "$RAW_SLICE" "$REMAPPED_UNSORTED"

# ---- Step 3: sort + index -------------------------------------------------
echo "[3/4] sort + index -> $FINAL_OUT"
samtools sort -@ "$THREADS" -o "$FINAL_OUT" "$REMAPPED_UNSORTED"
samtools index "$FINAL_OUT"

# ---- Step 4: validate header ----------------------------------------------
echo "[4/4] validating header against resources/mini_genome.fa.fai"
if diff \
    <(samtools view -H "$FINAL_OUT" | awk -F'\t' '$1=="@SQ"' \
        | sed -E 's/.*SN:([^\t]*).*LN:([0-9]*).*/\1\t\2/' | sort) \
    <(cut -f1,2 resources/mini_genome.fa.fai | sort) \
    > "${OUT_DIR}/${SAMPLE}.header_diff.txt"
then
    echo "[$SAMPLE] OK: header matches mini_genome.fa.fai exactly."
    rm -f "${OUT_DIR}/${SAMPLE}.header_diff.txt"
else
    echo "[$SAMPLE] FAIL: header mismatch vs mini_genome.fa.fai." >&2
    echo "        See ${OUT_DIR}/${SAMPLE}.header_diff.txt -- do NOT" >&2
    echo "        trust this BAM until resolved." >&2
    exit 1
fi

echo "[$SAMPLE] read count / flag summary:"
samtools flagstat "$FINAL_OUT"

# Clean up intermediates on success; they're left in place if any
# earlier step exited non-zero (set -e stops the script before this line).
rm -f "$RAW_SLICE" "${RAW_SLICE}.bai" "$REMAPPED_UNSORTED"

echo "[$SAMPLE] done."