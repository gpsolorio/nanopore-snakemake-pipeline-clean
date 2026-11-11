#!/usr/bin/env bash
#SBATCH --job-name=modkit_pileup
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=400G
#SBATCH --output=logs/slurm_modkit_pileup.out
#SBATCH --error=logs/slurm_modkit_pileup.err

set -euo pipefail

# ============================================================================
# Modkit Pileup - Generate bedMethyl files from Nanopore BAM files
# ============================================================================
#
# Usage:
#   sbatch modkit_pileup_combined.sh --bam-dir PATH --reference PATH [--output-dir PATH]
#
# ============================================================================

BAM_DIR=""
REFERENCE=""
OUTPUT_DIR=""

# ============================================================================
# Help function
# ============================================================================

show_help() {
    cat << EOF
Modkit Pileup - Generate bedMethyl files from Nanopore BAM files

Required options:
  --bam-dir PATH      Directory containing input BAM files
  --reference PATH    Reference genome FASTA file

Optional:
  --output-dir PATH   Output directory for bedMethyl files (default: ./01_modkit_pileup_beds)
  --help              Show this help message

Example:
  sbatch modkit_pileup_combined.sh \
    --bam-dir /path/to/bams \
    --reference /path/to/hg38.fa \
    --output-dir /path/to/output
EOF
    exit 0
}

# ============================================================================
# Parse command-line arguments
# ============================================================================

while [[ $# -gt 0 ]]; do
    case $1 in
        --bam-dir)
            BAM_DIR="$2"
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
        --help|-h)
            show_help
            ;;
        *)
            echo "Unknown option: $1" >&2
            show_help
            ;;
    esac
done

# ============================================================================
# Validate required arguments
# ============================================================================

if [[ -z "$BAM_DIR" ]]; then
    echo "ERROR: --bam-dir is required" >&2
    exit 1
fi

if [[ -z "$REFERENCE" ]]; then
    echo "ERROR: --reference is required" >&2
    exit 1
fi

# ============================================================================
# Setup paths
#
# This script lives at analysis/dmr/01_prepare_beds/modkit_pileup_combined.sh
# -- 3 levels below repo root, not 4. The previous version walked up 4
# levels and pointed at workflow/scripts/analysis/differential_methylation/,
# a path that didn't exist anywhere in the repo.
# ============================================================================

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"   # analysis/dmr/01_prepare_beds -> repo root
DMR_DIR="$PROJECT_ROOT/analysis/dmr"
OUT_DIR="${OUTPUT_DIR:-$DMR_DIR/01_prepare_beds/output}"
LOG_DIR="$OUT_DIR/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# Redirect output to log files
exec >"$LOG_DIR/modkit_pileup.out" 2>"$LOG_DIR/modkit_pileup.err"

echo "Script started at $(date)"
echo "Log directory: $LOG_DIR"
echo "Output directory: $OUT_DIR"

# ============================================================================
# Activate conda environment
# ============================================================================

source "$(conda info --base)/etc/profile.d/conda.sh"
ENV_YAML="$PROJECT_ROOT/workflow/envs/modkit.yml"
ENV_NAME="tool-modkit"

if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "Creating conda environment from $ENV_YAML..."
    conda env create -n "$ENV_NAME" -f "$ENV_YAML"
fi

echo "Activating conda environment: $ENV_NAME"
conda activate "$ENV_NAME"

# ============================================================================
# Validation
# ============================================================================

if [[ ! -d "$BAM_DIR" ]]; then
    echo "ERROR: BAM directory does not exist: $BAM_DIR" >&2
    exit 1
fi

if [[ ! -f "$REFERENCE" ]]; then
    echo "ERROR: Reference genome does not exist: $REFERENCE" >&2
    exit 1
fi

if ! command -v modkit &>/dev/null; then
    echo "ERROR: modkit not found in PATH" >&2
    exit 1
fi

# ============================================================================
# Run modkit pileup on each BAM file
# ============================================================================

START_TIME=$(date +%s)
PROCESSED=0
SKIPPED=0
FAILED=0
THREADS="${SLURM_CPUS_PER_TASK}"

for INPUT_BAM in "$BAM_DIR"/*_5mCG-5hmCG.sorted.bam; do
    [[ -e "$INPUT_BAM" ]] || continue  # skip if no files

    BASENAME=$(basename "$INPUT_BAM" .sorted.bam)
    OUTPUT_BED="${OUT_DIR}/${BASENAME}_modkit_calls.bed"
    PILEUP_LOG="${LOG_DIR}/${BASENAME}_modkit.log"

    echo "Processing $BASENAME ..."

    if [[ ! -f "${INPUT_BAM}.bai" ]]; then
        echo "SKIPPED — BAM index not found: ${INPUT_BAM}.bai"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    if [[ -f "$OUTPUT_BED" ]]; then
        echo "SKIPPED — Output already exists: $OUTPUT_BED"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    if modkit pileup "$INPUT_BAM" "$OUTPUT_BED" \
        --log-filepath "$PILEUP_LOG" \
        --cpg \
        --combine-strands \
        --modified-bases 5mC \
        --ref "$REFERENCE" \
        --threads "$THREADS"; then

        if [[ -s "$OUTPUT_BED" ]]; then
            echo "SUCCESS — $(wc -l < "$OUTPUT_BED") records"
            PROCESSED=$((PROCESSED + 1))
        else
            echo "FAILED — Output file empty"
            FAILED=$((FAILED + 1))
        fi
    else
        echo "FAILED — modkit command failed"
        FAILED=$((FAILED + 1))
    fi
done

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo "Processed: $PROCESSED, Skipped: $SKIPPED, Failed: $FAILED, Elapsed: $ELAPSED seconds"

if [[ $FAILED -gt 0 ]]; then
    exit 1
fi

echo "Modkit pileup completed successfully!"
exit 0
