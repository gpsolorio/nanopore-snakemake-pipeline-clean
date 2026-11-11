#!/bin/bash
##############################################################################
# SLURM Submission Script for dmrseq DMR Analysis
##############################################################################
#
# Usage: sbatch submit_dmrseq.sh [DATASET] [CUTOFF] [MIN_CPGS] [MAX_PERMS] [CORES]
#
# Examples:
#   sbatch submit_dmrseq.sh mini 0.1 5 10 8
#   sbatch submit_dmrseq.sh full 0.1 5 500 32
#   sbatch submit_dmrseq.sh              # Uses defaults
#
##############################################################################

#SBATCH --job-name=dmrseq_dmr
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --output=logs/slurm/dmrseq_dmr_%j.out
#SBATCH --error=logs/slurm/dmrseq_dmr_%j.err
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

# ============================================================================
# Setup
# ============================================================================

echo "=========================================="
echo "dmrseq DMR Analysis - SLURM Submission"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Date: $(date)"
echo "=========================================="
echo ""

# Parse command-line arguments with defaults
DATASET=${1:-mini}
CUTOFF=${2:-0.1}
MIN_CPGS=${3:-5}
MAX_PERMS=${4:-10}
CORES=${5:-$SLURM_CPUS_PER_TASK}

echo "Parameters:"
echo "  Dataset: $DATASET"
echo "  Cutoff: $CUTOFF"
echo "  Min CpGs: $MIN_CPGS"
echo "  Max Permutations: $MAX_PERMS"
echo "  Cores: $CORES"
echo ""

# ============================================================================
# Conda Environment Setup
# ============================================================================

# Load conda if needed
if ! command -v conda &> /dev/null; then
    if [[ -n "${CONDA_SH:-}" && -f "$CONDA_SH" ]]; then
        source "$CONDA_SH"
    else
        echo "ERROR: conda not found. Load it (e.g. 'module load conda') or set" >&2
        echo "       CONDA_SH=/path/to/etc/profile.d/conda.sh" >&2
        exit 1
    fi
fi

# Activate the dmrseq-DSS environment
echo "Activating conda environment: dmrseq-DSS"
conda activate dmrseq-DSS

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to activate dmrseq-DSS environment"
    echo "Please ensure the environment exists:"
    echo "  conda env create -f dmrseq_DSS_env.yml"
    exit 1
fi

echo "Environment activated successfully"
echo "R version: $(Rscript --version 2>&1)"
echo ""

# ============================================================================
# Create Log Directory
# ============================================================================

LOG_DIR="../02_dmr_calling/logs/slurm"
mkdir -p "$LOG_DIR"

# ============================================================================
# Run dmrseq Analysis
# ============================================================================

echo "=========================================="
echo "Starting dmrseq Analysis"
echo "=========================================="
echo ""

# Navigate to DMR calling directory
SCRIPT_DIR="../02_dmr_calling"
cd "$SCRIPT_DIR" || { echo "ERROR: Cannot navigate to $SCRIPT_DIR"; exit 1; }
echo "Working directory: $(pwd)"
echo ""

# Run the R script
Rscript run_dmrseq.R \
    --dataset "$DATASET" \
    --cutoff "$CUTOFF" \
    --minCpGs "$MIN_CPGS" \
    --maxPerms "$MAX_PERMS" \
    --cores "$CORES"

EXIT_CODE=$?

echo ""
echo "=========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "dmrseq Analysis completed successfully!"
else
    echo "dmrseq Analysis failed with exit code: $EXIT_CODE"
fi
echo "=========================================="

# Deactivate conda environment
conda deactivate

exit $EXIT_CODE
