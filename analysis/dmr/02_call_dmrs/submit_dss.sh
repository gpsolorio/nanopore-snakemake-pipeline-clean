#!/bin/bash
##############################################################################
# SLURM Submission Script for DSS DMR Analysis
##############################################################################
#
# Usage: sbatch submit_dss.sh [DATASET] [DELTA] [P_THRESHOLD] [CORES]
#
# Examples:
#   sbatch submit_dss.sh mini 0.10 0.001 8
#   sbatch submit_dss.sh full 0.15 0.001 16
#   sbatch submit_dss.sh              # Uses defaults
#
##############################################################################

#SBATCH --job-name=dss_dmr
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --output=logs/slurm/dss_dmr_%j.out
#SBATCH --error=logs/slurm/dss_dmr_%j.err
#SBATCH --time=0-06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# ============================================================================
# Setup
# ============================================================================

echo "=========================================="
echo "DSS DMR Analysis - SLURM Submission"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Date: $(date)"
echo "=========================================="
echo ""

# Parse command-line arguments with defaults
DATASET=${1:-mini}
DELTA=${2:-0.10}
P_THRESHOLD=${3:-0.001}
CORES=${4:-$SLURM_CPUS_PER_TASK}

echo "Parameters:"
echo "  Dataset: $DATASET"
echo "  Delta: $DELTA"
echo "  P-value threshold: $P_THRESHOLD"
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
# Run DSS Analysis
# ============================================================================

echo "=========================================="
echo "Starting DSS Analysis"
echo "=========================================="
echo ""

# Navigate to DMR calling directory
SCRIPT_DIR="../02_dmr_calling"
cd "$SCRIPT_DIR" || { echo "ERROR: Cannot navigate to $SCRIPT_DIR"; exit 1; }
echo "Working directory: $(pwd)"
echo ""

# Run the R script
Rscript run_dss.R \
    --dataset "$DATASET" \
    --delta "$DELTA" \
    --p.threshold "$P_THRESHOLD" \
    --cores "$CORES"

EXIT_CODE=$?

echo ""
echo "=========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "DSS Analysis completed successfully!"
else
    echo "DSS Analysis failed with exit code: $EXIT_CODE"
fi
echo "=========================================="

# Deactivate conda environment
conda deactivate

exit $EXIT_CODE
