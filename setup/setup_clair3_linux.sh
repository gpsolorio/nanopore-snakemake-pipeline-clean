#!/bin/bash
#SBATCH --job-name=setup_clair3_singularity
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --output=./logs/%x.out
#SBATCH --error=./logs/%x.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

# ============================================================================
# Setup Clair3 Singularity/Apptainer Container
# ============================================================================
# Pulls the Clair3 Docker image (CPU or GPU) and converts it to a .sif,
# at the path common.smk expects (workflow/tools/singularity/clair3_*.sif).
#
# Only needs to be run once per machine (or when bumping the Clair3 version),
# and once per variant (CPU / GPU) you actually use.
#
# If you don't have Singularity/Apptainer available, or you're on Apple
# Silicon, use setup_clair3_source.sh instead.
# ============================================================================

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_PATH="$(scontrol show job "${SLURM_JOB_ID}" | grep -oP 'Command=\K[^\s]+')"
    SCRIPT_DIR="$(dirname "${SCRIPT_PATH}")"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "${LOG_DIR}"

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    exec 1>"${LOG_DIR}/${SLURM_JOB_NAME}.out" 2>"${LOG_DIR}/${SLURM_JOB_NAME}.err"
    if [[ -f "${SLURM_SUBMIT_DIR}/${SLURM_JOB_NAME}.out" ]]; then
        cat "${SLURM_SUBMIT_DIR}/${SLURM_JOB_NAME}.out" 2>/dev/null || true
    fi
fi

# script lives at setup/setup_clair3_singularity.sh -- one level below repo root
REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel 2>/dev/null || echo "$SCRIPT_DIR/..")"
REPO_ROOT="$(cd "$REPO_ROOT" && pwd)"

CLAIR3_BASE_VERSION="v2.0.2"
CONTAINER_DIR="${REPO_ROOT}/workflow/tools/singularity"
MODE="auto"   # auto | cpu | gpu | both

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

OPTIONS:
    --cpu                   Pull the CPU image only (${CLAIR3_BASE_VERSION})
    --gpu                   Pull the GPU image only (${CLAIR3_BASE_VERSION}_gpu)
    --both                  Pull both CPU and GPU images
    -o, --output-dir DIR    Directory to store the .sif file(s)
                             (default: workflow/tools/singularity/)
    -h, --help              Display this help message

    Default (no flag): auto-detect via nvidia-smi. GPU present -> pulls
    both CPU and GPU (CPU is still needed for e.g. login-node testing /
    CPU-only rules); no GPU -> pulls CPU only.

NOTES:
    - Pulls hkubal/clair3:${CLAIR3_BASE_VERSION}[_gpu] from Docker Hub
    - Creates .sif file(s) (multi-GB) at the paths workflow/snakemake/common.smk
      expects (_CLAIR3_IMAGE_CPU_SIF / _CLAIR3_IMAGE_GPU_SIF)
    - Only needs to be run once per image (or when updating Clair3)
EOF
    exit "$1"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cpu)            MODE="cpu"; shift ;;
        --gpu)             MODE="gpu"; shift ;;
        --both)            MODE="both"; shift ;;
        -o|--output-dir)   CONTAINER_DIR="$2"; shift 2 ;;
        -h|--help)         usage 0 ;;
        *)                 echo "Error: Unknown option '$1'" >&2; usage 1 ;;
    esac
done

# Prefer singularity; fall back to apptainer if that's what's installed.
# (Matches common.smk's SIF_ENGINES = ("singularity", "apptainer") and its
# PATH-based fallback in _resolve_sif_engine_bin().)
if command -v singularity >/dev/null 2>&1; then
    ENGINE="singularity"
elif command -v apptainer >/dev/null 2>&1; then
    ENGINE="apptainer"
else
    echo "[ERROR] Neither 'singularity' nor 'apptainer' found in PATH." >&2
    echo "        Activate the conda env that provides it (e.g. 'conda activate smk-env')," >&2
    echo "        or run setup_apptainer.sh for a repo-local install." >&2
    exit 1
fi

if [[ "$MODE" == "auto" ]]; then
    if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
        echo "[setup_clair3_singularity] NVIDIA GPU detected -- pulling both CPU and GPU images"
        MODE="both"
    else
        echo "[setup_clair3_singularity] No NVIDIA GPU detected -- pulling CPU image only"
        MODE="cpu"
    fi
fi

mkdir -p "$CONTAINER_DIR"
cd "$CONTAINER_DIR"

pull_and_test() {
    local variant="$1"        # "cpu" or "gpu"
    local image_tag="$2"      # e.g. hkubal/clair3:v2.0.2 or v2.0.2_gpu
    local sif_name="$3"       # e.g. clair3_v2.0.2.sif

    echo ""
    echo "=========================================="
    echo " Setting up Clair3 (${variant}) Singularity Container"
    echo "=========================================="
    echo "Engine:           $ENGINE"
    echo "Image:            $image_tag"
    echo "Output directory: $CONTAINER_DIR"
    echo ""

    local start_time
    start_time=$(date +%s)

    if [[ -f "${CONTAINER_DIR}/${sif_name}" ]]; then
        echo "[setup_clair3_singularity] ${sif_name} already exists -- skipping pull"
        echo "                            (delete it first if you want to re-pull)"
        return 0
    fi

    echo "Pulling Clair3 ${variant} image (this may take several minutes)..."
    "$ENGINE" pull "$sif_name" "docker://${image_tag}"

    if [[ ! -f "${CONTAINER_DIR}/${sif_name}" ]]; then
        echo "[ERROR] Failed to create ${sif_name}" >&2
        exit 1
    fi

    local file_size
    file_size=$(du -h "${CONTAINER_DIR}/${sif_name}" | cut -f1)

    echo ""
    echo "Testing Clair3 container..."
    "$ENGINE" exec "${CONTAINER_DIR}/${sif_name}" /opt/bin/run_clair3.sh --help > /dev/null 2>&1
    echo "[SUCCESS] Container test successful"

    if [[ "$variant" == "gpu" ]]; then
        echo ""
        echo "Verifying GPU passthrough (requires --nv and a GPU on this machine)..."
        if "$ENGINE" exec --nv "${CONTAINER_DIR}/${sif_name}" nvidia-smi > /dev/null 2>&1; then
            echo "[SUCCESS] GPU visible inside container"
        else
            echo "[WARN] GPU not visible inside container. This is expected if you're"
            echo "       building this on a CPU-only login node -- verify manually on"
            echo "       a GPU node with:"
            echo "         $ENGINE exec --nv ${CONTAINER_DIR}/${sif_name} nvidia-smi"
        fi
        echo ""
        echo "Reminder: GPU runs of run_clair3.sh also need --cleanenv --env TMPDIR=/tmp"
        echo "(otherwise 'parallel' can fail if the host TMPDIR isn't visible in-container)."
    fi

    local end_time elapsed elapsed_min elapsed_sec
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))
    elapsed_min=$((elapsed / 60))
    elapsed_sec=$((elapsed % 60))

    echo ""
    echo "Container location: ${CONTAINER_DIR}/${sif_name}"
    echo "Container size:     ${file_size}"
    echo "Elapsed time:        ${elapsed_min}m ${elapsed_sec}s"
}

case "$MODE" in
    cpu)
        pull_and_test "cpu" "hkubal/clair3:${CLAIR3_BASE_VERSION}" "clair3_${CLAIR3_BASE_VERSION}.sif"
        ;;
    gpu)
        pull_and_test "gpu" "hkubal/clair3:${CLAIR3_BASE_VERSION}_gpu" "clair3_${CLAIR3_BASE_VERSION}_gpu.sif"
        ;;
    both)
        pull_and_test "cpu" "hkubal/clair3:${CLAIR3_BASE_VERSION}" "clair3_${CLAIR3_BASE_VERSION}.sif"
        pull_and_test "gpu" "hkubal/clair3:${CLAIR3_BASE_VERSION}_gpu" "clair3_${CLAIR3_BASE_VERSION}_gpu.sif"
        ;;
esac

echo ""
echo "=========================================="
echo " SETUP COMPLETE"
echo "=========================================="
echo "Set cluster.container_engine: \"${ENGINE}\" in your config, and"
echo "resources.call_variants.<profile>.use_gpu accordingly."
echo ""
echo "Next: see docs/clair3_gpu_setup.md to verify PyTorch CUDA access and"
echo "set up squashfuse/GPU concurrency caps before running the pipeline (GPU only)."
echo ""

exit 0