#!/usr/bin/env bash

# ============================================================================
# NOTE: This is the FALLBACK Clair3 install path -- source build with conda
# + PyTorch, no container runtime required. Use this only if:
#   - you're on Apple Silicon (no GPU container support), or
#   - Singularity/Docker isn't available on your machine
#
# On Linux with Singularity available, use setup_clair3_singularity.sh
# instead -- it's simpler and matches what the pipeline's GPU path expects
# (see docs/clair3_gpu_setup.md).
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ============================================================================
# Clair3 v2 Source Installation Helper
# ============================================================================
# This script automates a source install of Clair3 v2 with mamba/conda,
# PyTorch, compiled C components, mandatory PyPy, and optional model download.
#
# Supports:
# - Linux CPU installs
# - Linux NVIDIA installs (CUDA wheel index configurable)
# - Apple Silicon installs (uses default PyTorch wheels)
# ============================================================================

usage() {
    cat << 'EOF'
Usage:
  setup_clair3_v2.sh [options]

Options:
  --env-name NAME           Conda env name (default: clair3_v2)
    --install-dir DIR         Clair3 clone directory (default: workflow/tools/clair3)
    --torch-build BUILD       PyTorch build: auto|cpu|cu126|custom (default: auto)
  --torch-index-url URL     Custom PyTorch index URL (required if --torch-build custom)
  --download-models         Download pretrained Clair3 models
  --model-name NAME         Default model name for run command example
                            (default: r1041_e82_400bps_sup_v500)
  -h, --help                Show help

Examples:
  ./setup_clair3_v2.sh --torch-build cpu
  ./setup_clair3_v2.sh --torch-build cu126 --download-models
  ./setup_clair3_v2.sh --torch-build custom \
      --torch-index-url https://download.pytorch.org/whl/cu126

Notes:
  - Requires conda or mamba (mamba preferred for speed).
  - On Apple Silicon, --torch-build auto uses default PyTorch wheels.
    - PyPy installation is mandatory and downloaded into the Clair3 directory.
EOF
}

log() {
    printf '[INFO] %s\n' "$*"
}

warn() {
    printf '[WARN] %s\n' "$*" >&2
}

die() {
    printf '[ERROR] %s\n' "$*" >&2
    exit 1
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

ENV_NAME="clair3_v2"
INSTALL_DIR="${REPO_ROOT}/workflow/tools/clair3"
TORCH_BUILD="auto"
TORCH_INDEX_URL=""
DOWNLOAD_MODELS="true"
MODEL_NAME="r1041_e82_400bps_sup_v500"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --env-name)
            ENV_NAME="$2"
            shift 2
            ;;
        --install-dir)
            INSTALL_DIR="$2"
            shift 2
            ;;
        --torch-build)
            TORCH_BUILD="$2"
            shift 2
            ;;
        --torch-index-url)
            TORCH_INDEX_URL="$2"
            shift 2
            ;;
        --download-models)
            DOWNLOAD_MODELS="true"
            shift
            ;;
        --model-name)
            MODEL_NAME="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown argument: $1"
            ;;
    esac
done

case "$TORCH_BUILD" in
    auto|cpu|cu126|custom)
        ;;
    *)
        die "Invalid --torch-build value '$TORCH_BUILD'. Use: auto|cpu|cu126|custom"
        ;;
esac

if [[ "$TORCH_BUILD" == "custom" && -z "$TORCH_INDEX_URL" ]]; then
    die "--torch-index-url is required when --torch-build custom"
fi

log "Starting Clair3 v2 setup"
log "Environment: $ENV_NAME"
log "Install dir: $INSTALL_DIR"
log "Torch build: $TORCH_BUILD"

# Step 0: dependencies and conda shell integration
require_cmd git
require_cmd make
require_cmd wget
require_cmd tar

if command -v mamba >/dev/null 2>&1; then
    CONDA_EXE="mamba"
elif command -v conda >/dev/null 2>&1; then
    CONDA_EXE="conda"
else
    die "Neither mamba nor conda found in PATH. Install Miniforge first."
fi

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

# Step 1: create and activate environment
log "Step 1/7: Creating conda environment (if needed)"
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    log "Conda env '$ENV_NAME' already exists; reusing it"
else
    "$CONDA_EXE" create -n "$ENV_NAME" -c conda-forge -c bioconda -y \
        python=3.11 samtools whatshap parallel \
        zstd xz zlib bzip2 automake make gcc gxx curl pigz
fi

# Some conda activate.d hooks reference optional unset vars on macOS.
# Temporarily relax nounset for activation, then restore strict mode.
set +u
conda activate "$ENV_NAME"
set -u

# Force use of the activated environment's Python toolchain.
export PATH="${CONDA_PREFIX}/bin:${PATH}"
hash -r
PYTHON_BIN="${CONDA_PREFIX}/bin/python"

"$PYTHON_BIN" -m pip install --upgrade pip
"$PYTHON_BIN" -m pip install --upgrade uv

# Clair3's Darwin Makefile expects ${BUILD_PREFIX}/share/gnuconfig/config.*
# and skips configure if missing, which can leave stale/broken build settings.
if [[ "$(uname -s)" == "Darwin" ]]; then
    "$CONDA_EXE" install -n "$ENV_NAME" -c conda-forge -y gnuconfig
fi

# Step 2: install PyTorch
log "Step 2/7: Installing PyTorch"
if [[ -n "$TORCH_INDEX_URL" ]]; then
    "$PYTHON_BIN" -m uv pip install torch torchvision --index-url "$TORCH_INDEX_URL"
else
    ARCH="$(uname -m)"
    OS="$(uname -s)"

    if [[ "$TORCH_BUILD" == "auto" ]]; then
        if [[ "$OS" == "Darwin" && "$ARCH" == "arm64" ]]; then
            log "Detected Apple Silicon; using default PyTorch wheels"
            "$PYTHON_BIN" -m uv pip install torch torchvision
        elif command -v nvidia-smi >/dev/null 2>&1; then
            warn "Detected NVIDIA GPU; defaulting to CUDA 12.6 wheels. Override with --torch-build or --torch-index-url if needed."
            "$PYTHON_BIN" -m uv pip install torch torchvision --index-url https://download.pytorch.org/whl/cu126
        else
            log "No NVIDIA GPU detected; using CPU wheels"
            "$PYTHON_BIN" -m uv pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu
        fi
    elif [[ "$TORCH_BUILD" == "cpu" ]]; then
        "$PYTHON_BIN" -m uv pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu
    elif [[ "$TORCH_BUILD" == "cu126" ]]; then
        "$PYTHON_BIN" -m uv pip install torch torchvision --index-url https://download.pytorch.org/whl/cu126
    else
        die "Unexpected torch build state: $TORCH_BUILD"
    fi
fi

# Step 3: clone Clair3
log "Step 3/7: Cloning Clair3"
if [[ -d "$INSTALL_DIR/.git" ]]; then
    log "Existing git repo found at $INSTALL_DIR; pulling latest changes"
    git -C "$INSTALL_DIR" pull --ff-only
else
    if [[ -d "$INSTALL_DIR" && -f "$INSTALL_DIR/run_clair3.sh" ]]; then
        log "Existing Clair3 checkout found at $INSTALL_DIR; reusing it"
    elif [[ -e "$INSTALL_DIR" ]]; then
        die "Install path exists but is not a Clair3 checkout or git repo: $INSTALL_DIR"
    else
        mkdir -p "$(dirname "$INSTALL_DIR")"
        git clone https://github.com/HKU-BAL/Clair3.git "$INSTALL_DIR"
    fi
fi

cd "$INSTALL_DIR"
CLAIR3_PATH="$PWD"
export CLAIR3_PATH

# Step 4: python deps + build
log "Step 4/7: Installing Clair3 Python dependencies and building C sources"
"$PYTHON_BIN" -m uv pip install numpy h5py hdf5plugin numexpr tqdm cffi torchmetrics

# On macOS, force Apple toolchain to avoid conda gcc linker/runtime mismatches.
if [[ "$(uname -s)" == "Darwin" ]]; then
    PATH="${CONDA_PREFIX}/bin:${PATH}" make clean >/dev/null 2>&1 || true
    PATH="${CONDA_PREFIX}/bin:${PATH}" \
        make PREFIX="${CONDA_PREFIX}" PYTHON="$PYTHON_BIN" \
        BUILD_PREFIX="${CONDA_PREFIX}" CC_PATH=/usr/bin/clang \
        CC=/usr/bin/clang CXX=/usr/bin/clang++ AR=/usr/bin/ar RANLIB=/usr/bin/ranlib
else
    PATH="${CONDA_PREFIX}/bin:${PATH}" make PREFIX="${CONDA_PREFIX}" PYTHON="$PYTHON_BIN"
fi

# Step 5: PyPy (mandatory)
log "Step 5/7: Installing PyPy3.11 (mandatory)"
OS_NAME="$(uname -s)"
ARCH_NAME="$(uname -m)"
PYPY_VERSION="pypy3.11-v7.3.20"

case "${OS_NAME}-${ARCH_NAME}" in
    Linux-x86_64)
        PYPY_PLATFORM="linux64"
        ;;
    Darwin-arm64)
        PYPY_PLATFORM="macos_arm64"
        ;;
    Darwin-x86_64)
        PYPY_PLATFORM="macos_x86_64"
        ;;
    *)
        die "Unsupported OS/arch for PyPy download: ${OS_NAME}-${ARCH_NAME}"
        ;;
esac

PYPY_DIR="${CLAIR3_PATH}/${PYPY_VERSION}-${PYPY_PLATFORM}"
PYPY_TARBALL="${PYPY_VERSION}-${PYPY_PLATFORM}.tar.bz2"
PYPY_URL="https://downloads.python.org/pypy/${PYPY_TARBALL}"

if [[ ! -d "$PYPY_DIR" ]]; then
    (cd "$CLAIR3_PATH" && wget "$PYPY_URL" && tar -xjf "$PYPY_TARBALL" && rm -f "$PYPY_TARBALL")
fi

ln -sf "${PYPY_DIR}/bin/pypy3" "${CONDA_PREFIX}/bin/pypy3"
ln -sf "${PYPY_DIR}/bin/pypy3" "${CONDA_PREFIX}/bin/pypy"

pypy3 -m ensurepip
pypy3 -m pip install mpmath==1.2.1

# Step 6: optional model download
log "Step 6/7: Model download"
if [[ "$DOWNLOAD_MODELS" == "true" ]]; then
    mkdir -p "${CLAIR3_PATH}/models"
    wget -r -np -nH --cut-dirs=2 -R "index.html*" -P "${CLAIR3_PATH}/models" \
        https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/
    log "Model download complete"
else
    log "Skipping model download (use --download-models to enable)"
fi

# Step 7: usage summary
log "Step 7/7: Setup summary"
cat << EOF

====================================================================
Clair3 v2 setup complete.

Environment:
  conda activate ${ENV_NAME}

Clair3 path:
  ${CLAIR3_PATH}

Example run:
  THREADS=8
  OUTPUT_DIR=./clair3_output
  ${CLAIR3_PATH}/run_clair3.sh \\
    --bam_fn=input.bam \\
    --ref_fn=ref.fa \\
    --threads=\${THREADS} \\
    --platform=ont \\
    --model_path=${CLAIR3_PATH}/models/${MODEL_NAME} \\
    --output=\${OUTPUT_DIR}

Equivalent Python entrypoint:
  python3 ${CLAIR3_PATH}/run_clair3.py [same arguments]
====================================================================
EOF
