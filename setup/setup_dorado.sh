#!/bin/bash
#SBATCH --job-name=setup_dorado
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --output=./logs/%x.out
#SBATCH --error=./logs/%x.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

# ============================================================================
# Dorado Installation Script
# ============================================================================
# This script downloads and sets up Dorado basecaller in workflow/tools.
# Works in both HPC (SLURM) and local environments (macOS/Linux).
# ============================================================================

set -euo pipefail

# Get the actual script path (not the SLURM copy in /var/spool)
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_PATH="$(scontrol show job ${SLURM_JOB_ID} | grep -oP 'Command=\K[^\s]+')"
    SCRIPT_DIR="$(dirname "${SCRIPT_PATH}")"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# Get repository root (script is in setup/)
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "${LOG_DIR}"

# Move SLURM logs to proper location
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    exec 1>"${LOG_DIR}/${SLURM_JOB_NAME}.out" 2>"${LOG_DIR}/${SLURM_JOB_NAME}.err"
fi

# Configuration
DORADO_INSTALL_DIR="${REPO_ROOT}/workflow/tools/dorado"
DORADO_VERSION="1.4.0"

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    RUN_CONTEXT="HPC/SLURM"
else
    RUN_CONTEXT="LOCAL"
fi

echo ""
echo "=============================================================================="
echo "Dorado Installation"
echo "=============================================================================="
echo ""
echo "Context: $RUN_CONTEXT"
echo "Installation directory: $DORADO_INSTALL_DIR"
echo "Version: $DORADO_VERSION"
echo ""

# Step 0: Detect platform
echo "[STEP 0] Detecting platform..."
OS="$(uname -s)"
ARCH="$(uname -m)"

echo "OS: $OS, Architecture: $ARCH"

if [[ "$OS" == "Linux" && "$ARCH" == "x86_64" ]]; then
    FILE="dorado-${DORADO_VERSION}-linux-x64.tar.gz"
    ARCHIVE_TYPE="tar.gz"
    
elif [[ "$OS" == "Darwin" && "$ARCH" == "arm64" ]]; then
    FILE="dorado-${DORADO_VERSION}-osx-arm64.zip"
    ARCHIVE_TYPE="zip"
    
elif [[ "$OS" == "Darwin" && "$ARCH" == "x86_64" ]]; then
    FILE="dorado-${DORADO_VERSION}-osx-x64.zip"
    ARCHIVE_TYPE="zip"
    
else
    echo "[ERROR] Unsupported platform: $OS $ARCH" >&2
    exit 1
fi

echo "✓ Detected platform: $FILE"
echo ""

# Step 1: Create install directory
echo "[STEP 1] Setting up installation directory..."
mkdir -p "$(dirname "$DORADO_INSTALL_DIR")"
echo "✓ Directory ready: $(dirname "$DORADO_INSTALL_DIR")"
echo ""

# Step 2: Download Dorado
echo "[STEP 2] Downloading Dorado..."
if [[ -d "$DORADO_INSTALL_DIR" ]] && [[ -f "$DORADO_INSTALL_DIR/bin/dorado" ]]; then
    echo "✓ Dorado already installed at $DORADO_INSTALL_DIR"
else
    URL="https://cdn.oxfordnanoportal.com/software/analysis/${FILE}"
    
    # Create temporary directory for download
    TEMP_DIR=$(mktemp -d)
    trap "rm -rf $TEMP_DIR" EXIT
    
    cd "$TEMP_DIR"
    echo "Downloading from: $URL"
    if ! curl -L "$URL" -o "$FILE"; then
        echo "[ERROR] Failed to download Dorado from $URL" >&2
        exit 1
    fi
    echo "✓ Download complete"

    FILE_SIZE=$(stat -f%z "$FILE" 2>/dev/null || stat -c%s "$FILE")

    if [[ "$FILE_SIZE" -lt 1000000 ]]; then
        echo "[ERROR] Downloaded file is suspiciously small (${FILE_SIZE} bytes)" >&2
        echo "[ERROR] URL may be invalid: $URL" >&2
        exit 1
    fi
    
    echo ""

    echo "[STEP 3] Extracting Dorado..."

    if [[ "$ARCHIVE_TYPE" == "tar.gz" ]]; then
        tar -xzf "$FILE"
    elif [[ "$ARCHIVE_TYPE" == "zip" ]]; then
        unzip -q "$FILE"
    else
        echo "[ERROR] Unknown archive type: $ARCHIVE_TYPE" >&2
        exit 1
    fi
    
    # Move extracted folder to final location
    EXTRACTED_DIR=$(ls -d dorado-* 2>/dev/null | head -1)
    if [[ -z "$EXTRACTED_DIR" ]]; then
        echo "[ERROR] Could not find extracted Dorado directory" >&2
        exit 1
    fi
    
    mkdir -p "$(dirname "$DORADO_INSTALL_DIR")"
    mv "$EXTRACTED_DIR" "$DORADO_INSTALL_DIR"
    cd - > /dev/null
    
    echo "✓ Dorado extracted to $DORADO_INSTALL_DIR"
fi
echo ""

# Step 4: Verify installation
echo "[STEP 4] Verifying Dorado installation..."
if [[ ! -f "$DORADO_INSTALL_DIR/bin/dorado" ]]; then
    echo "[ERROR] Dorado binary not found at $DORADO_INSTALL_DIR/bin/dorado" >&2
    exit 1
fi

"$DORADO_INSTALL_DIR/bin/dorado" --version
echo "✓ Dorado command is available"
echo ""

# Summary
echo "=============================================================================="
echo "✓ DORADO INSTALLATION COMPLETE"
echo "=============================================================================="
echo ""
echo "Installation Details:"
echo "  Installation path: $DORADO_INSTALL_DIR"
echo "  Binary path: $DORADO_INSTALL_DIR/bin/dorado"
echo ""
if [[ "$RUN_CONTEXT" = "HPC/SLURM" ]]; then
    echo "To use Dorado on HPC:"
    echo "  - Add to PATH: export PATH=\"$DORADO_INSTALL_DIR/bin:\$PATH\""
    echo "  - Or call directly: $DORADO_INSTALL_DIR/bin/dorado basecaller ..."
else
    echo "To use Dorado locally:"
    echo "  - Add to PATH: export PATH=\"$DORADO_INSTALL_DIR/bin:\$PATH\""
    echo "  - Or call directly: $DORADO_INSTALL_DIR/bin/dorado basecaller ..."
fi
echo ""
echo "Tool installed within repository for reproducibility:"
echo "  - workflow/tools/dorado/"
echo "=============================================================================="
