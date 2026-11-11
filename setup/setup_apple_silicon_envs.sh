#!/usr/bin/env bash
# setup/setup_apple_silicon_envs.sh
#
# Apple Silicon (M1/M2/M3/M4) has no native osx-arm64 build for
# bioconductor-bsseq / bioconductor-dss. Three of this pipeline's conda
# envs depend on them:
#
#   - dss-env        (dss_call_dmrs, dss_smooth_betas)
#   - wgcna-env       (wgcna_modules)
#   - pheatmap-env    (plotting rules)
#
# Snakemake's normal `--use-conda` auto-creation will always fail for these
# three on Apple Silicon, because it solves natively for osx-arm64. The fix:
# build them once, manually, under Rosetta (osx-64 subdir) at the EXACT path
# Snakemake itself computes from each env yaml's content. Once built there,
# Snakemake sees them as already-created and never attempts to solve them
# natively again.
#
# This version routes every Snakemake conda-create-envs-only call through a
# log file instead of letting it stream to the terminal, then classifies any
# CreateCondaEnvironmentException it finds:
#   - failures for envs in ROSETTA_ENVS are EXPECTED -> printed as a short,
#     calm note (no raw conda traceback shown)
#   - failures for any other env are UNEXPECTED -> full log context is
#     printed and the script exits non-zero immediately
#
# Run this ONCE per fresh clone / per .snakemake/ wipe, before your first
# real pipeline run. Safe to re-run — it's idempotent (skips envs that
# already exist).
#
# Usage:
#   ./setup/setup_apple_silicon_envs.sh [-s SNAKEFILE] [-c CONFIGFILE]
#
# Requires: conda (or mamba as CONDA_EXE frontend), snakemake, this repo's
# workflow/envs/*.yml files.

set -euo pipefail

# ===========================================================================
# Config (override via flags or env vars)
# ===========================================================================
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SNAKEFILE="${SNAKEFILE:-workflow/snakemake/nanopore_adaptive_pipeline.smk}"
CONFIGFILE="${CONFIGFILE:-configs/mini_data_public_2023_mac_local.yml}"

# Envs known to need the Rosetta workaround. If you add a new R/Bioconductor
# env that pulls in bsseq/dss, add its `name:` field (== its yml basename)
# here.
ROSETTA_ENVS=(dss-env wgcna-env pheatmap-env)

LOG_DIR="$(mktemp -d /tmp/apple_silicon_env_setup.XXXXXX)"

while getopts "s:c:" opt; do
  case "$opt" in
    s) SNAKEFILE="$OPTARG" ;;
    c) CONFIGFILE="$OPTARG" ;;
    *) echo "Usage: $0 [-s SNAKEFILE] [-c CONFIGFILE]" >&2; exit 1 ;;
  esac
done

cd "$REPO_ROOT"

# ===========================================================================
# Guardrails
# ===========================================================================
ARCH="$(uname -m)"
OS="$(uname -s)"

if [[ "$OS" != "Darwin" || "$ARCH" != "arm64" ]]; then
    echo "This script is only needed on Apple Silicon Macs (Darwin/arm64)."
    echo "Detected: ${OS}/${ARCH} -- nothing to do here. Exiting."
    exit 0
fi

if ! command -v conda &>/dev/null; then
    echo "ERROR: conda not found on PATH. Install Miniforge/Miniconda first." >&2
    exit 1
fi

if ! command -v snakemake &>/dev/null; then
    echo "ERROR: snakemake not found on PATH. Activate your pipeline conda env first." >&2
    exit 1
fi

safe_rm() {
    local target="$1"
    # Refuse to touch anything that isn't a real hash-named env dir under
    # .snakemake/conda/ -- this is the guardrail that would have prevented
    # an earlier real incident where an empty variable turned into
    # `rm -rf .snakemake/conda/` and wiped every built env.
    if [[ -z "$target" ]]; then
        echo "REFUSING to delete: empty path" >&2
        exit 1
    fi
    case "$target" in
        .snakemake/conda|.snakemake/conda/|.snakemake|.)
            echo "REFUSING to delete suspicious path: '$target'" >&2
            exit 1
            ;;
    esac
    rm -rf -- "$target"
}

# ===========================================================================
# Helper: run a conda-create-envs-only pass, fully logged (never streamed
# raw to the terminal), and echo back only the lines worth seeing.
# ===========================================================================
run_conda_create_pass() {
    local logfile="$1"

    # IMPORTANT: stream through `tee` + a line-buffered grep filter rather
    # than fully buffering to $logfile and only echoing after the whole
    # pass exits. A single conda-create-envs-only pass can take minutes
    # (downloads/solves for every env in the DAG); if nothing prints until
    # it's entirely done, the script looks hung even when it's working
    # fine. This way you still see each "Creating..."/"Environment ...
    # created" line live, in real time, while the full raw output is also
    # saved to $logfile for later inspection.
    set +o pipefail
    snakemake -s "$SNAKEFILE" --configfile "$CONFIGFILE" \
        --use-conda --conda-create-envs-only --forceall \
        --cores 1 --keep-going --conda-frontend conda \
        2>&1 | tee "$logfile" \
        | grep --line-buffered -E "^(Creating conda environment|Environment for .* created)" || true
    local result="${PIPESTATUS[0]}"
    set -o pipefail

    return "$result"
}

# ===========================================================================
# Helper: inspect a pass's log for CreateCondaEnvironmentException failures.
# Prints a calm note for expected (Rosetta-list) envs. For any other env,
# prints full context and exits non-zero -- that's a real problem, not the
# known arm64/bsseq gap.
# ===========================================================================
check_for_unexpected_failures() {
    local logfile="$1"

    local failing_ymls
    failing_ymls=$(grep "Could not create conda environment from" "$logfile" \
        | grep -oE "workflow/envs/[A-Za-z0-9_.-]+\.ya?ml" \
        | xargs -n1 basename 2>/dev/null \
        | sed -E 's/\.ya?ml$//' \
        | sort -u || true)

    [[ -z "$failing_ymls" ]] && return 0

    while IFS= read -r name; do
        [[ -z "$name" ]] && continue
        local is_known=0
        for known in "${ROSETTA_ENVS[@]}"; do
            [[ "$name" == "$known" ]] && is_known=1 && break
        done

        if [[ "$is_known" -eq 1 ]]; then
            echo "  (expected) ${name}: no osx-arm64 build available -- will build under Rosetta (osx-64) below."
        else
            echo ""
            echo "ERROR: unexpected conda env failure for '${name}' (not in ROSETTA_ENVS)." >&2
            echo "This is NOT the known bsseq/dss arm64 gap -- treating it as a real problem." >&2
            echo "Relevant log context (${logfile}):" >&2
            echo "---" >&2
            grep -A 25 "Could not create conda environment from.*${name}\.ya\?ml" "$logfile" >&2 || cat "$logfile" >&2
            echo "---" >&2
            exit 1
        fi
    done <<< "$failing_ymls"
}

# ===========================================================================
# Step 1: Iteratively generate env specs + build any Rosetta-only envs.
#
# Each round: run a create-envs-only pass (logged, not streamed). Classify
# any failures. Build any newly-appeared Rosetta env yamls under osx-64.
# Repeat until a pass produces zero new Rosetta builds and exits cleanly.
# ===========================================================================
echo "=== Building conda environments (Rosetta workaround applied silently where expected) ==="
mkdir -p .snakemake/conda

MAX_ROUNDS=10
for round in $(seq 1 $MAX_ROUNDS); do
    echo ""
    echo "--- Round ${round} ---"
    logfile="${LOG_DIR}/round_${round}.log"

    set +e
    run_conda_create_pass "$logfile"
    pass_result=$?
    set -e

    check_for_unexpected_failures "$logfile"

    built_something=0
    for name in "${ROSETTA_ENVS[@]}"; do
        yaml=$(grep -l "^name: ${name}\$" .snakemake/conda/*.yaml 2>/dev/null || true)
        if [[ -z "$yaml" ]]; then
            continue
        fi
        hash=$(basename "$yaml" .yaml)
        prefix=".snakemake/conda/${hash}"

        if [[ -d "$prefix" ]] && conda run --prefix "$prefix" true &>/dev/null; then
            continue  # already built and functional
        fi

        echo "  Building ${name} under Rosetta (osx-64) at ${prefix} ..."
        safe_rm "$prefix"
        CONDA_SUBDIR=osx-64 conda env create -f "$yaml" --prefix "$prefix" --quiet
        conda config --set subdir osx-64 --prefix "$prefix"
        echo "  -> ${name} built successfully."
        built_something=1
    done

    if [[ $built_something -eq 0 && $pass_result -eq 0 ]]; then
        echo ""
        echo "=== All envs built successfully after ${round} round(s). ==="
        break
    fi

    if [[ $round -eq $MAX_ROUNDS ]]; then
        echo "ERROR: reached MAX_ROUNDS (${MAX_ROUNDS}) without converging." >&2
        echo "See logs under ${LOG_DIR} for details." >&2
        exit 1
    fi
done

# ===========================================================================
# Step 2: Final verification -- dry run must show zero env creation
# ===========================================================================
echo ""
echo "=== Final verification: dry run should show no env creation ==="
verify_log="${LOG_DIR}/verify.log"
snakemake -s "$SNAKEFILE" --configfile "$CONFIGFILE" \
    --use-conda -n --conda-create-envs-only \
    > "$verify_log" 2>&1

if grep -q "Creating conda environment" "$verify_log"; then
    echo "WARNING: some envs still appear to need creation. See:" >&2
    grep "Creating conda environment" "$verify_log" >&2
    exit 1
fi

echo "All conda environments are built and recognized by Snakemake."

# ===========================================================================
# Step 3: plink2 binary (NOT a conda env)
#
# plink2 has no osx-arm64 bioconda package. The pipeline resolves it through
# tools.plink2_path in config (see common.smk::_resolve_plink2()).
#
# Install a native Apple Silicon binary here rather than forcing Rosetta or
# contaminating the conda environments.
#
# The URL is intentionally pinned to a dated upstream release for
# reproducibility. Update PLINK2_VERSION and PLINK2_DATE when upgrading.
# ===========================================================================

PLINK2_DEST_DIR="${REPO_ROOT}/workflow/tools/plink2"
PLINK2_BIN="${PLINK2_DEST_DIR}/plink2"

PLINK2_VERSION="alpha7"
PLINK2_DATE="20260504"
PLINK2_URL="https://s3.amazonaws.com/plink2-assets/${PLINK2_VERSION}/plink2_mac_arm64_${PLINK2_DATE}.zip"

echo ""
echo "=== Checking plink2 (external binary, not conda-managed) ==="

if [[ -x "$PLINK2_BIN" ]]; then
    echo "plink2 already present:"
    "$PLINK2_BIN" --version || true
else
    echo "Installing native Apple Silicon plink2..."
    echo "Source:"
    echo "  ${PLINK2_URL}"

    mkdir -p "$PLINK2_DEST_DIR"

    TMP_DIR="$(mktemp -d)"
    trap 'rm -rf "$TMP_DIR"' EXIT

    curl -fsSL "$PLINK2_URL" \
        -o "${TMP_DIR}/plink2.zip"

    unzip -q "${TMP_DIR}/plink2.zip" -d "$TMP_DIR"

    if [[ ! -f "${TMP_DIR}/plink2" ]]; then
        echo "ERROR: downloaded archive does not contain plink2 binary." >&2
        exit 1
    fi

    mv "${TMP_DIR}/plink2" "$PLINK2_BIN"
    chmod +x "$PLINK2_BIN"

    echo ""
    echo "Installed:"
    echo "  ${PLINK2_BIN}"
fi

echo ""
echo "Checking executable architecture:"
file "$PLINK2_BIN"

echo ""
echo "plink2 version:"
"$PLINK2_BIN" --version || true

echo ""
echo "Add this to your config YAML:"
echo ""
echo "  tools:"
echo "    plink2_path: ${PLINK2_BIN}"

echo ""
echo "(Full per-round logs kept at: ${LOG_DIR} -- safe to delete once you're satisfied.)"