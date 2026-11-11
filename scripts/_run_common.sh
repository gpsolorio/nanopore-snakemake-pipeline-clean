#!/usr/bin/env bash
##############################################################################
# _run_common.sh — Shared setup for scripts/run_mac.sh / scripts/run_hpc.sh
#
# This file is SOURCED by the two entrypoint scripts, never run directly.
# It exists so the two environment-specific scripts only contain what's
# actually different between a Mac laptop and an HPC cluster; everything
# else (arg parsing, preflight checks, the reproducibility manifest) lives
# here once instead of being copy-pasted twice.
#
# Before sourcing this file, the caller must set:
#   RUN_LABEL              short name, e.g. "run_mac" (used in log/manifest names)
#   RUN_DEFAULT_CONFIGFILE default --configfile if the caller doesn't pass one
#   RUN_HELP_FN            name of a function to call on -h/--help
#
# Provides back to the caller:
#   REPO_ROOT               absolute path to the repo (cwd is chdir'd here)
#   SNAKEFILE                path to the main Snakefile
#   CONFIGFILE               resolved --configfile
#   SNAKEMAKE_ARGS            array of any other args, passed through untouched
##############################################################################

set -euo pipefail

# --- Make behavior independent of the caller's working directory ----------
# This file lives in scripts/, next to scripts/run_mac.sh / scripts/run_hpc.sh, so the repo
# root is its parent directory -- this works whether someone runs
# `bash scripts/run_mac.sh` or `bash /full/path/to/scripts/run_mac.sh`
# from anywhere on disk.
REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." &>/dev/null && pwd)"
cd "$REPO_ROOT"

SNAKEFILE="workflow/snakemake/nanopore_adaptive_pipeline.smk"

# --- Preflight: fail with a clear message instead of a snakemake traceback -
for _bin in snakemake conda; do
    if ! command -v "$_bin" >/dev/null 2>&1; then
        echo "ERROR: '$_bin' not found on PATH." >&2
        echo "Activate the environment that provides it (e.g. 'conda activate <base-env>') before running this script." >&2
        exit 1
    fi
done
unset _bin

# --- Shared CLI parsing ------------------------------------------------------
# --configfile overrides RUN_DEFAULT_CONFIGFILE; -h/--help calls RUN_HELP_FN;
# everything else is collected and passed straight through to snakemake.
: "${RUN_DEFAULT_CONFIGFILE:?RUN_DEFAULT_CONFIGFILE must be set before sourcing _run_common.sh}"
: "${RUN_LABEL:?RUN_LABEL must be set before sourcing _run_common.sh}"
: "${RUN_HELP_FN:?RUN_HELP_FN must be set before sourcing _run_common.sh}"

CONFIGFILE="$RUN_DEFAULT_CONFIGFILE"
SNAKEMAKE_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --configfile)
            CONFIGFILE="$2"
            shift 2
            ;;
        -h|--help)
            "$RUN_HELP_FN"
            exit 0
            ;;
        *)
            SNAKEMAKE_ARGS+=("$1")
            shift
            ;;
    esac
done

if [[ ! -f "$CONFIGFILE" ]]; then
    echo "ERROR: config file not found: $CONFIGFILE" >&2
    echo "Pass a valid path with --configfile, e.g.:" >&2
    echo "  bash scripts/${RUN_LABEL}.sh --configfile configs/your_config.yml" >&2
    exit 1
fi

# --- Reproducibility manifest ------------------------------------------------
# One small record per invocation: exactly what ran, against what code and
# config, from where. Independent of anything Snakemake itself logs, and
# cheap enough to write unconditionally (including on dry runs).
_MANIFEST_DIR="logs/run_manifest"
mkdir -p "$_MANIFEST_DIR"
_MANIFEST_FILE="${_MANIFEST_DIR}/$(date -u +%Y%m%dT%H%M%SZ)_${RUN_LABEL}.log"

{
    echo "entrypoint      : $RUN_LABEL"
    echo "timestamp_utc   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "host            : $(hostname)"
    echo "user            : $(whoami)"
    echo "repo_root       : $REPO_ROOT"
    echo "git_commit      : $(git rev-parse HEAD 2>/dev/null || echo 'not a git repo')"
    echo "git_dirty       : $(git diff --quiet HEAD 2>/dev/null && echo clean || echo 'uncommitted changes present')"
    echo "snakemake       : $(snakemake --version 2>/dev/null || echo 'unknown')"
    echo "conda           : $(conda --version 2>/dev/null || echo 'unknown')"
    echo "configfile      : $CONFIGFILE"
    echo "extra_args      : ${SNAKEMAKE_ARGS[*]:-<none>}"
} | tee "$_MANIFEST_FILE"
echo ">>> Run manifest: $_MANIFEST_FILE"
echo ""
