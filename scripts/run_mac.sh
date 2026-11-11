#!/usr/bin/env bash
##############################################################################
# scripts/run_mac.sh — Run the pipeline on a Mac laptop (local development / testing)
#
# On Apple Silicon, CONDA_SUBDIR=osx-64 is exported so that conda envs
# Snakemake creates via `conda: <path-to-yml>` (which Snakemake manages
# itself, hashed by file content -- NOT by matching a pre-existing named
# env) solve against Intel/Rosetta builds. This is needed because some
# bioconda/Bioconductor packages used by this pipeline have no native
# osx-arm64 build yet.
#
# Usage:
#   bash scripts/run_mac.sh -n                                        # dry-run, default config
#   bash scripts/run_mac.sh                                            # full run, default config
#   bash scripts/run_mac.sh --configfile configs/mini_data_public_2023_mac_local.yml  # different config
#   bash scripts/run_mac.sh --until dss_call_dmrs -n                   # dry-run up to a rule
#   bash scripts/run_mac.sh -h                                         # help
##############################################################################

set -euo pipefail

run_mac_help() {
    cat <<EOF
scripts/run_mac.sh — run the nanopore adaptive pipeline on a Mac (local)

Usage:
  bash scripts/run_mac.sh [-n] [--configfile PATH] [SNAKEMAKE_ARGS...]

Options:
  --configfile PATH   Path to a Snakemake config YAML.
                       Default: ${RUN_DEFAULT_CONFIGFILE:-<unset>}
  -h, --help           Show this help and exit.

Any other arguments (e.g. -n, --until RULE) are passed through to snakemake
untouched. Uses the 'profiles/local' Snakemake profile.
EOF
}

RUN_LABEL="run_mac"
RUN_DEFAULT_CONFIGFILE="configs/mini_data_public_2023_mac_local.yml"
RUN_HELP_FN="run_mac_help"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
# shellcheck source=./_run_common.sh
source "${SCRIPT_DIR}/_run_common.sh"

# --- Mac-specific: Apple Silicon needs osx-64 conda envs (Rosetta) ---------
if [[ "$(uname)" == "Darwin" && "$(uname -m)" == "arm64" ]]; then
    echo ">>> Apple Silicon detected: exporting CONDA_SUBDIR=osx-64 for conda env solves"
    export CONDA_SUBDIR=osx-64
fi

exec snakemake \
    -s "$SNAKEFILE" \
    --configfile "$CONFIGFILE" \
    --profile profiles/local \
    --resources mem_gb=18 \
    -j 2 \
    --use-conda \
    "${SNAKEMAKE_ARGS[@]+"${SNAKEMAKE_ARGS[@]}"}"
