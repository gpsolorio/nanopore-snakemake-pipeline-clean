#!/usr/bin/env bash
##############################################################################
# scripts/run_hpc.sh — Run the pipeline on an HPC cluster (SLURM)
#
# Unlike the Mac entrypoint, no CONDA_SUBDIR/Rosetta workaround is needed
# here: HPC nodes are linux-64, and all bioconda/Bioconductor packages used
# by this pipeline ship native linux-64 builds. Job submission is handled
# by the Snakemake 'profiles/slurm' profile (executor plugin config,
# per-rule resources, etc.) rather than by this script.
#
# Usage:
#   bash scripts/run_hpc.sh -n                                        # dry-run, default config
#   bash scripts/run_hpc.sh                                            # full run, default config
#   bash scripts/run_hpc.sh --configfile configs/full_data.yml         # different config
#   bash scripts/run_hpc.sh --until dss_call_dmrs -n                    # dry-run up to a rule
#   bash scripts/run_hpc.sh -h                                          # help
#
# Env overrides (set before invoking, e.g. HPC_MAX_JOBS=50 bash scripts/run_hpc.sh):
#   HPC_MAX_JOBS   max concurrent jobs Snakemake will submit (default: 50)
#   HPC_MEM_GB     mem_gb resource ceiling passed to --resources (default: 256)
##############################################################################

set -euo pipefail

run_hpc_help() {
    cat <<EOF
scripts/run_hpc.sh — run the nanopore adaptive pipeline on an HPC cluster (SLURM)

Usage:
  bash scripts/run_hpc.sh [-n] [--configfile PATH] [SNAKEMAKE_ARGS...]

Options:
  --configfile PATH   Path to a Snakemake config YAML.
                       Default: ${RUN_DEFAULT_CONFIGFILE:-<unset>}
  -h, --help           Show this help and exit.

Env overrides:
  HPC_MAX_JOBS   max concurrent SLURM jobs (default: 50)
  HPC_MEM_GB     mem_gb resource ceiling (default: 256)

Any other arguments (e.g. -n, --until RULE) are passed through to snakemake
untouched. Uses the 'profiles/slurm' Snakemake profile — make sure that
profile exists and is configured for your cluster's account/partition
before running a full job.
EOF
}

RUN_LABEL="run_hpc"
RUN_DEFAULT_CONFIGFILE="configs/mini_data_public_2023_hpc.yml"
RUN_HELP_FN="run_hpc_help"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
# shellcheck source=./_run_common.sh
source "${SCRIPT_DIR}/_run_common.sh"

# --- HPC-specific: submission profile must exist ----------------------------
if [[ ! -d "${REPO_ROOT}/profiles/slurm" ]]; then
    echo "ERROR: profiles/slurm not found." >&2
    echo "Create a Snakemake profile for your cluster (account, partition," >&2
    echo "per-rule resources, executor plugin config) before running on HPC." >&2
    exit 1
fi

HPC_MAX_JOBS="${HPC_MAX_JOBS:-50}"
HPC_MEM_GB="${HPC_MEM_GB:-256}"

exec snakemake \
    -s "$SNAKEFILE" \
    --configfile "$CONFIGFILE" \
    --profile profiles/slurm \
    --resources "mem_gb=${HPC_MEM_GB}" \
    -j "$HPC_MAX_JOBS" \
    --use-conda \
    --rerun-incomplete \
    --keep-going \
    "${SNAKEMAKE_ARGS[@]+"${SNAKEMAKE_ARGS[@]}"}"
