#!/usr/bin/env bash
#SBATCH --job-name=mini_bam_orchestrator
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=00:15:00
#SBATCH --output=./logs/orchestrator_%j.out
#SBATCH --error=./logs/orchestrator_%j.err
#
# sbatch_submit_mini_bam_batch.sh
#
# Thin orchestrator job. Does the same pre-flight directory/BAM-count
# checks as submit_mini_bam_batch.sh, but runs them on a COMPUTE NODE
# instead of the login node -- scanning 23 directories (some with tens
# of thousands of small BAM chunks) over network storage is exactly
# the kind of I/O that's polite to keep off a shared login node.
#
# It then submits the 23 per-sample extraction jobs via sbatch from
# inside this job -- nested sbatch calls (a job submitting more jobs)
# are supported by SLURM and are a normal pattern, not a hack.
#
# Usage:
#   sbatch workflow/scripts/mini_genomes/slurm/sbatch_submit_mini_bam_batch.sh
#   sbatch workflow/scripts/mini_genomes/slurm/sbatch_submit_mini_bam_batch.sh --dry-run

set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$REPO_ROOT"

echo "Orchestrator running on: $(hostname)"
echo "Mode: ${1:-submit}"
echo ""

# Delegates to the existing submit_mini_bam_batch.sh -- same logic,
# just executed here instead of on the login node.
exec workflow/scripts/mini_genomes/slurm/submit_mini_bam_batch.sh "${1:-}"