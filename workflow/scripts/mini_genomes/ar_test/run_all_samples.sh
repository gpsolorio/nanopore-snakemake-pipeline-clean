#!/bin/bash
set -euo pipefail

##############################################################################
# run_all_samples.sh
#
# Runs extract_ar_reads.sh + realign_to_mini.sh across every the private cohort proband
# under a given base directory, matching the naming convention:
#   <base>/<run_dir>/<ID>_5mC-5hmC/pass
#
# Usage: run_all_samples.sh <region> [batch_size]
#   region     : ORIGINAL hg38 coordinates, e.g. chrX:67544032-67730619
#   batch_size : optional, passed through to extract_ar_reads.sh (default 500)
##############################################################################

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <region> [batch_size]" >&2
    exit 1
fi

REGION="$1"
BATCH_SIZE="${2:-500}"

# Ensure the conda env with minimap2/winnowmap2/meryl is active in THIS shell.
# Running on a different node/session than where you last `conda activate`d
# means these tools won't be on PATH even if samtools happens to be reachable
# some other way -- that's exactly the failure mode this guards against.
if ! command -v minimap2 >/dev/null 2>&1 || ! command -v winnowmap >/dev/null 2>&1 || ! command -v meryl >/dev/null 2>&1; then
    echo "One or more required tools not on PATH — attempting to activate ngs-utils-env..."
    CONDA_BASE="$(conda info --base 2>/dev/null || true)"
    if [[ -n "$CONDA_BASE" && -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
        # shellcheck disable=SC1091
        source "${CONDA_BASE}/etc/profile.d/conda.sh"
        conda activate ngs-utils-env
    fi
    # Re-check after attempting activation
    if ! command -v minimap2 >/dev/null 2>&1 || ! command -v winnowmap >/dev/null 2>&1 || ! command -v meryl >/dev/null 2>&1; then
        echo "ERROR: minimap2/winnowmap/meryl still not found after attempting to activate ngs-utils-env." >&2
        echo "Run 'conda activate ngs-utils-env' manually in this shell first, then rerun this script." >&2
        exit 1
    fi
    echo "ngs-utils-env activated successfully."
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${1:?Usage: $0 samples.tsv  (columns: sample_id<TAB>pass_dir)}"
[[ -f "$MANIFEST" ]] || { echo "Manifest not found: $MANIFEST" >&2; exit 1; }

while IFS=$'\t' read -r sample_id pass_dir; do
    [[ -z "$sample_id" || "$sample_id" == \#* ]] && continue
    if [[ ! -d "$pass_dir" ]]; then
        echo "  [skip] no pass/ dir for ${sample_id}: ${pass_dir}" >&2
        continue
    fi
    label="${sample_id}_AR"
    echo ""
    echo "=== Sample ${sample_id} -> label ${label} ==="
    "${SCRIPT_DIR}/extract_ar_reads.sh" "$pass_dir" "$REGION" "$label" "$BATCH_SIZE"
    "${SCRIPT_DIR}/realign_to_mini.sh" "$label"
done < "$MANIFEST"

echo ""
echo "All samples processed. Results under ../mini-data/aligned/<sample>_AR/"