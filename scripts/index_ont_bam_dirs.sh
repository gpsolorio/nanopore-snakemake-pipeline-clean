#!/usr/bin/env bash
##############################################################################
# scripts/index_ont_bam_dirs.sh — Build a flat BAM_INPUT_DIR symlink index
#
# common.smk's get_bam_paths() has a fast path:
#   BAM_INPUT_DIR/{barcode}/*.bam
# and only falls back to a recursive `BAM_INPUT_DIR/**/*.bam` walk (over
# the WHOLE tree, matching barcode-as-substring) if that's empty. Real
# methylation-calling delivery folders are named "<sample_id>_5mC-5hmC"
# (layout A/, layout C/) or "barcode0N_<sample_id>_5mC-5hmC" (layout B/),
# never a bare sample ID -- so every sample falls through to the slow
# recursive walk over the entire data tree on EVERY snakemake invocation
# (dry-runs included), which only gets worse as more pools land.
#
# IMPORTANT: only "*_5mC-5hmC" directories are indexed. Per genome center
# (Ruta), BAM files anywhere else in a run directory (bam_pass/bam_fail,
# pod5/basecalling/pass/barcodeNN/) are fast-mode, alignment-only output
# from in-run adaptive-sampling basecalling and do NOT contain methylation
# calls. Only the "*_5mC-5hmC" delivery folders (re-basecalled offline in
# accurate mode with 5mCG_5hmCG modbase calling) are valid pipeline input.
# This script will not find/index anything outside those folders.
#
# CHANGED: some "*_5mC-5hmC" delivery folders (confirmed so far: all
# checked layout A/ dirs) split their BAMs one level deeper into a pass/
# subdirectory (reads that passed the basecall quality filter) alongside
# fail/ (discarded reads -- never indexed). Others (confirmed: layout B
# barcode dirs) have BAMs directly inside the "*_5mC-5hmC" dir with no
# pass/fail split at all. Since this has already differed between two
# folders in the SAME data tree, and per the person running this the
# project has moving parts, this script now resolves each delivery folder
# dynamically per-directory instead of assuming one convention:
#   1. If a pass/ subdirectory exists AND contains *.bam files directly
#      inside it -> index that pass/ subdirectory (skip fail/ always).
#   2. Else if *.bam files exist directly inside the "*_5mC-5hmC" dir
#      itself -> index that directly (old flat-layout behavior).
#   3. Else -> nothing usable found; loudly skip and report, don't guess.
# This auto-detects per sample, so it keeps working regardless of which
# layout any given batch/pool uses, without needing a manual flag.
#
# This does a recursive scan of the real data root(s) and creates
# symlinks:
#   <index_dir>/<sample_id>  ->  <resolved BAM-containing directory>
# so subsequent Snakemake runs hit get_bam_paths()'s fast, non-recursive
# path. No BAM data is copied or moved -- these are symlinks only, and
# nothing under the real data root(s) is ever written to, renamed, or
# deleted by this script.
#
# Sample ID extraction handles both conventions seen in this data tree:
#   "sample_001_5mC-5hmC"          -> sample_001   (layout A/, layout C/)
#   "barcode07_sample_006_5mC-5hmC" -> sample_006  (layout B/, barcode prefix stripped)
#
# Usage:
#   bash scripts/index_ont_bam_dirs.sh [--verify <barcode_reference.tsv>] \
#       <index_dir> <data_root> [<data_root> ...]
#
# --verify <tsv>  Optional. A 2-column TSV (header: barcode<TAB>sample_id)
#   -- e.g. pool_barcode_reference.tsv, built from the genome center
#   spreadsheet. For every "barcode0N_<sample_id>_5mC-5hmC" directory
#   found, cross-checks that the embedded sample_id matches what the TSV
#   says barcode0N should be. Mismatches are hard errors -- this is the
#   patient-data-correctness safety net: it's the only file being trusted
#   as a check, not as the source of truth (the directory name is), so a
#   future stale/wrong spreadsheet can only cause a loud failure here, not
#   a silent misattribution.
#
# Example:
#   bash scripts/index_ont_bam_dirs.sh \
#       --verify pool_barcode_reference.tsv \
#       /path/to/bam_index \
#       /path/to/data \
#       /path/to/data
#
# Then in your config:
#   paths:
#     inputs:
#       bam: /path/to/bam_index
##############################################################################

set -euo pipefail

VERIFY_TSV=""
if [[ "${1:-}" == "--verify" ]]; then
    VERIFY_TSV="$2"
    shift 2
fi

if [[ $# -lt 2 ]]; then
    echo "Usage: bash scripts/index_ont_bam_dirs.sh [--verify <barcode_reference.tsv>] <index_dir> <data_root> [<data_root> ...]" >&2
    exit 1
fi

INDEX_DIR="$1"; shift
DATA_ROOTS=("$@")

mkdir -p "$INDEX_DIR"

# Load the verification reference, if given: barcode -> expected sample_id
declare -A EXPECTED_SAMPLE
if [[ -n "$VERIFY_TSV" ]]; then
    if [[ ! -f "$VERIFY_TSV" ]]; then
        echo "ERROR: --verify file not found: $VERIFY_TSV" >&2
        exit 1
    fi
    while IFS=$'\t' read -r barcode sample_id; do
        [[ "$barcode" == "barcode" ]] && continue   # skip header
        [[ -z "$barcode" ]] && continue
        EXPECTED_SAMPLE["$barcode"]="$sample_id"
    done < "$VERIFY_TSV"
    echo ">>> Loaded ${#EXPECTED_SAMPLE[@]} barcode->sample_id reference entries from $VERIFY_TSV"
fi

# "sample_001_5mC-5hmC"           -> sample_id=sample_001, barcode=(none)
# "barcode07_sample_006_5mC-5hmC" -> sample_id=sample_006, barcode=barcode07
_extract_sample_id() {
    basename "$1" | sed -E 's/^(barcode[0-9]+_)?([0-9]+)_.*/\2/'
}
_extract_barcode() {
    local base; base="$(basename "$1")"
    if [[ "$base" =~ ^(barcode[0-9]+)_ ]]; then
        echo "${BASH_REMATCH[1]}"
    else
        echo ""
    fi
}

# Read-only. Given a "*_5mC-5hmC" delivery directory, decides where the
# real BAM files actually live and prints that path on stdout. Never
# creates, deletes, or modifies anything -- only `find`/test operations.
# Returns non-zero (and prints nothing) if no usable BAMs were found.
_resolve_bam_dir() {
    local delivery_dir="$1"

    # Prefer a pass/ subdirectory (case-insensitive) if one exists AND
    # actually has *.bam files directly inside it. fail/ is intentionally
    # never considered -- those are discarded reads regardless of layout.
    local pass_dir
    pass_dir="$(find "$delivery_dir" -mindepth 1 -maxdepth 1 -type d -iname "pass" 2>/dev/null | head -1)"
    if [[ -n "$pass_dir" ]]; then
        if find "$pass_dir" -maxdepth 1 -iname "*.bam" -print -quit 2>/dev/null | grep -q .; then
            echo "$pass_dir"
            return 0
        else
            echo "  NOTE: $delivery_dir has a pass/ subdir but no *.bam directly inside it -- checking flat layout instead" >&2
        fi
    fi

    # No pass/ dir, or it existed but was empty of BAMs: check the old
    # flat layout (BAMs directly inside the "*_5mC-5hmC" dir, as in the
    # layout B barcode folders).
    if find "$delivery_dir" -maxdepth 1 -iname "*.bam" -print -quit 2>/dev/null | grep -q .; then
        echo "$delivery_dir"
        return 0
    fi

    return 1
}

declare -A seen   # sample_id -> source delivery dir, to detect duplicates loudly

n_found=0
n_linked=0
n_conflicts=0
n_verify_fail=0
n_no_bam=0

for root in "${DATA_ROOTS[@]}"; do
    if [[ ! -d "$root" ]]; then
        echo "WARNING: data root not found, skipping: $root" >&2
        continue
    fi
    echo ">>> Scanning $root for *_5mC-5hmC directories..."

    while IFS= read -r -d '' dir; do
        n_found=$((n_found + 1))
        sample_id="$(_extract_sample_id "$dir")"
        barcode="$(_extract_barcode "$dir")"

        if [[ ! "$sample_id" =~ ^[0-9]+$ ]]; then
            echo "  SKIP (couldn't extract a numeric sample ID): $dir" >&2
            continue
        fi

        # Cross-check against the genome center reference, if provided.
        if [[ -n "$barcode" && -n "$VERIFY_TSV" ]]; then
            expected="${EXPECTED_SAMPLE[$barcode]:-}"
            if [[ -z "$expected" ]]; then
                echo "  VERIFY WARNING: $barcode not found in $VERIFY_TSV (dir: $dir)" >&2
                n_verify_fail=$((n_verify_fail + 1))
            elif [[ "$expected" != "$sample_id" ]]; then
                echo "  VERIFY FAIL: $dir" >&2
                echo "    Directory name says sample_id=$sample_id for $barcode" >&2
                echo "    $VERIFY_TSV says $barcode -> sample_id=$expected" >&2
                echo "    NOT indexing this directory -- resolve the discrepancy first." >&2
                n_verify_fail=$((n_verify_fail + 1))
                continue
            fi
        fi

        if [[ -n "${seen[$sample_id]:-}" && "${seen[$sample_id]}" != "$dir" ]]; then
            echo "  CONFLICT: sample $sample_id maps to two different directories:" >&2
            echo "      ${seen[$sample_id]}" >&2
            echo "      $dir" >&2
            echo "    Not creating/overwriting a symlink for $sample_id -- resolve manually" >&2
            echo "    (duplicate run? re-sequenced sample? pooled + proband both present?)." >&2
            n_conflicts=$((n_conflicts + 1))
            continue
        fi
        seen[$sample_id]="$dir"

        # CHANGED: resolve pass/ vs flat layout per-directory instead of
        # always linking straight to the "*_5mC-5hmC" dir.
        bam_dir="$(_resolve_bam_dir "$dir" || true)"
        if [[ -z "$bam_dir" ]]; then
            echo "  SKIP (no *.bam found in pass/ or directly in the delivery dir): $dir" >&2
            n_no_bam=$((n_no_bam + 1))
            continue
        fi

        link_path="${INDEX_DIR}/${sample_id}"
        if [[ -L "$link_path" ]]; then
            current_target="$(readlink -f "$link_path")"
            real_target="$(readlink -f "$bam_dir")"
            if [[ "$current_target" == "$real_target" ]]; then
                continue   # already correctly linked
            fi
            echo "  UPDATING stale symlink: $sample_id ($current_target -> $real_target)"
            rm "$link_path"
        elif [[ -e "$link_path" ]]; then
            echo "  SKIP: $link_path exists and is not a symlink -- not touching it" >&2
            n_conflicts=$((n_conflicts + 1))
            continue
        fi

        ln -s "$bam_dir" "$link_path"
        n_linked=$((n_linked + 1))
    done < <(find "$root" -type d -iname '*_5mC-5hmC' -print0)
done

echo ""
echo ">>> Found $n_found candidate director(ies), linked $n_linked, $n_conflicts conflict(s)/skip(s), $n_verify_fail verify issue(s), $n_no_bam with no BAMs found (pass/ or flat)."
echo ">>> Index ready at: $INDEX_DIR"
if [[ "$n_conflicts" -gt 0 || "$n_verify_fail" -gt 0 ]]; then
    echo ">>> Review the CONFLICT/VERIFY lines above before running the pipeline -- those" >&2
    echo ">>> samples are NOT indexed until resolved." >&2
    exit 1
fi
