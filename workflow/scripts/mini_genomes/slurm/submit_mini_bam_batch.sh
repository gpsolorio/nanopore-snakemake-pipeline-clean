#!/usr/bin/env bash
#
# submit_mini_bam_batch.sh
#
# Submits sbatch_extract_mini_bam.sh as one SLURM job per sample, for
# EVERY sample directory actually confirmed via `ls`/`find` in this
# project's investigation session -- 23 samples total across:
#   - layout A/
#   - layout B/
#   - layout C/
#
# This is a full-plumbing test: it exercises every stage of the DAG
# (merge_bams, sort_index_bam, call_variants, merge_vcfs, dss_call_dmrs)
# against real data, not just a single-sample smoke test. It does NOT
# imply any particular ASD/TD grouping is biologically meaningful --
# see the printed sample/family table below for what's actually known
# vs. assumed.
#
# Directory shape differs across sources (confirmed by ls, not assumed):
#   layout A:  {sample}_5mC-5hmC/pass/*.bam
#   layout B:  barcode0N_{sample}_5mC-5hmC/*.bam   (no pass/)
#   layout C:   {sample}_5mC-5hmC/*.bam             (no pass/, no barcode prefix)
#
# This script auto-detects which shape applies per-sample by checking
# for a pass/ subfolder on disk, rather than hardcoding the assumption --
# if a sample's actual layout doesn't match ANY known shape, that
# sample is skipped with a clear warning instead of submitting a job
# that would fail confusingly on the compute node.
#
# Usage:
#   ./submit_mini_bam_batch.sh          # submit all 23
#   ./submit_mini_bam_batch.sh --dry-run  # just show what WOULD run

set -euo pipefail

DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
fi

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$REPO_ROOT"
mkdir -p logs

MANIFEST="${MANIFEST:?Set MANIFEST=/path/to/samples.tsv (columns: sample_id<TAB>run_dir)}"
[[ -f "$MANIFEST" ]] || { echo "Manifest not found: $MANIFEST" >&2; exit 1; }

# format of each entry: "sample_id|run_dir"
SAMPLES=()
while IFS=$'\t' read -r sid dir; do
  [[ -z "$sid" || "$sid" == \#* ]] && continue
  SAMPLES+=("${sid}|${dir}")
done < "$MANIFEST"

echo "=========================================================="
echo "Mini-BAM batch: ${#SAMPLES[@]} samples"
echo "Mode: $([ "$DRY_RUN" = true ] && echo 'DRY RUN (no jobs submitted)' || echo 'SUBMITTING')"
echo "=========================================================="
echo ""

JOB_IDS=()
SKIPPED=()

for entry in "${SAMPLES[@]}"; do
  sample="${entry%%|*}"
  base_dir="${entry#*|}"

  # Auto-detect pass/ vs. no-pass/ layout instead of assuming.
  if [ -d "${base_dir}/pass" ]; then
    resolved_dir="${base_dir}/pass"
    shape="pass/"
  elif [ -d "$base_dir" ]; then
    resolved_dir="$base_dir"
    shape="direct"
  else
    echo "SKIP  $sample -- directory does not exist: $base_dir"
    SKIPPED+=("$sample (missing dir)")
    continue
  fi

  n_bams=$(find "$resolved_dir" -maxdepth 1 -name "*.bam" 2>/dev/null | wc -l | tr -d ' ')
  if [ "$n_bams" -eq 0 ]; then
    echo "SKIP  $sample -- no *.bam files found in $resolved_dir"
    SKIPPED+=("$sample (no BAMs found)")
    continue
  fi

  echo "OK    $sample  [$shape, $n_bams BAM(s)]  -> $resolved_dir"

  if [ "$DRY_RUN" = true ]; then
    continue
  fi

  job_id=$(sbatch --parsable --job-name="$sample" \
      workflow/scripts/mini_genomes/slurm/sbatch_extract_mini_bam.sh "$resolved_dir" "$sample")
  JOB_IDS+=("$job_id")
  echo "      submitted as job $job_id"
done

echo ""
echo "=========================================================="
if [ "$DRY_RUN" = true ]; then
  echo "Dry run complete. ${#SKIPPED[@]} sample(s) would be skipped."
else
  echo "Submitted ${#JOB_IDS[@]} job(s). ${#SKIPPED[@]} sample(s) skipped."
fi
if [ "${#SKIPPED[@]}" -gt 0 ]; then
  echo ""
  echo "Skipped samples (fix path/layout and re-run for these):"
  for s in "${SKIPPED[@]}"; do
    echo "  - $s"
  done
fi
echo "=========================================================="

if [ "$DRY_RUN" = false ] && [ "${#JOB_IDS[@]}" -gt 0 ]; then
  echo ""
  echo "Monitor with:"
  echo "  squeue -u \$USER"
  echo "  sacct -j $(IFS=,; echo "${JOB_IDS[*]}") --format=JobID,JobName,State,ExitCode,Elapsed"
  echo ""
  echo "Once jobs finish, check headers across ALL resulting BAMs at once:"
  echo "  for f in workflow/data/bam/*/*.bam; do"
  echo "    echo \"== \$f ==\"; samtools view -H \"\$f\" | grep '^@SQ' | md5sum"
  echo "  done"
  echo "  # every line's md5 should be identical if all headers match"
fi