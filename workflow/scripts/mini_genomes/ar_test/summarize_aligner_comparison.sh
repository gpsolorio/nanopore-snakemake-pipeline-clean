#!/bin/bash
set -euo pipefail

##############################################################################
# summarize_aligner_comparison.sh
#
# Prints a compact table comparing minimap2 vs winnowmap2 mapping stats
# across every sample in ../mini-data/aligned/.
##############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALIGNED_DIR="${SCRIPT_DIR}/../mini-data/aligned"

printf "%-10s %-12s %10s %10s %10s %10s\n" "sample" "aligner" "primary" "mapped" "pct_mapped" "supp"
printf "%-10s %-12s %10s %10s %10s %10s\n" "------" "-------" "-------" "------" "----------" "----"

for sample_dir in "$ALIGNED_DIR"/*_AR; do
    sample=$(basename "$sample_dir" | sed 's/_AR$//')
    for aligner in minimap2 winnowmap2; do
        bam="${sample_dir}/${aligner}.bam"
        [[ -f "$bam" ]] || continue

        stats=$(samtools flagstat "$bam")
        primary=$(echo "$stats" | grep -m1 " primary$" | awk '{print $1}')
        primary_mapped_line=$(echo "$stats" | grep "primary mapped")
        primary_mapped=$(echo "$primary_mapped_line" | awk '{print $1}')
        pct=$(echo "$primary_mapped_line" | grep -oP '\(\K[0-9.]+(?=%)')
        supp=$(echo "$stats" | grep -m1 " supplementary$" | awk '{print $1}')

        printf "%-10s %-12s %10s %10s %9s%% %10s\n" "$sample" "$aligner" "$primary" "$primary_mapped" "$pct" "$supp"
    done
done