#!/usr/bin/env bash
# reorganize_mini_resources.sh
# Run from the repo root: nanopore-snakemake-pipeline/
# Uses `git mv` so history is preserved. Safe to run once; re-running
# will just fail on files that no longer exist at the old path.

set -euo pipefail

mkdir -p resources/mini/scripts
mkdir -p resources/mini/beds

# --- existing files that need to move ------------------------------------
[ -f resources/make_mini_genome.py ] && git mv resources/make_mini_genome.py resources/mini/scripts/make_mini_genome.py
[ -f resources/mini_genome.fa ]      && git mv resources/mini_genome.fa      resources/mini/mini_genome.fa
[ -f resources/mini_genome.fa.fai ]  && git mv resources/mini_genome.fa.fai  resources/mini/mini_genome.fa.fai
[ -f resources/mini_genes.gtf ]      && git mv resources/mini_genes.gtf     resources/mini/mini_genes.gtf

echo "Moved existing mini-genome files into resources/mini/"
echo ""
echo "Now copy in the new scripts (not yet in the repo):"
echo "  cp window_spec.py make_mini_gtf.py translate_bed.py build_locus_windows.py resources/mini/scripts/"
echo ""
echo "Then regenerate everything from the new scripts so make_mini_genome.py's"
echo "output matches what's tracked in git (old resources/mini_genome.fa was"
echo "built by the old script, before locus windows existed):"
echo "  cd resources/mini/scripts"
echo "  python make_mini_genome.py ../../hg38.fa .."
echo "  python make_mini_gtf.py ../../hg38_genes.gtf ../mini_genes.gtf"
echo "  python translate_bed.py ../../beds/cpg_islands.bed ../beds/cpg_islands.mini.bed"
echo "  python translate_bed.py ../../beds/methylartist_genes.bed ../beds/methylartist_genes.mini.bed"