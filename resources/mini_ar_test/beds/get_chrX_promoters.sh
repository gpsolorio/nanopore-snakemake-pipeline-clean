#!/bin/bash
set -euo pipefail

##############################################################################
# get_chrX_promoters.sh
#
# Builds a BED file of promoter regions for every chrX RefSeq transcript,
# defined as a strand-aware fixed window around the TSS (transcription
# start site): --upstream bases before TSS, --downstream bases after.
#
# Source: UCSC refGene table dump (hgdownload.soe.ucsc.edu) -- this is
# UCSC's own bulk-download endpoint, built for exactly this kind of use.
# NOT GeneCards -- GeneCards' terms of use explicitly prohibit scraping
# and bulk automated extraction, and it isn't a primary coordinate source
# anyway (it aggregates from RefSeq/Ensembl/GeneHancer for human reading,
# not bulk programmatic export).
#
# Usage:
#   ./get_chrX_promoters.sh [--genome hg38] [--upstream 2000] [--downstream 1000] [--out FILE]
#
# Output: BED6 (chrom, start, end, name, score, strand), sorted, one row
# per transcript. Multiple transcripts of the same gene with different
# TSS positions each get their own row -- see notes at the bottom for how
# to collapse to one row per gene, or merge overlapping windows.
##############################################################################

GENOME="hg38"
UPSTREAM=2000
DOWNSTREAM=1000
OUT_BED=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome)     GENOME="$2"; shift 2 ;;
        --upstream)   UPSTREAM="$2"; shift 2 ;;
        --downstream) DOWNSTREAM="$2"; shift 2 ;;
        --out)        OUT_BED="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$OUT_BED" ]]; then
    OUT_BED="${GENOME}_chrX_promoters_up${UPSTREAM}_down${DOWNSTREAM}.bed"
fi

TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

echo "Genome:      $GENOME"
echo "Window:      -${UPSTREAM} / +${DOWNSTREAM} around TSS (strand-aware)"
echo "Output:      $OUT_BED"
echo ""

REFGENE_URL="https://hgdownload.soe.ucsc.edu/goldenPath/${GENOME}/database/refGene.txt.gz"
echo "Downloading refGene table from UCSC..."
curl -fSL "$REFGENE_URL" -o "${TMP_DIR}/refGene.txt.gz"

if [[ ! -s "${TMP_DIR}/refGene.txt.gz" ]]; then
    echo "Download appears to have failed: ${TMP_DIR}/refGene.txt.gz is missing or empty" >&2
    exit 1
fi

echo "Filtering to chrX and computing promoter windows..."
# refGene.txt columns (0-indexed in UCSC docs, 1-indexed here for awk):
#   $2  name    (RefSeq transcript accession, e.g. NM_...)
#   $3  chrom
#   $4  strand
#   $5  txStart
#   $6  txEnd
#   $13 name2   (gene symbol)
#
# Use gunzip -c, not zcat: macOS's built-in /usr/bin/zcat is the BSD
# variant for .Z (compress) files, not .gz (gzip) -- it will fail with
# a confusing "No such file" error on any real .gz file. gunzip -c
# handles .gz correctly on both macOS and Linux.
gunzip -c "${TMP_DIR}/refGene.txt.gz" | awk -v up="$UPSTREAM" -v down="$DOWNSTREAM" '
BEGIN { OFS = "\t" }
$3 == "chrX" {
    chrom  = $3
    strand = $4
    txStart = $5
    txEnd   = $6
    name2   = $13
    name    = $2

    if (strand == "+") {
        tss   = txStart
        start = tss - up
        end   = tss + down
    } else {
        tss   = txEnd
        start = tss - down
        end   = tss + up
    }
    if (start < 0) start = 0

    print chrom, start, end, name2"_"name, 0, strand
}' > "${TMP_DIR}/promoters.unsorted.bed"

sort -k1,1 -k2,2n "${TMP_DIR}/promoters.unsorted.bed" | uniq > "$OUT_BED"

echo ""
echo "Wrote $(wc -l < "$OUT_BED" | tr -d ' ') promoter regions to $OUT_BED"
echo ""
echo "Notes:"
echo "  - One row per transcript. Genes with multiple RefSeq transcripts"
echo "    (different TSS) get multiple rows -- this is usually what you"
echo "    want for methylation work, since alternative promoters can"
echo "    have independently regulated methylation."
echo "  - To collapse to one row per gene (using the union of all its"
echo "    transcript promoter windows), pipe through:"
echo "      sort -k1,1 -k2,2n $OUT_BED | bedtools merge -i - -c 4 -o distinct"
echo "  - Verify chromosome naming matches your reference FASTA."
echo "    UCSC uses 'chrX'; Ensembl-style references often use bare 'X'."
echo "    Your pipeline already has workflow/scripts/utils/chr_name_conv.txt"
echo "    for this -- reuse it here if needed, e.g.:"
echo "      bcftools annotate --rename-chrs chr_name_conv.txt ..."
echo "    (or a simple sed 's/^chr//' if going the other direction)."