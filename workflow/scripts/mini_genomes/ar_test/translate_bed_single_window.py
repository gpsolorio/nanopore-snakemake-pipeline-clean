#!/usr/bin/env python3
"""
Translate hg38 BED features into mini-genome coordinates for a single
extraction window (e.g. the AR-only mini test genome).

Usage:
    translate_bed_single_window.py \
        <in_bed> <win_chrom> <win_start> <win_end> <new_chrom> <out_bed>

    in_bed    : source BED file in ORIGINAL hg38 coordinates
                (e.g. resources/mini_production/beds/hg38_CpG_island_coords.bed)
    win_chrom : original chromosome of the window, e.g. chrX
    win_start : 0-based window start -- MUST exactly match the value passed
                to make_mini_genome.sh when the genome was built
    win_end   : window end -- same requirement
    new_chrom : mini-genome contig name, e.g. AR_locus
    out_bed   : full output path, e.g.
                resources/mini_ar_test/beds/mini_CpG_island_coords.bed

Features outside the window are dropped. No path auto-resolution -- both
in_bed and out_bed are used exactly as given, since input (production
resources) and output (AR-test resources) now live in separate trees.
"""
import sys


def translate(in_bed, win_chrom, win_start, win_end, new_chrom, out_bed):
    n_kept = 0
    with open(in_bed) as fin, open(out_bed, "w") as fout:
        for line in fin:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            fields = line.rstrip("\n").split()
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom != win_chrom:
                continue
            if start < win_start or end > win_end:
                continue
            fields[0] = new_chrom
            fields[1] = str(start - win_start)
            fields[2] = str(end - win_start)
            fout.write("\t".join(fields) + "\n")
            n_kept += 1
    return n_kept


def main():
    if len(sys.argv) != 7:
        print(
            "Usage: translate_bed_single_window.py "
            "<in_bed> <win_chrom> <win_start> <win_end> <new_chrom> <out_bed>",
            file=sys.stderr,
        )
        sys.exit(1)

    in_bed, win_chrom, win_start, win_end, new_chrom, out_bed = sys.argv[1:7]
    win_start, win_end = int(win_start), int(win_end)

    n = translate(in_bed, win_chrom, win_start, win_end, new_chrom, out_bed)
    print(f"Kept {n} feature(s) inside {win_chrom}:{win_start}-{win_end}")
    print(f"Wrote: {out_bed}")


if __name__ == "__main__":
    main()