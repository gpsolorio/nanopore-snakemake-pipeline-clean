#!/usr/bin/env python3
"""
make_mini_genome.py

Reads hg38 and extracts each window in window_spec.ALL_WINDOWS as a
contig in a new mini FASTA. Also writes:
  - mini_coverage_spans.bed  (the hi/lo spans, in MINI coordinates --
    directly usable as your CpG-histogram/coverage-QC input)
  - coordinate_map.tsv       (mini_name, mini_start, mini_end,
    hg38_chrom, hg38_start, hg38_end -- one row per window, for
    translating any other BED/GTF through translate_bed.py)

Requires: samtools (uses `samtools faidx` for indexed region extraction
so it never loads the full hg38 FASTA into memory).

Usage:
    python make_mini_genome.py <hg38.fa> <out_dir>
"""
import subprocess
import sys
import os

from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "common"))

from window_spec import ALL_WINDOWS, validate_windows


def extract_region(hg38_fa, chrom, start, end):
    """samtools faidx region extraction; start/end are 1-based inclusive."""
    region = f"{chrom}:{start}-{end}"
    result = subprocess.run(
        ["samtools", "faidx", hg38_fa, region],
        capture_output=True, text=True, check=True,
    )
    lines = result.stdout.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise RuntimeError(f"samtools faidx returned unexpected output for {region}")
    seq = "".join(lines[1:])
    return seq


def write_fasta_record(fh, name, seq, width=60):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    hg38_fa, out_dir = sys.argv[1], sys.argv[2]

    validate_windows(ALL_WINDOWS)
    os.makedirs(out_dir, exist_ok=True)

    fasta_path = os.path.join(out_dir, "mini_genome.fa")
    beds_dir = os.path.join(out_dir, "beds")
    os.makedirs(beds_dir, exist_ok=True)
    spans_bed_path = os.path.join(beds_dir, "mini_spans.bed")
    map_path = os.path.join(out_dir, "coordinate_map.tsv")

    with open(fasta_path, "w") as fasta_fh, \
         open(spans_bed_path, "w") as bed_fh, \
         open(map_path, "w") as map_fh:

        map_fh.write("mini_name\tmini_start\tmini_end\thg38_chrom\thg38_start\thg38_end\n")

        for w in ALL_WINDOWS:
            print(f"extracting {w.mini_name}: {w.chrom}:{w.hg38_start}-{w.hg38_end} "
                  f"({w.seq:,} bp)", file=sys.stderr)
            seq = extract_region(hg38_fa, w.chrom, w.hg38_start, w.hg38_end)
            if len(seq) != w.seq:
                raise RuntimeError(
                    f"{w.mini_name}: expected {w.seq} bp, got {len(seq)} bp "
                    f"-- check hg38 contig length / window bounds"
                )
            write_fasta_record(fasta_fh, w.mini_name, seq)

            map_fh.write(f"{w.mini_name}\t1\t{w.seq}\t{w.chrom}\t{w.hg38_start}\t{w.hg38_end}\n")

            for (start, end, label) in w.spans:
                # BED is 0-based half-open
                bed_fh.write(f"{w.mini_name}\t{start - 1}\t{end}\t{w.mini_name}_{label}\n")

    print(f"\nWrote:\n  {fasta_path}\n  {spans_bed_path}\n  {map_path}", file=sys.stderr)
    print("Next: samtools faidx mini_genome.fa   (to index the new reference)", file=sys.stderr)


if __name__ == "__main__":
    main()