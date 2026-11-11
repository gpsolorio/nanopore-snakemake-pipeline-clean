#!/usr/bin/env python3
"""
translate_bed.py

Takes any hg38-coordinate BED file (CpG islands, methylartist gene
regions, whatever) and translates it into mini-genome coordinates using
window_spec.ALL_WINDOWS. Entries with no overlap in any window are
dropped and counted -- run this BEFORE assuming an annotation file will
be usable against the mini-genome.

Usage:
    python translate_bed.py <input.bed> <output.bed>

BED is 0-based half-open: [start, end).
"""
import sys
from window_spec import ALL_WINDOWS, validate_windows


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    in_path, out_path = sys.argv[1], sys.argv[2]

    validate_windows(ALL_WINDOWS)

    windows_by_chrom = {}
    for w in ALL_WINDOWS:
        windows_by_chrom.setdefault(w.chrom, []).append(w)

    kept, dropped, malformed = 0, 0, 0
    dropped_examples = []
    malformed_examples = []

    with open(in_path) as fin, open(out_path, "w") as fout:
        for lineno, line in enumerate(fin, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.split()
            if len(fields) < 3:
                malformed += 1
                if len(malformed_examples) < 5:
                    malformed_examples.append(f"line {lineno}: {line!r}")
                continue
            try:
                chrom, start0, end0 = fields[0], int(fields[1]), int(fields[2])
            except ValueError:
                malformed += 1
                if len(malformed_examples) < 5:
                    malformed_examples.append(f"line {lineno}: {line!r}")
                continue
            rest = fields[3:]

            # BED half-open -> 1-based inclusive for window math
            start1, end1 = start0 + 1, end0

            candidates = windows_by_chrom.get(chrom, [])
            wrote = False
            for w in candidates:
                ov_start = max(start1, w.hg38_start)
                ov_end = min(end1, w.hg38_end)
                if ov_start > ov_end:
                    continue

                mini_start = w.from_hg38(ov_start)
                mini_end = w.from_hg38(ov_end)
                # back to 0-based half-open for BED output
                fout.write("\t".join(
                    [w.mini_name, str(mini_start - 1), str(mini_end)] + rest
                ) + "\n")
                wrote = True

            if wrote:
                kept += 1
            else:
                dropped += 1
                if len(dropped_examples) < 5:
                    dropped_examples.append(line)

    print(f"kept {kept}, dropped {dropped} (no overlap with any window), "
          f"skipped {malformed} malformed line(s)", file=sys.stderr)
    if malformed_examples:
        print("malformed lines (skipped, not translated):", file=sys.stderr)
        for ex in malformed_examples:
            print(f"  {ex}", file=sys.stderr)
    if dropped_examples:
        print("example dropped entries:", file=sys.stderr)
        for ex in dropped_examples:
            print(f"  {ex}", file=sys.stderr)


if __name__ == "__main__":
    main()