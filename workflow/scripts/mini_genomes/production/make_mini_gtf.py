#!/usr/bin/env python3
"""
make_mini_gtf.py

Filters a full hg38 GTF/GFF3 down to features that overlap any window
in window_spec.ALL_WINDOWS, translates their coordinates into mini-genome
space, and clips features that cross a window boundary.

Works for both truncation windows (hg38_start=1, translation is a no-op)
and locus windows (hg38_start=arbitrary, translation matters) using the
same code path -- that's the point of keeping them in one Window list.

Usage:
    python make_mini_gtf.py <hg38.gtf> <out.gtf>

GTF/GFF3 are both tab-delimited with columns:
    seqname source feature start end score strand frame attribute
start/end are 1-based inclusive in both formats, so no offset fix needed
between formats.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "common"))

from window_spec import ALL_WINDOWS, validate_windows


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    in_path, out_path = sys.argv[1], sys.argv[2]

    validate_windows(ALL_WINDOWS)

    # index windows by hg38 chrom for fast lookup; a chrom can have
    # multiple windows (e.g. a truncation window + several locus windows)
    windows_by_chrom = {}
    for w in ALL_WINDOWS:
        windows_by_chrom.setdefault(w.chrom, []).append(w)

    kept, dropped, clipped = 0, 0, 0

    with open(in_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqname, source, feature, start, end, score, strand, frame, attr = fields
            start, end = int(start), int(end)

            candidates = windows_by_chrom.get(seqname, [])
            wrote_this_feature = False
            for w in candidates:
                # does [start,end] overlap this window at all?
                ov_start = max(start, w.hg38_start)
                ov_end = min(end, w.hg38_end)
                if ov_start > ov_end:
                    continue  # no overlap with this window

                if ov_start != start or ov_end != end:
                    clipped += 1

                mini_start = w.from_hg38(ov_start)
                mini_end = w.from_hg38(ov_end)

                fout.write("\t".join([
                    w.mini_name, source, feature,
                    str(mini_start), str(mini_end),
                    score, strand, frame, attr,
                ]) + "\n")
                wrote_this_feature = True

            if wrote_this_feature:
                kept += 1
            else:
                dropped += 1

    print(f"kept {kept} features ({clipped} clipped at a window boundary), "
          f"dropped {dropped} (no overlap with any window)", file=sys.stderr)


if __name__ == "__main__":
    main()
