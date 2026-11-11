#!/usr/bin/env python3
"""
print_hg38_regions.py

Prints one samtools-style region string per window in window_spec.ALL_WINDOWS
(hg38 coordinates -- these are the regions to slice OUT of the original,
hg38-aligned BAMs). Single source of truth: whenever you add a locus window
in window_spec.py, this picks it up automatically -- no need to hand-edit
extract_mini_bam.sh's MINI_REGIONS separately.

Usage:
    python print_hg38_regions.py
    # or, to use directly in bash:
    MINI_REGIONS=($(python print_hg38_regions.py))
"""
from window_spec import ALL_WINDOWS, validate_windows

if __name__ == "__main__":
    validate_windows(ALL_WINDOWS)
    for w in ALL_WINDOWS:
        print(f"{w.chrom}:{w.hg38_start}-{w.hg38_end}")
