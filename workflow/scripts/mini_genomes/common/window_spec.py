"""
window_spec.py

Defines the set of hg38 regions used to build the mini-genome, and
validates them. Two kinds of windows live in the same structure:

  1. "Truncation" windows (the professor's spec) -- always start at
     hg38 position 1, sized to exercise coverage QC / pipeline plumbing.
     hi/lo spans are expected-coverage labels for testing, not biology.

  2. "Locus" windows -- small slices centered on a real gene/CpG-island
     of interest, extracted from wherever that gene actually sits in
     hg38. These exist so that annotation-heavy code paths (methylartist,
     CpG-island histograms) have real, non-empty, biologically meaningful
     data to plot.

Both kinds share the same shape so make_mini_genome.py / make_mini_gtf.py
/ translate_bed.py don't need to know which is which -- they just need
`hg38_start` (defaults to 1) and `seq` (the mini-contig length).
"""

from dataclasses import dataclass, field


@dataclass
class Window:
    chrom: str          # hg38 contig name, e.g. "chr11"
    seq: int             # length of the mini contig
    hg38_start: int = 1   # 1-based hg38 coordinate the mini contig starts at
    spans: list = field(default_factory=list)   # [(start, end, label), ...] in MINI coords (1-based)
    mini_name: str = None  # name to give this contig in the mini genome; defaults to chrom

    def __post_init__(self):
        if self.mini_name is None:
            self.mini_name = self.chrom

    @property
    def hg38_end(self):
        # inclusive
        return self.hg38_start + self.seq - 1

    def to_hg38(self, mini_pos):
        """Translate a 1-based mini-contig coordinate to 1-based hg38 coordinate."""
        return self.hg38_start + (mini_pos - 1)

    def from_hg38(self, hg38_pos):
        """Translate a 1-based hg38 coordinate to 1-based mini-contig coordinate.
        Returns None if hg38_pos falls outside this window."""
        if not (self.hg38_start <= hg38_pos <= self.hg38_end):
            return None
        return hg38_pos - self.hg38_start + 1


def validate_windows(windows):
    """Raise a clear error on span overlaps/gaps/out-of-range issues.
    Fail loud at startup, per the pipeline's existing convention --
    do NOT let a malformed spec silently produce truncated/empty output."""
    errors = []
    for w in windows:
        if w.chrom != w.chrom.strip():
            errors.append(f"{w.chrom}: whitespace in chrom name")
        if w.chrom.startswith("chr") and w.chrom[3:].islower() and w.chrom[3:] not in ("x", "y", "m"):
            pass  # numeric chroms are fine either way
        if w.chrom in ("chrx", "chry", "chrm"):
            errors.append(
                f"{w.chrom}: hg38 uses '{w.chrom[:3]}{w.chrom[3:].upper()}' "
                f"(case-sensitive contig names) -- fix before extraction."
            )

        prev_end = 0
        for (start, end, label) in sorted(w.spans):
            if start <= prev_end:
                errors.append(
                    f"{w.mini_name}: span ({start},{end},'{label}') overlaps "
                    f"or duplicates previous span ending at {prev_end}"
                )
            if end > w.seq:
                errors.append(
                    f"{w.mini_name}: span end {end} exceeds contig length {w.seq}"
                )
            prev_end = end
        if w.spans and prev_end != w.seq:
            errors.append(
                f"{w.mini_name}: spans cover up to {prev_end} but contig length is "
                f"{w.seq} (gap of {w.seq - prev_end} bp uncovered)"
            )

    if errors:
        raise ValueError(
            "window_spec validation failed:\n  " + "\n  ".join(errors)
        )


# ---------------------------------------------------------------------
# Job A: the professor's literal truncation spec.
# hg38_start defaults to 1 for all of these -- no translation needed,
# mini coordinate == hg38 coordinate.
# ---------------------------------------------------------------------
TRUNCATION_WINDOWS = [
    Window("chr11", 200_000, spans=[(1, 11_881, "lo"), (11_882, 115_217, "hi"), (115_218, 200_000, "lo")]),
    Window("chr5",  200_000, spans=[(1, 160_168, "hi"), (160_169, 200_000, "lo")]),
    Window("chr22", 200_000, spans=[(1, 20_202, "lo"), (20_203, 175_539, "hi"), (175_540, 200_000, "lo")]),
    Window("chr10", 300_000, spans=[(1, 39_448, "lo"), (39_449, 237_719, "hi"), (237_720, 300_000, "lo")]),  # NOTE: fixed 29_448 -> 39_448, see below
    Window("chr12", 300_000, spans=[(1, 138_570, "lo"), (138_571, 238_837, "hi"), (238_838, 300_000, "lo")]),
    Window("chr7",  300_000, spans=[(1, 254_591, "hi"), (254_592, 300_000, "lo")]),
    Window("chr19", 300_000, spans=[(1, 62_589, "lo"), (62_590, 259_769, "hi"), (259_770, 300_000, "lo")]),
    Window("chrX", 4_200_000, spans=[(1, 4_200_000, "hi")]),
]

# ---------------------------------------------------------------------
# Job B: locus windows for real genes/CpG islands, so annotation-driven
# plots (methylartist, CpG histograms) have real, non-empty data.
# Fill this in via build_locus_windows.py (below) or by hand.
# Example entry, once populated:
#   Window("chrX", 100_000, hg38_start=67_500_000, mini_name="AR_locus",
#          spans=[(1, 100_000, "hi")]),
# ---------------------------------------------------------------------
LOCUS_WINDOWS = [
    Window("chrX", 196_598, hg38_start=67_539_022, mini_name="AR_locus", spans=[(1, 196_598, "hi")]),
    Window("chrX", 56_707, hg38_start=40_046_252, mini_name="BCOR_locus", spans=[(1, 56_707, "hi")]),
    # add more via: python build_locus_windows.py methylartist_genes.bed [flank_bp]
]

ALL_WINDOWS = TRUNCATION_WINDOWS + LOCUS_WINDOWS


if __name__ == "__main__":
    validate_windows(ALL_WINDOWS)
    print(f"OK: {len(ALL_WINDOWS)} windows validated "
          f"({len(TRUNCATION_WINDOWS)} truncation, {len(LOCUS_WINDOWS)} locus)")
