#!/usr/bin/env python3
"""
remap_bam_to_mini.py

REQUIRED post-processing step after slicing reads out of hg38 by region.
extract_mini_bam.sh's samtools view -L slicing extracts the right READS
but does NOT touch their coordinates -- the resulting BAM still has
hg38 RNAME/POS values and an hg38 header. This script rewrites both so
the BAM actually matches resources/mini_genome.fa:

  - Truncation windows (chr11, chr5, ... chrX first-4.2Mb): RNAME/POS
    don't move (hg38_start=1), but the HEADER still needs fixing --
    without this step the BAM header claims chrX is ~156Mb long while
    mini_genome.fa's chrX is 4.2Mb, and samtools/modkit will reject the
    header/reference length mismatch.

  - Locus windows (AR_locus, BCOR_locus): RNAME/POS genuinely need
    translating -- e.g. a read at hg38 chrX:67,540,000 must become
    AR_locus:963 to land correctly on the mini contig, which starts
    numbering at 1 from hg38 position 67,539,022.

Because the mini FASTA sequence at each window is an EXACT copy of the
corresponding hg38 substring (built via samtools faidx), CIGAR and MD
stay valid across this shift -- this is coordinate translation, not
realignment. Reads that are unmapped, fall outside every window, or
whose alignment would hang off the edge of a mini contig are dropped
(counted and reported, not silently lost).

Usage:
    python remap_bam_to_mini.py <window_spec_dir> <input.bam> <output.bam>

Requires pysam (pip install pysam / conda install -c bioconda pysam).
Run BEFORE the final samtools sort + index step -- output here is
coordinate-shifted but not yet guaranteed sorted.
"""
import sys
import pysam


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    window_spec_dir, in_bam, out_bam = sys.argv[1:4]
    sys.path.insert(0, window_spec_dir)
    from window_spec import ALL_WINDOWS, validate_windows  # noqa: E402

    validate_windows(ALL_WINDOWS)

    windows_by_chrom = {}
    for w in ALL_WINDOWS:
        windows_by_chrom.setdefault(w.chrom, []).append(w)

    # New header: one @SQ per mini contig, same order make_mini_genome.py
    # wrote them in (== ALL_WINDOWS order), so it matches mini_genome.fa.fai.
    new_header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": w.mini_name, "LN": w.seq} for w in ALL_WINDOWS],
    }
    mini_index = {w.mini_name: i for i, w in enumerate(ALL_WINDOWS)}

    kept = 0
    dropped_no_window = 0
    dropped_edge = 0
    dropped_unmapped = 0

    with pysam.AlignmentFile(in_bam, "rb") as fin:
        with pysam.AlignmentFile(out_bam, "wb", header=new_header) as fout:
            for read in fin.fetch(until_eof=True):
                if read.is_unmapped:
                    dropped_unmapped += 1
                    continue

                chrom = read.reference_name
                start1 = read.reference_start + 1  # to 1-based hg38 coord

                target = None
                for w in windows_by_chrom.get(chrom, []):
                    if w.hg38_start <= start1 <= w.hg38_end:
                        target = w
                        break
                if target is None:
                    dropped_no_window += 1
                    continue

                offset = target.hg38_start - 1  # 0-based shift
                new_ref_start = read.reference_start - offset

                aln_len = read.reference_length  # None if no CIGAR
                if new_ref_start < 0 or aln_len is None or \
                   new_ref_start + aln_len > target.seq:
                    dropped_edge += 1
                    continue

                new_read = pysam.AlignedSegment(fout.header)
                new_read.query_name = read.query_name
                new_read.query_sequence = read.query_sequence
                new_read.flag = read.flag & ~0x1 & ~0x8  # clear paired/mate-unmapped bits; ONT reads are single-end
                new_read.reference_id = mini_index[target.mini_name]
                new_read.reference_start = new_ref_start
                new_read.mapping_quality = read.mapping_quality
                new_read.cigartuples = read.cigartuples
                new_read.next_reference_id = -1
                new_read.next_reference_start = -1
                new_read.template_length = 0
                new_read.query_qualities = read.query_qualities
                new_read.set_tags(read.get_tags())

                fout.write(new_read)
                kept += 1

    print(
        f"kept {kept}, dropped {dropped_no_window} (outside any window), "
        f"dropped {dropped_edge} (would hang off contig edge), "
        f"dropped {dropped_unmapped} (unmapped)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()