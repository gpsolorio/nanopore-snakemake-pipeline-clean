# Building mini-test data from giab_2023.05 (direct region-streaming — fast path)

> **Why this doc exists separately from `build-mini-test-data-giab.md`:**
> that doc builds from `giab_2025.01` (HG001-HG007, CRAM+CRAI) and had to
> use local-download-then-extract because direct remote region streaming
> reliably fails against that bucket's CRAMs (`Container header CRC32
> failure`). **`giab_2023.05` is a different bucket with a different
> format — per-sample merged BAM+BAI instead of per-flowcell CRAM+CRAI —
> and direct region streaming against it was tested and confirmed
> working** (34.5s, correct MM/ML + haplotag presence, on `hg001_sup_60x`).
> That makes this the faster route, but it only covers **HG001-HG004**,
> not all seven samples the other doc covers.
>
> **Don't run both builds into the same fixture directory or barcode
> names.** This doc uses `{sample}_AR` output under
> `resources/mini_ar_test/samples_public_2023/`, distinct from the
> `giab_2025.01`-sourced `resources/mini_ar_test/samples_public/` — same
> barcode naming convention would otherwise silently overwrite fixtures
> built from different source data (this bucket's `sup_60x` is
> 60x-downsampled; `giab_2025.01`'s CRAMs are full per-flowcell depth).

## Why HG001-HG004 (and not more)

- This bucket's `variant_calling/` directory only has HG001-HG004 —
  confirmed via bucket listing, no HG005/HG006/HG007 here (those only
  exist in `giab_2025.01`). If you need all seven samples, you need both
  buckets, not just this one.
- **HG004 has no `sup_60x` tier — confirmed.** Only `hac_all` and
  `sup_all` exist for HG004; it uses `sup_all` as the fallback (see
  step 1/2), while HG001-HG003 use `sup_60x`. This means HG004's fixture
  is built from full (not 60x-downsampled) depth — a real, intentional
  asymmetry in this build, not an oversight.
- **HG001 and HG004 are female (XX)**, HG002 and HG003 are male (XY) —
  same as noted in the other build doc; still relevant here since it's
  the same underlying samples.
- `sup_60x` (SUP-accuracy basecalls, downsampled to 60x depth) is the
  smallest tier that still matches the SUP model class
  (`r1041_e82_400bps_sup_v410`) used elsewhere in this repo's configs —
  prefer it over `hac_*` or `*_all` unless there's a specific reason to
  need full depth or HAC-tier data.
- Every file here already has per-read MM/ML tags and haplotag (`HP`)
  info in one place, same as `giab_2025.01`'s CRAMs — confirmed directly
  in the streaming test, not just assumed from file naming.

## 0. One-time setup

Same repo-root assumption and directory layout as the other build doc.
If you've already done this for `giab_2025.01`, nothing new to set up:

```bash
cd /path/to/nanopore-snakemake-pipeline   # repo root

mkdir -p results/public-mini-test-2023/bam
```

No reference download needed for this step — unlike CRAM, BAM doesn't
need an external reference to decode. (The reference is still needed
later, in step 2, for `realign_to_mini.sh`'s fastq-conversion +
minimap2 remapping — that part is unchanged from the other doc.)

## 1. Confirm the tier exists for every sample before building

`hg004_sup_60x` does **not** exist — confirmed, not just suspected. HG004
only has `hac_all` and `sup_all` tiers, so it needs the `sup_all` fallback
while the other three use `sup_60x`. Check per-sample tier like this:

```bash
for SAMPLE in hg001 hg002 hg003 hg004; do
    if [[ "$SAMPLE" == "hg004" ]]; then
        RUN="sup_all"
    else
        RUN="sup_60x"
    fi
    echo "== ${SAMPLE}_${RUN} =="
    aws s3 ls --no-sign-request \
        "s3://ont-open-data/giab_2023.05/analysis/variant_calling/${SAMPLE}_${RUN}/" \
        | grep "haplotagged.bam$"
done
```

If a future bucket refresh adds `hg004_sup_60x`, this loop will still work
unmodified for the other three samples — just drop the `hg004` special
case once confirmed present.

## 2. Direct region-streaming (confirmed working — no full download)

```bash
REGION="chrX:67544031-67730619"

build_one_2023 () {
    local SAMPLE_LOWER=$1   # e.g. hg001
    local TIER=$2           # e.g. sup_60x -- or sup_all as fallback for HG004 if needed

    local BAM_URL="https://ont-open-data.s3.eu-west-1.amazonaws.com/giab_2023.05/analysis/variant_calling/${SAMPLE_LOWER}_${TIER}/${SAMPLE_LOWER}.haplotagged.bam"
    local OUT_BAM="results/public-mini-test-2023/bam/${SAMPLE_LOWER}_AR_sorted.bam"

    echo "== ${SAMPLE_LOWER} (${TIER}) =="
    time samtools view -b "$BAM_URL" "$REGION" -o "$OUT_BAM"
    samtools sort -o "${OUT_BAM}.sorted" "$OUT_BAM" && mv "${OUT_BAM}.sorted" "$OUT_BAM"
    samtools index "$OUT_BAM"

    echo "-- sanity check --"
    samtools flagstat "$OUT_BAM"
    samtools view "$OUT_BAM" | head -1 | grep -qo 'MM:Z:' && echo "MM/ML tags: present" || echo "MM/ML tags: MISSING -- do not use this file"
    samtools view "$OUT_BAM" | head -1 | grep -qo 'HP:i:' && echo "haplotag: present" || echo "haplotag: MISSING -- do not use this file"
}

build_one_2023 hg001 sup_60x
build_one_2023 hg002 sup_60x
build_one_2023 hg003 sup_60x
build_one_2023 hg004 sup_all   # hg004 has no sup_60x tier, sup_all is the only SUP option
```

No temp-file cleanup step here, unlike the `giab_2025.01` doc's
`build_one` — there's no ~50-300GB intermediate to delete, since the
region-only fetch is the whole download.

**Extra sort step, why it's here:** `samtools view -b URL region` on a
remote BAM returns reads in the order they appear in the queried blocks,
which is normally already coordinate-sorted for a single contiguous
region — but re-sorting locally is cheap at this size (a few thousand
reads) and removes any doubt before this feeds into `realign_to_mini.sh`,
which expects sorted input.

## 3. Feed into the existing mini-build pipeline

Same script as the other doc, pointed at the new output root:

```bash
set -e

for SAMPLE in hg001 hg002 hg003 hg004; do
    echo "== realigning ${SAMPLE} =="
    MINI_OUT_ROOT="$(pwd)/results/public-mini-test-2023" \
        bash workflow/scripts/mini_genomes/ar_test/realign_to_mini.sh "${SAMPLE}_AR"

    OUT_BAM="results/public-mini-test-2023/aligned/${SAMPLE}_AR/minimap2.bam"
    if [[ ! -s "$OUT_BAM" ]]; then
        echo "ERROR: ${SAMPLE} realign did not produce ${OUT_BAM}" >&2
        exit 1
    fi
done

echo "All samples realigned successfully:"
ls -lh results/public-mini-test-2023/aligned/*/minimap2.bam
```

## 4. Promote into a separate tracked fixture directory

**Uses `samples_public_2023/`, not `samples_public/`** — keeps this build
from colliding with (or silently replacing) any `giab_2025.01`-sourced
fixtures that already exist under the same barcode names:

```bash
set -e

for SAMPLE in hg001 hg002 hg003 hg004; do
    SRC_BAM="results/public-mini-test-2023/aligned/${SAMPLE}_AR/minimap2.bam"
    if [[ ! -s "$SRC_BAM" ]]; then
        echo "ERROR: missing ${SRC_BAM} -- rerun step 3 for ${SAMPLE} first" >&2
        exit 1
    fi

    mkdir -p "resources/mini_ar_test/samples_public_2023/${SAMPLE}_AR"
    cp "$SRC_BAM" "resources/mini_ar_test/samples_public_2023/${SAMPLE}_AR/"
    cp "${SRC_BAM}.bai" "resources/mini_ar_test/samples_public_2023/${SAMPLE}_AR/"
done

git lfs track "resources/mini_ar_test/samples_public_2023/**/*.bam"
git lfs track "resources/mini_ar_test/samples_public_2023/**/*.bam.bai"
git add resources/mini_ar_test/samples_public_2023/ .gitattributes
git commit -m "Add public GIAB mini-test fixtures (giab_2023.05, HG001-HG004, direct-stream method)"
git push origin main
```

Verify before trusting the push:

```bash
git lfs ls-files   # every {sample}_AR path under samples_public_2023 should list, marked with *
file resources/mini_ar_test/samples_public_2023/hg001_AR/minimap2.bam   # should report BGZF, not ASCII text
```

## 5. Config

Use `configs/mini_data_public_2023.yml` (see companion file) rather than
editing the existing `configs/mini_data_public.yml` in place — that file
is wired to the `giab_2025.01`/the private cohort-analog 3-sample build under
`samples_public/`. Keeping them separate means you can smoke-test either
build independently without one config accidentally pointing at the
other's fixture directory.

## Caveats

- **Only HG001-HG004 available in this bucket** — for HG005-HG007 you
  need `giab_2025.01` and the other build doc, full stop.
- **HG004 uses `sup_all`, not `sup_60x`** (confirmed absent for HG004) —
  HG001-HG003 are 60x-downsampled, HG004 is full depth. Don't treat
  AR-locus coverage as uniform across all four samples in this build; it
  isn't, by necessity rather than choice.
- **Coverage depth also isn't comparable across buckets** — this build's
  depth (60x for three samples, full for HG004) differs from both the
  the private cohort adaptive-sampling data and `giab_2025.01`'s full-depth CRAMs.
- **Direct region-streaming being confirmed working here is bucket/format
  specific** — don't assume it'll work against `giab_2025.01`'s CRAMs
  just because it worked here; that's a different code path (BAM+BAI vs
  CRAM+CRAI) and was independently tested and found to fail there.
- **License is CC BY-NC 4.0** — non-commercial use only, same as the
  other buckets.
- **Truth-set VCFs available but not used by this build**:
  `giab_2023.05/analysis/benchmarking/truthset/HG00{1,2,3,4}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`
  (+ `_noinconsistent.bed` confident regions) are real GIAB gold-standard
  calls, not pipeline output — worth pulling separately (they're small,
  no region-extraction needed) if you ever want to validate what
  `call_variants` produces at the AR locus against a real truth set,
  rather than just checking the pipeline runs without crashing.