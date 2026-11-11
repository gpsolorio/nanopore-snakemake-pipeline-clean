# Building mini-test data from public GIAB/ONT data (no DUA required)

> **Revision note:** this replaces the HG002/HG003/HG004/HG006-only draft.
> Key differences:
> - Expanded to all six samples you have flowcell listings for:
>   **HG001, HG002, HG003, HG004, HG005, HG006**.
> - Every sample in the current listing has **two flowcells**, not just
>   HG002. Don't assume any sample is single-flowcell without checking —
>   this listing already shows six-for-six with two apiece:
>
>   | Sample | Flowcell 1 | Flowcell 2 |
>   |--------|-----------|-----------|
>   | HG001  | PAW79146  | PAW81754  |
>   | HG002  | PAW70337  | PAW71238  |
>   | HG003  | PAY87794  | PAY87954  |
>   | HG004  | PAY87778  | PAY88428  |
>   | HG005  | PAW87816  | PAW88001  |
>   | HG006  | PAY77227  | PBA16846  |
>
>   This is just the top-level `PRE` listing (directory prefixes) — it does
>   **not** yet confirm `haplotagged.cram` exists under each. Run the
>   deeper `ls` from step 1 before trusting any of these.
>   Note HG006 previously showed as "already done in an earlier session"
>   using only `PAY77227` — now that `PBA16846` is visible too, that
>   earlier single-flowcell build is stale with respect to this listing.
>   Decide fresh whether HG006 stays single-flowcell or gets the merge
>   treatment, same as everything else.
> - Uses `giab_2025.01` (confirmed real bucket layout), not the older
>   `giab_lsk114_2022.12`.
> - The step-1 CRAM-confirmation loop had a copy-paste bug the first time it
>   was run: it looped every sample against **HG001's** flowcells
>   (`PAW79146`, `PAW81754`) instead of each sample's own pair, so
>   HG002-HG006 came back empty. Fixed below with a per-sample associative
>   array so the same mistake can't happen silently again.
> - CRAM files are large — confirmed for HG001: `PAW79146` is ~66.8GB,
>   `PAW81754` is ~76.9GB (**~144GB for HG001 alone across both
>   flowcells**). That's already above the doc's earlier "90-100GB per
>   flowcell" estimate — treat that number as a floor, not a ceiling, until
>   the rest of the sizes are confirmed.
> - **Direct remote region-streaming (`samtools view -b -T ref.fa CRAM_URL region`)
>   does not work reliably against this bucket in this environment** — it
>   fails with `[E::cram_read_container] Container header CRC32 failure`
>   after a partial range read, consistently, even with a modern
>   htslib build (1.23.1, full libcurl/S3/GCS support confirmed). Root cause
>   wasn't nailed down (network path is the leading suspect — proxy was ruled
>   out, but something between here and S3 is mangling multi-range HTTP
>   responses). **The working method is local-download-then-extract**,
>   documented below as the primary path, not a fallback.

## Why HG001 / HG002 / HG003 / HG004 / HG005 / HG006

- **Sex / allele complexity for the AR region (chrX):**
  - Male (XY) — clean single AR allele, no phasing ambiguity: **HG002,
    HG003, HG005, HG006**.
  - Female (XX) — two AR alleles, useful as a heterozygous test case but
    not "simple" the way the male samples are: **HG001, HG004**.
  - Decide deliberately which set you actually want per test purpose;
    don't assume "all six are simple" just because most are male.
- SUP-accuracy basecalls (`dna_r10.4.1_e8.2_400bps@v5.0.0`), matching the
  model class used elsewhere in this repo's configs.
- The haplotagged CRAM retains per-read MM/ML modification tags (needed for
  `plot_methylation`/`methylartist`) and haplotag info (needed for
  `call_variants --enable_phasing`). These are two separate pieces of
  information that happen to live in the same file — MM/ML comes from
  modkit, haplotagging comes from whatshap acting on variant calls.

## 0. One-time setup

```bash
cd /path/to/nanopore-snakemake-pipeline   # repo root

mkdir -p resources/mini_ar_test/_public_build_cache/ref
mkdir -p results/public-mini-test/bam
```

Confirm the cache dir is gitignored:

```bash
grep -q "_public_build_cache" .gitignore || \
    echo "resources/mini_ar_test/_public_build_cache/" >> .gitignore
```

Reference used by all samples (already present in this env as
`GRCh38_no_alt_analysis_set.fna` — confirm it's indexed):

```bash
samtools faidx GRCh38_no_alt_analysis_set.fna   # no-op if .fai already exists
```

## 1. Find each sample's flowcell(s) — don't assume one

```bash
for SAMPLE in HG001 HG002 HG003 HG004 HG005 HG006; do
    echo "== ${SAMPLE} =="
    aws s3 ls --no-sign-request \
        "s3://ont-open-data/giab_2025.01/analysis/wf-human-variation/sup/${SAMPLE}/"
done
```

You already ran this — see the table above. Every sample currently shows
two flowcells. For each flowcell, confirm the CRAM actually exists *and*
get its size before using it. Use a per-sample associative array rather
than a shared flowcell list for the loop — a shared list is exactly what
caused HG002-HG006 to silently come back empty the first time this was
run (the loop was checking HG001's flowcells against every sample):

```bash
declare -A FLOWCELLS=(
  [HG001]="PAW79146 PAW81754"
  [HG002]="PAW70337 PAW71238"
  [HG003]="PAY87794 PAY87954"
  [HG004]="PAY87778 PAY88428"
  [HG005]="PAW87816 PAW88001"
  [HG006]="PAY77227 PBA16846"
)

for SAMPLE in "${!FLOWCELLS[@]}"; do
    for FLOWCELL in ${FLOWCELLS[$SAMPLE]}; do
        echo "== ${SAMPLE} / ${FLOWCELL} =="
        LINE=$(aws s3 ls --no-sign-request \
            "s3://ont-open-data/giab_2025.01/analysis/wf-human-variation/sup/${SAMPLE}/${FLOWCELL}/output/" \
            | grep "haplotagged.cram$")
        echo "$LINE"
        BYTES=$(echo "$LINE" | awk '{print $3}')
        if [[ -n "$BYTES" ]]; then
            awk -v b="$BYTES" 'BEGIN{printf "  -> %.1f GB\n", b/1e9}'
        fi
    done
done
```

Confirmed so far (fill in the rest as you run it):

| Sample | Flowcell | CRAM size (GB) |
|--------|----------|-----------------|
| HG001  | PAW79146 | ~66.8 GB |
| HG001  | PAW81754 | ~76.9 GB |
| HG002  | PAW70337 | TBD |
| HG002  | PAW71238 | TBD |
| HG003  | PAY87794 | TBD |
| HG003  | PAY87954 | TBD |
| HG004  | PAY87778 | TBD |
| HG004  | PAY88428 | TBD |
| HG005  | PAW87816 | TBD |
| HG005  | PAW88001 | TBD |
| HG006  | PAY77227 | TBD |
| HG006  | PBA16846 | TBD |

**Decide once, up front, whether you're using one flowcell per sample or
merging all of a sample's flowcells** — and apply that decision consistently
across all six. Mixing "single flowcell for this one, merged for that one"
makes the fixture set inconsistent in ways that are easy to forget later.
(This doc assumes single-flowcell-per-sample for simplicity; see the end of
step 3 for the merge variant if you want fuller coverage.)

If you go single-flowcell, pick one flowcell ID per sample from the table
above now, e.g.:

```
HG001 -> PAW79146
HG002 -> PAW70337
HG003 -> PAY87794
HG004 -> PAY87778
HG005 -> PAW87816
HG006 -> PAY77227
```

(arbitrary "first listed" choice — swap in whichever you actually confirmed
has the CRAM, if the first one is missing it).

## 2. Local-download-then-extract (the method that actually works here)

For each sample/flowcell, this streams the full CRAM through a local sort
(so it never touches disk uncompressed/unsorted), then extracts the AR
region from the resulting local BAM, then deletes the large intermediate
immediately. Expect ~90-100GB of download and disk churn *per flowcell* —
budget time and bandwidth accordingly, and don't run multiple samples in
parallel unless your disk has room for more than one ~90GB intermediate at
once.

```bash
REGION="chrX:67544031-67730619"
REF="$(pwd)/GRCh38_no_alt_analysis_set.fna"

build_one () {
    local SAMPLE_LOWER=$1   # e.g. hg005
    local SAMPLE_UPPER=$2   # e.g. HG005
    local FLOWCELL=$3       # e.g. PAW87816

    local CRAM_URL="https://ont-open-data.s3.eu-west-1.amazonaws.com/giab_2025.01/analysis/wf-human-variation/sup/${SAMPLE_UPPER}/${FLOWCELL}/output/SAMPLE.haplotagged.cram"
    local TMP_FULL="/tmp/${SAMPLE_LOWER}_${FLOWCELL}_full.bam"
    local OUT_BAM="results/public-mini-test/bam/${SAMPLE_LOWER}_AR_sorted.bam"

    echo "== ${SAMPLE_UPPER} (${FLOWCELL}) =="
    samtools view -u --reference "$REF" "$CRAM_URL" \
        | samtools sort -@4 --write-index -o "$TMP_FULL" -

    samtools view -b "$TMP_FULL" "$REGION" -o "$OUT_BAM"
    samtools index "$OUT_BAM"

    echo "-- sanity check --"
    samtools flagstat "$OUT_BAM"
    samtools view -H "$OUT_BAM" | grep "^@RG"
    samtools view "$OUT_BAM" | head -1 | grep -qo 'MM:Z:' && echo "MM/ML tags: present" || echo "MM/ML tags: MISSING -- do not use this file"

    rm -f "$TMP_FULL" "${TMP_FULL}.csi"
}

# fill in real, CRAM-confirmed flowcell IDs from step 1
build_one hg001 HG001 PAW79146
build_one hg002 HG002 PAW70337
build_one hg003 HG003 PAY87794
build_one hg004 HG004 PAY87778
build_one hg005 HG005 PAW87816
build_one hg006 HG006 PAY77227   # re-verify: earlier "already done" build predates PBA16846 showing up
```

Run these one at a time, not backgrounded in parallel, unless you've
confirmed you have >200GB free — each `$TMP_FULL` is a full-genome BAM
before it gets deleted. Six samples at ~90-100GB each is on the order of
half a terabyte of transient download if you go through all of them
sequentially in one sitting — plan the session accordingly.

**Merge variant**, if a sample has multiple flowcells and you want combined
coverage rather than picking one (all six qualify now, going by the
current listing):

```bash
# after extracting each flowcell's region individually to its own file, e.g.:
#   /tmp/hg001_PAW79146_AR.bam
#   /tmp/hg001_PAW81754_AR.bam
samtools merge -f results/public-mini-test/bam/hg001_AR_sorted.bam \
    /tmp/hg001_PAW79146_AR.bam /tmp/hg001_PAW81754_AR.bam
samtools sort -o /tmp/hg001_merged_sorted.bam results/public-mini-test/bam/hg001_AR_sorted.bam
mv /tmp/hg001_merged_sorted.bam results/public-mini-test/bam/hg001_AR_sorted.bam
samtools index results/public-mini-test/bam/hg001_AR_sorted.bam
```

(repeat per sample, substituting the sample's own two flowcell-extracted
files)

## 3. Feed into the existing mini-build pipeline unchanged

Once `results/public-mini-test/bam/{sample}_AR_sorted.bam` exists for every
sample you're building:

```bash
set -e

for SAMPLE in hg001 hg002 hg003 hg004 hg005 hg006; do
    echo "== realigning ${SAMPLE} =="
    MINI_OUT_ROOT="$(pwd)/results/public-mini-test" \
        bash workflow/scripts/mini_genomes/ar_test/realign_to_mini.sh "${SAMPLE}_AR"

    OUT_BAM="results/public-mini-test/aligned/${SAMPLE}_AR/minimap2.bam"
    if [[ ! -s "$OUT_BAM" ]]; then
        echo "ERROR: ${SAMPLE} realign did not produce ${OUT_BAM}" >&2
        exit 1
    fi
done

echo "All samples realigned successfully:"
ls -lh results/public-mini-test/aligned/*/minimap2.bam
```

## 4. Promote into the tracked fixture directory

```bash
set -e

for SAMPLE in hg001 hg002 hg003 hg004 hg005 hg006; do
    SRC_BAM="results/public-mini-test/aligned/${SAMPLE}_AR/minimap2.bam"
    if [[ ! -s "$SRC_BAM" ]]; then
        echo "ERROR: missing ${SRC_BAM} -- rerun step 3 for ${SAMPLE} first" >&2
        exit 1
    fi

    mkdir -p "resources/mini_ar_test/samples_public/${SAMPLE}_AR"
    cp "$SRC_BAM" "resources/mini_ar_test/samples_public/${SAMPLE}_AR/"
    cp "${SRC_BAM}.bai" "resources/mini_ar_test/samples_public/${SAMPLE}_AR/"
done

git lfs track "resources/mini_ar_test/samples_public/**/*.bam"
git lfs track "resources/mini_ar_test/samples_public/**/*.bam.bai"
git add resources/mini_ar_test/samples_public/ .gitattributes
git commit -m "Refresh public GIAB mini-test fixtures (giab_2025.01, HG001-HG006, local-extract method)"
git push origin main
```

Verify before trusting the push:

```bash
git lfs ls-files   # every {sample}_AR path should list, marked with *
file resources/mini_ar_test/samples_public/hg004_AR/minimap2.bam   # should report BGZF, not ASCII text
```

## Caveats

- **Direct remote CRAM region streaming doesn't work here** — don't retry
  it as a "cleaner" alternative without first figuring out the actual root
  cause (suspected multi-range HTTP mangling somewhere in the network path,
  unconfirmed). Wastes ~2 seconds failing per attempt, which is cheap, but
  don't build automation that depends on it working.
  - In principle this *should* work: every CRAM here has a `.crai`
    sidecar, and htslib's normal path uses that index to issue HTTP range
    requests and only pull the blocks covering the target region — no
    full download needed. That's the thing that's failing
    (`Container header CRC32 failure` after a partial range read).
  - If it's ever worth re-investigating (e.g. bandwidth becomes a hard
    blocker across six samples × ~700GB+), these are different angles
    than what's already been ruled out (proxy, htslib version):
    - Manually compute byte offsets from the `.crai` and issue a plain
      `curl -r start-end` against the CRAM URL, bypassing htslib's own
      range logic — isolates whether it's htslib's multi-range handling
      or the network path itself.
    - Check whether the failure is specific to **multi-range** requests
      (some proxies/CDNs handle single-range fine but choke on the
      multi-part ranges CRAI-driven region queries typically issue).
    - Look for an **htsget** endpoint for this bucket, if one exists —
      that protocol is purpose-built for "give me this region from this
      CRAM" and avoids raw multi-range HTTP entirely.
  - Until/unless one of those pans out, local-download-then-extract is
    the working method — treat the above as a future optimization, not a
    blocker for building the fixtures now.
- **HG001 and HG004 are female (XX)**; HG002, HG003, HG005, HG006 are male
  (XY) — see "Why HG001 / HG002 / ... / HG006" above.
- **Every sample currently listed has two flowcells** (HG001, HG002,
  HG003, HG004, HG005, HG006 all show a `PRE .../ PRE .../` pair). Still
  worth re-running the step 1 listing before a build in case that changes,
  but as of this listing there's no single-flowcell sample among the six.
- **HG006's earlier "already done" build is stale**: it only used
  `PAY77227`; `PBA16846` wasn't visible at that time. Rebuild or explicitly
  decide to keep it single-flowcell with `PAY77227` before promoting.
- **Bandwidth/time cost is real**: ~90-100GB downloaded and immediately
  discarded per flowcell, since this bucket doesn't support the region-only
  streaming that was supposed to avoid this. Six samples at two flowcells
  each, if you go the merge route, is well over half a terabyte of
  transient download. Budget for it, and don't run samples in parallel
  unless disk has room for multiple ~90GB intermediates at once.
- **Coverage depth** reflects standard WGS, not the the private cohort data's
  adaptive-sampling enrichment — fine for pipeline smoke-testing, not a
  substitute for adaptive-sampling-specific validation.
- **License is CC BY-NC 4.0** — non-commercial use only.
- **`_public_build_cache/` and `/tmp/*_full.bam` are scratch, never
  committed.** The `build_one` function above deletes its own temp file at
  the end of each call — if a run is interrupted mid-sample, check `/tmp`
  for orphaned `*_full.bam` files before your disk fills up.