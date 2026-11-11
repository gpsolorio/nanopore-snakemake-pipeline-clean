# Building mini-test data from public GIAB/ONT data (no DUA required)

The pipeline's `mini_ar_test` smoke-test data was originally built from real
the private cohort patient samples (see `workflow/scripts/mini_genomes/ar_test/`).
That data can only be used on machines covered by the the private cohort data use
agreement.

This doc walks through building an equivalent, fully public mini-test dataset
from Oxford Nanopore's own open reference data (GIAB HG002/HG003/HG004), so
the pipeline can be smoke-tested on any machine with no data agreement
required.

Source: [ont-open-data](https://registry.opendata.aws/ont-open-data/)
(license: CC BY-NC 4.0 — non-commercial use).

All commands below assume you're running from the **repo root**
(`nanopore-snakemake-pipeline/`). Everything is written into
`resources/mini_ar_test/_public_build_cache/` and `results/public-mini-test/`
— nothing goes under `$HOME`. `_public_build_cache/` holds large,
freely-redownloadable intermediate files (the ~3.1GB GRCh38 reference,
per-sample subset BAMs) and is gitignored — it's a scratch cache, not
something to commit.

## Why HG002 / HG003 / HG004

- These are the GIAB trio (Coriell NA24385/NA24149/NA24143). HG002 is male
  (XY) — same relevance for the X-linked AR locus as our actual subjects.
- `giab_lsk114_2022.12/analysis/wf-human-var-output/{sample}_sup_v4/` was
  basecalled with SUP accuracy, matching the model class
  (`r1041_e82_400bps_sup_v410`) used by `tools.clair3.model_path` in our
  own configs.
- The CRAM in each folder retains per-read MM/ML modification tags, which
  is what `plot_methylation` / `methylartist` and `call_variants`
  (`--enable_phasing`) both require.

## 0. One-time setup

```bash
cd /path/to/nanopore-snakemake-pipeline   # repo root

mkdir -p resources/mini_ar_test/_public_build_cache/ref
mkdir -p resources/mini_ar_test/_public_build_cache/bam
mkdir -p results/public-mini-test/bam
```

Confirm `resources/mini_ar_test/_public_build_cache/` is gitignored before
downloading anything into it:

```bash
grep -q "_public_build_cache" .gitignore || \
    echo "resources/mini_ar_test/_public_build_cache/" >> .gitignore
```

Install the AWS CLI (no AWS account needed — the bucket allows anonymous
access via `--no-sign-request`):

```bash
brew install awscli          # macOS
# or: pip install awscli
```

Install `git-lfs` now too, even though it's only needed in step 5 — you'll
need it eventually to promote/commit the fixtures, and installing it via a
static binary rather than a conda-env package means it's on `PATH`
regardless of which env is active later (a conda-env-scoped install, e.g.
`mamba install -n some-env git-lfs`, silently stops working the moment you
switch envs, which is a confusing failure to hit mid-commit):

```bash
# macOS
brew install git-lfs

# Linux without sudo — static binary into ~/.local/bin (works in any env)
mkdir -p ~/.local/bin
curl -fSL https://github.com/git-lfs/git-lfs/releases/download/v3.5.1/git-lfs-linux-amd64-v3.5.1.tar.gz -o /tmp/git-lfs.tar.gz
tar xzf /tmp/git-lfs.tar.gz -C /tmp
mv /tmp/git-lfs-3.5.1/git-lfs ~/.local/bin/
chmod +x ~/.local/bin/git-lfs

git lfs version   # confirm it resolves before moving on
git lfs install   # user-level, no sudo needed
```

## 1. Browse the bucket (for reference — already done, paths below are confirmed)

```bash
aws s3 ls --no-sign-request s3://ont-open-data/giab_lsk114_2022.12/
#   PRE analysis/
#   PRE flowcells/

aws s3 ls --no-sign-request s3://ont-open-data/giab_lsk114_2022.12/analysis/wf-human-var-output/
#   PRE hg002_sup_v4/   <- SUP model, matches our own config
#   PRE hg003_sup_v4/
#   PRE hg004_sup_v4/
#   ... *_hac_v4/ also available, not used here

aws s3 ls --no-sign-request s3://ont-open-data/giab_lsk114_2022.12/analysis/wf-human-var-output/hg002_sup_v4/
```

Relevant files found (same layout for hg003/hg004, substitute the sample name):

| File                                                   | Size    | Notes                                                                                                                               |
| ------------------------------------------------------ | ------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `{sample}.pass.cram`                                  | ~120 GB | Aligned, haplotagged, per-read MM/ML tags. **Never download in full** — stream just the region we need (step 3).             |
| `{sample}.pass.cram.crai`                             | ~1 MB   | Index — required for region-based streaming.                                                                                       |
| `GCA_000001405.15_GRCh38_no_alt_analysis_set.fa`     | ~3.1 GB | Reference used to align the CRAM. CRAM (unlike BAM) requires the *exact* reference to decode — download once, shared across all three samples. |
| `GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai` | 8 KB    | Index for the reference above.                                                                                                      |
| `{sample}.methyl.cpg.bed.gz`                          | ~800+MB | Aggregated site-level methylation (bedMethyl). **Not used here** — we need per-read tags from the CRAM, not aggregated calls. |

## 2. Download the reference (one-time, ~3.1 GB, shared across all samples)

```bash
REF_DIR="resources/mini_ar_test/_public_build_cache/ref"

aws s3 cp --no-sign-request \
    s3://ont-open-data/giab_lsk114_2022.12/analysis/wf-human-var-output/hg002_sup_v4/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    "${REF_DIR}/"
aws s3 cp --no-sign-request \
    s3://ont-open-data/giab_lsk114_2022.12/analysis/wf-human-var-output/hg002_sup_v4/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai \
    "${REF_DIR}/"
```

## 3. Stream the AR-locus region from each remote CRAM, straight to its final destination

`samtools` issues HTTP range requests against the public S3 URL directly —
this pulls only the bytes covering our region of interest, not the full
~120GB file per sample. Streams directly into `results/public-mini-test/bam/`
using the `{barcode}_AR_sorted.bam` naming `realign_to_mini.sh` expects, so
there's no separate staging copy to manage.

Uses the same window as the the private cohort-derived mini test
(`chrX:67544031-67730619`, matching `make_mini_genome.sh`).

```bash
REF="$(pwd)/resources/mini_ar_test/_public_build_cache/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
REGION="chrX:67544031-67730619"

for SAMPLE in hg002 hg003 hg004; do
    CRAM_URL="https://ont-open-data.s3.eu-west-1.amazonaws.com/giab_lsk114_2022.12/analysis/wf-human-var-output/${SAMPLE}_sup_v4/${SAMPLE}.pass.cram"
    OUT_BAM="results/public-mini-test/bam/${SAMPLE}_AR_sorted.bam"

    echo "== ${SAMPLE} =="
    samtools view -b -T "$REF" "$CRAM_URL" "$REGION" > "$OUT_BAM"
    samtools index "$OUT_BAM"
    samtools flagstat "$OUT_BAM"
done
```

> **Note:** this GRCh38 reference uses NCBI-style contig naming
> (`GCA_000001405.15_...`). Confirm the chromosome name for X matches
> what's in each CRAM header before running the above — check with:
> `samtools view -H "$CRAM_URL" | grep -i "SN:chrX\|SN:X"` and adjust
> `REGION` accordingly if it's named `X` instead of `chrX`. Do this once;
> all three samples in this trio share the same alignment/reference setup.

## 4. Feed into the existing mini-build pipeline unchanged

Each subset BAM from step 3 is already sitting exactly where
`realign_to_mini.sh` expects it (`results/public-mini-test/bam/${LABEL}_sorted.bam`)
— it converts to fastq and realigns to the tracked mini genome:

```bash
set -e   # stop immediately if any sample's realign step fails, rather than
         # silently continuing to the next sample with a partial result set

for SAMPLE in hg002 hg003 hg004; do
    echo "== realigning ${SAMPLE} =="
    MINI_OUT_ROOT="$(pwd)/results/public-mini-test" \
        bash workflow/scripts/mini_genomes/ar_test/realign_to_mini.sh "${SAMPLE}_AR"

    OUT_BAM="results/public-mini-test/aligned/${SAMPLE}_AR/minimap2.bam"
    if [[ ! -s "$OUT_BAM" ]]; then
        echo "ERROR: ${SAMPLE} realign did not produce ${OUT_BAM}" >&2
        exit 1
    fi
done

echo "All 3 samples realigned successfully:"
ls -lh results/public-mini-test/aligned/*/minimap2.bam
```

This produces `results/public-mini-test/aligned/{sample}_AR/minimap2.bam` for
each sample, directly analogous to the the private cohort-derived
`results/mini_test_data/aligned/<barcode>_AR/minimap2.bam` outputs, but
built entirely from public data — usable on any machine, no DUA required.

## 5. Promote into the tracked fixture directory (if replacing/updating the committed fixtures)

If this run is meant to refresh the fixtures already checked into
`resources/mini_ar_test/samples_public/` (tracked via Git LFS — see repo
history), copy the final aligned BAMs there. Same fail-fast pattern as
step 4 — this confirms every source file actually exists before touching
the tracked fixture directory, rather than leaving it partially updated if
one sample's `cp` fails partway through:

```bash
set -e

for SAMPLE in hg002 hg003 hg004; do
    SRC_BAM="results/public-mini-test/aligned/${SAMPLE}_AR/minimap2.bam"
    if [[ ! -s "$SRC_BAM" ]]; then
        echo "ERROR: missing ${SRC_BAM} -- rerun step 4 for ${SAMPLE} first" >&2
        exit 1
    fi

    mkdir -p "resources/mini_ar_test/samples_public/${SAMPLE}_AR"
    cp "$SRC_BAM" "resources/mini_ar_test/samples_public/${SAMPLE}_AR/"
    cp "${SRC_BAM}.bai" "resources/mini_ar_test/samples_public/${SAMPLE}_AR/"
done

git lfs track "resources/mini_ar_test/samples_public/**/*.bam"
git lfs track "resources/mini_ar_test/samples_public/**/*.bam.bai"
git add resources/mini_ar_test/samples_public/ .gitattributes
git commit -m "Refresh public GIAB trio mini-test fixtures"
git push origin main
```

Verify before trusting the push — confirm these are real LFS-tracked binary
content, not pointer files that failed to smudge:

```bash
git lfs ls-files   # all six paths (hg002/3/4 x .bam/.bam.bai) should list, marked with *
file resources/mini_ar_test/samples_public/hg004_AR/minimap2.bam   # should report BGZF, not ASCII text
```

## 6. Config

`configs/mini_data_public.yml` already points at these `{sample}_AR` barcodes
the same way `configs/mini_data.yml` points at the the private cohort barcodes, with
`experiment.use_prebuilt_alignments: true` so `call_variants` pulls directly
from `resources/mini_ar_test/samples_public/` instead of re-running
alignment. No changes needed here if you're just refreshing existing
fixtures with this doc's steps.

## Caveats

- **Coverage depth may differ** from the the private cohort adaptive-sampling data —
  these CRAMs are whole-genome sequencing, not locus-targeted adaptive
  sampling, so AR-locus depth here reflects standard WGS coverage rather
  than the enriched depth adaptive sampling produces. Fine for smoke-testing
  pipeline mechanics; not a substitute for validating adaptive-sampling-
  specific behavior.
- **License is CC BY-NC 4.0** — non-commercial use only. Fine for internal
  pipeline testing/development; check before any commercial or for-profit
  use case.
- **`_public_build_cache/` is a scratch cache, not a deliverable.** Safe to
  `rm -rf resources/mini_ar_test/_public_build_cache/` at any time — rerun
  steps 2-3 to rebuild it. Never commit its contents.