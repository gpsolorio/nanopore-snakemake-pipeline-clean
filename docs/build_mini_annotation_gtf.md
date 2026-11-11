# Building the Mini-Genome Annotation GTF (AR locus)

This documents how `analysis.gtf_file` for the `mini` genome-mode configs
(e.g. `mini_data_public_2023_hpc_a100.yml`) was built, and how to rebuild
it or extend it to a new locus. Written for the `annotate_dmr_context`
Snakemake rule, which needs a GTF with real `transcript`-type rows and a
`gene_name` attribute — most GTFs on the cluster don't satisfy both, so
this isn't a one-line fix.

## Background: why this isn't just "point at a GTF"

`mini_genome.fa` is a synthetic single-contig reference (contig name
`AR_locus`) built from a real hg38 region, not real hg38 itself. Every
downstream file — modkit BEDs, DMR results, this annotation GTF — has to
agree on that same `AR_locus` contig name and its own 1-based coordinate
system (starting at 1, not the real hg38 position). A GTF subset directly
from an hg38-coordinate file will not overlap anything, silently, because
`findOverlaps()` just returns zero matches rather than erroring.

Separately, not every "GTF" file has what this script needs:
- **`gtf_data[gtf_data$type == "transcript"]`** — some locally-built GTFs
  only carry `CDS`/`exon`/`start_codon`/`stop_codon` rows and no
  `transcript` rows at all. Check for this before trusting a source file.
- **`mcols(gtf_transcripts)$gene_name`** — some GTFs use `gene_id` or a
  bare `gene` attribute instead of `gene_name`. NCBI RefSeq GTFs in
  particular use accession-style contig names (`NC_000023.11`, not
  `chrX`) *and* a bare `gene` attribute — both wrong for this script.

`refdata-gex-GRCh38-2024-A/genes/genes.gtf` (the 10x Cell Ranger reference,
GENCODE-derived) has both. That's the source used below.

## Step 1 — Confirm the mini-genome's source window

The mini genome's contig was built from a specific hg38 region. For the
AR locus, that region is defined in `configs/old_configs/mini_data.yml`:

```bash
grep -A3 "mini_data:" configs/old_configs/mini_data.yml
```

```yaml
mini_data:
  region: "chrX:67544031-67730619"
```

Confirm this window's length matches the mini FASTA's actual sequence
length (they should match, modulo any padding):

```bash
awk '/^>/{if(seq)print length(seq); seq=""} !/^>/{seq=seq$0} END{print length(seq)}' \
  resources/mini_ar_test/genome/mini_genome.fa
# -> 186588
# 67730619 - 67544031 = 186588  ✓ matches (approximately -- see the
# resolved 1bp note near the end of Step 3: this check compares an
# exclusive difference against an inclusive length and happens to land
# on the same number, which masked a real 1bp start-side error in the
# declared window. Don't treat this check alone as sufficient proof.)
```

Confirm the contig name used throughout the pipeline:

```bash
head -1 resources/mini_ar_test/genome/mini_genome.fa
cut -f1 resources/mini_ar_test/beds/mini_AR_genes.bed | sort -u
cut -f1 results/<STUDY>/samples/<any_barcode>/methylation/<barcode>_modkit.bed | sort -u
# all three -> AR_locus
```

If any of these disagree, stop — the window or contig name assumed below
is wrong for your data and needs re-deriving before continuing.

## Step 2 — Find a GTF that actually has what's needed

Check candidate GTFs on the cluster for both requirements before picking
one. Do not assume a file works because it's named `*.gtf`.

```bash
# Does it have transcript-type rows overlapping the locus?
awk -F'\t' '$1=="chrX" && $4<=67730619 && $5>=67544031 {print $3}' \
  <candidate.gtf> | sort -u

# Does it have a real gene_name attribute (not just gene_id/gene)?
awk -F'\t' '$1=="chrX" && $4<=67730619 && $5>=67544031 && $3=="gene"' \
  <candidate.gtf>
```

For the AR locus, `refdata-gex-GRCh38-2024-A/genes/genes.gtf` passed both
checks:

```
chrX  HAVANA  gene  67544021  67730619  .  +  .  gene_id "ENSG00000169083"; ... gene_name "AR"; ...
```

(Coordinates here — `67544021`–`67730619` — are the gene's own span, not
the mini-genome window; they don't need to match exactly, only overlap.)

Rejected in this case:
- `/path/to/hg38_genes.gtf` — no `transcript`
  rows over this region (only `CDS`/`exon`/`start_codon`/`stop_codon`),
  and its `gene_id`/`transcript_id` attributes are RefSeq accessions
  (e.g. `NM_000044`) with no `gene_name` field at all.
- `/path/to/GCF_000001405.40_GRCh38.p14_genomic.gtf`
  (NCBI RefSeq) — uses accession-style contig names (`NC_000023.11`, not
  `chrX`) and a bare `gene` attribute, not `gene_name`.

## Step 3 — Subset and translate coordinates into mini-genome space

GTF is 1-based and fully-closed on both ends, so translation is a
constant per-row offset using the confirmed window start (`ws =
67544032` -- verified by direct sequence comparison against the mini
FASTA, see the resolved 1bp note below; the naively-declared
`67544031` from `configs/old_configs/mini_data.yml` is off by one).

```bash
awk -F'\t' -v OFS='\t' -v ws=67544032 -v we=67730619 '
  $1=="chrX" && $5>=ws && $4<=we {
    start = ($4 < ws) ? ws : $4
    end   = ($5 > we) ? we : $5
    $1 = "AR_locus"
    $4 = start - ws + 1
    $5 = end   - ws + 1
    print
  }' /path/to/genes.gtf \
  > /tmp/AR_locus_mini.gtf
```

This preserves every column and attribute (`gene_id`, `gene_name`,
`gene_type`, etc.) and only shifts the two coordinate columns.

### Confirmed: strict containment silently drops whole transcripts, not just the gene envelope

An earlier version of this command used strict containment
(`$4>=ws && $5<=we`) instead of overlap+clip. Run against the real AR
locus, that version silently dropped **the `gene`-type row entirely, and
two full transcripts (`AR-202`, `AR-207`)** — both start at the same
`67544021` boundary as the gene, 10bp before `ws=67544031`, so all three
failed `$4>=ws` and were dropped whole rather than clipped. This wasn't
cosmetic: `annotate_dmr_context.R` filters on `type == "transcript"`, so
those two transcripts would have been silently invisible to the
gene-overlap step. (Practical impact was likely small here, since several
other AR transcripts share the same exon-1 span and would still resolve
`gene_name = "AR"` for most regions — but any region unique to those two
transcripts' structure would have come back unannotated with no error.)

Always use the overlap+clip version below, and always check the
feature-type list after generating a translated GTF — a clean row/
gene_name count can still be hiding a silently dropped feature type or
transcript:

```bash
awk -F'\t' -v OFS='\t' -v ws=67544032 -v we=67730619 '
  $1=="chrX" && $5>=ws && $4<=we {
    start = ($4 < ws) ? ws : $4
    end   = ($5 > we) ? we : $5
    $1 = "AR_locus"
    $4 = start - ws + 1
    $5 = end   - ws + 1
    print
  }' /path/to/genes.gtf \
  > /tmp/AR_locus_mini.gtf

awk -F'\t' '{print $3}' /tmp/AR_locus_mini.gtf | sort -u
# -> CDS, UTR, exon, gene, start_codon, stop_codon, transcript  ✓
```

### 1bp precision gap — resolved

The declared window in `configs/old_configs/mini_data.yml`
(`chrX:67544031-67730619`) is off by one base on the **start** side.
Confirmed by direct sequence comparison: `mini_genome.fa`'s actual first
bases (`CGAGATCCCGGGG...`, right after the `>AR_locus` header) match real
hg38 `chrX:67544032-67544130` exactly — one base later than declared.
The **end** of the window (`67730619`) was already correct, confirmed the
same way against the mini FASTA's last bases.

**The true window is `chrX:67544032-67730619`** (186588bp inclusive,
exactly matching the FASTA's real length). Use `ws=67544032` in the
translation command below, not `67544031`.

## Step 4 — Verify the translated file

```bash
wc -l /tmp/AR_locus_mini.gtf
grep -c "gene_name" /tmp/AR_locus_mini.gtf     # should be > 0
awk -F'\t' '{print $3}' /tmp/AR_locus_mini.gtf | sort -u   # should include "transcript"
grep "gene_name \"AR\"" /tmp/AR_locus_mini.gtf  # should show at least one row
```

## Step 5 — Move it to a permanent location and wire it into the config

```bash
mkdir -p resources/mini_ar_test/gtf
cp /tmp/AR_locus_mini.gtf resources/mini_ar_test/gtf/AR_locus_mini.gtf
```

```yaml
analysis:
  gtf_file: resources/mini_ar_test/gtf/AR_locus_mini.gtf
  genome: hg38
```

## Known limitation — read this before expecting full annotation output

`annotate_dmr_context.R` does two independent things:

1. **Gene-name lookup by overlap with the GTF above** — works correctly
   once the GTF is properly translated, as done here.
2. **CpG/genic-context annotation via `annotatr::build_annotations(genome
   = "hg38", ...)`** — this builds annotation objects for a real, named
   genome build. There is no such build for a synthetic single-locus
   genome. Setting `analysis.genome: hg38` here only makes this call
   succeed mechanically (it builds real hg38-wide annotation objects);
   comparing those against `AR_locus`-space regions will not produce
   meaningful CpG/genic-context columns. Treat those specific output
   columns as expected-empty for any `genome: mini` config. This is a
   structural limitation of running this rule against a synthetic
   reference, not something a better GTF or more translation fixes.

## Extending this to additional loci

The gene list below (curated NDD/X-linked genes, consistent with the private cohort's
ASD focus) is a candidate for a future multi-locus mini genome, but that's
a separate, larger effort — it needs the multi-window machinery in
`workflow/scripts/mini_genomes/production/` (`window_spec.py`), not the
single-window tools used here, and isn't required to unblock the
single-AR-locus annotation rule above.

```
chrX    67544021    67730619    AR
chrX    40051251    40097958    BCOR
chrX    41514934    41923554    CASK
chrX    120523858   120560962   CUL4B
chrX    149476988   149505306   IDS
chrX    28587446    29956718    IL1RAPL1
chrX    43656299    43746817    MAOA
chrX    136216814   136251400   MAP7D3
chrX    154088132   154104842   MECP2
chrX    10445556    10620585    MID1
chr22   49043919    49052549    NHIP
chrX    48831096    48835610    PCSK1N
chrX    78104248    78129295    PGK1
chrX    23334849    23404374    PTCHD1
chrX    70925579    70931096    SLC7A3
chrX    21940709    21994837    SMS
chrX    9465058     9719740     TBL1X
```

If/when this is pursued, the same two checks from Step 2 (transcript rows
present, `gene_name` attribute present) need to be re-run per gene, since
GTF quality/coverage can vary across the genome even within one source
file.