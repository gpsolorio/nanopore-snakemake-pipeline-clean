##############################################################################
# nanopore_adaptive_pipeline.smk
#
# Entry point for the Nanopore adaptive sequencing Snakemake pipeline.
#
# Execution modes (experiment.genome in config):
#   mini  →  mini_genome.fa + stub BEDs  — fast local dry-run
#   hg38  →  hg38.fa + full annotation BEDs — full cluster analysis
#
# Target presets (targets: in config):
#   full           →  complete pipeline through TensorQTL (default)
#   modkit_pileup  →  stop after per-sample methylation pileup
#   plots          →  methylartist locus plots + CpG histograms only
#
# Usage:
#   # Full run (cluster) — set genome: hg38 in config
#   snakemake --configfile configs/full_data/TEMPLATE_full_data.yml
#
#   # Mini dry-run (local) — set genome: mini, targets: [modkit_pileup] in config
#   snakemake --configfile configs/mini_data/mini_data_public_2023_mac_local.yml --dry-run
##############################################################################


# ============================================================================
# Shared globals — MUST be included first.
#
# Defines: BARCODES, GENOME_MODE, REFERENCE_GENOME, OUTPUT_BASE_DIR,
#          sample(), cohort(), sample_log(), cohort_log(), res(), ...
# ============================================================================

include: "common.smk"


# ============================================================================
# Mini-data AR-locus extraction / aligner comparison (opt-in side pipeline)
#
# Only pulled in when the config explicitly defines a `mini_data:` section
# — keeps the main full/plots targets unaffected for configs that don't
# use this (e.g. configs/full_data/TEMPLATE_full_data.yml).
# ============================================================================

if "mini_data" in config:
    include: "rules/mini_genomes/mini_data.smk"


# ============================================================================
# Target resolution
#
# Controls which outputs rule all requires. Use targets: [modkit_pileup]
# in config for mini dry-runs to avoid planning the full DMR + TensorQTL
# chain, which requires a real cohort to produce meaningful outputs.
# ============================================================================

_TARGETS = config.get("targets", ["full"])


def _all_inputs():
    full_outputs = [
        cohort("dmrs", "smoothed", "DMR_individual_smoothed_methylation.bed"),
        cohort("dmrs", "smoothed", "DMR_individual_smoothed_methylation.bed.gz"),
        cohort("dmrs", "smoothed", "DMR_individual_smoothed_methylation.bed.gz.tbi"),
        cohort("methylqtl", "tensorqtl", "DMR_cisQTL.cis_qtl_results.txt"),
        cohort("methylqtl", "tensorqtl", "DMR_cisQTL.significant_dmr_pairs.txt"),
        cohort("dmrs", "wgcna", "DMR_Modules.rds"),
        cohort("dmrs", "dmr_annotations"),
        cohort("dmrs", "enrichr"),
        cohort("summaries", "dmr_heatmap", "heatmap.pdf"),
        cohort("summaries", "xci_population_summary.xlsx"),
        cohort("summaries", "xci_population_plots"),
    ]

    out = []

    if "full" in _TARGETS:
        out += full_outputs

    if "modkit_pileup" in _TARGETS:
        out += expand(
            sample("{barcode}", "methylation", "{barcode}_modkit.bed"),
            barcode=BARCODES,
        )

    if "plots" in _TARGETS:
        out += expand(
            sample("{barcode}", "histograms"),
            barcode=BARCODES,
        )
        out += [cohort("summaries", "methylartist")]

    return out or full_outputs


rule all:
    input:
        _all_inputs()


# ============================================================================
# Per-sample processing
#
# Chain per barcode:
#   pod5 → basecall_pod5 → merge_bams (merge+sort+index, streamed -- see
#   rules/alignment/merge_bams.smk) → call_variants
# ============================================================================

include: "rules/basecalling/basecall_pod5.smk"  # Dorado basecalling + minimap2 alignment
include: "rules/alignment/merge_bams.smk"       # Merge per-chunk BAMs into one sorted+indexed BAM per sample (streamed merge|sort, no unsorted intermediate on disk)
include: "rules/variants/call_variants.smk"     # Clair3 SNV/indel calling + WhatsHap phasing


# ============================================================================
# Methylation analysis
# ============================================================================

include: "rules/methylation/modkit_pileup.smk"  # Per-CpG 5mC pileup (modkit)
include: "rules/plots/methylartist.smk"         # Locus-level methylation plots + reheader
include: "rules/plots/plot_cpg_hists.smk"       # Per-CpG-island methylation histograms
include: "rules/plots/aggregate_xci_stats.smk"


# ============================================================================
# DMR calling  (dss recommended for small cohorts; dmrseq for larger ones)
# ============================================================================

DMR_CALLER = config.get("dmr", {}).get("caller", "dss")

if DMR_CALLER == "dss":
    include: "rules/dmr/dss_call_dmrs.smk"
elif DMR_CALLER == "dmrseq":
    include: "rules/dmr/dmrseq_call_dmrs.smk"
else:
    raise ValueError(
        f"dmr.caller must be 'dss' or 'dmrseq', got '{DMR_CALLER}'. "
        f"Check the dmr: section of your config."
    )

include: "rules/dmr/dss_smooth_betas.smk"       # Per-sample beta smoothing for TensorQTL


# WGCNA module detection, context annotation, and heatmap visualization --
# all cohort-level, downstream of whichever DMR caller ran above.
# Order matters: wgcna_modules.smk defines _DMRS_CSV, which
# plot_dmr_heatmap.smk reuses.
include: "rules/annotations/wgcna_modules.smk"
include: "rules/annotations/annotate_dmr_context.smk"
include: "rules/annotations/pathway_enrichment.smk"
include: "rules/plots/plot_dmr_heatmap.smk"


# ============================================================================
# Methylation QTL (TensorQTL)
#
# Chain: merge_vcfs → reheader_vcf → vcf_to_pgen → compute_pca → tensorqtl
# ============================================================================

include: "rules/methylqtl/merge_vcfs.smk"
include: "rules/methylqtl/reheader_vcf.smk"
include: "rules/methylqtl/vcf_to_pgen.smk"
include: "rules/methylqtl/compute_pca_covariates.smk"
include: "rules/methylqtl/run_tensorqtl.smk"
