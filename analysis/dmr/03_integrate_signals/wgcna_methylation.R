##############################################################################
# wgcna_methylation.R — WGCNA Co-methylation Module Detection
##############################################################################
#
# Purpose: Identify co-methylation modules among DSS DMRs using WGCNA,
#          generating weighted network analysis of regional methylation
#          patterns.
#
# Note: this rebuilds its own BSseq object from the modkit BEDs via
# build_bsseq_from_metadata() (config/utils.R) rather than depending on a
# saved object from run_dss.R -- run_dss.R never wrote one, and this keeps
# 02_call_dmrs / 03_individual_betas / 03_integrate_signals independently
# runnable in parallel, matching how the rest of this pipeline is structured.
#
# Called by Snakemake rule wgcna_modules via:
#   Rscript wgcna_methylation.R --repo_root <path> --metadata <path> \
#       --dmr_file <path> --outdir <path> [--group_col diagnosis]
#       [--group_labels ASD,TD] [--min_cluster_size 30]
#
##############################################################################

suppressPackageStartupMessages({
  library(bsseq)
  library(DSS)
  library(data.table)
  library(GenomicRanges)
  library(WGCNA)
  library(conflicted)
  library(glue)
})

conflicted::conflicts_prefer(WGCNA::cor, .quiet = TRUE)
conflicted::conflicts_prefer(base::setdiff, .quiet = TRUE)
conflicted::conflicts_prefer(MatrixGenerics::rowVars, .quiet = TRUE)


# ============================================================================
# 1. CLI Arguments
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(args, flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}

repo_root    <- get_arg(args, "--repo_root")
metadata_tsv <- get_arg(args, "--metadata")
dmr_file     <- get_arg(args, "--dmr_file")
outdir       <- get_arg(args, "--outdir")
group_col    <- get_arg(args, "--group_col", "diagnosis")
group_labels <- get_arg(args, "--group_labels", NA_character_)
min_cluster_size_arg <- get_arg(args, "--min_cluster_size", NA_character_)
na_thresh_arg         <- get_arg(args, "--na_thresh", NA_character_)

if (is.null(repo_root))    stop("--repo_root is required")
if (is.null(metadata_tsv)) stop("--metadata is required")
if (is.null(dmr_file))     stop("--dmr_file is required")
if (is.null(outdir))       stop("--outdir is required")

config_dir <- file.path(repo_root, "analysis", "dmr", "config")
source(file.path(config_dir, "utils.R"))
source(file.path(config_dir, "params.R"))

min_cluster_size <- if (!is.na(min_cluster_size_arg)) as.integer(min_cluster_size_arg) else WGCNA_PARAMS$min_module_size
na_thresh         <- if (!is.na(na_thresh_arg)) as.numeric(na_thresh_arg) else 0.5

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

print_config("wgcna_modules", list(
  metadata         = metadata_tsv,
  dmr_file         = dmr_file,
  outdir           = outdir,
  group_col        = group_col,
  min_cluster_size = min_cluster_size,
  na_thresh        = na_thresh
))


# ============================================================================
# 2. Load sample metadata + DMRs
# ============================================================================

meta <- read.table(metadata_tsv, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, colClasses = "character")

if (!group_col %in% colnames(meta)) {
  stop(glue("Metadata TSV is missing group column: {group_col}"))
}

if (!is.na(group_labels)) {
  wanted <- trimws(strsplit(group_labels, ",")[[1]])
  meta   <- meta[meta[[group_col]] %in% wanted, ]
} else {
  meta <- meta[!is.na(meta[[group_col]]) & meta[[group_col]] != "", ]
}

if (nrow(meta) == 0) {
  stop("No samples found in metadata TSV after group filtering.")
}

ds <- list(
  sample_names = meta$sample_id,
  groups       = meta[[group_col]],
  bed_paths    = meta$bed_path
)

message("Step 1: Loading DMRs...")
dmrs <- read.csv(dmr_file, header = TRUE)
message("  Loaded ", nrow(dmrs), " DMRs")


# ============================================================================
# 3. Build BSseq object (shared helper -- no longer expects a pre-saved RDS)
# ============================================================================

message("Step 2: Building BSseq object from modkit BEDs...")

BSobj <- build_bsseq_from_metadata(
  bed_paths       = ds$bed_paths,
  sample_names    = ds$sample_names,
  groups          = ds$groups,
  coverage_params = COVERAGE,
  metadata        = meta
)

message("  Loaded BSseq object with ", ncol(BSobj), " samples")


# ============================================================================
# 4. Create Methylation Matrix
# ============================================================================

message("Step 3: Creating methylation matrix...")

dmr_gr <- GRanges(
  seqnames = dmrs$chr,
  ranges   = IRanges(start = dmrs$start, end = dmrs$end)
)

meth_mat <- getMeth(BSobj, regions = dmr_gr, type = "raw", what = "perRegion")

region_ids <- paste0("region_", seq_len(nrow(meth_mat)))
rownames(meth_mat) <- region_ids

region_metadata <- as.data.frame(dmr_gr)
region_metadata$RegionID <- region_ids

message("  Methylation matrix: ", nrow(meth_mat), " regions x ", ncol(meth_mat), " samples")


# ============================================================================
# 5. Filter Regions
# ============================================================================

message("Step 4: Filtering regions...")

filter_regions <- function(mat, na_thresh) {
  mat <- mat[rowSums(is.na(mat)) / ncol(mat) <= na_thresh, , drop = FALSE]
  mat <- mat[apply(mat, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]
  mat
}

meth_mat_clean <- filter_regions(meth_mat, na_thresh)
message("  Kept ", nrow(meth_mat_clean), " regions after filtering")

# ----------------------------------------------------------------------------
# Degenerate-input guards, added after a crash on a 1-region mini/GIAB smoke
# test: goodSamplesGenes() internally does
#   colSums(!is.na(datExpr[useSamples, useGenes]))
# and matrix `[` subsetting drops to a bare vector whenever exactly one row
# or column survives (R's default drop = TRUE) -- colSums() then fails with
# "'x' must be an array of at least two dimensions". That happens inside the
# WGCNA package itself, so it can't be patched there; it has to be caught
# here, before anything gets handed to WGCNA. Correlation networks are also
# not *meaningful* at n=2 even where they don't crash -- these minimums are
# about producing a clear error instead of a cryptic one, not a claim that
# n=2/n=3 is a statistically adequate cohort.
# ----------------------------------------------------------------------------

n_regions_for_wgcna <- nrow(meth_mat_clean)
n_samples_for_wgcna <- ncol(meth_mat_clean)

if (n_regions_for_wgcna < 2) {
  stop(sprintf(paste(
    "[wgcna_modules] Only %d region(s) survived filtering -- WGCNA requires",
    "at least 2 to build a co-methylation network (pairwise correlation is",
    "undefined for a single region). This is expected on small/synthetic",
    "test cohorts (e.g. a mini/GIAB smoke-test config) where very few DMRs",
    "get called upstream. Run against a full-size DMR set, relax the",
    "upstream DMR-calling thresholds so more regions are called, or skip",
    "wgcna_modules for this config.\n",
    "  DMRs loaded from dmr_file:           %d\n",
    "  Regions after NA/variance filtering: %d"
  ), n_regions_for_wgcna, nrow(dmrs), n_regions_for_wgcna))
}

if (n_samples_for_wgcna < 3) {
  stop(sprintf(paste(
    "[wgcna_modules] Only %d sample(s) available after group filtering --",
    "correlation-based WGCNA needs at least 3 to compute a meaningful",
    "adjacency matrix. Check --group_col/--group_labels and the metadata TSV."
  ), n_samples_for_wgcna))
}

if (n_regions_for_wgcna < min_cluster_size) {
  warning(sprintf(paste(
    "[wgcna_modules] Only %d region(s) survived filtering, below",
    "min_cluster_size=%d. cutreeDynamic() will very likely detect 0 real",
    "modules (everything falls into the unassigned/grey module) even though",
    "the run will complete without error. Consider lowering",
    "--min_cluster_size for small cohorts, or treat this run as a smoke",
    "test rather than real module detection."
  ), n_regions_for_wgcna, min_cluster_size))
}


# ============================================================================
# 6. WGCNA Analysis
# ============================================================================

message("Step 5: Running WGCNA...")

datExpr <- t(meth_mat_clean)

message("  Checking sample and region quality...")
gsg <- goodSamplesGenes(datExpr, verbose = 0)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
message("  Retained ", nrow(datExpr), " samples and ", ncol(datExpr), " regions")

# Second guard, same failure mode as above: goodSamplesGenes() can itself
# flag enough samples/regions as low-quality (too many NAs, near-zero
# variance) to collapse datExpr down to a single row/column even when the
# pre-filter counts above looked fine. Catch that here too, with a message
# that points at *why* (how many were flagged), rather than letting it hit
# pickSoftThreshold()/adjacency() and fail there instead.
if (ncol(datExpr) < 2 || nrow(datExpr) < 3) {
  stop(sprintf(paste(
    "[wgcna_modules] After goodSamplesGenes() quality filtering, only %d",
    "region(s) and %d sample(s) remain -- too few to proceed (need >=2",
    "regions, >=3 samples). goodSamplesGenes() flagged %d region(s) and",
    "%d sample(s) as low quality (excess NAs and/or near-zero variance);",
    "investigate those before rerunning."
  ), ncol(datExpr), nrow(datExpr), sum(!gsg$goodGenes), sum(!gsg$goodSamples)))
}

message("  Selecting soft power threshold...")
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
softPower <- if (is.na(sft$powerEstimate)) WGCNA_PARAMS$soft_power else sft$powerEstimate
message("  Selected soft power: ", softPower,
        if (is.na(sft$powerEstimate)) " (auto-selection failed, using params.R WGCNA_PARAMS$soft_power fallback)" else "")

# ----------------------------------------------------------------------------
# WGCNA correctness check: scale-free topology fit R^2 at the selected power.
# WGCNA assumes the co-methylation network is approximately scale-free --
# if this fit is poor, the modules downstream are not well-supported by the
# data, regardless of how clean the dendrogram/module plot looks. There was
# previously no visibility into this number at all. WGCNA's own guidance
# targets fit R^2 >= ~0.8-0.9; below that, treat modules as exploratory only.
# Sample count matters too -- WGCNA's authors recommend >=15-20 samples for
# a correlation network to be stable; below that, small cohorts (a lot of
# nanopore studies) will show unstable modules run-to-run.
fit_row <- sft$fitIndices[sft$fitIndices$Power == softPower, ]
fit_r2  <- if (nrow(fit_row) == 1) fit_row$SFT.R.sq else NA_real_
message(sprintf(
  "  Scale-free topology fit R^2 at power %s: %s%s",
  softPower,
  if (is.na(fit_r2)) "NA" else round(fit_r2, 3),
  if (!is.na(fit_r2) && fit_r2 < 0.8) "  <-- WARNING: below 0.8, network may not be scale-free; treat modules as exploratory" else ""
))
message(sprintf(
  "  Samples going into WGCNA: %d%s",
  nrow(datExpr),
  if (nrow(datExpr) < 15) "  <-- WARNING: WGCNA typically needs >=15-20 samples for stable modules" else ""
))
write.csv(sft$fitIndices, file.path(outdir, "soft_threshold_fit.csv"), row.names = FALSE)

message("  Building adjacency and TOM matrices...")
adjacency_mat <- adjacency(datExpr, power = softPower)
TOM <- TOMsimilarity(adjacency_mat)
dissTOM <- 1 - TOM

message("  Detecting modules...")
geneTree <- hclust(as.dist(dissTOM), method = "average")
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = min_cluster_size
)
moduleColors <- labels2colors(dynamicMods)

message("  Identified ", length(unique(moduleColors)), " modules")


# ============================================================================
# 7. Prepare Module Results
#
# Saved under $regions so annotate_dmr_context.R's
# `modules <- readRDS(...); regions <- modules$regions` reads it directly.
# ============================================================================

message("Step 6: Preparing results...")

region_ids_wgcna <- colnames(datExpr)
region_metadata_filtered <- region_metadata[
  match(region_ids_wgcna, region_metadata$RegionID),
]
stopifnot(nrow(region_metadata_filtered) == length(moduleColors))
region_metadata_filtered$module <- moduleColors

message("  Module distribution:")
print(table(region_metadata_filtered$module))

modules <- list(regions = region_metadata_filtered)


# ============================================================================
# 8. Save Results
# ============================================================================

message("Step 7: Saving results...")

# Restrict to the samples goodSamplesGenes() actually kept for WGCNA --
# BSmooth(BSobj) on the full, unfiltered object would let low-quality
# samples (flagged and excluded from module detection above) still show up
# in bsseq_smoothed.rds, which dmr_heatmap.R reads directly. That mismatch
# meant the heatmap could visually include samples WGCNA itself excluded.
good_sample_names <- rownames(datExpr)
BSobj_good <- BSobj[, colnames(BSobj) %in% good_sample_names]
bs_filtered_smoothed <- BSmooth(BSobj_good)

saveRDS(modules, file = file.path(outdir, "DMR_Modules.rds"))
saveRDS(bs_filtered_smoothed, file = file.path(outdir, "bsseq_smoothed.rds"))

write.table(
  region_metadata_filtered,
  file = file.path(outdir, "wgcna_modules.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Results saved to: ", outdir)
message("  - DMR_Modules.rds        (consumed by annotate_dmr_context.R)")
message("  - bsseq_smoothed.rds     (consumed by dmr_heatmap.R; now matches WGCNA's good-sample subset)")
message("  - wgcna_modules.tsv")
message("  - soft_threshold_fit.csv (scale-free topology fit per power -- check this before trusting modules)")
message("Analysis complete.")