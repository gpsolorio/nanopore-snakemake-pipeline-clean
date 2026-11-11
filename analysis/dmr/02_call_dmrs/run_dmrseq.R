##############################################################################
# run_dmrseq_pipeline.R — dmrseq DMR calling for Snakemake
#
# Called by rule dmrseq_call_dmrs via:
#   Rscript run_dmrseq_pipeline.R \
#       --outfile <path> --repo_root <path> --barcodes <b1,b2,...> \
#       [--cutoff 0.1] [--min_cpgs 5] [--max_perms 10] [--cores 8]
##############################################################################

suppressPackageStartupMessages({
  library(dmrseq)
  library(data.table)
  library(bsseq)
  library(BiocParallel)
  library(GenomicRanges)
  library(glue)
  library(plyranges)
  library(dplyr)
})

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
outfile      <- get_arg(args, "--outfile")
metadata_tsv <- get_arg(args, "--metadata")
group_column <- get_arg(args, "--group_column", "diagnosis")
group_values_raw <- get_arg(args, "--group_values", "ASD,TD")
group_values <- strsplit(group_values_raw, ",")[[1]]
restrict_role <- get_arg(args, "--restrict_role", "")
cutoff       <- as.numeric(get_arg(args, "--cutoff",   "0.1"))
min_cpgs     <- as.numeric(get_arg(args, "--min_cpgs", "5"))
max_perms    <- as.numeric(get_arg(args, "--max_perms","10"))
cores        <- as.numeric(get_arg(args, "--cores",    "4"))
genome       <- get_arg(args, "--genome", "hg38")

if (is.null(repo_root))    stop("--repo_root is required")
if (is.null(outfile))      stop("--outfile is required")
if (is.null(metadata_tsv)) stop("--metadata is required")
if (length(group_values) != 2) {
  stop(sprintf("--group_values must specify exactly 2 values, got: %s", group_values_raw))
}

config_dir <- file.path(repo_root, "analysis", "dmr", "config")
source(file.path(config_dir, "utils.R"))
source(file.path(config_dir, "params.R"))
# NOTE: no longer sources datasets.R -- that file does not exist anywhere in
# the repo (only params.R and utils.R live in analysis/dmr/config/), so
# dmr.caller: dmrseq could never actually run before this fix. Group info
# now comes from the metadata TSV via --group_column/--group_values, the
# same mechanism run_dss.R already uses -- so both callers are driven by
# the same config.dmr.group_column/group_values, instead of dmrseq being
# hardcoded to a TD/ASD diagnosis split regardless of what's configured.

cat(glue("[dmrseq] cutoff={cutoff}, min_cpgs={min_cpgs}, max_perms={max_perms}, cores={cores}\n\n"))

# ============================================================================
# 2. Load sample metadata (shared helper -- same logic run_dss.R uses)
# ============================================================================

ds <- load_dmr_metadata(
  metadata_tsv  = metadata_tsv,
  group_column  = group_column,
  group_values  = group_values,
  restrict_role = restrict_role
)
meta <- ds$meta

cat(sprintf("[dmrseq] %d samples loaded\n", length(ds$sample_names)))
cat(sprintf("[dmrseq] %s: %s\n", group_values[1],
            paste(ds$sample_names[ds$groups == group_values[1]], collapse=", ")))
cat(sprintf("[dmrseq] %s: %s\n\n", group_values[2],
            paste(ds$sample_names[ds$groups == group_values[2]], collapse=", ")))

# ============================================================================
# 3. Build coverage-filtered BSseq object
#
# Uses the shared helper (config/utils.R) instead of the hand-rolled
# union-of-sites + per-group coverage filter that used to live here --
# that copy had drifted to hardcode perGroup <- 0.75 rather than reading
# COVERAGE$per_group/$min_cov from params.R, the way run_dss.R and
# wgcna_methylation.R already do. The shared helper no longer depends on
# the DSS package either (see utils.R), so this env can stay DSS-free.
# ============================================================================

BSobj <- build_bsseq_from_metadata(
  bed_paths       = ds$bed_paths,
  sample_names    = ds$sample_names,
  groups          = ds$groups,
  coverage_params = COVERAGE,
  metadata        = meta
)

kept_samples <- colnames(BSobj)
kept_mask    <- ds$sample_names %in% kept_samples
ds$sample_names <- ds$sample_names[kept_mask]
ds$groups       <- ds$groups[kept_mask]

# Prints a per-group sample summary and enforces a real minimum-per-group
# size (not just "both groups exist") -- see validate_group_sizes() in
# utils.R.
validate_group_sizes(ds$sample_names, ds$groups, min_per_group = 2, label = "dmrseq")

# build_bsseq_from_metadata() attaches every metadata column (incl.
# group_column) as pData(), but as character -- dmrseq's testCovariate and
# dmr_heatmap.R's annotation coloring both expect a factor.
pData(BSobj)[[group_column]] <- factor(pData(BSobj)[[group_column]], levels = group_values)

cat(glue("[dmrseq] {ncol(BSobj)} samples retained after coverage filter, {nrow(BSobj)} CpG sites\n\n"))

# ============================================================================
# 4. Run dmrseq
# ============================================================================

register(MulticoreParam(workers = cores))

cat(glue("[dmrseq] Running dmrseq (cutoff={cutoff}, minNumRegion={min_cpgs}, maxPerms={max_perms})...\n"))
start_time <- Sys.time()

regions <- dmrseq(
  bs             = BSobj,
  testCovariate  = group_column,
  cutoff         = cutoff,
  minNumRegion   = min_cpgs,
  maxPerms       = max_perms,
  BPPARAM        = MulticoreParam(workers = cores, progressbar = TRUE)
)

cat(glue("[dmrseq] Completed in {round(difftime(Sys.time(), start_time, units='mins'), 1)} min\n\n"))

# ============================================================================
# 5. Select significant regions
# ============================================================================

if (sum(regions$qval < 0.05, na.rm = TRUE) >= 100) {
  sigRegions <- regions %>% plyranges::filter(qval < 0.05)
  cat(glue("[dmrseq] Using q-value < 0.05: {length(sigRegions)} DMRs\n"))
} else if (sum(regions$pval < 0.05, na.rm = TRUE) > 0) {
  sigRegions <- regions %>% plyranges::filter(pval < 0.05)
  cat(glue("[dmrseq] Using p-value < 0.05: {length(sigRegions)} DMRs\n"))
} else {
  stop("[dmrseq] No significant DMRs detected")
}

# ============================================================================
# 6. Write normalized CSV (shared cross-method schema)
#
# Previously this wrote its own hand-built column set that happened to
# resemble but not exactly match run_dss.R's. Both callers now go through
# normalize_dmr_schema() in utils.R, so downstream comparison code sees one
# consistent schema regardless of which caller produced a given file.
# ============================================================================

df <- as.data.frame(sigRegions) %>%
  dplyr::transmute(
    chr       = as.character(seqnames),
    start     = start,
    end       = end,
    nCG       = L,           # number of CpGs in region
    stat      = stat,        # dmrseq test statistic
    pval      = pval,
    qval      = qval,
    direction = dplyr::case_when(stat > 0 ~ "Hypermethylated", TRUE ~ "Hypomethylated")
  )

normalized <- normalize_dmr_schema(df, method = "dmrseq")

dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
write.csv(normalized, outfile, row.names = FALSE)
cat(glue("\n[dmrseq] {nrow(normalized)} DMRs written to {outfile}\n"))