##############################################################################
# run_dss.R — DSS DML test + DMR calling
#
# Called by Snakemake rule dss_call_dmrs via:
#   Rscript run_dss.R --outfile <path> --repo_root <path>
#
# Steps:
#   1. Parse CLI args and load config
#   2. Load + filter sample metadata (shared helper, config/utils.R)
#   3. Read modkit BED files → BSseq object (shared helper, config/utils.R)
#   4. Coverage filter (also in the shared helper)
#   5. DMLtest (with smoothing)
#   6. callDMR with DSS_PARAMS from params.R
#   7. Write DMR_results.csv (raw DSS columns, unchanged downstream contract)
#      + DMR_results_normalized.csv (shared cross-method schema -- see
#        normalize_dmr_schema() in utils.R for why this is a second file
#        rather than replacing DMR_results.csv for now)
##############################################################################

suppressPackageStartupMessages({
  library(DSS)
  library(data.table)
  library(glue)
  library(GenomicRanges)   # for summarize_region_evidence()'s GRanges/findOverlaps
                           # (DSS imports this transitively, but don't rely on that)
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


if (is.null(repo_root))    stop("--repo_root is required")
if (is.null(outfile))      stop("--outfile is required")
if (is.null(metadata_tsv)) stop("--metadata is required")
if (length(group_values) != 2) {
  stop(sprintf("--group_values must specify exactly 2 values, got: %s", group_values_raw))
}


config_dir <- file.path(repo_root, "analysis", "dmr", "config")
source(file.path(config_dir, "utils.R"))
source(file.path(config_dir, "params.R"))
# NOTE: no longer sourcing datasets.R — group info comes from metadata TSV


# ============================================================================
# 2. Load + filter sample metadata (shared helper -- same logic dmrseq uses)
# ============================================================================

ds <- load_dmr_metadata(
  metadata_tsv  = metadata_tsv,
  group_column  = group_column,
  group_values  = group_values,
  restrict_role = restrict_role
)

cat(sprintf("\n[DSS] %d samples loaded\n", length(ds$sample_names)))
cat(sprintf("[DSS] %s: %s\n", group_values[1],
            paste(ds$sample_names[ds$groups == group_values[1]], collapse=", ")))
cat(sprintf("[DSS] %s: %s\n\n", group_values[2],
            paste(ds$sample_names[ds$groups == group_values[2]], collapse=", ")))


# ============================================================================
# 3-4. Read modkit BEDs -> coverage-filtered BSseq object
# ============================================================================

BSobj <- build_bsseq_from_metadata(
  bed_paths       = ds$bed_paths,
  sample_names    = ds$sample_names,
  groups          = ds$groups,
  coverage_params = COVERAGE
)

# build_bsseq_from_metadata() drops missing BED files internally; re-derive
# sample_names/groups from the columns of the returned object so downstream
# group assignment for DMLtest stays in sync even if files were dropped.
kept_samples <- colnames(BSobj)
kept_mask    <- ds$sample_names %in% kept_samples
ds$sample_names <- ds$sample_names[kept_mask]
ds$groups       <- ds$groups[kept_mask]

# Prints a per-group sample summary and enforces a real minimum-per-group
# size (not just "both groups exist") -- see validate_group_sizes() in
# utils.R.
validate_group_sizes(ds$sample_names, ds$groups, min_per_group = 2, label = "DSS")


# ============================================================================
# 5. DML Test
# ============================================================================

group1_samples <- ds$sample_names[ds$groups == group_values[1]]
group2_samples <- ds$sample_names[ds$groups == group_values[2]]

cat(glue("[DSS] DMLtest: {group_values[1]}={paste(group1_samples, collapse=',')} ",
         "vs {group_values[2]}={paste(group2_samples, collapse=',')}\n"))
cat("[DSS] smoothing = TRUE (recommended for Nanopore whole-genome data)\n\n")

dmlTest <- DMLtest(
  BSobj,
  group1    = group1_samples,
  group2    = group2_samples,
  smoothing = TRUE
)

cat(glue("[DSS] DMLtest complete: {nrow(dmlTest)} sites tested\n\n"))


# ============================================================================
# 6. Call DMRs
# ============================================================================

dmrs <- callDMR(
  dmlTest,
  delta       = DSS_PARAMS$delta,
  p.threshold = DSS_PARAMS$p.threshold,
  minlen      = DSS_PARAMS$minlen,
  minCG       = DSS_PARAMS$minCpG,
  dis.merge   = DSS_PARAMS$dis.merge,
  pct.sig     = DSS_PARAMS$pct.sig
)


# ============================================================================
# 6b. Per-region evidence from the underlying per-CpG DML test
#
# callDMR() never computes a region-level p-value -- see the docstring on
# normalize_dmr_schema() in utils.R for why (summing correlated per-CpG
# test statistics without correcting for spatial correlation would be
# anti-conservative). What DSS DOES compute correctly, with proper per-site
# dispersion shrinkage, is a per-CpG p-value/FDR in dmlTest. This section
# summarizes those per-CpG values within each merged region -- min p-value,
# mean FDR, and fraction of member CpGs individually significant at
# fdr < 0.05 -- as transparent supporting evidence, NOT as a substitute
# region-level p/q-value. Use these to eyeball how strong the underlying
# per-CpG signal is within a region; don't treat them as directly
# comparable to dmrseq's permutation-based region q-value, which is a
# different kind of statistic answering a different question.
# ============================================================================

summarize_region_evidence <- function(dmrs, dmlTest, sig_fdr = 0.05) {
  empty <- data.frame(
    min_dml_pval     = numeric(0),
    mean_dml_fdr     = numeric(0),
    frac_cpg_fdr_sig = numeric(0)
  )
  if (is.null(dmrs) || nrow(dmrs) == 0) return(empty)

  dml_gr <- GenomicRanges::GRanges(
    seqnames = dmlTest$chr,
    ranges   = IRanges::IRanges(start = dmlTest$pos, width = 1L)
  )
  region_gr <- GenomicRanges::GRanges(
    seqnames = dmrs$chr,
    ranges   = IRanges::IRanges(start = dmrs$start, end = dmrs$end)
  )

  hits      <- GenomicRanges::findOverlaps(region_gr, dml_gr)
  by_region <- split(S4Vectors::subjectHits(hits), S4Vectors::queryHits(hits))

  rows <- lapply(seq_len(nrow(dmrs)), function(i) {
    idx <- by_region[[as.character(i)]]
    if (is.null(idx)) {
      return(c(min_dml_pval = NA_real_, mean_dml_fdr = NA_real_, frac_cpg_fdr_sig = NA_real_))
    }
    c(
      min_dml_pval     = min(dmlTest$pval[idx], na.rm = TRUE),
      mean_dml_fdr     = mean(dmlTest$fdr[idx],  na.rm = TRUE),
      frac_cpg_fdr_sig = mean(dmlTest$fdr[idx] < sig_fdr, na.rm = TRUE)
    )
  })

  as.data.frame(do.call(rbind, rows))
}

dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)

if (is.null(dmrs) || nrow(dmrs) == 0) {

  warning("No DMRs detected by callDMR.")

  empty <- data.frame(
    chr = character(),
    start = integer(),
    end = integer()
  )

  write.csv(empty, outfile, row.names = FALSE)
  normalized <- normalize_dmr_schema(NULL, method = "dss")

} else {

  region_evidence <- summarize_region_evidence(dmrs, dmlTest)
  dmrs_with_evidence <- cbind(dmrs, region_evidence)

  # Additive columns only (min_dml_pval, mean_dml_fdr, frac_cpg_fdr_sig) --
  # every column callDMR() originally produced is still present, in the
  # same order, so anything downstream reading this CSV by column name
  # (e.g. run_bsmooth.R) is unaffected. If anything downstream reads it by
  # column POSITION instead, that would break -- worth confirming before a
  # full-cohort run if you're not sure.
  write.csv(dmrs_with_evidence, outfile, row.names = FALSE)
  normalized <- normalize_dmr_schema(dmrs, method = "dss")

}

cat(glue("\n[DSS] Results written to {outfile}\n"))

# --------------------------------------------------------------------------
# Additive normalized-schema output. NOT yet the primary --outfile: DSS's
# raw columns (areaStat, diff.Methy, meanMethy1/2, length) are what
# run_bsmooth.R currently reads from DMR_results.csv via dss_smooth_betas.smk,
# and that script hasn't been checked against this schema yet. Also note
# DSS produces no region-level pval/qval (see normalize_dmr_schema()'s
# docstring in utils.R) -- those columns are NA here by design, not a bug.
# --------------------------------------------------------------------------
normalized_outfile <- file.path(dirname(outfile), "DMR_results_normalized.csv")
write.csv(normalized, normalized_outfile, row.names = FALSE)
cat(glue("[DSS] Normalized (cross-method) results written to {normalized_outfile}\n"))