##############################################################################
# params.R — Analysis parameters for DMR calling and downstream steps
#
# Usage: source("config/params.R")
#
# These are the values you change between runs. All downstream scripts read
# from this file — never hardcode parameters inside analysis scripts.
#
# Sections:
#   COVERAGE    — minimum read depth filters
#   DSS_PARAMS  — DSS DMR calling thresholds
#   DMRSEQ      — dmrseq DMR calling thresholds
#   INTEGRATION — WGCNA co-methylation network settings
#   ANNOTATION  — genomic context and pathway enrichment settings
##############################################################################


# ============================================================================
# Coverage Filters
# ============================================================================

COVERAGE <- list(
  min_cov   = 1,      # minimum read depth per CpG
  per_group = 0.75    # fraction of samples per group that must meet min_cov
                      # 1.0 = all samples, 0.75 = 75% of each group
)


# ============================================================================
# DSS Parameters
#
# Named DSS_PARAMS, not DSS -- the DSS package is also called DSS, and
# run_dss.R calls DSS::DMLtest()/DSS::callDMR() directly. Reusing the
# package's own name for this config list was legal R (DSS$delta doesn't
# collide with the DSS namespace) but reads ambiguously at every call site
# and invites exactly that collision if the code ever does library(DSS)
# without namespacing calls.
# ============================================================================

DSS_PARAMS <- list(
  delta       = 0.02,   # minimum methylation difference between groups
  p.threshold = 0.1,   # DML p-value cutoff
  minlen      = 20,     # minimum DMR length (bp)
  minCpG      = 2,      # minimum CpGs per DMR
  dis.merge   = 100,    # merge DMRs within this distance (bp)
  pct.sig     = 0.5     # fraction of CpGs in region that must be significant
)


# ============================================================================
# dmrseq Parameters
# ============================================================================

DMRSEQ <- list(
  cutoff  = 0.1,   # per-CpG methylation difference cutoff for candidate regions
  minCpGs = 5      # minimum CpGs per candidate region
)


# ============================================================================
# WGCNA / Signal Integration Parameters
#
# Named WGCNA_PARAMS, not WGCNA -- the WGCNA package is also called WGCNA,
# and wgcna_methylation.R does library(WGCNA). Same rationale as
# DSS_PARAMS above: WGCNA$min_module_size happened to work (package
# attachment doesn't create a global "WGCNA" variable), but it read
# ambiguously at the call site and is one accidental `WGCNA <- ...`
# reassignment away from actually breaking.
# ============================================================================

WGCNA_PARAMS <- list(
  soft_power      = 6,      # soft-thresholding power for adjacency
  min_module_size = 30,     # minimum genes per module
  merge_cut_height = 0.25   # dendrogram cut height for module merging
)


# ============================================================================
# Annotation Parameters
# ============================================================================

ANNOTATION <- list(
  genome          = "hg38",
  p.adjust.method = "BH",     # multiple testing correction (pathway enrichment)
  p.cutoff        = 0.05,     # adjusted p-value cutoff for enrichment results
  q.cutoff        = 0.2       # q-value cutoff for enrichment results
)