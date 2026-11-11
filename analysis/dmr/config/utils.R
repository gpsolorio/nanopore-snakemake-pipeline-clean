##############################################################################
# utils.R — Generic helper functions for DMR analysis scripts
#
# Usage: source("config/utils.R")
#
# Functions:
#   get_project_root()          — auto-detect repo root
#   create_output_dirs()        — create results/ subdirectory structure
#   generate_filename()         — build parameter-encoded output filenames
#   parse_arguments()           — parse CLI args with defaults
#   print_config()              — pretty-print analysis configuration
#   load_dmr_metadata()         — read + filter sample_metadata.tsv by group
#   validate_group_sizes()      — print + enforce a minimum samples/group
#   build_bsseq_from_metadata() — read modkit BEDs -> coverage-filtered BSseq
#   normalize_dmr_schema()      — map DSS/dmrseq output to one shared schema
##############################################################################


# ============================================================================
# Project Root Detection
# ============================================================================

#' Auto-detect the nanopore-snakemake-pipeline project root
#'
#' Walks up from the current working directory looking for a `workflow/`
#' directory. Falls back to $PROJECT_ROOT and standard install locations.
#'
#' @return Absolute path to project root (character)
get_project_root <- function() {
  current <- getwd()

  while (current != "/" && !file.exists(file.path(current, "workflow"))) {
    current <- dirname(current)
  }

  if (file.exists(file.path(current, "workflow"))) {
    return(current)
  }

  fallback_paths <- c(
    Sys.getenv("PROJECT_ROOT"),
    file.path(Sys.getenv("HOME"), "nanopore-snakemake-pipeline"),
    "/opt/nanopore-snakemake-pipeline"
  )

  for (path in fallback_paths) {
    if (nzchar(path) && file.exists(file.path(path, "workflow"))) {
      return(path)
    }
  }

  stop(
    "Could not detect project root. ",
    "Run from within the repo, or set the PROJECT_ROOT environment variable."
  )
}


# ============================================================================
# Output Directory Management
# ============================================================================

#' Create and return the standard output directory structure for an analysis step
#'
#' Creates: <results_root>/<step>/{output,logs,rdata,plots}
#'
#' @param step         Step directory name, e.g. "02_call_dmrs"
#' @param results_root Path to results root (default: <project_root>/results)
#'
#' @return Named list with paths: root, output, logs, rdata, plots
create_output_dirs <- function(step,
                               results_root = file.path(get_project_root(), "results")) {
  outdir  <- file.path(results_root, step)
  subdirs <- c("output", "logs", "rdata", "plots")

  for (s in subdirs) {
    dir.create(file.path(outdir, s), showWarnings = FALSE, recursive = TRUE)
  }

  as.list(setNames(
    c(outdir, file.path(outdir, subdirs)),
    c("root",  subdirs)
  ))
}


# ============================================================================
# Filename Generation
# ============================================================================

#' Build a parameter-encoded output filename
#'
#' Produces filenames like:
#'   GIAB_2023_05_dss_delta0.1_DMRs.tsv
#'   STUDY_dmrseq_cutoff0.1_minCpGs5_DMRs.bed
#'
#' Parameters are appended in sorted order so filenames are deterministic.
#'
#' @param study_label Cohort/study label, e.g. "MY_STUDY" or "GIAB_2023_05" --
#'   no default is provided. Earlier versions hardcoded "modkit_AS_TD_vs_ASD"
#'   here, which silently mislabeled outputs from any non-adaptive-sampling cohort
#'   (this pipeline now also runs the private cohort and GIAB data). Pass
#'   config$experiment$study (already present in every hardware config) from
#'   the calling script.
#' @param base_name  Descriptive suffix, e.g. "DMRs" or "BSobj"
#' @param method     "dss" or "dmrseq"
#' @param params     Named list of analysis parameters to encode
#' @param extension  File extension without dot, e.g. "tsv" (omit for none)
#'
#' @return Filename string
generate_filename <- function(study_label, base_name, method, params = list(), extension = "") {
  if (missing(study_label) || !nzchar(study_label)) {
    stop("generate_filename(): study_label is required (e.g. config$experiment$study) -- ",
         "no default is provided since a hardcoded default previously mislabeled ",
         "outputs from cohorts other than the one it was written for.")
  }

  valid_methods <- c("dss", "dmrseq")
  if (!method %in% valid_methods) {
    stop(sprintf("Invalid method '%s'. Must be one of: %s",
                 method, paste(valid_methods, collapse = ", ")))
  }

  parts <- c(study_label, method)

  if (length(params) > 0) {
    param_strs <- vapply(sort(names(params)), function(nm) {
      v <- params[[nm]]
      if (is.numeric(v)) {
        v <- sub("0+$", "", sub("\\.$", "", as.character(round(v, 4))))
      }
      paste0(nm, v)
    }, character(1))
    parts <- c(parts, param_strs)
  }

  fname <- paste(c(parts, base_name), collapse = "_")
  if (nzchar(extension)) fname <- paste0(fname, ".", extension)
  fname
}


# ============================================================================
# CLI Argument Parsing
# ============================================================================

#' Parse --key value style CLI arguments with fallback defaults
#'
#' Numeric-looking values are coerced to numeric. Bare flags (no following
#' value) are set to TRUE.
#'
#' @param args     Character vector, typically commandArgs(trailingOnly = TRUE)
#' @param defaults Named list of default values
#'
#' @return Named list merging defaults with any overrides from args
#'
#' @examples
#' parse_arguments(c("--delta", "0.1", "--p.threshold", "0.001"),
#'                 list(delta = 0.05, p.threshold = 0.01))
parse_arguments <- function(args, defaults) {
  parsed <- defaults
  i      <- 1L

  while (i <= length(args)) {
    if (startsWith(args[i], "--")) {
      key <- sub("^--", "", args[i])

      if (i < length(args) && !startsWith(args[i + 1L], "--")) {
        val <- args[i + 1L]
        if (grepl("^-?[0-9]+(\\.[0-9]+)?$", val)) val <- as.numeric(val)
        parsed[[key]] <- val
        i <- i + 2L
      } else {
        parsed[[key]] <- TRUE
        i <- i + 1L
      }
    } else {
      i <- i + 1L
    }
  }

  parsed
}


# ============================================================================
# Logging
# ============================================================================

#' Pretty-print an analysis configuration block
#'
#' @param analysis_name  Label shown in the header
#' @param params         Named list of scalar parameters to display
print_config <- function(analysis_name, params) {
  cat(sprintf("\n[%s] Configuration  %s\n", analysis_name,
              format(Sys.time(), "%d-%m-%Y %X")))

  for (nm in names(params)) {
    if (!is.list(params[[nm]])) {
      cat(sprintf("  %s: %s\n", nm, params[[nm]]))
    }
  }
  cat("\n")
}


# ============================================================================
# Shared Metadata Loading
#
# Extracted from what were near-identical duplicated blocks in run_dss.R
# and run_dmrseq_pipeline.R. Both callers now see exactly the same filtered
# sample set built from exactly the same logic -- previously a hand-edit to
# one script's filtering block could silently drift out of sync with the
# other's.
# ============================================================================

#' Read sample_metadata.tsv and filter to the two configured group values
#'
#' @param metadata_tsv Path to the metadata TSV (sample_id, bed_path, role,
#'   and a group column at minimum)
#' @param group_column  Column name to group samples on, e.g. "diagnosis"
#' @param group_values  Character vector of exactly 2 group values to keep
#' @param restrict_role Optional. If non-empty, restricts to rows whose
#'   `role` column equals this value (e.g. "PRI") before applying the
#'   group_column/group_values filter. Use this when the sample sheet's
#'   `roles:` list includes relatives alongside probands (so everyone gets
#'   processed through basecalling/alignment/pileup) but the group_column
#'   value on a relative's row is their FAMILY's label, not their own --
#'   comparing relatives as if diagnosed themselves would confound
#'   diagnosis with generation/relationship. Empty string (the default)
#'   means no role restriction, preserving old behavior for sample sheets
#'   with no role concept (id_col mode) or with roles already limited to
#'   probands only.
#'
#' @return A list:
#'   meta         — the filtered metadata data.frame (all columns preserved)
#'   sample_names — meta$sample_id
#'   groups       — meta[[group_column]]
#'   bed_paths    — meta$bed_path
load_dmr_metadata <- function(metadata_tsv, group_column, group_values, restrict_role = "") {
  if (length(group_values) != 2) {
    stop(sprintf(
      "load_dmr_metadata(): group_values must specify exactly 2 values, got: %s",
      paste(group_values, collapse = ",")
    ))
  }

  meta <- read.table(metadata_tsv, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, colClasses = "character")

  if (!group_column %in% colnames(meta)) {
    stop(sprintf(
      "Group column '%s' not found in metadata TSV.\n  Available columns: %s",
      group_column, paste(colnames(meta), collapse = ", ")
    ))
  }

  if (nzchar(restrict_role)) {
    if (!"role" %in% colnames(meta)) {
      stop("load_dmr_metadata(): --restrict_role was given but the metadata ",
           "TSV has no 'role' column. Regenerate sample_metadata.tsv with a ",
           "current common.smk (older runs predate the role column).")
    }
    n_before <- nrow(meta)
    meta <- meta[meta$role == restrict_role, ]
    cat(glue::glue(
      "[load_dmr_metadata] --restrict_role {restrict_role}: kept {nrow(meta)}/{n_before} ",
      "metadata rows (excluded relatives whose group_column value is their ",
      "family's label, not their own)\n\n"
    ))
    if (nrow(meta) == 0) {
      stop(sprintf(
        "restrict_role='%s' matched zero rows in the metadata TSV. Check the ",
        "'role' column's actual values (case-sensitive) against this setting.",
        restrict_role
      ))
    }
  }

  meta <- meta[meta[[group_column]] %in% group_values, ]

  if (nrow(meta) == 0) {
    stop(sprintf(
      "No samples found with '%s' in {%s} in metadata TSV.\n  Check your config's dmr.group_column/dmr.group_values%s and the metadata TSV contents.",
      group_column, paste(group_values, collapse = ", "),
      if (nzchar(restrict_role)) sprintf("/dmr.restrict_role='%s'", restrict_role) else ""
    ))
  }

  list(
    meta         = meta,
    sample_names = meta$sample_id,
    groups       = meta[[group_column]],
    bed_paths    = meta$bed_path
  )
}


# ============================================================================
# Shared Group-Size Validation
# ============================================================================

#' Print a per-group sample summary and enforce a minimum group size
#'
#' Replaces the old `length(unique(groups)) < 2` check, which only verified
#' that both groups were represented at all -- not that either had enough
#' samples for the statistical test to be meaningful.
#'
#' @param sample_names  Character vector of sample IDs
#' @param groups        Character vector, same length, group label per sample
#' @param min_per_group Minimum required samples per group (default 2)
#' @param label         Short caller label for printed output, e.g. "DSS" or "dmrseq"
#'
#' @return (invisibly) the table() of group counts
validate_group_sizes <- function(sample_names, groups, min_per_group = 2, label = "") {
  stopifnot(length(sample_names) == length(groups))

  prefix       <- if (nzchar(label)) sprintf("[%s] ", label) else ""
  group_counts <- table(groups)

  cat(sprintf("\n%sGroup summary:\n", prefix))
  for (g in names(group_counts)) {
    cat(sprintf("  %s (%d)\n", g, group_counts[[g]]))
    cat(paste0("    ", sample_names[groups == g], collapse = "\n"))
    cat("\n")
  }
  cat("\n")

  too_small <- group_counts[group_counts < min_per_group]
  if (length(too_small) > 0 || length(group_counts) < 2) {
    stop(sprintf(
      "%sEach of the 2 groups must contain at least %d sample(s). Got: %s",
      prefix, min_per_group,
      paste(names(group_counts), group_counts, sep = "=", collapse = ", ")
    ))
  }

  invisible(group_counts)
}


# ============================================================================
# Shared BSseq Construction
#
# Extracted from what were previously near-duplicate copies of this exact
# logic in run_dss.R and run_bsmooth.R. wgcna_methylation.R now uses it too,
# instead of expecting a saved BSobj.rds that nothing in the pipeline ever
# wrote (that mismatch was the reason WGCNA module detection couldn't run).
#
# build_bsseq_from_metadata() previously called DSS::makeBSseqData() to
# assemble the BSseq object, which meant the "shared" helper -- and every
# caller of it, including the supposedly DSS-independent dmrseq pipeline --
# had a hard dependency on the DSS package. It now builds the BSseq object
# directly via bsseq::BSseq(), so dmrseq-env.yml (which doesn't and
# shouldn't ship DSS) can run this unmodified. Only bsseq, data.table, and
# S4Vectors need to be loaded by the caller now.
# ============================================================================

#' Merge per-sample modkit read tables into BSseq-ready M/Cov matrices
#'
#' Union of all chr/pos sites seen across samples; sites a given sample has
#' no coverage at get 0/0 (methylated count / total count) for that sample,
#' matching DSS::makeBSseqData()'s prior behavior.
#'
#' @param site_list List of data.tables, one per sample, each with columns
#'   chr, pos, N (total count), X (methylated count)
#' @param sample_names Character vector, same length/order as site_list
#'
#' @return A bsseq::BSseq object
build_bsseq_native <- function(site_list, sample_names) {
  stopifnot(length(site_list) == length(sample_names))

  dt_all <- data.table::rbindlist(
    Map(function(dt, nm) { dt <- data.table::copy(dt); dt$sample <- nm; dt },
        site_list, sample_names)
  )

  # fun.aggregate = sum guards against duplicate chr/pos rows within a
  # single sample's BED (shouldn't happen for a pileup, but dcast() errors
  # out on ambiguous duplicates by default rather than silently picking one).
  wide_X <- data.table::dcast(dt_all, chr + pos ~ sample, value.var = "X",
                              fill = 0L, fun.aggregate = sum)
  wide_N <- data.table::dcast(dt_all, chr + pos ~ sample, value.var = "N",
                              fill = 0L, fun.aggregate = sum)

  data.table::setorder(wide_X, chr, pos)
  data.table::setorder(wide_N, chr, pos)

  # dcast's generated sample columns follow sort order of the "sample"
  # values, not the order sample_names was passed in -- reindex explicitly
  # so M/Cov columns line up with sample_names position-for-position.
  M   <- as.matrix(wide_X[, ..sample_names])
  Cov <- as.matrix(wide_N[, ..sample_names])

  bsseq::BSseq(
    chr         = wide_X$chr,
    pos         = wide_X$pos,
    M           = M,
    Cov         = Cov,
    sampleNames = sample_names
  )
}

#' Read modkit pileup BED files into a coverage-filtered BSseq object
#'
#' @param bed_paths Character vector of modkit BED file paths
#' @param sample_names Character vector, same length/order as bed_paths
#' @param groups Character/factor vector, same length/order as bed_paths --
#'   used only for the per-group coverage filter
#' @param coverage_params A COVERAGE-shaped list (see params.R): must have
#'   $min_cov and $per_group
#' @param metadata Optional full sample metadata data.frame (e.g. the parsed
#'   sample_metadata.tsv) with a `sample_id` column. If provided, every
#'   column is attached to the returned BSseq object's pData(), so
#'   downstream scripts (e.g. dmr_heatmap.R --test-covariate) can reference
#'   any of them, not just the group column used for filtering here.
#'
#' @return A coverage-filtered BSseq object, with pData() populated from
#'   `metadata` if supplied
build_bsseq_from_metadata <- function(bed_paths, sample_names, groups, coverage_params, metadata = NULL) {
  stopifnot(length(bed_paths) == length(sample_names),
            length(bed_paths) == length(groups))

  exists_mask <- file.exists(bed_paths)
  if (any(!exists_mask)) {
    warning(sprintf(
      "%d BED file(s) not found and will be excluded:\n  %s",
      sum(!exists_mask),
      paste(bed_paths[!exists_mask], collapse = "\n  ")
    ))
    bed_paths    <- bed_paths[exists_mask]
    sample_names <- sample_names[exists_mask]
    groups       <- groups[exists_mask]
  }

  read_modkit_to_dss <- function(file_path) {
    dt <- data.table::fread(file_path, select = c(1L, 2L, 10L, 12L))
    data.table::setnames(dt, c("chr", "start", "N", "X"))
    dt[, pos := start + 1L]   # 0-based -> 1-based
    dt[, start := NULL]
    dt
  }

  cat(glue::glue("[build_bsseq_from_metadata] Reading {length(bed_paths)} BED files...\n\n"))
  site_list <- lapply(bed_paths, read_modkit_to_dss)
  BSobj     <- build_bsseq_native(site_list, sample_names)
  cat(glue::glue("[build_bsseq_from_metadata] BSseq object: {nrow(BSobj)} CpG sites, {ncol(BSobj)} samples\n\n"))

    if (!is.null(metadata)) {
    if (!"sample_id" %in% colnames(metadata)) {
      stop("[build_bsseq_from_metadata] `metadata` must have a `sample_id` column.")
    }
    meta_ordered <- metadata[match(sample_names, metadata$sample_id), , drop = FALSE]
    if (any(is.na(meta_ordered$sample_id))) {
      missing_samples <- sample_names[is.na(meta_ordered$sample_id)]
      stop("[build_bsseq_from_metadata] No metadata row found for sample(s): ",
           paste(missing_samples, collapse = ", "))
    }
    meta_ordered <- meta_ordered[, !colnames(meta_ordered) %in% c("bed_path"), drop = FALSE]
    pData(BSobj) <- S4Vectors::DataFrame(meta_ordered, row.names = sample_names)
    cat(glue::glue("[build_bsseq_from_metadata] Attached pData() columns: ",
                    "{paste(colnames(meta_ordered), collapse = ', ')}\n\n"))
  }

  groups_f    <- factor(groups)
  group_sizes <- table(groups_f)
  cov_mat     <- bsseq::getCoverage(BSobj, type = "Cov")

  keep <- Reduce(`&`, lapply(levels(groups_f), function(g) {
    idx       <- which(groups == g)
    threshold <- ceiling(group_sizes[g] * coverage_params$per_group)
    rowSums(cov_mat[, idx, drop = FALSE] >= coverage_params$min_cov) >= threshold
  }))

  BSobj <- BSobj[keep, ]
  cat(glue::glue("[build_bsseq_from_metadata] Retained {nrow(BSobj)} CpG sites after coverage filter\n\n"))

  BSobj
}


# ============================================================================
# Shared Output Schema
#
# DSS::callDMR() and dmrseq's region output use different column sets and
# different statistical outputs (see normalize_dmr_schema() below for the
# specifics). This maps both onto one shared schema so a downstream
# comparison step doesn't need to know which caller produced a given file.
#
# NOTE: this is used unconditionally by run_dmrseq_pipeline.R now. It is
# NOT yet wired into run_dss.R's primary --outfile write, because
# dss_smooth_betas.smk feeds that same DMR_results.csv into run_bsmooth.R,
# and it's not yet confirmed which columns that script expects -- swapping
# DSS's output schema without checking that first risks silently breaking
# smoothing. run_dss.R instead writes the normalized schema to a *second*,
# clearly-labeled file alongside its existing raw output until that's
# verified.
# ============================================================================

#' Map a caller's raw DMR output onto the shared comparison schema
#'
#' Shared columns: chr, start, end, nCG, stat, pval, qval, direction, method
#'
#' IMPORTANT ASYMMETRY: DSS::callDMR() does not report a region-level
#' p-value or q-value (DSS only reports per-CpG DML p-values; region calls
#' come from areaStat/diff.Methy/pct.sig thresholds, not a region p-value).
#' pval/qval are therefore NA for method = "dss" -- this is not a bug in
#' this function, it's a genuine difference in what the two methods report.
#' Don't treat a DSS row's NA pval/qval as "not significant"; DSS already
#' applied its own significance criteria before a region reached this table.
#'
#' @param df     Raw result data.frame from callDMR() (DSS) or the dmrseq
#'   post-processing step in run_dmrseq_pipeline.R
#' @param method "dss" or "dmrseq"
#'
#' @return data.frame with columns: chr, start, end, nCG, stat, pval, qval,
#'   direction, method
normalize_dmr_schema <- function(df, method) {
  common_cols <- c("chr", "start", "end", "nCG", "stat", "pval", "qval", "direction", "method")

  if (!method %in% c("dss", "dmrseq")) {
    stop(sprintf("normalize_dmr_schema(): unknown method '%s'", method))
  }

  if (is.null(df) || nrow(df) == 0) {
    empty <- as.data.frame(matrix(ncol = length(common_cols), nrow = 0))
    colnames(empty) <- common_cols
    return(empty)
  }

  if (method == "dss") {
    out <- data.frame(
      chr       = df$chr,
      start     = df$start,
      end       = df$end,
      nCG       = df$nCG,
      stat      = df$areaStat,
      pval      = NA_real_,
      qval      = NA_real_,
      direction = ifelse(df$diff.Methy > 0, "Hypermethylated", "Hypomethylated"),
      method    = "dss",
      stringsAsFactors = FALSE
    )
  } else {
    out <- data.frame(
      chr       = df$chr,
      start     = df$start,
      end       = df$end,
      nCG       = df$nCG,
      stat      = df$stat,
      pval      = df$pval,
      qval      = df$qval,
      direction = df$direction,
      method    = "dmrseq",
      stringsAsFactors = FALSE
    )
  }

  out[, common_cols]
}