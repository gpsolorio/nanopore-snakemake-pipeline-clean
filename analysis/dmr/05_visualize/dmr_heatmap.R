#!/usr/bin/env Rscript
##############################################################################
# DMR Z-Score Heatmap
##############################################################################
# Wraps the smoothPheatmap plotting function as a standalone CLI script.
##############################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(bsseq)
  library(pheatmap)
  library(dplyr)
  library(RColorBrewer)
  library(glue)
  library(purrr)
  library(magrittr)
})

option_list <- list(
  make_option(c("--bsseq-file"), type = "character",
              help = "Smoothed BSseq object (.rds), e.g. bsseq_smoothed.rds from wgcna_methylation.R"),
  make_option(c("--dmrs-file"), type = "character",
              help = "DMR results table (chr/start/end columns), CSV -- defines sigRegions"),
  make_option(c("--test-covariate"), type = "character",
              help = "Column name in the BSseq pData() to color-annotate by (e.g. 'sex', 'diagnosis')"),
  make_option(c("--output-file"), type = "character", default = "./heatmap.pdf")
)

args <- parse_args(OptionParser(option_list = option_list))

if (is.null(args$`bsseq-file`) || is.null(args$`dmrs-file`) || is.null(args$`test-covariate`)) {
  stop("--bsseq-file, --dmrs-file, and --test-covariate are all required.")
}

bsseq_file <- normalizePath(args$`bsseq-file`, mustWork = TRUE)
dmrs_file  <- normalizePath(args$`dmrs-file`, mustWork = TRUE)
testCovariate <- args$`test-covariate`
output_file <- args$`output-file`

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

message("Loading smoothed BSseq object...")
bs.filtered.bsseq <- readRDS(bsseq_file)

message("Loading DMRs to build sigRegions...")
dmrs <- read.csv(dmrs_file, stringsAsFactors = FALSE)

message("nrow: ", nrow(dmrs), " | chr: ", length(dmrs$chr),
        " start: ", length(dmrs$start), " end: ", length(dmrs$end))
str(dmrs)

sigRegions <- GRanges(
  seqnames = dmrs$chr,
  ranges   = IRanges(start = dmrs$start, end = dmrs$end)
)

if (!testCovariate %in% colnames(pData(bs.filtered.bsseq))) {
  stop("--test-covariate '", testCovariate, "' not found in BSseq pData() columns: ",
       paste(colnames(pData(bs.filtered.bsseq)), collapse = ", "))
}

#' smoothPheatmap
#' @title DMR heatmap
#' @description Plot a heatmap of normalized individual smoothed methylation
#'   value z-scores for selected regions (i.e. significant DMRs)
#' @param bs.filtered.bsseq Smoothed \code{bsseq} object
#' @param sigRegions \code{GRanges} object of regions to plot a heatmap for
#' @param testCovariate The factor tested for differences between groups
#' @param filename Output PDF path
#' @param ... Additional arguments passed onto \link[pheatmap]{pheatmap}
#' @return Saves a pdf image of the heatmap
#' @references \url{https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/}
smoothPheatmap <- function(bs.filtered.bsseq = bs.filtered.bsseq,
                            sigRegions = sigRegions,
                            testCovariate = testCovariate,
                            filename = "./DMRs/heatmap.pdf",
                            ...) {
  cat("\n[DMR heatmap] \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

  meth_mat <- bsseq::getMeth(BSseq = bs.filtered.bsseq,
                              regions = sigRegions,
                              type = "smooth",
                              what = "perRegion") %>%
    na.omit() %>%
    as.matrix()

  # scale = "row" below z-scores each region by dividing by its across-sample
  # SD. A region with zero variance (constant methylation across every
  # sample -- easy to get with a loose DSS p.threshold/delta, see params.R)
  # divides by zero and produces NaN/Inf, corrupting that row and often the
  # whole heatmap's color scale. Nothing upstream filters for this -- the
  # WGCNA script has an equivalent filter (filter_regions() in
  # wgcna_methylation.R) but it's local to that script and never applied
  # here, since this script builds sigRegions straight from the DMR CSV.
  row_sd <- apply(meth_mat, 1, sd, na.rm = TRUE)
  n_dropped <- sum(row_sd == 0 | is.na(row_sd))
  if (n_dropped > 0) {
    message(sprintf(
      "[DMR heatmap] Dropping %d/%d region(s) with zero variance across samples (would produce NaN/Inf under scale='row')",
      n_dropped, nrow(meth_mat)
    ))
  }
  meth_mat <- meth_mat[row_sd > 0 & !is.na(row_sd), , drop = FALSE]

  if (nrow(meth_mat) == 0) {
    stop("[DMR heatmap] No regions remain after filtering NA/zero-variance rows. ",
         "Check --dmrs-file -- if most DMRs have near-constant methylation, ",
         "the DSS/dmrseq calling thresholds (params.R) may be too loose.")
  }

  meth_mat %>%
    pheatmap::pheatmap(.,
                        scale = "row",
                        annotation_col = bs.filtered.bsseq %>%
                          pData() %>%
                          as.data.frame() %>%
                          dplyr::mutate(dplyr::across(dplyr::everything(), as.factor)) %>%
                          dplyr::select_if(~ nlevels(.) > 1),
                        color = RColorBrewer::brewer.pal(11, name = "RdBu") %>%
                          rev(),
                        show_colnames = FALSE,
                        border_color = "grey",
                        main = glue::glue("Z-Scores of {nrow(meth_mat)} Differentially Methylated Regions"),
                        fontsize = 16,
                        filename = filename,
                        width = 11,
                        height = 8.5,
                        annotation_colors = list2env(setNames(
                          list(bs.filtered.bsseq %>%
                                 pData() %>%
                                 as.data.frame() %>%
                                 purrr::pluck(testCovariate) %>%
                                 unique() %>%
                                 sort() %>%
                                 rev() %>%
                                 { setNames(RColorBrewer::brewer.pal(max(3, length(.)), "Set1")[seq_along(.)], .) }),
                          testCovariate
                        )) %>% as.list(),
                        ...
    ) %>%
    return()
}

message("Generating heatmap...")
smoothPheatmap(
  bs.filtered.bsseq = bs.filtered.bsseq,
  sigRegions = sigRegions,
  testCovariate = testCovariate,
  filename = output_file
)
message("Saved: ", output_file)