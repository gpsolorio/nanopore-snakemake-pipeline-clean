#!/usr/bin/env Rscript
##############################################################################
# annotate_dmr_context.R — Annotate DMR modules with gene, CpG, and genic
# context, then generate distribution plots.
#
# Called by Snakemake rule annotate_dmr_context.
##############################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
  library(tidyr)
  library(forcats)
  library(stringr)
  library(annotatr)
  library(AnnotationHub)
})

option_list <- list(
  make_option(c("-i", "--input-dir"), type = "character", default = NULL,
              help = "Directory containing a *_DMR_Modules.rds file (e.g. wgcna_modules output dir)"),
  make_option(c("-g", "--gtf-file"), type = "character", default = NULL,
              help = "GTF/GFF3 annotation file"),
  make_option(c("-o", "--output-dir"), type = "character", default = "."),
  make_option(c("--genome"), type = "character", default = "hg38"),
  make_option(c("--module-pattern"), type = "character", default = "DMR_Modules\\.rds$"),
  make_option(c("--modules_of_interest"), type = "character", default = "")
)

parse_modules <- function(value) {
  if (is.null(value) || identical(value, "")) {
    return(character(0))
  }
  unique(trimws(unlist(strsplit(value, ",", fixed = TRUE))))
}

args <- parse_args(OptionParser(option_list = option_list))

if (is.null(args$`input-dir`) || is.null(args$`gtf-file`)) {
  stop("--input-dir and --gtf-file are both required -- no default paths are ",
       "provided (earlier defaults pointed into another project's directory, ",
       "which risked silently annotating the wrong data).")
}

input_dir <- normalizePath(args$`input-dir`, winslash = "/", mustWork = TRUE)
gtf_file <- normalizePath(args$`gtf-file`, winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args$`output-dir`, winslash = "/", mustWork = FALSE)
genome <- args$genome
module_pattern <- args$`module-pattern`
modules_of_interest <- parse_modules(args$modules_of_interest)

# Offset used to translate AR_locus-local mini-test coordinates back to real
# chrX coordinates. Applies to BOTH the GTF gene-name overlap below and the
# annotatr CpG/genic-context call further down -- single source of truth so
# the two steps can't silently drift out of sync with each other.
# See mini-genome build doc: chrX:67544031-67730619 (make_mini_genome.sh).
AR_LOCUS_HG38_OFFSET <- 67544032

# The one and only synthetic contig name used by the AR-locus mini-test
# data. Translation is gated on this exact match -- NOT on "not chrX" --
# because real genome-wide data legitimately contains many chromosome
# names (chr1, chr2, ..., chrX, etc.) and must never be relabeled/shifted.
# Checking only regions$seqnames[1] is safe because a given DMR-modules
# run is either entirely mini-test (single synthetic contig) or entirely
# real genome-wide data -- never a mix of both.
MINI_TEST_SEQNAME <- "AR_locus"

tables_dir <- file.path(output_dir, "tables")
plots_dir <- file.path(output_dir, "plots")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

message("Libraries loaded successfully\n")
message("Input directory: ", input_dir)
message("Output directory: ", output_dir)

message("Loading DMR modules...")
module_file <- list.files(path = input_dir, pattern = module_pattern, full.names = TRUE)
if (length(module_file) != 1) {
  stop("Expected exactly one module RDS file matching pattern '", module_pattern, "', but found: ", length(module_file))
}

message("Using module file: ", module_file)
modules <- readRDS(module_file)
regions <- modules$regions %>%
  filter(module != "grey")

if (!"RegionID" %in% names(regions)) {
  regions$RegionID <- paste0(seqnames(regions), ":", start(regions), "-", end(regions))
}

message("Loaded ", nrow(regions), " regions from ", length(unique(regions$module)), " modules\n")

is_mini_test_data <- as.character(regions$seqnames[1]) == MINI_TEST_SEQNAME
message("Detected ", if (is_mini_test_data) "mini-test (AR_locus)" else "real genome-wide", " coordinate space\n")

message("Annotating regions with gene names from GTF...")
gtf_data <- import(gtf_file, format = "gtf")
gtf_transcripts <- gtf_data[gtf_data$type == "transcript"]

gr_regions <- GRanges(
  seqnames = regions$seqnames,
  ranges = IRanges(start = regions$start, end = regions$end)
)

# GTF gene coordinates are real genome coordinates. Our regions may be in
# AR_locus-local space (1-based from the mini genome window start) for
# mini-test data, so translate a *copy* to real chrX coordinates just for
# this overlap call -- ONLY when the data is actually mini-test. Real
# genome-wide regions (any chromosome, including chrX) are left untouched.
gr_regions_real <- gr_regions
if (is_mini_test_data) {
  seqlevels(gr_regions_real) <- "chrX"
  seqnames(gr_regions_real) <- "chrX"
  ranges(gr_regions_real) <- IRanges(
    start = start(gr_regions) + AR_LOCUS_HG38_OFFSET - 1,
    end   = end(gr_regions) + AR_LOCUS_HG38_OFFSET - 1
  )
}

gr_genes <- GRanges(
  seqnames = seqnames(gtf_transcripts),
  ranges = IRanges(start = start(gtf_transcripts), end = end(gtf_transcripts)),
  gene_name = mcols(gtf_transcripts)$gene_name
)

overlaps <- findOverlaps(gr_regions_real, gr_genes)
regions$gene_name <- NA_character_
regions$gene_name[queryHits(overlaps)] <- gr_genes$gene_name[subjectHits(overlaps)]

message("Found gene annotations for ", sum(!is.na(regions$gene_name)), " regions\n")
write.csv(regions, file = file.path(tables_dir, "Modules_Regions_Genes.csv"), row.names = FALSE)
message("Saved: Modules_Regions_Genes.csv\n")

module_genes <- regions %>%
  group_by(module) %>%
  summarize(
    unique_genes = list(unique(gene_name[!is.na(gene_name)])),
    num_genes = length(unique(gene_name[!is.na(gene_name)])),
    .groups = "drop"
  )

module_genes_out <- module_genes %>%
  mutate(unique_genes = vapply(unique_genes, paste, collapse = ";", FUN.VALUE = character(1)))

write.csv(module_genes_out, file = file.path(tables_dir, "Module_Gene_Summary.csv"), row.names = FALSE)
print(module_genes %>% dplyr::select(module, num_genes))

annotations <- c(
  paste0(genome, "_cpgs"),
  paste0(genome, "_basicgenes"),
  paste0(genome, "_genes_intergenic"),
  paste0(genome, "_enhancers_fantom")
)

message("\nAnnotating regions with CpG and genic context...")

# annotatr's hg38 annotation sets use real genome coordinates. Same
# is_mini_test_data gate as above -- translate only for AR_locus mini-test
# data, leave real genome-wide regions (any chromosome) untouched.

gr <- GRanges(
  seqnames = regions$seqnames,
  ranges = IRanges(start = regions$start, end = regions$end),
  RegionID = regions$RegionID
)

gr_real <- gr
if (is_mini_test_data) {
  seqlevels(gr_real) <- "chrX"
  seqnames(gr_real) <- "chrX"
  ranges(gr_real) <- IRanges(
    start = start(gr) + AR_LOCUS_HG38_OFFSET - 1,
    end   = end(gr) + AR_LOCUS_HG38_OFFSET - 1
  )
}

annots <- build_annotations(genome = genome, annotations = annotations)
annotated <- annotate_regions(
  regions = gr_real,
  annotations = annots,
  ignore.strand = TRUE,
  quiet = TRUE
)

annotated_df <- as.data.frame(annotated)

gene_info <- annotated_df %>%
  filter(!is.na(annot.symbol)) %>%
  dplyr::select(RegionID, gene_symbol = annot.symbol, gene_id = annot.gene_id, tx_id = annot.tx_id, annot.type) %>%
  distinct()

cpg_context <- annotated_df %>%
  filter(grepl("cpg", annot.type, ignore.case = TRUE)) %>%
  dplyr::select(RegionID, CpG_context = annot.type) %>%
  distinct()

gene_context <- annotated_df %>%
  filter(grepl("gene", annot.type, ignore.case = TRUE)) %>%
  dplyr::select(RegionID, gene_context = annot.type) %>%
  distinct()

regions_annotated <- regions %>%
  left_join(gene_info, by = "RegionID") %>%
  left_join(gene_context, by = "RegionID") %>%
  left_join(cpg_context, by = "RegionID")

write.csv(
  regions_annotated,
  file = file.path(tables_dir, "Regions_with_Gene_and_CpG_Annotations.csv"),
  row.names = FALSE
)
message("Saved: Regions_with_Gene_and_CpG_Annotations.csv\n")

message("Generating plots...")
regions_clean <- regions_annotated %>%
  mutate(
    CpG_context = str_replace_all(CpG_context, paste0(genome, "_cpg_"), ""),
    gene_context = str_replace_all(gene_context, paste0(genome, "_genes_"), "")
  )

plot_context_distribution <- function(data, context_col, title, filename) {
  plot_data <- data %>%
    count(module, .data[[context_col]]) %>%
    group_by(module) %>%
    mutate(percentage = 100 * n / sum(n)) %>%
    ungroup()

  p <- ggplot(plot_data, aes(x = fct_reorder(module, -percentage), y = percentage, fill = .data[[context_col]])) +
    geom_bar(stat = "identity", position = "fill") +
    labs(title = title, x = "Module", y = "Proportion") +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(filename, plot = p, width = 10, height = 6)
  message("Saved: ", filename)
}

plot_context_distribution(
  regions_clean, "CpG_context",
  "CpG Context Distribution by Module",
  file.path(plots_dir, "CpG_Context_Distribution_by_Module.pdf")
)

plot_context_distribution(
  regions_clean, "gene_context",
  "Genic Context Distribution by Module",
  file.path(plots_dir, "Genic_Context_Distribution_by_Module.pdf")
)

if (length(modules_of_interest) > 0) {
  filtered_regions <- regions_clean %>%
    filter(module %in% modules_of_interest)

  plot_context_distribution(
    filtered_regions, "CpG_context",
    "CpG Context Distribution (Selected Modules)",
    file.path(plots_dir, "CpG_Context_Distribution_Selected_Modules.pdf")
  )

  plot_context_distribution(
    filtered_regions, "gene_context",
    "Genic Context Distribution (Selected Modules)",
    file.path(plots_dir, "Genic_Context_Distribution_Selected_Modules.pdf")
  )
}

message("\n=== Annotation complete ===")