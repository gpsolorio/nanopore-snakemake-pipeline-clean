#!/usr/bin/env Rscript
##############################################################################
# pathway_enrichment.R — Pathway Enrichment Analysis for DMR-Associated Genes
#
# Performs enrichment analysis on genes associated with DMRs using EnrichR
# and generates plots for top enriched terms.
#
# Usage:
#   Rscript pathway_enrichment.R --input-file Modules_Regions_Genes.csv \
#       --output-dir EnrichR [--top-n-terms 25]
##############################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(enrichR)
  library(dplyr)
  library(openxlsx)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-i", "--input-file"), type = "character", default = NULL,
              help = "Modules_Regions_Genes.csv (or similar) with a gene_name column -- typically annotate_dmr_context.R's output"),
  make_option(c("-o", "--output-dir"), type = "character", default = "EnrichR"),
  make_option(c("--top-n-terms"), type = "integer", default = 25)
)

args <- parse_args(OptionParser(option_list = option_list))

if (is.null(args$`input-file`)) {
  stop("--input-file is required (typically tables/Modules_Regions_Genes.csv from annotate_dmr_context.R).")
}

input_file <- normalizePath(args$`input-file`, mustWork = TRUE)
output_dir <- args$`output-dir`
top_n_terms <- args$`top-n-terms`

message("Libraries loaded successfully\n")

ENRICHR_LIBRARIES <- c(
  "GO_Biological_Process_2025",
  "GO_Cellular_Component_2025",
  "GO_Molecular_Function_2025",
  "KEGG_2021_Human",
  "Panther_2016",
  "Reactome_2022",
  "SynGO_2022"
)

setEnrichrSite("Enrichr")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("Loading DMR data from: ", input_file)
dmrs <- read.csv(input_file)

genes <- unique(na.omit(dmrs$gene_name))
message("Number of unique genes for enrichment: ", length(genes), "\n")

if (length(genes) == 0) {
  stop("No genes found in input file. Check that 'gene_name' column exists and has values.")
}

message("Running EnrichR analysis...")
results <- enrichr(genes, ENRICHR_LIBRARIES)
message("Enrichment analysis complete\n")

message("Saving results to Excel...")
wb <- createWorkbook()

for (lib in names(results)) {
  df <- results[[lib]]
  addWorksheet(wb, lib)
  writeData(wb, lib, df)
}

excel_file <- file.path(output_dir, "AllDMRs_Enrichr_Results.xlsx")
saveWorkbook(wb, file = excel_file, overwrite = TRUE)
message("Saved: ", excel_file, "\n")

message("Generating plots for top ", top_n_terms, " terms per library...\n")

for (lib in names(results)) {
  df <- results[[lib]] %>%
    arrange(P.value) %>%
    head(top_n_terms)

  if (nrow(df) > 0) {
    p <- ggplot(df, aes(x = reorder(Term, -log10(P.value)),
                        y = -log10(P.value))) +
      geom_col(fill = "steelblue") +
      coord_flip() +
      labs(
        title = paste("Top", min(top_n_terms, nrow(df)), "Terms -", lib),
        x = "",
        y = "-log10(P-value)"
      ) +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(hjust = 0.5))

    safe_lib_name <- gsub("[^A-Za-z0-9]", "_", lib)
    plot_file <- file.path(output_dir, paste0("AllDMRs_", safe_lib_name, ".pdf"))

    ggsave(filename = plot_file, plot = p, width = 10, height = 6)
    message("  Saved: ", basename(plot_file))
  } else {
    message("  No significant terms found for ", lib)
  }
}

message("\n=== Enrichment analysis complete ===")
message("Results saved in: ", normalizePath(output_dir))
