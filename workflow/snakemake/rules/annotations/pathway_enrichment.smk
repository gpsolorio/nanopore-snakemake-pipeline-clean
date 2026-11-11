##############################################################################
# Pathway Enrichment (EnrichR)
##############################################################################

PATHWAY_ENRICHMENT_SCRIPT = (REPO_ROOT / "analysis" / "dmr" / "04_annotate_dmrs" / "pathway_enrichment.R").resolve()


rule pathway_enrichment:
    input:
        # annotate_dmr_context.R writes this file inside its directory()
        # output at tables/Modules_Regions_Genes.csv -- referencing the
        # concrete file (not just the parent directory) keeps Snakemake's
        # dependency tracking correct.
        genes_csv = cohort("dmrs", "dmr_annotations", "tables", "Modules_Regions_Genes.csv")
    output:
        outdir = directory(cohort("dmrs", "enrichr"))
    log:
        cohort_log("pathway_enrichment.log")
    params:
        top_n_terms = config.get("analysis", {}).get("enrichr_top_n_terms", 25),
        script = str(PATHWAY_ENRICHMENT_SCRIPT)
    threads:
        res("pathway_enrichment", "threads")
    resources:
        mem_gb    = res("pathway_enrichment", "mem_gb"),
        time      = res("pathway_enrichment", "time"),
        account   = res("pathway_enrichment", "account"),
        partition = res("pathway_enrichment", "partition")
    conda:
        str(ENVS_DIR / "enrichment-env.yml")
    shell:
        SHELL_PRELUDE + """
        mkdir -p {output.outdir}

        "$CONDA_PREFIX/bin/Rscript" {params.script} \
            --input-file {input.genes_csv} \
            --output-dir {output.outdir} \
            --top-n-terms {params.top_n_terms}
        """ + SHELL_EPILOG
