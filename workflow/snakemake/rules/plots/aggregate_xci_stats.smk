##############################################################################
# Cohort-level XCI Ratio Summary
##############################################################################

_XCI_AGGREGATE_SCRIPT = (WORKFLOW_ROOT / "scripts" / "analysis" / "aggregate_xci_summary.py").resolve()

rule aggregate_xci_summary:
    input:
        # plot_cpg_hists declares directory(<sample>/histograms) as its
        # output, so Snakemake can't see files INSIDE that directory.
        # Depend on the directories themselves (which a rule does produce)
        # and pass the per-sample TSV paths via params below.
        hist_dirs = expand(
            sample("{barcode}", "histograms"),
            barcode=BARCODES,
        ),
        sample_metadata = SAMPLE_METADATA_TSV,
    output:
        xlsx      = cohort("summaries", "xci_population_summary.xlsx"),
        plots_dir = directory(cohort("summaries", "xci_population_plots")),
    log:
        cohort_log("aggregate_xci_summary.log")
    conda:
        str(ENVS_DIR / "ngs-utils-env.yml")   # confirm openpyxl is in this env -- see note below
    params:
        script        = str(_XCI_AGGREGATE_SCRIPT),
        diagnosis_col = DMR_GROUP_COLUMN,
        xci_tsvs = expand(
            sample("{barcode}", "histograms", "5mC", "tsvs", "xci_summary_5mC.tsv"),
            barcode=BARCODES,
        ),
    threads:
        res("aggregate_xci_summary", "threads", 1)
    resources:
        mem_gb    = res("aggregate_xci_summary", "mem_gb", 8),
        time      = res("aggregate_xci_summary", "time", "01:00:00"),
        account   = res("aggregate_xci_summary", "account"),
        partition = res("aggregate_xci_summary", "partition"),
    shell:
        SHELL_PRELUDE + """
        # Fail with a clear message if plot_cpg_hists finished but didn't
        # write the expected TSV for some sample.
        for f in {params.xci_tsvs}; do
            [ -s "$f" ] || {{ echo "ERROR: missing or empty XCI summary TSV: $f" >&2; exit 1; }}
        done

        mkdir -p {output.plots_dir}
        python {params.script} \\
            --xci-summary-tsvs {params.xci_tsvs} \\
            --sample-metadata {input.sample_metadata} \\
            --diagnosis-col {params.diagnosis_col:q} \\
            --output-xlsx {output.xlsx} \\
            --output-plots-dir {output.plots_dir}
        """ + SHELL_EPILOG
