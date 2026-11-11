##############################################################################
# DMR Context Annotation
##############################################################################

DMR_ANNOTATION_SCRIPT = (REPO_ROOT / "analysis" / "dmr" / "04_annotate_dmrs" / "annotate_dmr_context.R").resolve()


def _require_analysis_key(key, label):
    val = config.get("analysis", {}).get(key)
    if not val:
        raise ValueError(
            f"Missing required config key: analysis.{key}\n"
            f"({label} -- no default is provided for this key. Earlier "
            f"defaults pointed at other researchers' project directories "
            f"(other cohorts), which risked silently analyzing the "
            f"wrong data. Set analysis.{key} explicitly in your config.)"
        )
    return val


rule annotate_dmr_context:
    input:
        module_file = rules.wgcna_modules.output.modules_rds
    output:
        genes_csv        = cohort("dmrs", "dmr_annotations", "tables", "Modules_Regions_Genes.csv"),
        gene_summary_csv = cohort("dmrs", "dmr_annotations", "tables", "Module_Gene_Summary.csv"),
        annotated_csv    = cohort("dmrs", "dmr_annotations", "tables", "Regions_with_Gene_and_CpG_Annotations.csv"),
        cpg_plot         = cohort("dmrs", "dmr_annotations", "plots", "CpG_Context_Distribution_by_Module.pdf"),
        gene_plot        = cohort("dmrs", "dmr_annotations", "plots", "Genic_Context_Distribution_by_Module.pdf"),
        outdir           = directory(cohort("dmrs", "dmr_annotations"))
    log:
        cohort_log("annotate_dmr_context.log")
    params:
        gtf_file = lambda wc: _require_analysis_key("gtf_file", "GTF annotation file"),
        genome   = lambda wc: config.get("analysis", {}).get("genome", "hg38"),
        modules_of_interest_arg = lambda wc: (
            f"--modules_of_interest {config.get('analysis', {}).get('modules_of_interest')}"
            if config.get("analysis", {}).get("modules_of_interest")
            else ""
        ),
        script = str(DMR_ANNOTATION_SCRIPT)
    threads:
        res("annotate_dmr_context", "threads")
    resources:
        mem_gb    = res("annotate_dmr_context", "mem_gb"),
        time      = res("annotate_dmr_context", "time"),
        account   = res("annotate_dmr_context", "account"),
        partition = res("annotate_dmr_context", "partition")
    conda:
        str(ENVS_DIR / "dmr-annotation-env.yml")
    shell:
        SHELL_PRELUDE + """
        mkdir -p {output.outdir}

        "$CONDA_PREFIX/bin/Rscript" {params.script} \
            --input-dir $(dirname {input.module_file}) \
            --gtf-file {params.gtf_file} \
            --output-dir {output.outdir} \
            --genome {params.genome} \
            {params.modules_of_interest_arg}
        """ + SHELL_EPILOG
        