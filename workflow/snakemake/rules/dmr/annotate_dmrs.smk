##############################################################################
# DMR Context Annotation
##############################################################################

DMR_ANNOTATION_SCRIPT = (SCRIPTS_DIR / "analysis" / "differential_methylation" / "05_annotations" / "annotate_dmr_context.R").resolve()


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


def _dmr_input_dir(wc):
    raw = _require_analysis_key("dmr_input_dir", "DMR input directory")
    return raw.format(barcode=wc.barcode) if "{barcode}" in raw else raw


rule annotate_dmr_context:
    output:
        outdir = directory(results("{barcode}", "dmr_annotations"))
    log:
        get_log_path("{barcode}", "08_annotate_dmr_context.log")
    params:
        input_dir = _dmr_input_dir,
        gtf_file  = lambda wc: _require_analysis_key("gtf_file", "GTF annotation file"),
        genome    = lambda wc: config.get("analysis", {}).get("genome", "hg38"),
        modules_of_interest_arg = lambda wc: (
            f"--modules_of_interest {config.get('analysis', {}).get('modules_of_interest')}"
            if config.get("analysis", {}).get("modules_of_interest")
            else ""
        ),
        script = str(DMR_ANNOTATION_SCRIPT)
    threads:
        config["resources"]["annotate_dmrs"]["threads"]
    resources:
        mem_gb    = config["resources"]["annotate_dmrs"]["mem_gb"],
        time      = config["resources"]["annotate_dmrs"]["time"],
        account   = config["resources"]["annotate_dmrs"]["account"],
        partition = config["resources"]["annotate_dmrs"]["partition"],
        gpus      = config["resources"]["annotate_dmrs"].get("gpus", 0),
        slurm_log_dir = lambda wc: get_slurm_log_dir(wc)
    conda:
        str(ENVS_DIR / "dss-env.yml")
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1

        START=$(date +%s)
        START_FMT=$(date "+%Y-%m-%d %H:%M:%S")
        echo "Starting annotate_dmr_context for {wildcards.barcode} - $START_FMT"

        mkdir -p {output.outdir}

        "$CONDA_PREFIX/bin/Rscript" {params.script} \
            --input-dir {params.input_dir} \
            --gtf-file {params.gtf_file} \
            --output-dir {output.outdir} \
            --genome {params.genome} \
            {params.modules_of_interest_arg}

        END=$(date +%s)
        END_FMT=$(date "+%Y-%m-%d %H:%M:%S")
        ELAPSED=$((END - START))
        echo "Finished annotate_dmr_context for {wildcards.barcode} - $END_FMT"
        echo "Elapsed time: $(((ELAPSED%604800)/86400))d $(((ELAPSED%86400)/3600))h $(((ELAPSED%3600)/60))m $((ELAPSED%60))s"
        """