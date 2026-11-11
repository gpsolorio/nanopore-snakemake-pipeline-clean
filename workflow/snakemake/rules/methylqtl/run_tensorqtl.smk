##############################################################################
# cis-mQTL mapping with tensorQTL (Python API, two-pass)
#
# Pass 1 — permutation: one row per DMR, empirical p-value + q-value (FDR)
# Pass 2 — nominal:     all SNP–DMR pairs for FDR-significant DMRs only
#
# Calls: workflow/scripts/methylqtl/run_tensorqtl.py
#
# NOTE: run_tensorqtl.py computes q-values via tensorqtl.post.calculate_qvalues,
# which uses rpy2 to call R's `qvalue` package. tensorqtl-env.yml MUST include
# r-base, rpy2, and bioconductor-qvalue, or this rule will fail at runtime.
##############################################################################

rule run_tensorqtl:
    input:
        pgen       = cohort("methylqtl", "pgen", "genotypes_qc.pgen"),
        pvar       = cohort("methylqtl", "pgen", "genotypes_qc.pvar"),
        psam       = cohort("methylqtl", "pgen", "genotypes_qc.psam"),
        bed        = cohort("dmrs", "smoothed", "DMR_individual_smoothed_methylation.bed.gz"),
        covariates = cohort("methylqtl", "pgen", "tensorqtl_covariates.txt"),
    output:
        perm_results = cohort("methylqtl", "tensorqtl", "DMR_cisQTL.cis_qtl_results.txt"),
        nom_pairs    = cohort("methylqtl", "tensorqtl", "DMR_cisQTL.significant_dmr_pairs.txt"),
    log:
        cohort_log("tensorqtl.log")
    params:
        plink_prefix = cohort("methylqtl", "pgen", "genotypes_qc"),
        out_dir      = cohort("methylqtl", "tensorqtl"),
        prefix       = "DMR_cisQTL",
        cis_window   = config["tensorqtl"].get("cis_window", 1_000_000),
        fdr          = config["tensorqtl"].get("fdr", 0.05),
        seed         = config["tensorqtl"].get("seed", 123456),
        scripts_dir  = str(SCRIPTS_DIR),
    threads:
        res("run_tensorqtl", "threads")
    resources:
        mem_gb    = res("run_tensorqtl", "mem_gb"),
        time      = res("run_tensorqtl", "time"),
        account   = res("run_tensorqtl", "account"),
        partition = res("run_tensorqtl", "partition"),
        gpus      = res("run_tensorqtl", "gpus", 0),
    conda:
        str(ENVS_DIR / "tensorqtl-env.yml")
    shell:
        SHELL_PRELUDE + """
        export KMP_DUPLICATE_LIB_OK=TRUE
        export R_HOME=$("$CONDA_PREFIX/bin/R" RHOME)

        mkdir -p {params.out_dir}

        $CONDA_PREFIX/bin/python {params.scripts_dir}/methylqtl/run_tensorqtl.py \\
            --plink-prefix  {params.plink_prefix} \\
            --bed           {input.bed} \\
            --covariates    {input.covariates} \\
            --out-dir       {params.out_dir} \\
            --prefix        {params.prefix} \\
            --cis-window    {params.cis_window} \\
            --fdr           {params.fdr} \\
            --seed          {params.seed} \\
            --perm-out      {output.perm_results} \\
            --nom-out       {output.nom_pairs}
        """ + SHELL_EPILOG