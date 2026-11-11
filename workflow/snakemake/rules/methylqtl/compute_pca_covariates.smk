##############################################################################
# PCA on merged genotypes → covariates file for tensorQTL
#
# Output covariates TSV format (required by tensorQTL):
#   - Rows: SEX, PC1 … PC{n_pcs}  (and ChildAgeMos if sample_metadata_csv provided)
#   - Cols: sample IDs
#   - Top-left cell is blank (leading tab)
#
# n_pcs and the approx/exact PCA algorithm are both resolved from the actual
# sample count at runtime (see _resolve_n_pcs / _resolve_pca_flag below)
# rather than a fixed config value -- tensorqtl.n_pcs is a *ceiling*, not a
# guarantee, since you can never compute more than N_samples - 1 PCs.
#
# If you already have a hand-crafted covariates file, point
# tensorqtl.covariates_tsv at it in config and skip this rule.
##############################################################################


def _count_samples(psam_path):
    """Count non-header samples in a .psam file."""
    n = 0
    with open(psam_path) as fh:
        for line in fh:
            if not line.startswith("#"):
                n += 1
    return n


def _resolve_n_pcs(psam_path, configured_n_pcs):
    """Clamp requested PC count to what's actually possible for N samples."""
    n_samples = _count_samples(psam_path)
    max_pcs = max(n_samples - 1, 1)
    return min(configured_n_pcs, max_pcs)


def _resolve_pca_flag(psam_path, approx_min_n):
    """
    plink2's randomized/approximate PCA is designed for large-N biobank data
    and errors out (or is simply the wrong algorithm) for small cohorts.
    Only request 'approx' once N crosses the configured threshold.
    """
    n_samples = _count_samples(psam_path)
    return "approx" if n_samples >= approx_min_n else ""


rule compute_pca_covariates:
    input:
        pgen     = cohort("methylqtl", "pgen", "genotypes_qc.pgen"),
        pvar     = cohort("methylqtl", "pgen", "genotypes_qc.pvar"),
        psam     = cohort("methylqtl", "pgen", "genotypes_qc.psam"),
        metadata = SAMPLE_METADATA_TSV,
    output:
        eigenvec   = cohort("methylqtl", "pgen", "genotypes_qc.eigenvec"),
        eigenval   = cohort("methylqtl", "pgen", "genotypes_qc.eigenval"),
        covariates = cohort("methylqtl", "pgen", "tensorqtl_covariates.txt"),
    log:
        cohort_log("pca_covariates.log")
    params:
        plink2 = PLINK2_BIN,
        plink_prefix = lambda wc: cohort("methylqtl", "pgen", "genotypes_qc"),
        n_pcs        = lambda wc, input: _resolve_n_pcs(
            input.psam, config["tensorqtl"].get("n_pcs", 5)
        ),
        pca_flag     = lambda wc, input: _resolve_pca_flag(
            input.psam, config["tensorqtl"].get("pca_approx_min_n", 50)
        ),
        sample_meta  = config["tensorqtl"].get("sample_metadata_csv", ""),
        scripts_dir  = str(SCRIPTS_DIR),
    threads:
        res("compute_pca_covariates", "threads")
    resources:
        mem_gb    = res("compute_pca_covariates", "mem_gb"),
        time      = res("compute_pca_covariates", "time"),
        account   = res("compute_pca_covariates", "account"),
        partition = res("compute_pca_covariates", "partition"),
    conda:
        plink2_env_yml()
    shell:
        SHELL_PRELUDE + """

        # --new-id-max-allele-len: allow long REF/ALT alleles in variant IDs;
        # falls back to '.' if an allele exceeds 100 bp
        
        {params.plink2} \\
            --pfile {params.plink_prefix} \\
            --set-all-var-ids '@:#:$r:$a' \\
            --new-id-max-allele-len 100 missing \\
            --make-pgen \\
            --out {params.plink_prefix}_vid \\
            --threads {threads}
        
        {params.plink2} \\
            --pfile {params.plink_prefix}_vid \\
            --freq \\
            --out {params.plink_prefix}_vid \\
            --threads {threads}
        
        {params.plink2} \\
            --pfile {params.plink_prefix}_vid \\
            --read-freq {params.plink_prefix}_vid.afreq \\
            --pca {params.n_pcs} {params.pca_flag} \\
            --out {params.plink_prefix} \\
            --threads {threads}
        
        $CONDA_PREFIX/bin/python {params.scripts_dir}/methylqtl/format_pca_covariates.py \\
            --eigenvec {output.eigenvec} \\
            --n-pcs    {params.n_pcs} \\
            --meta     "{params.sample_meta}" \\
            --metadata {input.metadata} \\
            --out      {output.covariates}
        """ + SHELL_EPILOG
