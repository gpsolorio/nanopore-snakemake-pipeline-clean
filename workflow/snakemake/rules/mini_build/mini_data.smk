##############################################################################
# Mini-Data AR-Locus Extraction & Aligner Comparison
#
# Wraps extract_ar_reads.sh / realign_to_mini.sh as real Snakemake rules so
# this fits into the normal DAG (dry-runnable, resumable, resource-aware)
# instead of the standalone run_all_samples.sh loop.
#
# Output split, intentionally:
#   _MINI_STATIC_DIR   -- workflow/new-mini-data/mini-data/{genome,beds}
#                         Manually built once (make_mini_genome.sh etc.),
#                         includes the meryl k-mer index. Expensive to
#                         rebuild, NOT wiped between test runs.
#   _MINI_RESULTS_DIR  -- results/mini-test/{bam,aligned,logs}
#                         Fully Snakemake-managed, disposable. Safe to
#                         `rm -rf results/mini-test` for a clean pipeline
#                         smoke test without touching the genome/meryl index.
##############################################################################

if "mini_data" not in config:
    raise ValueError(
        "Missing required config key: mini_data\n"
        "Add a `mini_data:` section to your config YAML with base_dir, region, "
        "suffix, and samples (sample_id -> proband dirname) before including "
        "rules/mini_data.smk."
    )

_MINI_CFG       = config["mini_data"]
_MINI_BASE_DIR  = _MINI_CFG["base_dir"]
_MINI_REGION    = _MINI_CFG["region"]
_MINI_SUFFIX    = _MINI_CFG["suffix"]
_MINI_SAMPLES   = _MINI_CFG["samples"]   # dict: sample_id -> proband dirname

_MINI_SCRIPTS_DIR = (REPO_ROOT / "workflow" / "scripts" / "mini_genomes" / "ar_test").resolve()
_MINI_STATIC_DIR  = (REPO_ROOT / "resources" / "mini_ar_test").resolve()
_MINI_RESULTS_DIR = (REPO_ROOT / "results" / "mini-test").resolve()
_MINI_LOGS_DIR    = (_MINI_RESULTS_DIR / "logs").resolve()
_MINI_GENOME_DIR  = (_MINI_STATIC_DIR / "genome").resolve()


def _mini_pass_dir(wildcards):
    proband_dirname = _MINI_SAMPLES[wildcards.sample_id]
    return f"{_MINI_BASE_DIR}/{proband_dirname}/{wildcards.sample_id}_5mC-5hmC/pass"


rule extract_ar_reads:
    # NOTE: pass_dir is intentionally a `params`, not an `input` -- it's
    # thousands of read-only chunk BAMs on external storage. Declaring it as
    # a Snakemake input would force a timestamp check across every chunk on
    # every dry run, which is slow and pointless since this raw data is
    # immutable once sequenced.
    output:
        bam = str(_MINI_RESULTS_DIR / "bam") + "/{sample_id}_" + _MINI_SUFFIX + "_sorted.bam",
        bai = str(_MINI_RESULTS_DIR / "bam") + "/{sample_id}_" + _MINI_SUFFIX + "_sorted.bam.bai",
    params:
        pass_dir     = _mini_pass_dir,
        region       = _MINI_REGION,
        label        = lambda wc: f"{wc.sample_id}_{_MINI_SUFFIX}",
        script       = str(_MINI_SCRIPTS_DIR / "extract_ar_reads.sh"),
        logs_dir     = str(_MINI_LOGS_DIR),
        mini_out_root = str(_MINI_RESULTS_DIR),
        mini_genome_dir = str(_MINI_GENOME_DIR),
    log:
        str(_MINI_LOGS_DIR / "extract_{sample_id}.log")
    threads:
        res("extract_ar_reads", "threads")
    resources:
        mem_gb    = res("extract_ar_reads", "mem_gb"),
        time      = res("extract_ar_reads", "time"),
        account   = res("extract_ar_reads", "account"),
        partition = res("extract_ar_reads", "partition"),
    conda:
        str(ENVS_DIR / "ngs-utils-env.yml")
    shell:
        SHELL_PRELUDE + """
        export MINI_OUT_ROOT={params.mini_out_root:q}
        export MINI_GENOME_DIR={params.mini_genome_dir:q}
        mkdir -p {params.logs_dir}
        {params.script} {params.pass_dir:q} {params.region:q} {params.label:q}
        """ + SHELL_EPILOG


rule realign_to_mini:
    input:
        bam = rules.extract_ar_reads.output.bam,
        bai = rules.extract_ar_reads.output.bai,
        genome = str(_MINI_GENOME_DIR / "mini_genome.fa"),
    output:
        minimap2      = str(_MINI_RESULTS_DIR / "aligned") + "/{sample_id}_" + _MINI_SUFFIX + "/minimap2.bam",
        minimap2_bai  = str(_MINI_RESULTS_DIR / "aligned") + "/{sample_id}_" + _MINI_SUFFIX + "/minimap2.bam.bai",
        winnowmap     = str(_MINI_RESULTS_DIR / "aligned") + "/{sample_id}_" + _MINI_SUFFIX + "/winnowmap2.bam",
        winnowmap_bai = str(_MINI_RESULTS_DIR / "aligned") + "/{sample_id}_" + _MINI_SUFFIX + "/winnowmap2.bam.bai",
    params:
        label         = lambda wc: f"{wc.sample_id}_{_MINI_SUFFIX}",
        script        = str(_MINI_SCRIPTS_DIR / "realign_to_mini.sh"),
        logs_dir      = str(_MINI_LOGS_DIR),
        mini_out_root = str(_MINI_RESULTS_DIR),
        mini_genome_dir = str(_MINI_GENOME_DIR),
    log:
        str(_MINI_LOGS_DIR / "realign_{sample_id}.log")
    threads:
        res("realign_to_mini", "threads")
    resources:
        mem_gb    = res("realign_to_mini", "mem_gb"),
        time      = res("realign_to_mini", "time"),
        account   = res("realign_to_mini", "account"),
        partition = res("realign_to_mini", "partition"),
    conda:
        str(ENVS_DIR / "ngs-utils-env.yml")
    shell:
        SHELL_PRELUDE + """
        export MINI_OUT_ROOT={params.mini_out_root:q}
        export MINI_GENOME_DIR={params.mini_genome_dir:q}
        mkdir -p {params.logs_dir}
        {params.script} {params.label:q}
        """ + SHELL_EPILOG


rule all_mini_aligned:
    input:
        expand(
            str(_MINI_RESULTS_DIR / "aligned") + "/{sample_id}_" + _MINI_SUFFIX + "/winnowmap2.bam",
            sample_id=_MINI_SAMPLES.keys(),
        )
