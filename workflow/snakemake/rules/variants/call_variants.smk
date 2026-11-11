##############################################################################
# Variant Calling + Phasing (Clair3 via Docker, Singularity, or Apptainer)
##############################################################################

def _call_variants_bam(wc):
    if USE_PREBUILT_ALIGNMENTS:
        return str(mini_aligned_root() / f"{wc.barcode}_AR" / "minimap2.bam")
    return sample(wc.barcode, "bam", f"{wc.barcode}_merged_sorted.bam")

def _call_variants_bai(wc):
    if USE_PREBUILT_ALIGNMENTS:
        return str(mini_aligned_root() / f"{wc.barcode}_AR" / "minimap2.bam.bai")
    return sample(wc.barcode, "bam", f"{wc.barcode}_merged_sorted.bam.bai")

def _call_variants_outputs():
    """
    Build the output dict for call_variants. Phased outputs are only
    expected when tools.clair3.enable_phasing is true, since Clair3 only
    produces phased_merge_output.vcf.gz / phased_output.bam when
    --enable_phasing is actually passed on the command line (see
    phasing_flags below). Declaring them unconditionally causes a
    MissingOutputException whenever phasing is disabled (e.g. the mini
    smoke-test config).
    """
    outputs = {
        "vcf":       sample("{barcode}", "variants", "merge_output.vcf.gz"),
        "vcf_index": sample("{barcode}", "variants", "merge_output.vcf.gz.tbi"),
    }
    if CLAIR3_CFG.get("enable_phasing", True):
        outputs.update({
            "phased_vcf":       sample("{barcode}", "variants", "phased_merge_output.vcf.gz"),
            "phased_vcf_index": sample("{barcode}", "variants", "phased_merge_output.vcf.gz.tbi"),
            "phased_bam":       temp(sample("{barcode}", "variants", "phased_output.bam")),
            "phased_bam_index": temp(sample("{barcode}", "variants", "phased_output.bam.bai")),
        })
    return outputs


rule call_variants:
    input:
        bam = _call_variants_bam,
        bai = _call_variants_bai,
        ref = str(REFERENCE_GENOME)
    output:
        **_call_variants_outputs()
    log:
        sample_log("{barcode}", "call_variants.log")
    params:
        bam_dir                 = lambda wc, input: str(Path(input.bam).parent),
        out_dir                 = lambda wc: sample(wc.barcode, "variants"),
        model                   = CLAIR3_CFG["model_path"],
        platform                = CLAIR3_CFG["platform"],
        container_image         = _CLAIR3_CONTAINER,
        singularity_ref         = singularity_ref(_CLAIR3_CONTAINER),
        container_engine        = _CONTAINER_ENGINE,    # for logging only
        container_engine_bin    = CONTAINER_ENGINE_BIN, # actual command to run
        is_sif_engine           = _CONTAINER_ENGINE in SIF_ENGINES,
        workflow_root           = str(WORKFLOW_ROOT),
        sample_name             = lambda wc: wc.barcode,

        # Bind args are built once in common.smk from every config-driven
        # path (reference genome, annotations, GTF, BAM/POD5 index dirs,
        # repo dirs) so this rule doesn't need to know or guess which
        # filesystem the reference genome happens to live on -- fixes the
        # "Reference file not found" failure that occurred whenever the
        # reference lived outside the previously-hardcoded four bind
        # roots (workflow/results/resources/examples).
        bind_args = SINGULARITY_BIND_ARGS if _CONTAINER_ENGINE in SIF_ENGINES else DOCKER_BIND_ARGS,

        # Flags
        gpu_flag         = "--use_gpu" if RUN_CLAIR3_USE_GPU else "",
        nv_flag          = "--nv" if RUN_CLAIR3_USE_GPU else "",
        ctg_flag = (f"--ctg_name {CLAIR3_CFG['ctg_name']}"
            if GENOME_MODE == "mini" and CLAIR3_CFG.get("ctg_name")
            else "--include_all_ctgs"),
        phasing_flags = (
            "--enable_phasing --use_whatshap_for_final_output_phasing "
            "--use_whatshap_for_final_output_haplotagging"
            if CLAIR3_CFG.get("enable_phasing", True) else ""
        ),
    threads:
        res("call_variants", "threads")
    resources:
        mem_gb    = res("call_variants", "mem_gb"),
        time      = res("call_variants", "time"),
        account   = res("call_variants", "account"),
        partition = _clair3_profile["partition"],
        gpus      = _clair3_profile["gpus"]
    shell:
        SHELL_PRELUDE + """
        echo "Image:            {params.container_image}"
        echo "Container engine: {params.container_engine}"
        echo "Platform:         linux/amd64"
        echo "Bind targets:     {params.bind_args}"

        if [ "{params.is_sif_engine}" = "True" ]; then
            {params.container_engine_bin} exec {params.nv_flag} \\
                {params.bind_args} \\
                {params.singularity_ref} \\
                /opt/bin/run_clair3.sh \\
                -b {input.bam} \\
                -f {input.ref} \\
                -t {threads} \\
                -p {params.platform} \\
                -m /opt/models/{params.model} \\
                -o {params.out_dir} \\
                --sample_name {params.sample_name} \\
                {params.ctg_flag} \\
                {params.phasing_flags} \\
                {params.gpu_flag} 2>&1 | strip_ansi_pipe
        else
            docker run --rm \\
                --platform=linux/amd64 \\
                --user $(id -u):$(id -g) \\
                {params.bind_args} \\
                -w {params.workflow_root} \\
                {params.container_image} \\
                /opt/bin/run_clair3.sh \\
                -b {input.bam} \\
                -f {input.ref} \\
                -t {threads} \\
                -p {params.platform} \\
                -m /opt/models/{params.model} \\
                -o {params.out_dir} \\
                --sample_name {params.sample_name} \\
                {params.ctg_flag} \\
                {params.phasing_flags} \\
                {params.gpu_flag} 2>&1 | strip_ansi_pipe
        fi
        """ + SHELL_EPILOG
