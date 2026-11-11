##############################################################################
# POD5 to BAM Conversion (Dorado Basecaller)
##############################################################################

_DORADO_BIN = (TOOLS_DIR / "dorado" / "bin" / "dorado").resolve()


def get_pod5_dir(wildcards):
    pod5s = get_pod5_paths(wildcards)
    paths = [str(Path(p).resolve()) for p in pod5s]
    try:
        common = os.path.commonpath(paths)
    except Exception:
        common = str(Path(paths[0]).parent)
    return common


checkpoint basecall_pod5:
    input:
        pod5_dir = get_pod5_dir,
        ref      = str(REFERENCE_GENOME)
    output:
        outdir = directory(sample("{barcode}", "bam", "basecall"))
    log:
        sample_log("{barcode}", "basecall_pod5.log")
    conda:
        str(ENVS_DIR / "ngs-utils-env.yml")
    params:
        dorado_bin  = str(_DORADO_BIN),
        model       = RUN_DORADO_CFG["model"],
        model_dir   = str(DORADO_MODELS_DIR),
        pty_wrapper = str(_PTY_WRAPPER),
        device_flag = RUN_DEVICE_FLAG
    resources:
        mem_gb    = res("basecall_pod5", "mem_gb"),
        time      = res("basecall_pod5", "time"),
        account   = res("basecall_pod5", "account"),
        partition = res("basecall_pod5", "partition")
    threads:
        res("basecall_pod5", "threads")
    shell:
        SHELL_PRELUDE_ORDERED + """
        mkdir -p {output.outdir}

        PTY_LOG={log}.dorado.tmp python3 {params.pty_wrapper} \\
            {params.dorado_bin} basecaller \\
                {params.model} \\
                {input.pod5_dir} \\
                --recursive \\
                --output-dir {output.outdir} \\
                --reference {input.ref} \\
                --models-directory {params.model_dir} \\
                --modified-bases 5mCG_5hmCG \\
                {params.device_flag} \\
                --verbose > /dev/null 2>&1
        
        cat {log}.dorado.tmp >> {log}
        rm -f {log}.dorado.tmp
    
        BAM_COUNT=$(find {output.outdir} -name "*.bam" | wc -l)
        echo "[CHECK] Dorado produced $BAM_COUNT BAM file(s)" >> {log}
        
        if [ "$BAM_COUNT" -eq 0 ]; then
            echo "[ERROR] No BAM files found — dorado may have failed" >> {log}
            exit 1
        fi
        """ + SHELL_EPILOG_ORDERED