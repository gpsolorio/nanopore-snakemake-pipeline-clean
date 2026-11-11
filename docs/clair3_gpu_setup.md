## Clair3 GPU Setup (Singularity)

Clair3 variant calling can run on GPU via a pre-built Singularity image. This
must be built once per machine before running the pipeline with
`cluster.profile: gpu` + `cluster.container_engine: singularity`.

### 1. Pull and build the image

```bash
mkdir -p workflow/tools/singularity
singularity pull workflow/tools/singularity/clair3_v2.0.2_gpu.sif \
    docker://hkubal/clair3:v2.0.2_gpu
```

This is a multi-GB pull/conversion and only needs to be done once. The
pipeline will fail fast with a clear error at startup if this file is
missing (see `_validate_tools()` in `common.smk`).

> **Do not** point the pipeline at `docker://...` directly at run time --
> Singularity has to re-convert the OCI image to SIF on every single job
> invocation if there's no local `.sif`, which can add many minutes per job
> and gets worse the more jobs run concurrently.

### 2. Verify GPU passthrough works

```bash
singularity exec --nv workflow/tools/singularity/clair3_v2.0.2_gpu.sif nvidia-smi
singularity exec --nv workflow/tools/singularity/clair3_v2.0.2_gpu.sif \
    python3 -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"
```

Both should show your GPU. If not, this is a driver/`--nv`/CUDA problem
independent of Clair3 or Snakemake -- fix that first.

### 3. (Recommended) Install squashfuse to avoid per-job sandbox extraction

Without `squashfuse`/`fuse2fs`, Singularity extracts the entire `.sif` to a
temp sandbox on every `exec` call instead of mounting it directly -- slow
for a multi-GB image, repeated on every job. Install into your conda env
(no sudo required):

```bash
conda install -n <your-env-name> -c conda-forge squashfuse fuse2fs
```

### 4. Cap GPU concurrency

The pipeline may schedule multiple `call_variants` jobs concurrently. If
you only have one physical GPU, make sure your Snakemake profile enforces
this, e.g. in `profiles/local/config.yaml`:

```yaml
resources:
  gpus: 1
```

or pass `--resources gpus=1` on the command line. Without this cap,
multiple Clair3 GPU jobs will contend for the same device.

### Version notes

- CPU image: `hkubal/clair3:v2.0.2`
- GPU image: `hkubal/clair3:v2.0.2_gpu`

Keep these in sync -- `v2.0.0_gpu` was never an officially documented tag;
the GPU image was only introduced in the v2.0.1 release. If you bump one,
bump the other, both defined in `common.smk` (`_CLAIR3_IMAGE_CPU`,
`_CLAIR3_IMAGE_GPU_DOCKER`, `_CLAIR3_IMAGE_GPU_SIF`).