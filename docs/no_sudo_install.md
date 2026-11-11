# Installing a container engine without sudo access

The pipeline needs either **Docker** or **Singularity/Apptainer** to run
Clair3. On a shared machine or HPC login node you often won't have admin
rights to install either the normal way. This guide walks through a fully
unprivileged install of **Apptainer** (the modern name for Singularity) that
needs no `sudo` at any point.

If you *do* have sudo, don't bother with any of this — just ask your admin
to install Apptainer/Singularity or Docker system-wide, or check if it's
already available as an environment module (`module avail singularity
apptainer docker`).

## Docker (rootless, no sudo)

```bash
curl -fsSL https://get.docker.com/rootless -o get-docker.sh
sh get-docker.sh
```

See https://docs.docker.com/go/rootless/ for details.

## When you'll hit this

You'll see this error from the pipeline:

```
FileNotFoundError in file .../common.smk, line NNN:
Container engine 'docker' not found in PATH.
Install Docker or set cluster.container_engine: singularity
```

## Step 0 — Check whether this path is even possible

Unprivileged containers rely on Linux user namespaces. Check they're enabled:

```bash
sysctl kernel.unprivileged_userns_clone
```

- If it prints `kernel.unprivileged_userns_clone = 1` → keep going, you're
  in good shape.
- If it prints `0`, or the setting doesn't exist on your system, an
  unprivileged install is not possible here — this is a genuine dead end
  and you'll need to ask a system administrator to either enable this
  kernel setting or install a container engine for you.

## Step 1 — Try Docker rootless first (fast path)

```bash
curl -fsSL https://get.docker.com/rootless -o rootless-install.sh
sh rootless-install.sh
```

This commonly fails with a missing `uidmap` package, which itself requires
`sudo apt-get install uidmap` — if you don't have sudo, this is a dead end
too. Confirm with:

```bash
which newuidmap newgidmap
```

If both come back empty, skip to Step 2.

## Step 2 — Build Apptainer from source, unprivileged

This is the reliable no-sudo fallback. It only needs a C compiler and Go —
no root-owned packages.

### 2a. Check prerequisites

```bash
which gcc make git go
```

`gcc`, `make`, and `git` are present on almost every Linux system already.
`go` usually isn't — that's fine, we'll install it locally in the next step
if it's missing.

### 2b. Install Go locally (skip if `go` was already found above)

```bash
cd ~
curl -s https://go.dev/VERSION?m=text | head -1   # confirms current version, e.g. go1.23.4
wget https://go.dev/dl/go1.23.4.linux-amd64.tar.gz   # match the version above
mkdir -p $HOME/.local
tar -C $HOME/.local -xzf go1.23.4.linux-amd64.tar.gz
export PATH=$HOME/.local/go/bin:$PATH
echo 'export PATH=$HOME/.local/go/bin:$PATH' >> ~/.bashrc
go version   # should print go1.23.x
```

### 2c. Clone and configure Apptainer

```bash
cd ~
git clone https://github.com/apptainer/apptainer.git
cd apptainer
git checkout v1.3.6   # or a newer stable tag
./mconfig --without-suid --prefix=$HOME/.local
```

Read the output carefully:

- **"go not found!"** → Go isn't on `PATH`. Repeat step 2b, or double-check
  the `export PATH=...` line actually took effect (`echo $PATH`).
- **Missing a C library** (most commonly `libseccomp+headers`) → this is
  usually just an optional feature and safe to ignore; `mconfig` will note
  it and continue. If `mconfig` instead refuses to proceed because of a
  missing *required* library, that one specific package is the one thing
  worth asking a sysadmin to install — much smaller ask than "please set up
  Apptainer for me."
- **Ends with `Makefile ready, try: make`** → you're good, move on.

### 2d. Build and install

```bash
cd builddir
make
make install
```

This takes several minutes (compiling Go and C sources, downloading Go
module dependencies). Everything is installed under `$HOME/.local` — no
system directories are touched.

### 2e. Add Apptainer to your PATH

```bash
export PATH=$HOME/.local/bin:$PATH
echo 'export PATH=$HOME/.local/bin:$PATH' >> ~/.bashrc
apptainer --version
```

### 2f. Verify it actually works

```bash
apptainer run docker://hello-world
```

If you see `Hello from Docker!`, the install is fully functional —
Apptainer can pull and run Docker Hub images without any privileged
components.

## Step 3 — Point the pipeline at Apptainer

Apptainer installs a `singularity` binary alongside `apptainer` for backward
compatibility, so no code changes are needed — just tell the pipeline which
engine to use. In your config file (e.g. `configs/mini_data.yml`):

```yaml
cluster:
  container_engine: singularity
```

## Step 4 — Re-run the pipeline's tool check

```bash
snakemake all -n -p \
    -s workflow/snakemake/nanopore_adaptive_pipeline.smk \
    --configfile configs/mini_data.yml \
    --profile profiles/local
```

The container-engine validation step should now pass.

## Troubleshooting quick reference

| Symptom | Cause | Fix |
|---|---|---|
| `sudo: a password is required` / `No sudo access` | No admin rights on this machine | Use this guide's unprivileged path |
| Docker rootless asks for `apt-get install uidmap` | Rootless Docker needs a privileged package | Skip Docker, use Apptainer (Step 2) |
| `mconfig` says `go not found!` | Go isn't installed or not on `PATH` | Install Go locally (Step 2b) |
| `mconfig` fails on `unprivileged user namespaces` check | Kernel setting disabled | Ask a sysadmin — no user-level fix exists |
| `apptainer: command not found` after install | `~/.local/bin` not on `PATH` in a new shell | Confirm the `export PATH=...` line is in `~/.bashrc` and re-source it, or open a new terminal |