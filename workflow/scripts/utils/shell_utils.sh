#!/usr/bin/env bash

_SHELL_START_TS=""

shell_start() {
    local rule="${1:-unknown}"
    _SHELL_START_TS=$(date +%s)
    echo "[START] $(date '+%Y-%m-%d %H:%M:%S') — ${rule}"
    echo "[HOST]  $(hostname)"
    echo "[DIR]   $(pwd)"
}

shell_end() {
    local rule="${1:-unknown}"
    local end_ts elapsed
    end_ts=$(date +%s)
    elapsed=$(( end_ts - _SHELL_START_TS ))
    echo "[END]   $(date '+%Y-%m-%d %H:%M:%S') — ${rule}"
    printf "Elapsed: %dd %dh %dm %ds\n" \
        $(( elapsed / 86400 )) \
        $(( (elapsed % 86400) / 3600 )) \
        $(( (elapsed % 3600) / 60 )) \
        $(( elapsed % 60 ))
}

verify_conda_path() {
    # Fails loudly if PATH resolves one of this env's own binaries to a
    # DIFFERENT install than the one Snakemake just activated. Catches the
    # "conda env activated correctly, but a system dir (e.g. /usr/local/bin)
    # sits ahead of it on PATH" class of bug -- the Rscript/python3 shadowing
    # issue hit twice in this pipeline already.
    if [[ -z "${CONDA_PREFIX:-}" ]]; then
        return 0  # not a conda-managed rule -- nothing to check
    fi

    local bin_name expected resolved
    local mismatches=()

    for bin_name in "$@"; do
        expected="${CONDA_PREFIX}/bin/${bin_name}"
        [[ -x "$expected" ]] || continue   # this env doesn't ship it, skip
        resolved="$(command -v "$bin_name" 2>/dev/null || true)"
        if [[ "$resolved" != "$expected" ]]; then
            mismatches+=("  ${bin_name}: PATH resolves to '${resolved:-<not found>}', expected '${expected}'")
        fi
    done

    if (( ${#mismatches[@]} > 0 )); then
        echo "FATAL: conda env activated (${CONDA_PREFIX}) but PATH resolves elsewhere:" >&2
        printf '%s\n' "${mismatches[@]}" >&2
        echo "A system directory is likely ahead of the conda env's bin/ on PATH." >&2
        echo "Fix by calling \$CONDA_PREFIX/bin/<tool> explicitly in the rule." >&2
        exit 1
    fi
}

strip_ansi_from_log() {
    local log_path="$1"
    if [ -f "$log_path" ]; then
        $CONDA_PREFIX/bin/python -c "
import re, sys
ansi = re.compile(r'\x1b\[[0-9;]*[A-Za-z]|\x1b[^[]')
sys.stdout.write(ansi.sub('', open('$log_path').read()))
" > "${log_path}.tmp" && mv "${log_path}.tmp" "$log_path"
    fi
}

strip_ansi_pipe() {
    $CONDA_PREFIX/bin/python -c "
import re, sys
ansi = re.compile(r'\x1b\[[0-9;]*[A-Za-z]|\x1b[^[]')
for line in sys.stdin:
    sys.stdout.write(ansi.sub('', line))
sys.stdout.flush()
"
}