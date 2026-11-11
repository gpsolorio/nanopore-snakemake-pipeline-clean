#!/usr/bin/env python3
"""
run_with_pty.py — Run a command inside a pseudo-terminal to capture
output that tools write directly to /dev/tty (e.g. dorado).

Usage:
    PTY_LOG=/path/to/log python run_with_pty.py <command> [args...]
"""

import os
import pty
import re
import sys

# Match all ANSI escape sequences including color codes
_ANSI_ESCAPE = re.compile(rb"\x1b\[[\d;]*[A-Za-z]|\x1b[^[]|\r")


def strip_ansi(data: bytes) -> bytes:
    return _ANSI_ESCAPE.sub(b"", data)


def main():
    if len(sys.argv) < 2:
        print("Usage: run_with_pty.py <command> [args...]", file=sys.stderr)
        sys.exit(1)

    cmd      = sys.argv[1:]
    log_path = os.environ.get("PTY_LOG")
    log_fh   = open(log_path, "ab") if log_path else None

    seen = set()

    def read(fd):
        data = os.read(fd, 4096)
        clean = strip_ansi(data)

        if log_fh:
            for line in clean.splitlines(keepends=True):
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped not in seen:
                    seen.add(stripped)
                    log_fh.write(line)
            log_fh.flush()

        sys.stdout.buffer.write(data)
        sys.stdout.buffer.flush()
        return data

    try:
        pty.spawn(cmd, read)
    finally:
        if log_fh:
            log_fh.close()


if __name__ == "__main__":
    main()