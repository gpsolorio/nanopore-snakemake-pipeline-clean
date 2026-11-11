#!/usr/bin/env python3

import subprocess
import sys

jobid = sys.argv[-1]

if jobid == "Submitted":
    sys.stderr.write(
        "Invalid job ID. Did you forget --parsable in sbatch?\n"
    )
    sys.exit(1)

cmd = (
    f"sacct -j {jobid} --format State --noheader "
    "| head -1 | awk '{print $1}'"
)

state = subprocess.check_output(cmd, shell=True, text=True).strip()

running = [
    "PENDING",
    "CONFIGURING",
    "COMPLETING",
    "RUNNING",
    "SUSPENDED",
    "PREEMPTED",
]

if "COMPLETED" in state:
    print("success")
elif any(r in state for r in running):
    print("running")
else:
    print("failed")