#!/bin/bash
# Backward-compatible wrapper.
# The canonical location is now ccp4i2/scripts/run_job_safe.sh.
#
# This wrapper accepts either the old 3-arg format (python, manage.py, job_uuid)
# or the new 2-arg format (python, job_uuid) and forwards to the new script.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
NEW_SCRIPT="$SCRIPT_DIR/ccp4i2/scripts/run_job_safe.sh"

if [ "$#" -eq 3 ]; then
    # Old format: <python> <manage.py> <job_uuid> — drop the manage.py arg
    exec "$NEW_SCRIPT" "$1" "$3"
elif [ "$#" -eq 2 ]; then
    # New format: <python> <job_uuid>
    exec "$NEW_SCRIPT" "$1" "$2"
else
    echo "[run_job_safe] Usage: run_job_safe.sh <python> <job_uuid>" >&2
    echo "[run_job_safe]    or: run_job_safe.sh <python> <manage.py> <job_uuid> (legacy)" >&2
    exit 1
fi
