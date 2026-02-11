#!/bin/bash
# Crash-safe job runner wrapper.
#
# Catches process crashes (e.g., segfault in a C extension library) and marks
# the job as FAILED in the database, preventing stuck "running" jobs.
#
# When a C extension like iris_validation/mmdb2 segfaults, the Python interpreter
# is killed before any try/except or atexit handler can run. This wrapper detects
# the non-zero exit code from the outside and updates the job status via a
# separate management command.
#
# Used by context_run.py for local/Electron job execution.
# The Azure worker (worker.py) has its own subprocess supervision.
#
# Usage: run_job_safe.sh <python> <manage.py> <job_uuid>

PYTHON="$1"
MANAGE_PY="$2"
JOB_UUID="$3"

if [ -z "$PYTHON" ] || [ -z "$MANAGE_PY" ] || [ -z "$JOB_UUID" ]; then
    echo "[run_job_safe] Usage: run_job_safe.sh <python> <manage.py> <job_uuid>" >&2
    exit 1
fi

"$PYTHON" "$MANAGE_PY" run_job -ju "$JOB_UUID" -y
EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
    # Log the crash (output goes to Django server's stderr)
    if [ $EXIT_CODE -gt 128 ]; then
        SIGNAL=$((EXIT_CODE - 128))
        echo "[run_job_safe] Job $JOB_UUID: killed by signal $SIGNAL (exit code $EXIT_CODE)" >&2
    else
        echo "[run_job_safe] Job $JOB_UUID: exited with code $EXIT_CODE" >&2
    fi

    # Mark the job as FAILED so it doesn't appear stuck in "running" state
    echo "[run_job_safe] Marking job $JOB_UUID as FAILED" >&2
    "$PYTHON" "$MANAGE_PY" set_job_status -ju "$JOB_UUID" -s FAILED 2>/dev/null || \
        echo "[run_job_safe] WARNING: could not update job status" >&2
fi

exit $EXIT_CODE
