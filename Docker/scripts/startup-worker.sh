#!/bin/bash
set -e

echo "=== CCP4i2 Worker Container Startup ==="
echo "Time: $(date)"

# Note: Environment setup (CCP4, PYTHONPATH, DATABASE_URL) is done by entrypoint.sh
# This script focuses on worker-specific initialization

# Determine which Python to use
WORKER_PYTHON="${CCP4_PYTHON:-python3}"
echo "Using $WORKER_PYTHON for worker"

# Print configuration summary
echo ""
echo "=== Configuration Summary ==="
echo "CCP4_DATA_PATH: ${CCP4_DATA_PATH:-not set}"
echo "DJANGO_SETTINGS_MODULE: ${DJANGO_SETTINGS_MODULE:-not set}"
echo "CCP4 Available: $([ -n "$CCP4" ] && echo "Yes ($CCP4)" || echo 'No')"
echo "CCP4_PYTHON: ${CCP4_PYTHON:-'not set'}"
echo "PYTHONPATH: $PYTHONPATH"
echo "SERVICE_BUS_QUEUE_NAME: ${SERVICE_BUS_QUEUE_NAME:-'not set'}"
echo "============================="
echo ""

# Start the worker
echo "Starting worker process..."
exec $WORKER_PYTHON /usr/src/app/worker.py
