#!/bin/bash
set -e

echo "=== CCP4i2 Worker Container Startup ==="
echo "Time: $(date)"

# Configuration from environment
CCP4_DATA_PATH=${CCP4_DATA_PATH:-"/mnt/ccp4data"}
DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE:-"azure_extensions.settings"}

# Export for Django/worker
export CCP4_DATA_PATH
export DJANGO_SETTINGS_MODULE

# Database configuration (for status updates)
if [ -n "$DATABASE_URL" ]; then
    echo "Using DATABASE_URL from environment"
    export DATABASE_URL
elif [ -n "$DB_HOST" ]; then
    DB_USER=${DB_USER:-"ccp4i2"}
    DB_NAME=${DB_NAME:-"ccp4i2"}
    DB_PASSWORD=${DB_PASSWORD:-""}
    DB_PORT=${DB_PORT:-5432}

    ENCODED_DB_USER=$(python3 -c "import urllib.parse; print(urllib.parse.quote('$DB_USER', safe=''))")
    ENCODED_DB_PASSWORD=$(python3 -c "import urllib.parse; print(urllib.parse.quote('$DB_PASSWORD', safe=''))")

    export DATABASE_URL="postgresql://${ENCODED_DB_USER}:${ENCODED_DB_PASSWORD}@${DB_HOST}:${DB_PORT}/${DB_NAME}"
    echo "Constructed DATABASE_URL for PostgreSQL"
fi

# Setup CCP4 environment if available
CCP4_VERSION=${CCP4_VERSION:-""}
CCP4_DIR=""

if [ -n "$CCP4" ] && [ -d "$CCP4" ]; then
    CCP4_DIR="$CCP4"
elif [ -d "/opt/ccp4" ] && [ -f "/opt/ccp4/bin/ccp4.setup-sh" ]; then
    CCP4_DIR="/opt/ccp4"
elif [ -n "$CCP4_VERSION" ]; then
    CCP4_DIR="$CCP4_DATA_PATH/$CCP4_VERSION"
    if [ ! -d "$CCP4_DIR" ]; then
        echo "WARNING: Specified CCP4_VERSION=$CCP4_VERSION not found at $CCP4_DIR"
        CCP4_DIR=""
    fi
else
    for dir in "$CCP4_DATA_PATH"/ccp4-*; do
        if [ -d "$dir" ] && [ -f "$dir/bin/ccp4.setup-sh" ]; then
            CCP4_DIR="$dir"
            break
        fi
    done
fi

if [ -n "$CCP4_DIR" ] && [ -f "$CCP4_DIR/bin/ccp4.setup-sh" ]; then
    echo "Sourcing CCP4 environment from $CCP4_DIR/bin/ccp4.setup-sh"
    set +e
    . "$CCP4_DIR/bin/ccp4.setup-sh"
    set -e
    export CCP4_PYTHON="$CCP4_DIR/bin/ccp4-python"

    # Ensure container site-packages are at front of PYTHONPATH
    CONTAINER_SITE_PACKAGES=$(python3 -c "import site; print(site.getsitepackages()[0])")
    export PYTHONPATH="/usr/src/app:${CONTAINER_SITE_PACKAGES}:$PYTHONPATH"
    echo "CCP4 environment configured (CCP4=$CCP4_DIR)"
else
    echo "WARNING: CCP4 not found in $CCP4_DATA_PATH"
    echo "Job execution will not work without CCP4"
fi

# Determine which Python to use
WORKER_PYTHON="${CCP4_PYTHON:-python3}"
echo "Using $WORKER_PYTHON for worker"

# Change to app directory
cd /usr/src/app

# Print configuration summary
echo ""
echo "=== Configuration Summary ==="
echo "CCP4_DATA_PATH: $CCP4_DATA_PATH"
echo "DJANGO_SETTINGS_MODULE: $DJANGO_SETTINGS_MODULE"
echo "CCP4 Available: $([ -n "$CCP4_DIR" ] && echo "Yes ($CCP4_DIR)" || echo 'No')"
echo "CCP4_PYTHON: ${CCP4_PYTHON:-'not set'}"
echo "PYTHONPATH: $PYTHONPATH"
echo "SERVICE_BUS_QUEUE_NAME: ${SERVICE_BUS_QUEUE_NAME:-'not set'}"
echo "============================="
echo ""

# Start the worker
echo "Starting worker process..."
exec $WORKER_PYTHON /usr/src/app/worker.py
