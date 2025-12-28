#!/bin/bash
set -e

echo "=== CCP4i2 Server Container Startup ==="
echo "Time: $(date)"

# Configuration from environment
CCP4_DATA_PATH=${CCP4_DATA_PATH:-"/mnt/ccp4data"}
CCP4I2_PROJECTS_DIR=${CCP4I2_PROJECTS_DIR:-"/mnt/ccp4data/ccp4i2-projects"}
DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE:-"ccp4i2.config.settings"}
PORT=${PORT:-8000}

# Export for Django
export CCP4_DATA_PATH
export CCP4I2_PROJECTS_DIR
export DJANGO_SETTINGS_MODULE

# Validate required environment
if [ -z "$SECRET_KEY" ]; then
    echo "ERROR: SECRET_KEY environment variable is required"
    exit 1
fi
export SECRET_KEY

# Database configuration (SQLite for local, PostgreSQL for production)
if [ -n "$DATABASE_URL" ]; then
    echo "Using DATABASE_URL from environment"
elif [ -n "$DB_HOST" ]; then
    # Build PostgreSQL URL from components
    DB_USER=${DB_USER:-"ccp4i2"}
    DB_NAME=${DB_NAME:-"ccp4i2"}
    DB_PASSWORD=${DB_PASSWORD:-""}
    DB_PORT=${DB_PORT:-5432}

    # URL encode credentials
    ENCODED_DB_USER=$(python3 -c "import urllib.parse; print(urllib.parse.quote('$DB_USER', safe=''))")
    ENCODED_DB_PASSWORD=$(python3 -c "import urllib.parse; print(urllib.parse.quote('$DB_PASSWORD', safe=''))")

    export DATABASE_URL="postgresql://${ENCODED_DB_USER}:${ENCODED_DB_PASSWORD}@${DB_HOST}:${DB_PORT}/${DB_NAME}"
    echo "Constructed DATABASE_URL for PostgreSQL"
else
    # Default to SQLite in projects directory
    export DATABASE_URL="sqlite:///${CCP4I2_PROJECTS_DIR}/ccp4i2.sqlite"
    echo "Using SQLite database at ${CCP4I2_PROJECTS_DIR}/ccp4i2.sqlite"
fi

# Setup CCP4 environment if available
# Use CCP4_VERSION env var if set, otherwise auto-detect
CCP4_VERSION=${CCP4_VERSION:-""}
CCP4_DIR=""

if [ -n "$CCP4_VERSION" ]; then
    # Explicit version specified
    CCP4_DIR="$CCP4_DATA_PATH/$CCP4_VERSION"
    if [ ! -d "$CCP4_DIR" ]; then
        echo "WARNING: Specified CCP4_VERSION=$CCP4_VERSION not found at $CCP4_DIR"
        CCP4_DIR=""
    fi
else
    # Auto-detect CCP4 installation
    for dir in "$CCP4_DATA_PATH"/ccp4-*; do
        if [ -d "$dir" ] && [ -f "$dir/bin/ccp4.setup-sh" ]; then
            CCP4_DIR="$dir"
            break
        fi
    done
fi

if [ -n "$CCP4_DIR" ] && [ -f "$CCP4_DIR/bin/ccp4.setup-sh" ]; then
    echo "Sourcing CCP4 environment from $CCP4_DIR/bin/ccp4.setup-sh"
    # Source without set -e to handle missing optional components (like arp_warp)
    set +e
    . "$CCP4_DIR/bin/ccp4.setup-sh"
    set -e
    export CCP4_PYTHON="$CCP4_DIR/bin/ccp4-python"

    # Get the container's site-packages directory (where psycopg2-binary is installed)
    CONTAINER_SITE_PACKAGES=$(python3 -c "import site; print(site.getsitepackages()[0])")

    # Ensure container site-packages and app paths are at front of PYTHONPATH
    # This ensures container-installed packages (like psycopg2-binary) take precedence
    # over CCP4's py-packages which may have incompatible binary modules
    export PYTHONPATH="/usr/src/app:${CONTAINER_SITE_PACKAGES}:$PYTHONPATH"
    echo "CCP4 environment configured (CCP4=$CCP4_DIR)"
else
    echo "WARNING: CCP4 not found in $CCP4_DATA_PATH"
    echo "Job execution will not work without CCP4"
fi

# Determine which Python to use - prefer ccp4-python if available
DJANGO_PYTHON="${CCP4_PYTHON:-python}"
echo "Using $DJANGO_PYTHON for Django server"

# Ensure projects directory exists
mkdir -p "$CCP4I2_PROJECTS_DIR"

# Change to app directory
cd /usr/src/app

# Run migrations
echo "Running Django migrations..."
$DJANGO_PYTHON manage.py migrate --run-syncdb

# Collect static files (for admin interface)
echo "Collecting static files..."
$DJANGO_PYTHON manage.py collectstatic --noinput 2>/dev/null || true

# Print configuration summary
echo ""
echo "=== Configuration Summary ==="
echo "CCP4_DATA_PATH: $CCP4_DATA_PATH"
echo "CCP4I2_PROJECTS_DIR: $CCP4I2_PROJECTS_DIR"
echo "DJANGO_SETTINGS_MODULE: $DJANGO_SETTINGS_MODULE"
echo "PORT: $PORT"
echo "CCP4 Available: $([ -n "$CCP4_DIR" ] && echo "Yes ($CCP4_DIR)" || echo 'No')"
echo "============================="
echo ""

# Start the server
echo "Starting Django server on port $PORT..."

if [ "${USE_GUNICORN:-true}" = "true" ]; then
    # Production: gunicorn with uvicorn workers
    exec $DJANGO_PYTHON -m gunicorn asgi:application \
        -b 0.0.0.0:$PORT \
        --worker-class uvicorn.workers.UvicornWorker \
        --workers ${WORKERS:-2} \
        --timeout 120 \
        --keep-alive 10 \
        --graceful-timeout 30 \
        --access-logfile - \
        --error-logfile -
else
    # Development: uvicorn directly with reload
    exec $DJANGO_PYTHON -m uvicorn asgi:application \
        --host 0.0.0.0 \
        --port $PORT \
        --reload
fi
