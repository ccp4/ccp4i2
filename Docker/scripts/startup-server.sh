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
CCP4_SETUP_SCRIPT="$CCP4_DATA_PATH/ccp4-9/bin/ccp4.setup-sh"
if [ -f "$CCP4_SETUP_SCRIPT" ]; then
    echo "Sourcing CCP4 environment from $CCP4_SETUP_SCRIPT"
    . "$CCP4_SETUP_SCRIPT"
    export CCP4_PYTHON="$CCP4_DATA_PATH/ccp4-9/bin/ccp4-python"

    # Ensure app paths are at front of PYTHONPATH
    export PYTHONPATH="/usr/src/app:$PYTHONPATH"
    echo "CCP4 environment configured"
else
    echo "WARNING: CCP4 not found at $CCP4_SETUP_SCRIPT"
    echo "Job execution will not work without CCP4"
fi

# Ensure projects directory exists
mkdir -p "$CCP4I2_PROJECTS_DIR"

# Change to app directory
cd /usr/src/app

# Run migrations
echo "Running Django migrations..."
python manage.py migrate --run-syncdb

# Collect static files (for admin interface)
echo "Collecting static files..."
python manage.py collectstatic --noinput 2>/dev/null || true

# Print configuration summary
echo ""
echo "=== Configuration Summary ==="
echo "CCP4_DATA_PATH: $CCP4_DATA_PATH"
echo "CCP4I2_PROJECTS_DIR: $CCP4I2_PROJECTS_DIR"
echo "DJANGO_SETTINGS_MODULE: $DJANGO_SETTINGS_MODULE"
echo "PORT: $PORT"
echo "CCP4 Available: $([ -f "$CCP4_SETUP_SCRIPT" ] && echo 'Yes' || echo 'No')"
echo "============================="
echo ""

# Start the server
echo "Starting Django server on port $PORT..."

if [ "${USE_GUNICORN:-true}" = "true" ]; then
    # Production: gunicorn with uvicorn workers
    exec python -m gunicorn asgi:application \
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
    exec python -m uvicorn asgi:application \
        --host 0.0.0.0 \
        --port $PORT \
        --reload
fi
