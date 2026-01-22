#!/bin/bash
set -e

echo "=== CCP4i2 Server Container Startup ==="
echo "Time: $(date)"

# Note: Environment setup (CCP4, PYTHONPATH, DATABASE_URL) is done by entrypoint.sh
# This script focuses on server-specific initialization

PORT=${PORT:-8000}

# Validate required environment
if [ -z "$SECRET_KEY" ]; then
    echo "ERROR: SECRET_KEY environment variable is required"
    exit 1
fi

# Determine which Python to use - prefer ccp4-python if available
DJANGO_PYTHON="${CCP4_PYTHON:-python3}"
echo "Using $DJANGO_PYTHON for Django server"

# Debug: Show Python path configuration
echo "=== Python Environment Debug ==="
echo "PYTHONPATH: $PYTHONPATH"
echo "Checking azure_packages:"
ls -la /usr/src/app/azure_packages/ 2>/dev/null | head -5 || echo "azure_packages directory not found!"
echo "Python sys.path:"
$DJANGO_PYTHON -c "import sys; print('\\n'.join(sys.path))" 2>&1 || echo "Failed to get sys.path"
echo "================================"

# Install Pillow if not already present
# Required for Django ImageField in compounds app
# Use --user to install to writable user site-packages (works with read-only CCP4 mount)
echo "Ensuring Pillow is available..."
$DJANGO_PYTHON -m pip install --user --quiet Pillow 2>/dev/null || echo "Note: Pillow install skipped"

# Ensure projects directory exists
mkdir -p "${CCP4I2_PROJECTS_DIR:-/mnt/projects}"

# Run migrations
# --fake-initial: If tables already exist but migrations aren't recorded, mark them as applied
# This handles cases where tables were created before migrations were added to source control
echo "Running Django migrations..."
$DJANGO_PYTHON manage.py migrate --run-syncdb --fake-initial

# Collect static files (for admin interface)
echo "Collecting static files..."
$DJANGO_PYTHON manage.py collectstatic --noinput 2>/dev/null || true

# Print configuration summary
echo ""
echo "=== Configuration Summary ==="
echo "CCP4_DATA_PATH: ${CCP4_DATA_PATH:-not set}"
echo "CCP4I2_PROJECTS_DIR: ${CCP4I2_PROJECTS_DIR:-not set}"
echo "DJANGO_SETTINGS_MODULE: ${DJANGO_SETTINGS_MODULE:-not set}"
echo "PORT: $PORT"
echo "CCP4 Available: $([ -n "$CCP4" ] && echo "Yes ($CCP4)" || echo 'No')"
echo "PYTHONPATH: $PYTHONPATH"
echo "============================="
echo ""

# Test Django URL loading before starting server
echo "=== Testing Django URL Configuration ==="
$DJANGO_PYTHON -c "
import django
django.setup()
from django.urls import get_resolver
from django.conf import settings
print(f'ROOT_URLCONF: {settings.ROOT_URLCONF}')
print(f'INSTALLED_APPS: {settings.INSTALLED_APPS}')
resolver = get_resolver()
print('URL patterns:')
for pattern in resolver.url_patterns:
    print(f'  {pattern.pattern}')
" || echo "WARNING: Django URL test failed"
echo "========================================"
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
