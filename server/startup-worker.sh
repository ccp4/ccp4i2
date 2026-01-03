#!/bin/bash
# Worker startup script for CCP4i2 Django Azure
# This script launches the queue monitor instead of the Django server

set -e

# Specific variable checks
echo "DB_HOST: ${DB_HOST:-NOT_SET}"
echo "DB_USER: ${DB_USER:-NOT_SET}"
echo "DB_NAME: ${DB_NAME:-NOT_SET}"
echo "DB_PASSWORD: ${DB_PASSWORD:+SET}"  # Shows SET if variable has value, nothing if empty
echo "SECRET_KEY: ${SECRET_KEY:+SET}"
echo "DJANGO_SETTINGS_MODULE: ${DJANGO_SETTINGS_MODULE:-NOT_SET}"
echo "DB_SSL_MODE: ${DB_SSL_MODE:-NOT_SET}"
echo "DB_SSL_ROOT_CERT: ${DB_SSL_ROOT_CERT:-NOT_SET}"
echo "DB_SSL_REQUIRE_CERT: ${DB_SSL_REQUIRE_CERT:-NOT_SET}"

# Access environment variables (including secrets passed as env vars)
CCP4_DATA_PATH=${CCP4_DATA_PATH:-"/mnt/ccp4data"}
CCP4I2_PROJECTS_DIR=${CCP4I2_PROJECTS_DIR:-"/mnt/ccp4data/ccp4i2-projects"}
DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE:-"ccp4i2.config.settings"}
SECRET_KEY=${SECRET_KEY}
DB_HOST=${DB_HOST:-"db-host"}
DB_USER=${DB_USER:-"db-user"}
DB_NAME=${DB_NAME:-"db-name"}
DB_PASSWORD=${DB_PASSWORD:-"db-password"}

# SSL-related environment variables
DB_SSL_MODE=${DB_SSL_MODE}
DB_SSL_ROOT_CERT=${DB_SSL_ROOT_CERT}
DB_SSL_REQUIRE_CERT=${DB_SSL_REQUIRE_CERT}

# URL encode using Python (more reliable for complex characters)
ENCODED_DB_USER=$(python3 -c "import urllib.parse; print(urllib.parse.quote('$DB_USER', safe=''))")
ENCODED_DB_PASSWORD=$(python3 -c "import urllib.parse; print(urllib.parse.quote('$DB_PASSWORD', safe=''))")

# Build query parameters for SSL configuration
QUERY_PARAMS=""

# Add SSL mode if set
if [ -n "$DB_SSL_MODE" ]; then
    ENCODED_SSL_MODE=$(python3 -c "import urllib.parse; print(urllib.parse.quote('$DB_SSL_MODE', safe=''))")
    QUERY_PARAMS="${QUERY_PARAMS}sslmode=${ENCODED_SSL_MODE}&"
    echo "Adding SSL mode: $DB_SSL_MODE"
fi

# Add SSL root certificate if set
if [ -n "$DB_SSL_ROOT_CERT" ]; then
    ENCODED_SSL_ROOT_CERT=$(python3 -c "import urllib.parse; print(urllib.parse.quote('$DB_SSL_ROOT_CERT', safe=''))")
    QUERY_PARAMS="${QUERY_PARAMS}sslrootcert=${ENCODED_SSL_ROOT_CERT}&"
    echo "Adding SSL root certificate: $DB_SSL_ROOT_CERT"
fi

# Add SSL require certificate if set (maps to sslcert parameter)
if [ -n "$DB_SSL_REQUIRE_CERT" ]; then
    ENCODED_SSL_REQUIRE_CERT=$(python3 -c "import urllib.parse; print(urllib.parse.quote('$DB_SSL_REQUIRE_CERT', safe=''))")
    QUERY_PARAMS="${QUERY_PARAMS}sslcert=${ENCODED_SSL_REQUIRE_CERT}&"
    echo "Adding SSL certificate: $DB_SSL_REQUIRE_CERT"
fi

# Remove trailing & if query parameters exist
if [ -n "$QUERY_PARAMS" ]; then
    QUERY_PARAMS=$(echo "$QUERY_PARAMS" | sed 's/&$//')
    QUERY_PARAMS="?${QUERY_PARAMS}"
fi

# Construct DATABASE_URL from secure components with encoded credentials and SSL parameters
export DATABASE_URL="postgresql://${ENCODED_DB_USER}:${ENCODED_DB_PASSWORD}@${DB_HOST}:5432/${DB_NAME}${QUERY_PARAMS}"

# Validate required environment variables
if [ -z "$DATABASE_URL" ]; then
    echo "ERROR: Failed to construct DATABASE_URL"
    echo "DB_HOST: $DB_HOST"
    echo "DB_USER: $DB_USER" 
    echo "DB_NAME: $DB_NAME"
    echo "DB_PASSWORD: [REDACTED]"
    exit 1
fi

# Show constructed URL (with password masked for security)
MASKED_DATABASE_URL=$(echo "$DATABASE_URL" | sed 's/:\/\/[^:]*:[^@]*@/:\/\/[USER]:[PASSWORD]@/')
echo "Constructed DATABASE_URL: $MASKED_DATABASE_URL"

if [ -z "$SECRET_KEY" ]; then
    echo "ERROR: SECRET_KEY environment variable is required"
    exit 1
fi

# Export variables for subprocesses (e.g., uvicorn)
export CCP4_DATA_PATH
export CCP4I2_PROJECTS_DIR
export DATABASE_URL
export DJANGO_SETTINGS_MODULE
export SECRET_KEY

# Print environment variables for debugging (avoid printing sensitive info)
echo "=== CCP4i2 Container Startup ==="
echo "CCP4_DATA_PATH: $CCP4_DATA_PATH"
echo "CCP4I2_PROJECTS_DIR: $CCP4I2_PROJECTS_DIR"
echo "DATABASE_URL: [HIDDEN FOR SECURITY]"
echo "DJANGO_SETTINGS_MODULE: $DJANGO_SETTINGS_MODULE"
echo "SECRET_KEY: [HIDDEN FOR SECURITY]"

# Quick health check endpoint using Python 3
python3 -c "
import http.server
import socketserver
import threading

class HealthHandler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path == '/health/':
            self.send_response(200)
            self.end_headers()
            self.wfile.write(b'OK')
        else:
            self.send_response(404)
            self.end_headers()

# Start health server in background
server = socketserver.TCPServer(('', 8000), HealthHandler)
thread = threading.Thread(target=server.serve_forever)
thread.daemon = True
thread.start()
print('Health check server started on port 8000')
" &
HEALTH_PID=$!  # Store the PID of the Python process
echo "Health server PID: $HEALTH_PID"

# CCP4 is pre-transferred to the file share
# Use CCP4_VERSION env var if set, otherwise default to ccp4-9
echo "DEBUG: CCP4_VERSION env var before default: '${CCP4_VERSION}'"
CCP4_VERSION=${CCP4_VERSION:-"ccp4-9"}
CCP4_DIR="$CCP4_DATA_PATH/$CCP4_VERSION"
echo "CCP4 distribution is pre-transferred to $CCP4_DIR (version: $CCP4_VERSION)"

# Setup CCP4 environment with retry logic for mounting delays
CCP4_SETUP_SCRIPT="$CCP4_DIR/bin/ccp4.setup-sh"
echo "Waiting for CCP4 setup script ${CCP4_SETUP_SCRIPT} to become available..."
WAIT_COUNT=0
MAX_WAIT=120  # 2 minutes

while [ $WAIT_COUNT -lt $MAX_WAIT ]; do
  if [ -f "$CCP4_SETUP_SCRIPT" ]; then
    echo "CCP4 setup script found after ${WAIT_COUNT} seconds"
    . "$CCP4_SETUP_SCRIPT"
    export CCP4_PYTHON="$CCP4_DIR/bin/ccp4-python"
    echo "CCP4 environment configured successfully (version: $CCP4_VERSION)"
    break
  else
    echo "Waiting for CCP4 setup script... (${WAIT_COUNT}/${MAX_WAIT}s)"
    sleep 1
    WAIT_COUNT=$((WAIT_COUNT + 1))
  fi
done

# Change to app directory
cd /usr/src/app

# After sourcing CCP4 setup, ensure PYTHONPATH starts with version-specific py-packages and app directory
export PYTHONPATH="/mnt/ccp4data/py-packages-${CCP4_VERSION}:/usr/src/app:$PYTHONPATH"
echo "PYTHONPATH manually corrected: $PYTHONPATH"

echo "=========================================="
echo "CCP4i2 Worker Container Starting"
echo "=========================================="
echo "Date: $(date)"
echo "Container: Worker"
echo "=========================================="

# Wait for database to be ready (reuse from entrypoint.sh logic if needed)
echo "Checking database connectivity..."
# Add database health check here if needed

# Launch the worker queue monitor
echo "Starting queue monitor..."
exec $CCP4_PYTHON  /usr/src/app/worker.py
