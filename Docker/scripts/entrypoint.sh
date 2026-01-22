#!/bin/bash
#
# Container entrypoint script for CCP4i2 server
#
# This script sets up the environment for ALL commands run in the container,
# including interactive shells (az containerapp exec) and management commands.
#
# The environment setup includes:
#   - CCP4 environment (ccp4.setup-sh)
#   - PYTHONPATH configuration
#   - Database configuration
#   - Required environment variables
#
# Usage:
#   ENTRYPOINT ["/usr/src/app/entrypoint.sh"]
#   CMD ["startup.sh"]  # or any other command
#

# Configuration from environment
CCP4_DATA_PATH=${CCP4_DATA_PATH:-"/mnt/ccp4data"}
CCP4I2_PROJECTS_DIR=${CCP4I2_PROJECTS_DIR:-"/mnt/projects"}
DJANGO_SETTINGS_MODULE=${DJANGO_SETTINGS_MODULE:-"azure_extensions.settings"}

# Export for Django
export CCP4_DATA_PATH
export CCP4I2_PROJECTS_DIR
export DJANGO_SETTINGS_MODULE

# Database configuration (SQLite for local, PostgreSQL for production)
if [ -n "$DATABASE_URL" ]; then
    : # Use DATABASE_URL as-is
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
else
    # Default to SQLite in projects directory
    export DATABASE_URL="sqlite://${CCP4I2_PROJECTS_DIR}/ccp4i2.sqlite"
fi

# Setup CCP4 environment if available
# Priority:
#   1. CCP4 env var (set by bundled image or explicit config)
#   2. /opt/ccp4 (bundled in image)
#   3. CCP4_DATA_PATH/CCP4_VERSION (mounted)
#   4. Auto-detect in CCP4_DATA_PATH
CCP4_VERSION=${CCP4_VERSION:-""}
CCP4_DIR=""

if [ -n "$CCP4" ] && [ -d "$CCP4" ]; then
    CCP4_DIR="$CCP4"
elif [ -d "/opt/ccp4" ] && [ -f "/opt/ccp4/bin/ccp4.setup-sh" ]; then
    CCP4_DIR="/opt/ccp4"
elif [ -n "$CCP4_VERSION" ]; then
    CCP4_DIR="$CCP4_DATA_PATH/$CCP4_VERSION"
    if [ ! -d "$CCP4_DIR" ]; then
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
    # Source CCP4 environment (suppress errors for missing optional components)
    . "$CCP4_DIR/bin/ccp4.setup-sh" 2>/dev/null || true
    export CCP4_PYTHON="$CCP4_DIR/bin/ccp4-python"
fi

# Add azure_extensions and azure_packages to PYTHONPATH for management commands
# This ensures import_legacy_* commands can find all required modules
# azure_packages MUST be first so ccp4-python finds our packages before CCP4's
export PYTHONPATH="/usr/src/app/azure_packages:/usr/src/app/azure_extensions:/usr/src/app:$PYTHONPATH"

# Note: Django packages are pre-installed to azure_packages during Docker build
# The PYTHONPATH set above ensures ccp4-python finds them

# Convenience: make 'python' and 'ccp4-python' work correctly
# The default Python should be the one with all Django dependencies
alias python='python3'
alias ccp4-python='${CCP4_PYTHON:-python3}'

# Change to app directory
cd /usr/src/app

# Execute the command passed to the entrypoint
# If no command provided, start an interactive bash shell
if [ $# -eq 0 ]; then
    exec /bin/bash
else
    exec "$@"
fi
