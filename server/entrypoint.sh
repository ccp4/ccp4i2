#!/bin/bash
#
# CONTAINER ENTRYPOINT -- not used by the desktop application. Wraps one of the
# startup*.sh scripts in the containerised deployment; the Dockerfiles that use
# it live in the newcastleuniversity/materia repository (commit eb130b163).
#
set -e

# Debug: Show all environment variables that Container Apps should have injected
echo "=== ENTRYPOINT DEBUG ==="
echo "Current time: $(date)"
echo "=== END ENTRYPOINT DEBUG ==="

# Execute the startup script with proper environment variable inheritance
exec "$@"
