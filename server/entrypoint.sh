#!/bin/bash
set -e

# Debug: Show all environment variables that Container Apps should have injected
echo "=== ENTRYPOINT DEBUG ==="
echo "Current time: $(date)"
echo "=== END ENTRYPOINT DEBUG ==="

# Execute the startup script with proper environment variable inheritance
exec "$@"
