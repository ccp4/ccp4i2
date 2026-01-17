#!/usr/bin/env bash
set -euo pipefail

#############################################
# AzCopy wrapper for Azure Files
#
# Reads configuration from Docker/.env and uploads
# local files/directories to Azure Files share.
#
# Usage:
#   ./azcopy-files.sh <local_path> [remote_path]
#
# Examples:
#   # Upload untarred CCP4 to file share root
#   ./azcopy-files.sh /tmp/ccp4-9 ccp4-9
#
#   # Upload a single file
#   ./azcopy-files.sh ./myfile.txt
#
#   # Sync CCP4 directory (only upload changed files)
#   ./azcopy-files.sh --sync /tmp/ccp4-9 ccp4-9
#############################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${SCRIPT_DIR}/../.env"

# Parse arguments
SYNC_MODE=false
if [[ "${1:-}" == "--sync" ]]; then
    SYNC_MODE=true
    shift
fi

LOCAL_PATH="${1:-}"
REMOTE_PATH="${2:-}"

# Load .env file if it exists
if [[ -f "$ENV_FILE" ]]; then
    echo ">>> Loading configuration from $ENV_FILE"
    # Export variables from .env (skip comments and empty lines)
    set -a
    source <(grep -v '^#' "$ENV_FILE" | grep -v '^$' | grep '=')
    set +a
else
    echo "WARNING: No .env file found at $ENV_FILE"
    echo "You can create one from .env.example or set environment variables directly."
fi

# Required Azure configuration (from .env or environment)
AZURE_STORAGE_ACCOUNT="${AZURE_STORAGE_ACCOUNT:-}"
AZURE_RESOURCE_GROUP="${AZURE_RESOURCE_GROUP:-}"
AZURE_FILE_SHARE="${AZURE_FILE_SHARE:-ccp4data}"
AZURE_TENANT_ID="${AZURE_TENANT_ID:-}"
EXPIRY_HOURS="${AZURE_SAS_EXPIRY_HOURS:-8}"

# Validate required parameters
if [[ -z "$LOCAL_PATH" ]]; then
    echo "Usage: $0 [--sync] <local_path> [remote_path]"
    echo ""
    echo "Options:"
    echo "  --sync    Use azcopy sync instead of copy (only upload changes)"
    echo ""
    echo "Examples:"
    echo "  $0 /tmp/ccp4-9 ccp4-9           # Upload directory to share/ccp4-9/"
    echo "  $0 --sync /tmp/ccp4-9 ccp4-9    # Sync directory (faster for updates)"
    echo "  $0 ./file.txt                    # Upload single file to share root"
    echo ""
    echo "Required environment variables (set in Docker/.env):"
    echo "  AZURE_STORAGE_ACCOUNT   - Azure storage account name"
    echo "  AZURE_RESOURCE_GROUP    - Resource group containing the storage account"
    echo "  AZURE_FILE_SHARE        - File share name (default: ccp4data)"
    exit 1
fi

if [[ -z "$AZURE_STORAGE_ACCOUNT" ]]; then
    echo "ERROR: AZURE_STORAGE_ACCOUNT is not set."
    echo "Add it to Docker/.env or set it as an environment variable."
    exit 1
fi

if [[ -z "$AZURE_RESOURCE_GROUP" ]]; then
    echo "ERROR: AZURE_RESOURCE_GROUP is not set."
    echo "Add it to Docker/.env or set it as an environment variable."
    exit 1
fi

if [[ ! -e "$LOCAL_PATH" ]]; then
    echo "ERROR: Local path does not exist: $LOCAL_PATH"
    exit 1
fi

# Check dependencies
command -v az >/dev/null 2>&1 || { echo "ERROR: Azure CLI (az) is not installed or not in PATH."; exit 1; }
command -v azcopy >/dev/null 2>&1 || { echo "ERROR: AzCopy is not installed or not in PATH."; exit 1; }

# Check if already logged in
echo ">>> Checking Azure login status..."
if ! az account show >/dev/null 2>&1; then
    echo ">>> Not logged in. Logging in to Azure..."
    if [[ -n "$AZURE_TENANT_ID" ]]; then
        az login --tenant "$AZURE_TENANT_ID"
    else
        az login
    fi
fi

# Show current account
echo ">>> Using Azure subscription:"
az account show --query "{Name:name, ID:id}" --output table

# Compute expiry timestamp in ISO 8601 (UTC)
NOW_UTC="$(date -u +%Y-%m-%dT%H:%MZ)"
if command -v gdate >/dev/null 2>&1; then
    # macOS with coreutils
    EXPIRY_UTC="$(gdate -u -d "+${EXPIRY_HOURS} hours" +%Y-%m-%dT%H:%MZ)"
elif date -u -d "+${EXPIRY_HOURS} hours" +%Y-%m-%dT%H:%MZ >/dev/null 2>&1; then
    # GNU date
    EXPIRY_UTC="$(date -u -d "+${EXPIRY_HOURS} hours" +%Y-%m-%dT%H:%MZ)"
else
    # BSD date fallback (macOS without coreutils)
    EXPIRY_UTC="$(date -u -r $(( $(date +%s) + EXPIRY_HOURS*3600 )) +%Y-%m-%dT%H:%MZ)"
fi

echo ">>> SAS valid from: $NOW_UTC"
echo ">>> SAS expires at: $EXPIRY_UTC"

# Generate SAS token for the file share
echo ">>> Generating SAS token for file share: $AZURE_FILE_SHARE"
SAS_TOKEN="$(az storage share generate-sas \
    --account-name "$AZURE_STORAGE_ACCOUNT" \
    --name "$AZURE_FILE_SHARE" \
    --permissions rcwdl \
    --start "$NOW_UTC" \
    --expiry "$EXPIRY_UTC" \
    --auth-mode login \
    --as-user \
    --output tsv 2>/dev/null)" || {
    # Fallback: use account key if user delegation fails
    echo ">>> User delegation SAS failed, trying with account key..."
    ACCOUNT_KEY="$(az storage account keys list \
        --resource-group "$AZURE_RESOURCE_GROUP" \
        --account-name "$AZURE_STORAGE_ACCOUNT" \
        --query '[0].value' -o tsv)"

    SAS_TOKEN="$(az storage share generate-sas \
        --account-name "$AZURE_STORAGE_ACCOUNT" \
        --account-key "$ACCOUNT_KEY" \
        --name "$AZURE_FILE_SHARE" \
        --permissions rcwdl \
        --start "$NOW_UTC" \
        --expiry "$EXPIRY_UTC" \
        --output tsv)"
}

if [[ -z "$SAS_TOKEN" ]]; then
    echo "ERROR: Failed to generate SAS token."
    echo "Check your permissions and ensure the file share exists."
    exit 1
fi

# Build destination URL
# For Azure Files: https://<account>.file.core.windows.net/<share>/<path>?<sas>
if [[ -n "$REMOTE_PATH" ]]; then
    DEST_URL="https://${AZURE_STORAGE_ACCOUNT}.file.core.windows.net/${AZURE_FILE_SHARE}/${REMOTE_PATH}?${SAS_TOKEN}"
else
    DEST_URL="https://${AZURE_STORAGE_ACCOUNT}.file.core.windows.net/${AZURE_FILE_SHARE}?${SAS_TOKEN}"
fi

# Determine source path display
if [[ -d "$LOCAL_PATH" ]]; then
    FILE_COUNT=$(find "$LOCAL_PATH" -type f | wc -l | tr -d ' ')
    echo ">>> Source: $LOCAL_PATH ($FILE_COUNT files)"
else
    echo ">>> Source: $LOCAL_PATH"
fi
echo ">>> Destination: https://${AZURE_STORAGE_ACCOUNT}.file.core.windows.net/${AZURE_FILE_SHARE}/${REMOTE_PATH:-}"

# Perform the copy or sync
if [[ "$SYNC_MODE" == true ]]; then
    echo ">>> Starting AzCopy SYNC (only changed files)..."
    azcopy sync "$LOCAL_PATH" "$DEST_URL" --recursive=true
else
    echo ">>> Starting AzCopy COPY..."
    azcopy copy "$LOCAL_PATH" "$DEST_URL" --recursive=true
fi

echo ""
echo ">>> Upload complete!"
echo ""

# Show what was uploaded
echo ">>> Listing uploaded items:"
az storage file list \
    --account-name "$AZURE_STORAGE_ACCOUNT" \
    --share-name "$AZURE_FILE_SHARE" \
    --path "${REMOTE_PATH:-}" \
    --auth-mode login \
    --num-results 20 \
    --output table 2>/dev/null || echo "(Could not list files - check permissions)"

echo ""
echo ">>> Done."
