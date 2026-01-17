#!/usr/bin/env bash 
set -euo pipefail 
 
 
############################################# 
# Edit these parameters 
############################################# 
STORAGE_ACCOUNT="newcastlecryoem"      # <--- your storage account name 
RESOURCE_GROUP="CryoEM"               # <--- resource group containing the storage account 
CONTAINER_NAME="mycontainer"         # <--- destination container 
LOCAL_PATH="/path/to/local/folder"   # <--- local file/folder to upload (can be a file or directory) 
EXPIRY_HOURS="8"                     # <--- SAS expiry duration (hours) 
 
 
# Optional: tenant hint if you have multiple tenants 
# AZURE_TENANT_ID="xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx" 
 
 
############################################# 
# No edits needed below (unless customizing) 
############################################# 
 
 
# Check dependencies 
command -v az >/dev/null 2>&1 || { echo "ERROR: Azure CLI (az) is not installed or not in PATH."; exit 1; } 
command -v azcopy >/dev/null 2>&1 || { echo "ERROR: AzCopy is not installed or not in PATH."; exit 1; } 
 
 
# Login (interactive). Comment out if you’re already logged in. 
echo ">>> Logging in to Azure..." 
if [[ -n "${AZURE_TENANT_ID:-}" ]]; then 
  az login --tenant "$AZURE_TENANT_ID" >/dev/null 
else 
  az login >/dev/null 
fi 
 
 
# Verify account context 
echo ">>> Using subscription/account:" 
az account show --output table 
 
 
# Basic validation 
if [[ -z "$STORAGE_ACCOUNT" || -z "$RESOURCE_GROUP" || -z "$CONTAINER_NAME" || -z "$LOCAL_PATH" ]]; then 
  echo "ERROR: One or more required parameters are empty." 
  echo "Please set STORAGE_ACCOUNT, RESOURCE_GROUP, CONTAINER_NAME, and LOCAL_PATH at the top of the script." 
  exit 1 
fi 
 
 
# Compute expiry timestamp in ISO 8601 (UTC) 
# GNU date (Linux) path; macOS BSD 'date' fallback. 
NOW_UTC="$(date -u +%Y-%m-%dT%H:%MZ || true)" 
if command -v gdate >/dev/null 2>&1; then 
  # macOS with coreutils installed: use gdate 
  EXPIRY_UTC="$(gdate -u -d "+${EXPIRY_HOURS} hours" +%Y-%m-%dT%H:%MZ)" 
else 
  # Try GNU date; if BSD date, approximate using seconds 
  if date -u -d "+${EXPIRY_HOURS} hours" +%Y-%m-%dT%H:%MZ >/dev/null 2>&1; then 
    EXPIRY_UTC="$(date -u -d "+${EXPIRY_HOURS} hours" +%Y-%m-%dT%H:%MZ)" 
  else 
    # BSD fallback: add seconds 
    EXPIRY_UTC="$(date -u -r $(( $(date +%s) + EXPIRY_HOURS*3600 )) +%Y-%m-%dT%H:%MZ)" 
  fi 
fi 
 
 
echo ">>> SAS start: $NOW_UTC" 
echo ">>> SAS expiry: $EXPIRY_UTC" 
 
 
# Ensure container exists (idempotent) 
echo ">>> Ensuring container exists: $CONTAINER_NAME" 
az storage container create \ 
  --account-name "$STORAGE_ACCOUNT" \ 
  --name "$CONTAINER_NAME" \ 
  --auth-mode login \ 
  --public-access off \ 
  --only-show-errors >/dev/null || true 
 
 
# Generate a User Delegation SAS for the container using Azure AD 
# For uploads, we typically need: create (c), write (w), add (a), list (l), read (r) as needed. 
# 'racwl' grants read/add/create/write/list at the container scope. 
echo ">>> Generating User Delegation SAS for container..." 
SAS_TOKEN="$(az storage container generate-sas \ 
  --account-name "$STORAGE_ACCOUNT" \ 
  --name "$CONTAINER_NAME" \ 
  --permissions racwl \ 
  --start "$NOW_UTC" \ 
  --expiry "$EXPIRY_UTC" \ 
  --auth-mode login \ 
  --as-user \ 
  --output tsv)" 
 
 
if [[ -z "$SAS_TOKEN" ]]; then 
  echo "ERROR: Failed to generate SAS token. Check your role (Storage Blob Data Contributor) and parameters." 
  exit 1 
fi 
 
 
DEST_URL="https://${STORAGE_ACCOUNT}.blob.core.windows.net/${CONTAINER_NAME}?${SAS_TOKEN}" 
 
 
# Perform the copy. For directories, use --recursive; for a single file you can omit it. 
# AzCopy will create target virtual directories as needed. 
echo ">>> Starting AzCopy upload from: $LOCAL_PATH" 
echo ">>> Destination: $DEST_URL" 
# If LOCAL_PATH is a directory, keep --recursive. If it’s a single file, --recursive is ignored. 
azcopy copy "$LOCAL_PATH" "$DEST_URL" --recursive=true 
 
 
echo ">>> Upload complete." 
 
 
# (Optional) Verify by listing top-level blobs (requires 'l' permission) 
echo ">>> Listing uploaded items (top-level):" 
az storage blob list \ 
  --account-name "$STORAGE_ACCOUNT" \ 
  --container-name "$CONTAINER_NAME" \ 
  --num-results 20 \ 
  --auth-mode login \ 
  --output table || true 
 
 
echo ">>> Done." 
 