#!/bin/bash
set -e

# Change to the script's directory to ensure correct relative paths
cd "$(dirname "$0")"

# Source environment variables from parent directory (bicep/)
if [ -f ../.env.deployment ]; then
    source ../.env.deployment
fi

# Fallback defaults if not set
RESOURCE_GROUP="${RESOURCE_GROUP:-ccp4i2-bicep-rg-uksouth}"
LOCATION="${LOCATION:-uksouth}"
PREFIX="${PREFIX:-ccp4i2-bicep}"
KEYVAULT_NAME="${KEYVAULT_NAME:-}"  # Will be set by .env.deployment after infrastructure deployment

az deployment group create \
  --resource-group "$RESOURCE_GROUP" \
  --template-file ../infrastructure/servicebus.bicep \
  --parameters location="$LOCATION" prefix="$PREFIX" keyVaultName="$KEYVAULT_NAME"

echo "Service Bus queue deployed."
echo "Connection string stored in Key Vault as 'servicebus-connection'."