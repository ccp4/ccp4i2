#!/bin/bash
# Get Azure Service Bus connection string for this project
# Uses .env.deployment for resource group and infers namespace name from Bicep convention

set -e

# Load environment variables
ENV_FILE=".env.deployment"
if [ -f "$ENV_FILE" ]; then
  source "$ENV_FILE"
else
  echo "Environment file $ENV_FILE not found."
  exit 1
fi

# Infer namespace and queue names from Bicep convention
PREFIX="ccp4i2-bicep"
NAMESPACE="${PREFIX}-servicebus"
QUEUE="${PREFIX}-jobs"
AUTH_RULE="RootManageSharedAccessKey"

# Get connection string
CONN_STR=$(az servicebus namespace authorization-rule keys list \
  --resource-group "$RESOURCE_GROUP" \
  --namespace-name "$NAMESPACE" \
  --name "$AUTH_RULE" \
  --query primaryConnectionString -o tsv)

if [ -z "$CONN_STR" ]; then
  echo "❌ Could not retrieve Service Bus connection string."
  exit 1
fi

# Output connection string and queue name for use in test scripts
cat <<EOF
SERVICEBUS_CONNECTION_STRING="$CONN_STR"
SERVICEBUS_QUEUE="$QUEUE"
EOF

echo "✅ Service Bus connection string and queue name ready for use."
