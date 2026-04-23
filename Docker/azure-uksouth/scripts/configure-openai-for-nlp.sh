#!/bin/bash
#
# Configure the shared Azure OpenAI resource for NLP-query managed-identity auth.
#
# Context (see apps/compounds/docs/NLP_QUERY_PROPOSAL.md §16.2):
# - An Azure OpenAI resource `ddu-openai` already exists in the ccp4i2 RG and
#   is used by the aacr-abstracts container app via API-key auth.
# - The NLP query backend (slices 1-6) authenticates via managed identity (§11
#   decision 12). That requires:
#     1. `customSubDomainName` set on the resource (https://<name>.openai.azure.com/)
#     2. "Cognitive Services OpenAI User" role assignment on the resource,
#        scoped to the existing `containerAppsIdentity`.
# - aacr-abstracts keeps working because (a) the regional endpoint it uses
#   continues to work alongside a custom-subdomain endpoint, and (b) API-key
#   auth is still enabled (disableLocalAuth is not set).
#
# This script is idempotent — safe to run repeatedly. It performs the steps
# above via az-cli so they can be applied without a full infrastructure.bicep
# redeploy. The same changes are also mirrored in infrastructure.bicep so a
# fresh deployment picks them up without needing to run this script.
#
# Requirements: az CLI logged in, selected subscription matches the ccp4i2 RG.

set -e

# Resolve colours only when stdout is a terminal, so redirected output stays clean.
if [ -t 1 ]; then
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    YELLOW='\033[1;33m'
    BLUE='\033[0;34m'
    NC='\033[0m'
else
    RED=''; GREEN=''; YELLOW=''; BLUE=''; NC=''
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="$SCRIPT_DIR/../.env.deployment"
if [ -f "$ENV_FILE" ]; then
    # shellcheck disable=SC1090
    source "$ENV_FILE"
else
    echo -e "${RED}Cannot find $ENV_FILE${NC}"
    exit 1
fi

OPENAI_RESOURCE="${AZURE_OPENAI_RESOURCE_NAME:-ddu-openai}"
CUSTOM_SUBDOMAIN="${AZURE_OPENAI_CUSTOM_SUBDOMAIN:-$OPENAI_RESOURCE}"
ROLE_GUID='5e0bd9bd-7b93-4f28-af87-19fc36ad61bd'  # Cognitive Services OpenAI User

echo -e "${GREEN}=== Azure OpenAI → NLP query backend wiring ===${NC}"
echo "  Resource group:   $RESOURCE_GROUP"
echo "  OpenAI resource:  $OPENAI_RESOURCE"
echo "  Managed identity: $(basename "$CONTAINER_APPS_IDENTITY_ID")"
echo ""

# ---------------------------------------------------------------------------
# Verify the resource exists and we can reach it.
# ---------------------------------------------------------------------------
if ! az cognitiveservices account show \
        --name "$OPENAI_RESOURCE" \
        --resource-group "$RESOURCE_GROUP" \
        --query name -o tsv >/dev/null 2>&1; then
    echo -e "${RED}Cannot find Azure OpenAI resource '$OPENAI_RESOURCE' in resource group '$RESOURCE_GROUP'.${NC}"
    echo "Check AZURE_OPENAI_RESOURCE_NAME in .env.deployment, or verify with:"
    echo "  az cognitiveservices account list --resource-group $RESOURCE_GROUP -o table"
    exit 1
fi

# ---------------------------------------------------------------------------
# Step 1: ensure customSubDomainName is set. Required for managed-identity auth.
# Once set, cannot be removed — chosen conservatively (matches resource name).
# ---------------------------------------------------------------------------
CURRENT_SUBDOMAIN=$(az cognitiveservices account show \
    --name "$OPENAI_RESOURCE" \
    --resource-group "$RESOURCE_GROUP" \
    --query "properties.customSubDomainName" -o tsv 2>/dev/null)

if [ -z "$CURRENT_SUBDOMAIN" ] || [ "$CURRENT_SUBDOMAIN" = "null" ]; then
    echo -e "${BLUE}Step 1/2: Setting customSubDomainName=$CUSTOM_SUBDOMAIN${NC}"
    echo -e "${YELLOW}  Note: this is a one-way change (cannot be removed). aacr-abstracts${NC}"
    echo -e "${YELLOW}  will continue to work against the regional endpoint.${NC}"
    az cognitiveservices account update \
        --name "$OPENAI_RESOURCE" \
        --resource-group "$RESOURCE_GROUP" \
        --custom-domain "$CUSTOM_SUBDOMAIN" \
        --output none
    CURRENT_SUBDOMAIN="$CUSTOM_SUBDOMAIN"
    echo -e "${GREEN}  ✓ customSubDomainName set${NC}"
else
    echo -e "${GREEN}Step 1/2: customSubDomainName already set to '$CURRENT_SUBDOMAIN' — skipping${NC}"
fi

MI_ENDPOINT="https://${CURRENT_SUBDOMAIN}.openai.azure.com/"
echo "  Managed-identity endpoint: $MI_ENDPOINT"
echo ""

# ---------------------------------------------------------------------------
# Step 2: ensure the "Cognitive Services OpenAI User" role is assigned to the
# containerAppsIdentity on this resource.
# ---------------------------------------------------------------------------
OPENAI_RESOURCE_ID=$(az cognitiveservices account show \
    --name "$OPENAI_RESOURCE" \
    --resource-group "$RESOURCE_GROUP" \
    --query id -o tsv)

IDENTITY_PRINCIPAL_ID=$(az identity show \
    --ids "$CONTAINER_APPS_IDENTITY_ID" \
    --query principalId -o tsv 2>/dev/null)

if [ -z "$IDENTITY_PRINCIPAL_ID" ]; then
    echo -e "${RED}Cannot resolve principalId of $CONTAINER_APPS_IDENTITY_ID${NC}"
    exit 1
fi

EXISTING_ASSIGNMENT=$(az role assignment list \
    --assignee "$IDENTITY_PRINCIPAL_ID" \
    --scope "$OPENAI_RESOURCE_ID" \
    --role "$ROLE_GUID" \
    --query "[].id" -o tsv 2>/dev/null)

if [ -z "$EXISTING_ASSIGNMENT" ]; then
    echo -e "${BLUE}Step 2/2: Assigning 'Cognitive Services OpenAI User' to containerAppsIdentity${NC}"
    az role assignment create \
        --assignee-object-id "$IDENTITY_PRINCIPAL_ID" \
        --assignee-principal-type ServicePrincipal \
        --role "$ROLE_GUID" \
        --scope "$OPENAI_RESOURCE_ID" \
        --output none
    echo -e "${GREEN}  ✓ role assignment created${NC}"
else
    echo -e "${GREEN}Step 2/2: Role assignment already present — skipping${NC}"
fi
echo ""

# ---------------------------------------------------------------------------
# Print the AZURE_OPENAI_ENDPOINT value that should go into .env.deployment*
# ---------------------------------------------------------------------------
echo -e "${GREEN}=== Done ===${NC}"
echo ""
echo "Set the following in each instance's .env file to enable the NLP endpoint:"
echo ""
echo "  AZURE_OPENAI_ENDPOINT=$MI_ENDPOINT"
echo "  AZURE_OPENAI_MODEL=gpt-4o"
echo "  COMPOUNDS_NLP_ENABLED=true   # only once you're ready to expose the feature"
echo ""
echo "Then redeploy the server container app:"
echo "  ./scripts/deploy-applications.sh [--env .env.<instance>]"
