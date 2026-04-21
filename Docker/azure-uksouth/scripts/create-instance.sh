#!/bin/bash
# create-instance.sh — bootstrap a new CCP4i2 deployment in UK South that shares
# the ccp4i2-demo app registration and the ccp4i2-bicep shared ACR.
#
# Usage:
#   ./create-instance.sh <instance-name> <team-name-or-group-id>
#
# Example:
#   ./create-instance.sh york 'York Structural Biology Lab'
#   ./create-instance.sh york 1234abcd-...-5678        # raw group object ID
#
# What it does:
#   1. Creates resource group ccp4i2-<instance>-rg-uksouth.
#   2. Deploys infrastructure.bicep (VNet, Postgres, KV, Storage, Service Bus,
#      MI, [unused] per-instance ACR).
#   3. Seeds the Key Vault with DB password, Django secret, registry password.
#   4. Wires up the shared ACR (AcrPull on the new MI, private endpoint,
#      DNS zone group) so containers can pull from ccp4acrukbwmx.
#   5. Writes .env.<instance> using the ccp4i2-demo AAD app + the Team's group
#      object ID.
#   6. Deploys applications.bicep via deploy-applications.sh.
#   7. Patches the ccp4i2-demo app registration: groupMembershipClaims=All +
#      assigns the Team's group to the enterprise app + adds the new web
#      redirect URI.
#
# Prerequisites:
#   az login as an owner of the ccp4i2-demo app registration and with rights
#   to create resource groups and role assignments in the subscription.

set -e

# ── Inputs ──────────────────────────────────────────────────────────────

INSTANCE="${1:?instance name required (e.g. york)}"
TEAM="${2:?team name or group object ID required}"

if ! [[ "$INSTANCE" =~ ^[a-z][a-z0-9-]*$ ]]; then
    echo "instance name must be lowercase alphanumeric (optionally hyphens), e.g. 'york' or 'crick-fbs'"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BICEP_DIR="$(dirname "$SCRIPT_DIR")"

# ── Fixed config ────────────────────────────────────────────────────────

LOCATION=uksouth
SUB_ID=$(az account show --query id -o tsv)
TENANT_ID=$(az account show --query tenantId -o tsv)

DEMO_APP_CLIENT_ID=cc780b24-ca44-4fec-b8e6-48d0c696a888
DEMO_APP_OBJECT_ID=$(az ad app show --id "$DEMO_APP_CLIENT_ID" --query id -o tsv)
DEMO_APP_SP_ID=$(az ad sp show --id "$DEMO_APP_CLIENT_ID" --query id -o tsv)

SHARED_ACR_NAME=ccp4acrshareduk14fb
SHARED_ACR_RG=ccp4i2-shared-rg-uksouth
SHARED_ACR_ID=$(az acr show -g "$SHARED_ACR_RG" -n "$SHARED_ACR_NAME" --query id -o tsv)
SHARED_ACR_LOGIN_SERVER=${SHARED_ACR_NAME}.azurecr.io

PLATFORM_ADMIN_EMAILS="nmemn@newcastle.ac.uk,martin.noble@newcastle.ac.uk"

PREFIX=ccp4i2-$INSTANCE
RG=${PREFIX}-rg-uksouth

# ── Colors ──────────────────────────────────────────────────────────────

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'
step() { echo -e "${GREEN}== $* ==${NC}"; }
note() { echo -e "${YELLOW}$*${NC}"; }
fail() { echo -e "${RED}$*${NC}"; exit 1; }

# ── 0. Resolve team → M365 group ID ─────────────────────────────────────

step "Resolving Team '$TEAM'"
if [[ "$TEAM" =~ ^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$ ]]; then
    GROUP_ID=$TEAM
    GROUP_NAME=$(az ad group show --group "$GROUP_ID" --query displayName -o tsv 2>/dev/null || echo "")
    [ -z "$GROUP_NAME" ] && fail "Group object ID $GROUP_ID not found in this tenant"
    note "Matched group: $GROUP_NAME ($GROUP_ID)"
else
    MATCHES=$(az ad group list --display-name "$TEAM" --query "[?contains(groupTypes, 'Unified')].{id:id, name:displayName}" -o tsv)
    COUNT=$(echo "$MATCHES" | grep -c . || true)
    if [ "$COUNT" -eq 0 ]; then
        fail "No M365 group (Team) matches display name '$TEAM'"
    elif [ "$COUNT" -gt 1 ]; then
        echo "Multiple M365 groups match '$TEAM':"
        echo "$MATCHES"
        fail "Pass the group object ID directly to disambiguate."
    fi
    GROUP_ID=$(echo "$MATCHES" | head -1 | awk -F'\t' '{print $1}')
    GROUP_NAME=$(echo "$MATCHES" | head -1 | awk -F'\t' '{print $2}')
    note "Matched group: $GROUP_NAME ($GROUP_ID)"
fi

# ── 1. Resource group ───────────────────────────────────────────────────

step "Creating resource group $RG"
az group create -n "$RG" -l "$LOCATION" -o none

# ── 2. Infrastructure ───────────────────────────────────────────────────

step "Deploying infrastructure.bicep into $RG (this takes ~10 minutes)"
DB_PASSWORD=$(openssl rand -base64 16)
INFRA_NAME="infrastructure-$(date +%Y%m%d-%H%M%S)"

az deployment group create \
  --resource-group "$RG" \
  --template-file "$BICEP_DIR/infrastructure/infrastructure.bicep" \
  --parameters location="$LOCATION" \
               prefix="$PREFIX" \
               environment=uk \
               postgresAdminPassword="$DB_PASSWORD" \
               skipPostgresDeployment=false \
               skipCcp4Storage=true \
  --name "$INFRA_NAME" \
  --mode Incremental \
  -o none

get_out() {
    az deployment group show -g "$RG" -n "$INFRA_NAME" \
        --query "properties.outputs.$1.value" -o tsv
}
ACR_NAME_LOCAL=$(get_out acrName)
KV_NAME=$(get_out keyVaultName)
MI_ID=$(get_out containerAppsIdentityId)
MI_PRINCIPAL=$(get_out containerAppsIdentityPrincipalId)
VNET_NAME=$(get_out vnetName)

note "Managed identity principal: $MI_PRINCIPAL"
note "Key Vault: $KV_NAME"
note "VNet: $VNET_NAME"

# ── 3. Seed Key Vault secrets ───────────────────────────────────────────

step "Seeding $KV_NAME with db/django/registry secrets"
# Force IPv4 — KV network rules reject IPv6. curl -4 isn't always honoured on macOS.
MY_IP=$(curl -s https://api.ipify.org || curl -s https://ipv4.icanhazip.com)
if ! [[ "$MY_IP" =~ ^[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    fail "Could not determine IPv4 public IP (got: '$MY_IP'). Key Vault firewall needs IPv4."
fi
az keyvault network-rule add --name "$KV_NAME" --ip-address "$MY_IP" -o none
az keyvault update --name "$KV_NAME" --public-network-access Enabled -o none

ME=$(az ad signed-in-user show --query id -o tsv)
az role assignment create \
    --assignee "$ME" \
    --role "Key Vault Secrets Officer" \
    --scope "/subscriptions/$SUB_ID/resourceGroups/$RG/providers/Microsoft.KeyVault/vaults/$KV_NAME" \
    -o none 2>/dev/null || true

# Wait for role propagation
sleep 30

az keyvault secret set --vault-name "$KV_NAME" --name database-admin-password --value "$DB_PASSWORD" -o none
az keyvault secret set --vault-name "$KV_NAME" --name django-secret-key --value "$(openssl rand -base64 32)" -o none
# Per-instance ACR is unused but the infra expects this secret to exist for symmetry
ACR_PW=$(az acr credential show --name "$ACR_NAME_LOCAL" --query passwords[0].value -o tsv 2>/dev/null || echo "")
if [ -n "$ACR_PW" ]; then
    az keyvault secret set --vault-name "$KV_NAME" --name registry-password --value "$ACR_PW" -o none
fi

# Service Bus connection string (needed by server + worker)
SB_NAMESPACE=$(get_out serviceBusNamespaceName)
SB_CONN=$(az servicebus namespace authorization-rule keys list \
    -g "$RG" --namespace-name "$SB_NAMESPACE" --name RootManageSharedAccessKey \
    --query primaryConnectionString -o tsv 2>/dev/null || echo "")
if [ -n "$SB_CONN" ]; then
    az keyvault secret set --vault-name "$KV_NAME" --name servicebus-connection --value "$SB_CONN" -o none
fi

az keyvault network-rule remove --name "$KV_NAME" --ip-address "$MY_IP" -o none
az keyvault update --name "$KV_NAME" --public-network-access Disabled -o none

# ── 4. Wire up shared ACR ───────────────────────────────────────────────

step "Wiring shared ACR $SHARED_ACR_NAME into $VNET_NAME"

# 4a. AcrPull role assignment (idempotent)
az role assignment create \
    --assignee "$MI_PRINCIPAL" \
    --role AcrPull \
    --scope "$SHARED_ACR_ID" \
    -o none 2>/dev/null || true

# 4b. Private endpoint
PE_NAME=${SHARED_ACR_NAME}-pe-from-${INSTANCE}
if ! az network private-endpoint show -g "$RG" -n "$PE_NAME" -o none 2>/dev/null; then
    az network private-endpoint create \
      --name "$PE_NAME" \
      --resource-group "$RG" \
      --vnet-name "$VNET_NAME" \
      --subnet private-endpoints-subnet \
      --private-connection-resource-id "$SHARED_ACR_ID" \
      --group-id registry \
      --connection-name "${SHARED_ACR_NAME}-conn-${INSTANCE}" \
      -o none
fi

# 4c. DNS zone group — links PE IPs into the privatelink.azurecr.io zone
# (the zone was already created by infrastructure.bicep for the per-instance ACR)
if ! az network private-endpoint dns-zone-group show -g "$RG" --endpoint-name "$PE_NAME" --name default -o none 2>/dev/null; then
    az network private-endpoint dns-zone-group create \
      --resource-group "$RG" \
      --endpoint-name "$PE_NAME" \
      --name default \
      --private-dns-zone privatelink.azurecr.io \
      --zone-name privatelink-azurecr-io \
      -o none
fi

# ── 5. Write .env.<instance> ────────────────────────────────────────────

ENV_FILE="$BICEP_DIR/.env.$INSTANCE"
step "Writing $ENV_FILE"

LATEST_TAG=$(az acr repository show-tags --name "$SHARED_ACR_NAME" --repository ccp4i2/web --orderby time_desc --top 1 -o tsv 2>/dev/null | head -1)
[ -z "$LATEST_TAG" ] && fail "Could not find any image tags in shared ACR; push an image first."
note "Latest shared image tag: $LATEST_TAG"

cat > "$ENV_FILE" <<EOF
# Generated by create-instance.sh on $(date -u +%Y-%m-%dT%H:%M:%SZ)
# Instance: $INSTANCE
# Team:     $GROUP_NAME ($GROUP_ID)

# Shared ACR (images built/pushed once, consumed by all instances)
ACR_NAME=$SHARED_ACR_NAME
ACR_LOGIN_SERVER=$SHARED_ACR_LOGIN_SERVER
SHARED_ACR_RESOURCE_GROUP=$SHARED_ACR_RG

# Instance resource group
RESOURCE_GROUP=$RG
CONTAINER_APP_PREFIX=$PREFIX
CUSTOM_DOMAIN=

# Managed identity provisioned by infrastructure.bicep
CONTAINER_APPS_IDENTITY_ID=$MI_ID
CONTAINER_APPS_IDENTITY_PRINCIPAL_ID=$MI_PRINCIPAL

# Pinned image tags — bump when promoting a new build for this instance
IMAGE_TAG_WEB=$LATEST_TAG
IMAGE_TAG_SERVER=$LATEST_TAG

CCP4_VERSION=ccp4-20251105

# Auth — reuses the ccp4i2-demo app registration
ENABLE_AUTHENTICATION=true
NEXT_PUBLIC_REQUIRE_AUTH=true
NEXT_PUBLIC_AAD_CLIENT_ID=$DEMO_APP_CLIENT_ID
NEXT_PUBLIC_AAD_TENANT_ID=$TENANT_ID

CCP4I2_REQUIRE_AUTH=true
AZURE_AD_TENANT_ID=$TENANT_ID
AZURE_AD_CLIENT_ID=$DEMO_APP_CLIENT_ID

# Access control: only members of this Team's M365 group
ALLOWED_AZURE_AD_GROUPS=$GROUP_ID

PLATFORM_ADMIN_EMAILS=$PLATFORM_ADMIN_EMAILS

# Compound registration IDs — override per instance if the group wants a
# distinct prefix (e.g. NCLP for Kawamura). Defaults match the main instance.
# COMPOUND_ID_PREFIX=NCL
# COMPOUND_ID_DIGITS=8
# First reg number on a fresh DB (default: 1; 26000 only preserves legacy NCL).
# COMPOUND_ID_START=1
EOF

# ── 6. Deploy applications ──────────────────────────────────────────────

step "Deploying applications"
"$SCRIPT_DIR/deploy-applications.sh" --env ".env.$INSTANCE"

# ── 7. App registration fixes ───────────────────────────────────────────

step "Patching ccp4i2-demo app registration"

# 7a. groupMembershipClaims = All (so M365 group IDs appear in JWTs)
CURRENT_CLAIMS=$(az ad app show --id "$DEMO_APP_CLIENT_ID" --query groupMembershipClaims -o tsv)
if [ "$CURRENT_CLAIMS" != "All" ]; then
    note "Setting groupMembershipClaims=All (was: ${CURRENT_CLAIMS:-unset})"
    az rest --method PATCH \
        --uri "https://graph.microsoft.com/v1.0/applications/$DEMO_APP_OBJECT_ID" \
        --headers "Content-Type=application/json" \
        --body '{"groupMembershipClaims": "All"}'
else
    note "groupMembershipClaims already set to All"
fi

# 7b. Assign the Team's group to the enterprise application (so members get tokens)
EXISTING_ASSIGN=$(az rest --method GET \
    --uri "https://graph.microsoft.com/v1.0/servicePrincipals/$DEMO_APP_SP_ID/appRoleAssignedTo" \
    --query "value[?principalId=='$GROUP_ID'].id" -o tsv 2>/dev/null || echo "")
if [ -z "$EXISTING_ASSIGN" ]; then
    note "Assigning group $GROUP_NAME to ccp4i2-demo enterprise app"
    az rest --method POST \
        --uri "https://graph.microsoft.com/v1.0/servicePrincipals/$DEMO_APP_SP_ID/appRoleAssignments" \
        --headers "Content-Type=application/json" \
        --body "{\"principalId\": \"$GROUP_ID\", \"resourceId\": \"$DEMO_APP_SP_ID\", \"appRoleId\": \"00000000-0000-0000-0000-000000000000\"}" \
        -o none
else
    note "Group already assigned to enterprise app"
fi

# 7c. Add the new web FQDN to the SPA redirect URIs
WEB_FQDN=$(az containerapp show -n "${PREFIX}-web" -g "$RG" --query properties.configuration.ingress.fqdn -o tsv)
NEW_URI="https://$WEB_FQDN/auth/callback"
EXISTING_URIS=$(az ad app show --id "$DEMO_APP_CLIENT_ID" --query "spa.redirectUris" -o json)
if ! echo "$EXISTING_URIS" | grep -q "$WEB_FQDN"; then
    note "Adding SPA redirect URI: $NEW_URI"
    PATCH_BODY=$(echo "$EXISTING_URIS" | python3 -c "
import json, sys
uris = json.load(sys.stdin)
uris.append('$NEW_URI')
print(json.dumps({'spa': {'redirectUris': uris}}))
")
    az rest --method PATCH \
        --uri "https://graph.microsoft.com/v1.0/applications/$DEMO_APP_OBJECT_ID" \
        --headers "Content-Type=application/json" \
        --body "$PATCH_BODY" \
        -o none
else
    note "SPA redirect URI already present"
fi

# ── Done ────────────────────────────────────────────────────────────────

echo ""
step "Done"
echo "Instance:  $INSTANCE"
echo "RG:        $RG"
echo "Team:      $GROUP_NAME ($GROUP_ID)"
echo "Web URL:   https://$WEB_FQDN"
echo "Env file:  $ENV_FILE"
echo ""
echo "Next:"
echo "  - Wait ~1 min for DNS + app-reg propagation."
echo "  - Sign in at https://$WEB_FQDN as any member of '$GROUP_NAME'."
echo "  - Bump IMAGE_TAG_* in $ENV_FILE and redeploy when the shared ACR gets a new image."
