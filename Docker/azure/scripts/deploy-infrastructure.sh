#!/bin/bash

# Azure Infrastructure Deployment Script using Bicep

# Ensure Homebrew paths are available
export PATH="/opt/homebrew/bin:$PATH"

# Get the directory where the script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BICEP_DIR="$(dirname "$SCRIPT_DIR")"

# Configuration
RESOURCE_GROUP="ccp4i2-bicep-rg-ne"
LOCATION="northeurope"
SUBSCRIPTION_ID=$(az account show --query id -o tsv)

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}🏗️ Deploying Infrastructure with Bicep${NC}"

# Create resource group first (needed to check for existing resources)
echo -e "${YELLOW}📁 Creating/verifying resource group...${NC}"
az group create --name $RESOURCE_GROUP --location $LOCATION

# Check if PostgreSQL server already exists
echo -e "${YELLOW}🔍 Checking for existing PostgreSQL server...${NC}"
EXISTING_POSTGRES=$(az postgres flexible-server list \
    --resource-group $RESOURCE_GROUP \
    --query "[?contains(name, 'ccp4i2-bicep-db')].name" \
    -o tsv 2>/dev/null)

if [ -n "$EXISTING_POSTGRES" ]; then
    echo -e "${GREEN}✅ Found existing PostgreSQL server: $EXISTING_POSTGRES${NC}"
    SKIP_POSTGRES="true"
    echo -e "${YELLOW}🔑 Retrieving existing password from Key Vault...${NC}"

    # Get existing Key Vault name
    EXISTING_KV=$(az keyvault list --resource-group $RESOURCE_GROUP --query "[0].name" -o tsv 2>/dev/null)

    if [ -n "$EXISTING_KV" ]; then
        # Temporarily enable access to retrieve secret
        CURRENT_IP=$(curl -4 -s ifconfig.me 2>/dev/null || curl -s ipv4.icanhazip.com 2>/dev/null)
        az keyvault network-rule add --name $EXISTING_KV --ip-address $CURRENT_IP --output none 2>/dev/null || true
        az keyvault update --name $EXISTING_KV --public-network-access Enabled --output none 2>/dev/null || true
        sleep 5

        # Try to get existing password
        DB_PASSWORD=$(az keyvault secret show --vault-name $EXISTING_KV --name db-password --query value -o tsv 2>/dev/null)

        # Restore Key Vault access
        az keyvault network-rule remove --name $EXISTING_KV --ip-address $CURRENT_IP --output none 2>/dev/null || true
        az keyvault update --name $EXISTING_KV --public-network-access Disabled --output none 2>/dev/null || true

        if [ -n "$DB_PASSWORD" ]; then
            echo -e "${GREEN}✅ Retrieved existing database password from Key Vault${NC}"
        else
            echo -e "${RED}⚠️  Could not retrieve password from Key Vault${NC}"
            echo -e "${RED}⚠️  Database exists but password unknown - deployment may fail${NC}"
            echo -e "${YELLOW}Using placeholder password - you may need to reset it manually${NC}"
            DB_PASSWORD="PLACEHOLDER_EXISTING_DB"
        fi
    else
        echo -e "${RED}⚠️  No Key Vault found - using placeholder password${NC}"
        DB_PASSWORD="PLACEHOLDER_EXISTING_DB"
    fi
else
    echo -e "${YELLOW}📝 No existing PostgreSQL server - generating new password${NC}"
    SKIP_POSTGRES="false"
    DB_PASSWORD=$(openssl rand -base64 16)
fi

# Check for and purge any existing soft-deleted Key Vault
echo -e "${YELLOW}🔍 Checking for soft-deleted Key Vaults...${NC}"
DELETED_VAULTS=$(az keyvault list-deleted --query "[?contains(name, 'ccp4i2-bicep-kv-ne')]" -o tsv)
if [ ! -z "$DELETED_VAULTS" ]; then
    echo -e "${YELLOW}🗑️ Purging soft-deleted Key Vaults...${NC}"
    az keyvault list-deleted --query "[?contains(name, 'ccp4i2-bicep-kv-ne')].name" -o tsv | while read vault_name; do
        echo "Purging vault: $vault_name"
        az keyvault purge --name "$vault_name" --location $LOCATION || true
    done
    sleep 10
fi

# Deploy infrastructure
echo -e "${YELLOW}🚀 Deploying infrastructure...${NC}"
INFRA_DEPLOYMENT_NAME="infrastructure-$(date +%Y%m%d-%H%M%S)"

# Capture deployment output and errors
DEPLOYMENT_OUTPUT=$(az deployment group create \
  --resource-group $RESOURCE_GROUP \
  --template-file "$BICEP_DIR/infrastructure/infrastructure.bicep" \
  --parameters location=$LOCATION \
               prefix=ccp4i2-bicep \
               environment=ne \
               postgresAdminPassword="$DB_PASSWORD" \
               skipPostgresDeployment=$SKIP_POSTGRES \
  --name $INFRA_DEPLOYMENT_NAME \
  --mode Incremental 2>&1)

DEPLOYMENT_EXIT_CODE=$?

if [ $DEPLOYMENT_EXIT_CODE -eq 0 ]; then
    echo -e "${GREEN}✅ Infrastructure deployment successful${NC}"
    
    # Get outputs
    ACR_NAME=$(az deployment group show \
      --resource-group $RESOURCE_GROUP \
      --name $INFRA_DEPLOYMENT_NAME \
      --query properties.outputs.acrName.value \
      --output tsv)
    
    ACR_LOGIN_SERVER=$(az deployment group show \
      --resource-group $RESOURCE_GROUP \
      --name $INFRA_DEPLOYMENT_NAME \
      --query properties.outputs.acrLoginServer.value \
      --output tsv)
      
    KEY_VAULT_NAME=$(az deployment group show \
      --resource-group $RESOURCE_GROUP \
      --name $INFRA_DEPLOYMENT_NAME \
      --query properties.outputs.keyVaultName.value \
      --output tsv)
    
    # Temporarily enable public access to Key Vault for secret storage
    echo -e "${YELLOW}🔑 Temporarily enabling Key Vault public access for initial setup...${NC}"
    
    # Get current public IP and add to Key Vault firewall
    CURRENT_IP=$(curl -4 -s ifconfig.me || curl -s ipv4.icanhazip.com)
    echo -e "${YELLOW}🌐 Adding current IP ($CURRENT_IP) to Key Vault firewall...${NC}"
    az keyvault network-rule add --name $KEY_VAULT_NAME --ip-address $CURRENT_IP
    az keyvault update --name $KEY_VAULT_NAME --public-network-access Enabled
    
    # Get current user Object ID and assign Key Vault Secrets Officer role
    CURRENT_USER_ID=$(az ad signed-in-user show --query id -o tsv)
    echo -e "${YELLOW}🔐 Assigning Key Vault Secrets Officer role to current user...${NC}"
    az role assignment create \
      --assignee $CURRENT_USER_ID \
      --role "Key Vault Secrets Officer" \
      --scope "/subscriptions/$SUBSCRIPTION_ID/resourceGroups/$RESOURCE_GROUP/providers/Microsoft.KeyVault/vaults/$KEY_VAULT_NAME"
    
    # Wait a moment for the changes to propagate
    sleep 30
    
    # Store PostgreSQL password in Key Vault (only if we have a real password)
    if [[ "$DB_PASSWORD" != "PLACEHOLDER_EXISTING_DB" ]]; then
        echo -e "${YELLOW}🔐 Storing PostgreSQL password in Key Vault...${NC}"
        az keyvault secret set \
          --vault-name $KEY_VAULT_NAME \
          --name db-password \
          --value "$DB_PASSWORD" \
          --output none
    else
        echo -e "${YELLOW}⚠️  Skipping password storage - using existing Key Vault secret${NC}"
    fi
    
    # Store ACR password in Key Vault
    ACR_PASSWORD=$(az acr credential show --name $ACR_NAME --query passwords[0].value -o tsv)
    az keyvault secret set \
      --vault-name $KEY_VAULT_NAME \
      --name registry-password \
      --value "$ACR_PASSWORD" \
      --output none
    
    # Generate and store Django secret key
    DJANGO_SECRET_KEY=$(openssl rand -base64 32)
    az keyvault secret set \
      --vault-name $KEY_VAULT_NAME \
      --name django-secret-key \
      --value "$DJANGO_SECRET_KEY" \
      --output none
    
    # Disable public access again and remove IP from firewall
    echo -e "${YELLOW}🔒 Disabling Key Vault public access...${NC}"
    az keyvault network-rule remove --name $KEY_VAULT_NAME --ip-address $CURRENT_IP
    az keyvault update --name $KEY_VAULT_NAME --public-network-access Disabled
    
    # Get shared identity outputs
    CONTAINER_APPS_IDENTITY_ID=$(az deployment group show \
      --resource-group $RESOURCE_GROUP \
      --name $INFRA_DEPLOYMENT_NAME \
      --query properties.outputs.containerAppsIdentityId.value \
      --output tsv)
    
    CONTAINER_APPS_IDENTITY_PRINCIPAL_ID=$(az deployment group show \
      --resource-group $RESOURCE_GROUP \
      --name $INFRA_DEPLOYMENT_NAME \
      --query properties.outputs.containerAppsIdentityPrincipalId.value \
      --output tsv)
    
    echo -e "${YELLOW}📝 Infrastructure Details:${NC}"
    echo "ACR Name: $ACR_NAME"
    echo "ACR Login Server: $ACR_LOGIN_SERVER"
    echo "Key Vault: $KEY_VAULT_NAME (private access only)"
    echo "Shared Identity ID: $CONTAINER_APPS_IDENTITY_ID"
    echo "Shared Identity Principal ID: $CONTAINER_APPS_IDENTITY_PRINCIPAL_ID"
    
    # Store outputs in environment file for application deployment
    cat > "$BICEP_DIR/.env.deployment" << EOF
ACR_NAME=$ACR_NAME
ACR_LOGIN_SERVER=$ACR_LOGIN_SERVER
RESOURCE_GROUP=$RESOURCE_GROUP
CONTAINER_APPS_IDENTITY_ID=$CONTAINER_APPS_IDENTITY_ID
CONTAINER_APPS_IDENTITY_PRINCIPAL_ID=$CONTAINER_APPS_IDENTITY_PRINCIPAL_ID
IMAGE_TAG_WEB=latest
IMAGE_TAG_SERVER=latest
EOF
    
    echo -e "${GREEN}✅ Infrastructure deployment completed. Environment saved to $BICEP_DIR/.env.deployment${NC}"
else
    echo -e "${RED}❌ Infrastructure deployment failed${NC}"
    echo -e "${RED}Error details:${NC}"
    echo "$DEPLOYMENT_OUTPUT"
    exit 1
fi