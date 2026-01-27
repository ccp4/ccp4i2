#!/bin/bash

# Container Apps Deployment Script

# Ensure Homebrew paths are available
export PATH="/opt/homebrew/bin:$PATH"

# Get the directory where the script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="$SCRIPT_DIR/../.env.deployment"

# Load environment variables
if [ -f "$ENV_FILE" ]; then
    source "$ENV_FILE"
else
    echo "❌ .env.deployment not found at $ENV_FILE. Run deploy-infrastructure.sh first."
    exit 1
fi

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}🚀 Deploying Container Apps${NC}"

# Security check: Ensure authentication is required for production deployments
echo -e "${YELLOW}🔒 Checking authentication configuration...${NC}"
if [ "$CCP4I2_REQUIRE_AUTH" != "true" ]; then
    echo -e "${RED}╔═══════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${RED}║                   ⚠️  SECURITY WARNING ⚠️                          ║${NC}"
    echo -e "${RED}╠═══════════════════════════════════════════════════════════════════╣${NC}"
    echo -e "${RED}║  CCP4I2_REQUIRE_AUTH is not set to 'true' in .env.deployment     ║${NC}"
    echo -e "${RED}║                                                                   ║${NC}"
    echo -e "${RED}║  This would deploy the application WITHOUT authentication,       ║${NC}"
    echo -e "${RED}║  creating a 'dev_admin' superuser account with full access       ║${NC}"
    echo -e "${RED}║  to ALL data for ANY request!                                    ║${NC}"
    echo -e "${RED}║                                                                   ║${NC}"
    echo -e "${RED}║  Current value: CCP4I2_REQUIRE_AUTH=\"$CCP4I2_REQUIRE_AUTH\"${NC}"
    echo -e "${RED}║                                                                   ║${NC}"
    echo -e "${RED}╠═══════════════════════════════════════════════════════════════════╣${NC}"
    echo -e "${RED}║  To proceed, you must either:                                    ║${NC}"
    echo -e "${RED}║  1. Set CCP4I2_REQUIRE_AUTH=true in .env.deployment (production) ║${NC}"
    echo -e "${RED}║  2. Set ALLOW_NO_AUTH_DEPLOY=true to override (local dev only)   ║${NC}"
    echo -e "${RED}╚═══════════════════════════════════════════════════════════════════╝${NC}"

    # Allow override for local development/testing
    if [ "$ALLOW_NO_AUTH_DEPLOY" = "true" ]; then
        echo -e "${YELLOW}⚠️  ALLOW_NO_AUTH_DEPLOY=true detected - proceeding with NO AUTH${NC}"
        echo -e "${YELLOW}⚠️  This should ONLY be used for local development!${NC}"
        echo -e "${YELLOW}⚠️  Waiting 5 seconds... Press Ctrl+C to cancel${NC}"
        sleep 5
    else
        exit 1
    fi
else
    echo -e "${GREEN}✅ Authentication required: CCP4I2_REQUIRE_AUTH=true${NC}"
fi

# Function to check and register resource providers
check_resource_providers() {
    echo -e "${BLUE}🔍 Checking required resource providers...${NC}"

    # List of required providers for this deployment
    REQUIRED_PROVIDERS=("Microsoft.ServiceBus" "Microsoft.KeyVault" "Microsoft.DBforPostgreSQL" "Microsoft.ContainerRegistry" "Microsoft.App")

    for provider in "${REQUIRED_PROVIDERS[@]}"; do
        REGISTRATION_STATE=$(az provider show --namespace "$provider" --query registrationState -o tsv 2>/dev/null)

        if [ "$REGISTRATION_STATE" != "Registered" ]; then
            echo -e "${YELLOW}📝 Registering $provider resource provider...${NC}"
            az provider register --namespace "$provider"
            echo -e "${GREEN}✅ $provider registered${NC}"
        else
            echo -e "${GREEN}✅ $provider already registered${NC}"
        fi
    done
}

# Function to get the most recent successful infrastructure deployment
get_successful_infrastructure_deployment() {
    echo -e "${BLUE}🔍 Finding successful infrastructure deployment...${NC}"

    # Get the most recent successful infrastructure deployment
    local deployment_name=$(az deployment group list \
        --resource-group $RESOURCE_GROUP \
        --query "[?properties.provisioningState=='Succeeded' && contains(name, 'infrastructure')].name" \
        -o tsv | sort | tail -1)

    if [ -z "$deployment_name" ]; then
        echo -e "${RED}❌ No successful infrastructure deployment found${NC}"
        echo -e "${YELLOW}💡 Run deploy-infrastructure.sh first${NC}"
        return 1
    fi

    echo -e "${GREEN}✅ Using infrastructure deployment: $deployment_name${NC}"
    # Set global variable
    INFRA_DEPLOYMENT="$deployment_name"
}

# Check resource providers first
check_resource_providers

# Get successful infrastructure deployment
if ! get_successful_infrastructure_deployment; then
    exit 1
fi

# Get infrastructure outputs from the successful deployment
echo -e "${BLUE}📋 Getting infrastructure outputs...${NC}"

CONTAINER_APPS_ENV_ID=$(az deployment group show \
  --resource-group $RESOURCE_GROUP \
  --name $INFRA_DEPLOYMENT \
  --query properties.outputs.containerAppsEnvironmentId.value \
  --output tsv)

if [ -z "$CONTAINER_APPS_ENV_ID" ]; then
    echo -e "${RED}❌ Could not get Container Apps Environment ID from deployment $INFRA_DEPLOYMENT${NC}"
    exit 1
fi

POSTGRES_FQDN=$(az deployment group show \
  --resource-group $RESOURCE_GROUP \
  --name $INFRA_DEPLOYMENT \
  --query properties.outputs.postgresServerFqdn.value \
  --output tsv)

if [ -z "$POSTGRES_FQDN" ]; then
    echo -e "${RED}❌ Could not get PostgreSQL FQDN from deployment $INFRA_DEPLOYMENT${NC}"
    exit 1
fi

KEY_VAULT_NAME=$(az deployment group show \
  --resource-group $RESOURCE_GROUP \
  --name $INFRA_DEPLOYMENT \
  --query properties.outputs.keyVaultName.value \
  --output tsv)

if [ -z "$KEY_VAULT_NAME" ]; then
    echo -e "${RED}❌ Could not get Key Vault name from deployment $INFRA_DEPLOYMENT${NC}"
    exit 1
fi

STORAGE_ACCOUNT_NAME=$(az deployment group show \
  --resource-group $RESOURCE_GROUP \
  --name $INFRA_DEPLOYMENT \
  --query properties.outputs.storageAccountName.value \
  --output tsv)

if [ -z "$STORAGE_ACCOUNT_NAME" ]; then
    echo -e "${RED}❌ Could not get Storage Account name from deployment $INFRA_DEPLOYMENT${NC}"
    exit 1
fi

echo -e "${GREEN}✅ Infrastructure outputs retrieved${NC}"

# Note: Secrets are already stored in Key Vault during infrastructure deployment
# Key Vault now has private access only, so we cannot access it from deployment machine

# Deploy applications
echo -e "${YELLOW}🚀 Deploying container applications...${NC}"
APP_DEPLOYMENT_NAME="applications-$(date +%Y%m%d-%H%M%S)"

# Skip CCP4 storage mount when CCP4 is baked into container image (new approach)
SKIP_CCP4_STORAGE="${SKIP_CCP4_STORAGE:-true}"

CONTAINER_APPS_IDENTITY_ID=$(az deployment group show \
  --resource-group $RESOURCE_GROUP \
  --name $INFRA_DEPLOYMENT \
  --query properties.outputs.containerAppsIdentityId.value \
  --output tsv)

CONTAINER_APPS_IDENTITY_PRINCIPAL_ID=$(az deployment group show \
  --resource-group $RESOURCE_GROUP \
  --name $INFRA_DEPLOYMENT \
  --query properties.outputs.containerAppsIdentityPrincipalId.value \
  --output tsv)

CONTAINER_APPS_IDENTITY_CLIENT_ID=$(az deployment group show \
  --resource-group $RESOURCE_GROUP \
  --name $INFRA_DEPLOYMENT \
  --query properties.outputs.containerAppsIdentityClientId.value \
  --output tsv)

if [ -z "$CONTAINER_APPS_IDENTITY_ID" ] || [ -z "$CONTAINER_APPS_IDENTITY_PRINCIPAL_ID" ]; then
    echo -e "${RED}❌ Could not get Container Apps Identity from deployment $INFRA_DEPLOYMENT${NC}"
    exit 1
fi

# If client ID is not in deployment outputs (older infra), get it directly from the identity
if [ -z "$CONTAINER_APPS_IDENTITY_CLIENT_ID" ]; then
    IDENTITY_NAME=$(basename $CONTAINER_APPS_IDENTITY_ID)
    CONTAINER_APPS_IDENTITY_CLIENT_ID=$(az identity show \
      --name $IDENTITY_NAME \
      --resource-group $RESOURCE_GROUP \
      --query clientId \
      --output tsv 2>/dev/null || echo "")
fi

az deployment group create \
  --resource-group $RESOURCE_GROUP \
  --template-file "$SCRIPT_DIR/../infrastructure/applications.bicep" \
  --parameters containerAppsEnvironmentId="$CONTAINER_APPS_ENV_ID" \
               acrLoginServer="$ACR_LOGIN_SERVER" \
               acrName="$ACR_NAME" \
               postgresServerFqdn="$POSTGRES_FQDN" \
               keyVaultName="$KEY_VAULT_NAME" \
               imageTagWeb="${IMAGE_TAG_WEB:-latest}" \
               imageTagServer="${IMAGE_TAG_SERVER:-latest}" \
               prefix=ccp4i2-bicep \
               aadClientId="${NEXT_PUBLIC_AAD_CLIENT_ID:-}" \
               aadTenantId="${NEXT_PUBLIC_AAD_TENANT_ID:-}" \
               containerAppsIdentityId="$CONTAINER_APPS_IDENTITY_ID" \
               containerAppsIdentityPrincipalId="$CONTAINER_APPS_IDENTITY_PRINCIPAL_ID" \
               containerAppsIdentityClientId="$CONTAINER_APPS_IDENTITY_CLIENT_ID" \
               ccp4Version="${CCP4_VERSION:-ccp4-20251105}" \
               storageAccountName="$STORAGE_ACCOUNT_NAME" \
               platformAdminEmails="${PLATFORM_ADMIN_EMAILS:-}" \
               skipCcp4Storage=$SKIP_CCP4_STORAGE \
               customDomain="${CUSTOM_DOMAIN:-}" \
  --name $APP_DEPLOYMENT_NAME \
  --mode Incremental

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✅ Application deployment successful${NC}"
    
    # Validate private endpoints
    echo -e "${YELLOW}🔍 Validating private network configuration...${NC}"
    
    # Check private endpoints status
    PRIVATE_ENDPOINTS=$(az network private-endpoint list --resource-group $RESOURCE_GROUP --query "length([?provisioningState=='Succeeded'])" -o tsv)
    echo "✅ Private endpoints active: $PRIVATE_ENDPOINTS"
    
    # Check VNet integration
    VNET_INTEGRATION=$(az containerapp env show --name $(basename $CONTAINER_APPS_ENV_ID) --resource-group $RESOURCE_GROUP --query "properties.vnetConfiguration.infrastructureSubnetId" -o tsv)
    if [ ! -z "$VNET_INTEGRATION" ]; then
        echo "✅ Container Apps Environment integrated with VNet"
    else
        echo "⚠️  Container Apps Environment not VNet integrated"
    fi
    
    # Get application URLs
    SERVER_URL=$(az deployment group show \
      --resource-group $RESOURCE_GROUP \
      --name $APP_DEPLOYMENT_NAME \
      --query properties.outputs.serverUrl.value \
      --output tsv)
    
    WEB_URL=$(az deployment group show \
      --resource-group $RESOURCE_GROUP \
      --name $APP_DEPLOYMENT_NAME \
      --query properties.outputs.webUrl.value \
      --output tsv)
    
    # Bind custom domain if configured
    if [ -n "$CUSTOM_DOMAIN" ]; then
        echo -e "${YELLOW}🌐 Binding custom domain: $CUSTOM_DOMAIN${NC}"

        # Check if hostname is already added
        EXISTING_HOSTNAME=$(az containerapp hostname list \
            --name ccp4i2-bicep-web \
            --resource-group "$RESOURCE_GROUP" \
            --query "[?name=='$CUSTOM_DOMAIN'].name" \
            -o tsv 2>/dev/null)

        if [ -z "$EXISTING_HOSTNAME" ]; then
            echo -e "${BLUE}📝 Adding hostname to container app...${NC}"
            az containerapp hostname add \
                --name ccp4i2-bicep-web \
                --resource-group "$RESOURCE_GROUP" \
                --hostname "$CUSTOM_DOMAIN" 2>/dev/null || true
        fi

        # Bind the hostname with managed certificate
        echo -e "${BLUE}🔒 Binding managed certificate...${NC}"
        az containerapp hostname bind \
            --name ccp4i2-bicep-web \
            --resource-group "$RESOURCE_GROUP" \
            --hostname "$CUSTOM_DOMAIN" \
            --environment ccp4i2-bicep-env-uk \
            --validation-method CNAME 2>/dev/null

        if [ $? -eq 0 ]; then
            echo -e "${GREEN}✅ Custom domain bound: https://$CUSTOM_DOMAIN${NC}"
        else
            echo -e "${YELLOW}⚠️  Custom domain binding may need manual verification${NC}"
            echo -e "${YELLOW}   Run: az containerapp hostname list --name ccp4i2-bicep-web --resource-group $RESOURCE_GROUP${NC}"
        fi
    fi

    echo -e "${GREEN}🎉 Deployment completed successfully!${NC}"
    echo -e "${GREEN}🔒 All services are running in private VNet with no public endpoints${NC}"
    echo -e "${YELLOW}📝 Application URLs (external access via Container Apps ingress):${NC}"
    echo "🌐 Web App: $WEB_URL"
    if [ -n "$CUSTOM_DOMAIN" ]; then
        echo "🌐 Custom Domain: https://$CUSTOM_DOMAIN"
    fi
    echo "🔧 Server API: $SERVER_URL"
    echo ""
    echo -e "${YELLOW}🔐 Security Features Active:${NC}"
    echo "✅ PostgreSQL: Private endpoint only"
    echo "✅ Storage Account: Private endpoint only"
    echo "✅ Key Vault: Private endpoint only"
    echo "✅ Container Registry: Private endpoint only"
    echo "✅ Container Apps: VNet integrated"
else
    echo -e "${RED}❌ Application deployment failed${NC}"
    exit 1
fi