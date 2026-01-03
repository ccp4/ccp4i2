#!/bin/bash
#
# Complete CCP4i2 Django Azure Deployment
# This script automates the entire deployment from infrastructure to running application
#

set -e  # Exit on any error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
RESOURCE_GROUP="ccp4i2-bicep-rg-ne"
LOCATION="northeurope"
PREFIX="ccp4i2-bicep"

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  CCP4i2 Django Azure Deployment${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Check prerequisites
echo -e "${YELLOW}Checking prerequisites...${NC}"
command -v az >/dev/null 2>&1 || { echo -e "${RED}Error: Azure CLI not found${NC}"; exit 1; }
command -v jq >/dev/null 2>&1 || { echo -e "${RED}Error: jq not found${NC}"; exit 1; }

# Ensure logged in
az account show >/dev/null 2>&1 || { echo -e "${RED}Error: Not logged in to Azure CLI${NC}"; exit 1; }

SUBSCRIPTION=$(az account show --query name -o tsv)
echo -e "${GREEN}✓${NC} Logged in to Azure subscription: ${SUBSCRIPTION}"
echo ""

# Step 1: Deploy Infrastructure
echo -e "${BLUE}Step 1: Deploying Infrastructure${NC}"
echo "This will create:"
echo "  - Container Apps Environment"
echo "  - PostgreSQL Flexible Server"
echo "  - Key Vault"
echo "  - Storage Account with File Shares"
echo "  - Service Bus"
echo "  - Shared User-Assigned Managed Identity"
echo ""

cd "$(dirname "$0")"
./deploy-infrastructure.sh

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓${NC} Infrastructure deployed successfully"
else
    echo -e "${RED}✗${NC} Infrastructure deployment failed"
    exit 1
fi

# Get infrastructure outputs
echo ""
echo -e "${YELLOW}Retrieving infrastructure details...${NC}"
INFRA_OUTPUT=$(az deployment group show \
    --resource-group $RESOURCE_GROUP \
    --name infrastructure \
    --query properties.outputs -o json)

CONTAINER_APPS_ENV_ID=$(echo $INFRA_OUTPUT | jq -r '.containerAppsEnvironmentId.value')
IDENTITY_ID=$(echo $INFRA_OUTPUT | jq -r '.containerAppsIdentityId.value')
POSTGRES_FQDN=$(echo $INFRA_OUTPUT | jq -r '.postgresServerFqdn.value')
KEYVAULT_NAME=$(echo $INFRA_OUTPUT | jq -r '.keyVaultName.value')
ACR_LOGIN_SERVER=$(echo $INFRA_OUTPUT | jq -r '.acrLoginServer.value')

echo -e "${GREEN}✓${NC} Infrastructure details retrieved"
echo ""

# Step 2: Check CCP4 Installation
echo -e "${BLUE}Step 2: Checking CCP4 Installation${NC}"
echo "Checking if CCP4 is already extracted on file share..."

# Deploy management container if not exists
MGMT_EXISTS=$(az containerapp show \
    --name ${PREFIX}-management \
    --resource-group $RESOURCE_GROUP \
    --query name -o tsv 2>/dev/null || echo "")

if [ -z "$MGMT_EXISTS" ]; then
    echo "Deploying management container for checks..."
    ./deploy-management.sh
    sleep 10
fi

# Check if CCP4 setup script exists
CCP4_EXISTS=$(az containerapp exec \
    --name ${PREFIX}-management \
    --resource-group $RESOURCE_GROUP \
    --command "test -f /mnt/ccp4data/ccp4-9/bin/ccp4.setup-sh && echo 'YES' || echo 'NO'" 2>/dev/null | grep -o "YES" || echo "NO")

if [ "$CCP4_EXISTS" = "YES" ]; then
    echo -e "${GREEN}✓${NC} CCP4 installation found"
else
    echo -e "${YELLOW}⚠${NC} CCP4 not found. You need to:"
    echo "  1. Upload ccp4-9.tar.gz to /mnt/ccp4data/"
    echo "  2. Extract it: tar -xzf /mnt/ccp4data/ccp4-9.tar.gz -C /mnt/ccp4data/"
    echo ""
    read -p "Have you done this? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo -e "${RED}Deployment cannot continue without CCP4 installation${NC}"
        exit 1
    fi
fi
echo ""

# Step 3: Install Python Packages
echo -e "${BLUE}Step 3: Installing Python Packages${NC}"
echo "This will install:"
echo "  - Django 3.2.25 with compatible dependencies"
echo "  - NumPy 1.26.4 (critical for CCP4 compatibility)"
echo "  - SciPy 1.15.3, Pandas 2.3.3"
echo "  - Azure integration packages"
echo "  - Gunicorn, Uvicorn, and other web server packages"
echo ""

# Check if packages already installed
PKG_EXISTS=$(az containerapp exec \
    --name ${PREFIX}-management \
    --resource-group $RESOURCE_GROUP \
    --command "test -d /mnt/ccp4data/py-packages/numpy && echo 'YES' || echo 'NO'" 2>/dev/null | grep -o "YES" || echo "NO")

if [ "$PKG_EXISTS" = "YES" ]; then
    echo -e "${YELLOW}⚠${NC} Packages already exist in /mnt/ccp4data/py-packages/"
    read -p "Reinstall packages? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Cleaning existing packages..."
        az containerapp exec \
            --name ${PREFIX}-management \
            --resource-group $RESOURCE_GROUP \
            --command "rm -rf /mnt/ccp4data/py-packages" 2>/dev/null || true
    else
        echo -e "${GREEN}✓${NC} Skipping package installation"
        SKIP_PACKAGES=true
    fi
fi

if [ "$SKIP_PACKAGES" != "true" ]; then
    # Deploy maintenance job if not exists
    JOB_EXISTS=$(az containerapp job show \
        --name ${PREFIX}-maintenance-job \
        --resource-group $RESOURCE_GROUP \
        --query name -o tsv 2>/dev/null || echo "")

    if [ -z "$JOB_EXISTS" ]; then
        echo "Deploying maintenance job..."
        ./deploy-maintenance-job.sh
    fi

    echo "Starting package installation job..."
    EXECUTION_NAME="${PREFIX}-maintenance-$(date +%s)"
    az containerapp job start \
        --name ${PREFIX}-maintenance-job \
        --resource-group $RESOURCE_GROUP \
        >/dev/null 2>&1

    echo "Waiting for package installation to complete (this may take 3-5 minutes)..."
    
    # Wait for job completion
    for i in {1..60}; do
        sleep 5
        STATUS=$(az containerapp job execution show \
            --name ${PREFIX}-maintenance-job \
            --resource-group $RESOURCE_GROUP \
            --query "properties.status" -o tsv 2>/dev/null || echo "Unknown")
        
        if [ "$STATUS" = "Succeeded" ]; then
            echo -e "${GREEN}✓${NC} Package installation completed successfully"
            break
        elif [ "$STATUS" = "Failed" ]; then
            echo -e "${RED}✗${NC} Package installation failed"
            echo "Check logs with: az containerapp job logs show --name ${PREFIX}-maintenance-job --resource-group $RESOURCE_GROUP"
            exit 1
        fi
        
        echo -n "."
    done
    echo ""
fi
echo ""

# Step 4: Deploy Applications
echo -e "${BLUE}Step 4: Deploying Applications${NC}"
echo "This will deploy:"
echo "  - Server (Django API)"
echo "  - Worker (Background job processor)"
echo "  - Web (Next.js frontend)"
echo ""

./deploy-applications.sh

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓${NC} Applications deployed successfully"
else
    echo -e "${RED}✗${NC} Application deployment failed"
    exit 1
fi
echo ""

# Step 5: Verify Deployment
echo -e "${BLUE}Step 5: Verifying Deployment${NC}"
echo "Waiting for server to start..."
sleep 30

SERVER_FQDN=$(az containerapp show \
    --name ${PREFIX}-server \
    --resource-group $RESOURCE_GROUP \
    --query properties.configuration.ingress.fqdn -o tsv)

echo "Testing health endpoint..."
HEALTH_STATUS=$(curl -s -o /dev/null -w "%{http_code}" "https://${SERVER_FQDN}/health/" 2>/dev/null || echo "000")

if [ "$HEALTH_STATUS" = "200" ]; then
    echo -e "${GREEN}✓${NC} Health check passed"
else
    echo -e "${RED}✗${NC} Health check failed (HTTP $HEALTH_STATUS)"
    echo "Check logs with: az containerapp logs show --name ${PREFIX}-server --resource-group $RESOURCE_GROUP --tail 50"
fi

echo ""
echo "Testing projects API..."
PROJECTS_STATUS=$(curl -s -o /dev/null -w "%{http_code}" "https://${SERVER_FQDN}/projects/" 2>/dev/null || echo "000")

if [ "$PROJECTS_STATUS" = "200" ]; then
    echo -e "${GREEN}✓${NC} Projects API responding"
    echo ""
    echo "Sample response:"
    curl -s "https://${SERVER_FQDN}/projects/" | jq '.[0] // empty' 2>/dev/null || echo "(No projects yet)"
else
    echo -e "${RED}✗${NC} Projects API failed (HTTP $PROJECTS_STATUS)"
    echo "Check logs with: az containerapp logs show --name ${PREFIX}-server --resource-group $RESOURCE_GROUP --tail 50"
fi

echo ""
echo -e "${BLUE}========================================${NC}"
echo -e "${GREEN}  Deployment Complete!${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo "Application URLs:"
echo "  Server API: https://${SERVER_FQDN}"
echo "  Web Frontend: https://$(az containerapp show --name ${PREFIX}-web --resource-group $RESOURCE_GROUP --query properties.configuration.ingress.fqdn -o tsv)"
echo ""
echo "Useful commands:"
echo "  View server logs:"
echo "    az containerapp logs show --name ${PREFIX}-server --resource-group $RESOURCE_GROUP --tail 50 --follow"
echo ""
echo "  Check application status:"
echo "    az containerapp show --name ${PREFIX}-server --resource-group $RESOURCE_GROUP --query properties.runningStatus"
echo ""
echo "  Restart server:"
echo "    az containerapp update --name ${PREFIX}-server --resource-group $RESOURCE_GROUP --set-env-vars RESTART=\$(date +%s)"
echo ""
