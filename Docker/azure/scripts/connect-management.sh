#!/bin/bash

# Quick access script for management container
# This provides an interactive shell with all environment variables,
# volume mounts, and Key Vault secrets accessible

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Load environment
if [ -f .env.deployment ]; then
    source .env.deployment
else
    echo "❌ .env.deployment not found. Run deploy-infrastructure.sh first."
    exit 1
fi

MANAGEMENT_APP_NAME="ccp4i2-bicep-management"

echo -e "${BLUE}🔧 Connecting to management container...${NC}"
echo -e "${YELLOW}💡 You will have access to:${NC}"
echo -e "   • /mnt/ccp4data - CCP4 data and installations"
echo -e "   • /mnt/staticfiles - Django static files"
echo -e "   • /mnt/mediafiles - Django media files"
echo -e "   • All environment variables and Key Vault secrets"
echo -e "   • Python, Django, and development tools"
echo ""
echo -e "${GREEN}📦 Container: $MANAGEMENT_APP_NAME${NC}"
echo ""

# Connect to the container
az containerapp exec \
  --resource-group $RESOURCE_GROUP \
  --name $MANAGEMENT_APP_NAME \
  --command bash
