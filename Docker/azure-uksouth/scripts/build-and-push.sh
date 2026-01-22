#!/bin/bash

# Docker Build and Push Script
# Usage: ./build-and-push.sh [server|web]
#   - No argument: Build both server and web images (default)
#   - server: Build only server image
#   - web: Build only web image

# Ensure Homebrew paths are available
export PATH="/opt/homebrew/bin:$PATH"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Error handling function
cleanup_on_error() {
    echo -e "${RED}❌ Build failed. Cleaning up...${NC}"
    if [ ! -z "$ACR_NAME" ] && [ ! -z "$ORIGINAL_DEFAULT_ACTION" ]; then
        echo -e "${YELLOW}🔧 Restoring ACR network security settings...${NC}"
        az acr update --name $ACR_NAME --default-action $ORIGINAL_DEFAULT_ACTION 2>/dev/null || true
    fi
    if [ ! -z "$ACR_NAME" ]; then
        echo -e "${YELLOW}📋 Current ACR configuration:${NC}"
        az acr show --name $ACR_NAME --query "{publicNetworkAccess: publicNetworkAccess, defaultAction: networkRuleSet.defaultAction, allowedIPs: networkRuleSet.ipRules[].ipAddressOrRange}" -o table 2>/dev/null || true
    fi
    exit 1
}

# Set up error handling (will be disabled for IP detection)
trap cleanup_on_error ERR

# Change to the correct working directory (relative to bicep directory)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BICEP_DIR="$(dirname "$SCRIPT_DIR")"

# The source code directory is the parent of Docker/azure
# i.e., Docker/azure/scripts -> Docker/azure -> Docker -> ccp4i2 (project root)
WORKING_DIR="$(cd "$BICEP_DIR/../.." && pwd)"

# Parse command line arguments
BUILD_SERVER=true
BUILD_WEB=true

if [ $# -gt 0 ]; then
    case "$1" in
        "server")
            BUILD_WEB=false
            echo -e "${YELLOW}🔧 Building server image only${NC}"
            ;;
        "web")
            BUILD_SERVER=false
            echo -e "${YELLOW}🔧 Building web image only${NC}"
            ;;
        *)
            echo -e "${RED}❌ Invalid argument. Use 'server', 'web', or no argument for both${NC}"
            exit 1
            ;;
    esac
else
    echo -e "${YELLOW}🔧 Building both server and web images${NC}"
fi

echo -e "${GREEN}🐳 Building and pushing Docker images${NC}"
echo -e "${YELLOW}📁 Working directory: $WORKING_DIR${NC}"

# Check if working directory exists and has server/Dockerfile
if [ -z "$WORKING_DIR" ] || [ ! -d "$WORKING_DIR" ]; then
    echo -e "${RED}❌ Working directory not found or invalid: $WORKING_DIR${NC}"
    exit 1
fi

if [ ! -f "$WORKING_DIR/Docker/server/Dockerfile" ]; then
    echo -e "${RED}❌ Docker/server/Dockerfile not found in $WORKING_DIR${NC}"
    exit 1
fi

cd "$WORKING_DIR"

# Load environment variables from bicep directory
if [ -f "$BICEP_DIR/.env.deployment" ]; then
    source "$BICEP_DIR/.env.deployment"
else
    echo "❌ .env.deployment not found. Run deploy-infrastructure.sh first."
    exit 1
fi

# Login to ACR
echo -e "${YELLOW}🔑 Logging into Azure Container Registry...${NC}"

# Get current public IP address
echo -e "${YELLOW}🌐 Getting current public IP address...${NC}"
set +e  # Disable exit on error for IP detection
CURRENT_IP=$(curl -s https://ipinfo.io/ip 2>/dev/null || curl -s https://api.ipify.org 2>/dev/null || echo "unknown")
set -e  # Re-enable exit on error

if [ "$CURRENT_IP" = "unknown" ]; then
    echo -e "${RED}❌ Could not determine public IP address${NC}"
    echo -e "${YELLOW}⚠️  You may need to manually add your IP to ACR firewall${NC}"
    echo -e "${YELLOW}⚠️  Continuing with ACR login attempt...${NC}"
else
    echo -e "${GREEN}📍 Current public IP: $CURRENT_IP${NC}"
    
    # Check if current IP is already in ACR firewall rules
    echo -e "${YELLOW}🔍 Checking ACR firewall rules...${NC}"
    set +e  # Disable exit on error for this check
    EXISTING_IP=$(az acr show --name $ACR_NAME --query "networkRuleSet.ipRules[?contains(ipAddressOrRange, '$CURRENT_IP')].ipAddressOrRange" -o tsv 2>/dev/null)
    set -e  # Re-enable exit on error
    
    if [ -z "$EXISTING_IP" ]; then
        echo -e "${YELLOW}➕ Adding current IP to ACR firewall...${NC}"
        az acr network-rule add --name $ACR_NAME --ip-address $CURRENT_IP
        echo -e "${GREEN}✅ IP address added to ACR firewall${NC}"
    else
        echo -e "${GREEN}✅ Current IP already allowed in ACR firewall${NC}"
    fi
fi

# Ensure public access is enabled for ACR builds
echo -e "${YELLOW}🌐 Ensuring ACR public access is enabled...${NC}"
az acr update --name $ACR_NAME --public-network-enabled true

# Temporarily disable network restrictions for Azure build agents
echo -e "${YELLOW}🔧 Temporarily allowing all access for ACR builds...${NC}"
ORIGINAL_DEFAULT_ACTION=$(az acr show --name $ACR_NAME --query "networkRuleSet.defaultAction" -o tsv)
echo -e "${YELLOW}📝 Original default action: $ORIGINAL_DEFAULT_ACTION${NC}"

# Set default action to Allow to bypass firewall during builds
az acr update --name $ACR_NAME --default-action Allow

echo -e "${YELLOW}🔐 Logging into ACR...${NC}"
set +e  # Disable exit on error for ACR login
az acr login --name $ACR_NAME
LOGIN_RESULT=$?
set -e  # Re-enable exit on error

if [ $LOGIN_RESULT -ne 0 ]; then
    echo -e "${RED}❌ ACR login failed${NC}"
    echo -e "${YELLOW}💡 Try running 'az login' to refresh authentication${NC}"
    echo -e "${YELLOW}💡 Or check if your IP is allowed in ACR firewall rules${NC}"
    exit 1
fi

echo -e "${GREEN}✅ ACR login successful${NC}"

# Generate image tags
TIMESTAMP=$(date +%Y%m%d-%H%M%S)
IMAGE_TAG_WEB=""
IMAGE_TAG_SERVER=""

if [ "$BUILD_WEB" = true ]; then
    IMAGE_TAG_WEB=$TIMESTAMP
fi

if [ "$BUILD_SERVER" = true ]; then
    IMAGE_TAG_SERVER=$TIMESTAMP
fi

# Helper function to create filtered context tarball
# az acr build ignores .dockerignore when creating its upload archive,
# so we create a pre-filtered tarball and upload it to blob storage
#
# Note: tar doesn't support Docker's ! negation syntax, so we use explicit
# --exclude patterns instead of --exclude-from=.dockerignore
create_filtered_context() {
    local tarball_path="$1"
    echo -e "${YELLOW}📦 Creating filtered build context (respecting .dockerignore)...${NC}" >&2

    # Create tarball with explicit exclude patterns
    # This mirrors .dockerignore but uses tar-compatible syntax
    # Key difference: Instead of "Docker/azure-uksouth/" with !exceptions,
    # we explicitly exclude subdirectories we DON'T want (scripts/, bicep/, etc.)
    tar -czf "$tarball_path" \
        --exclude='node_modules' \
        --exclude='.next' \
        --exclude='.git' \
        --exclude='.gitignore' \
        --exclude='.vscode' \
        --exclude='.idea' \
        --exclude='__pycache__' \
        --exclude='*.pyc' \
        --exclude='*.pyo' \
        --exclude='*.egg-info' \
        --exclude='.eggs' \
        --exclude='dist' \
        --exclude='build' \
        --exclude='.pytest_cache' \
        --exclude='.DS_Store' \
        --exclude='*.swp' \
        --exclude='*.swo' \
        --exclude='*~' \
        --exclude='*.tsbuildinfo' \
        --exclude='.env' \
        --exclude='.env.*' \
        --exclude='*.local' \
        --exclude='server/ccp4i2/demo_data' \
        --exclude='Docker/data' \
        --exclude='./docs' \
        --exclude='./tests' \
        --exclude='**/test_data' \
        --exclude='Docker/azure/scripts' \
        --exclude='Docker/azure/bicep' \
        --exclude='Docker/azure/.env*' \
        --exclude='Docker/azure-uksouth/scripts' \
        --exclude='Docker/azure-uksouth/bicep' \
        --exclude='Docker/azure-uksouth/.env*' \
        --exclude='client/renderer/public/baby-gru/rota500-arg.data' \
        --exclude='client/renderer/public/baby-gru/rota500-lys.data' \
        --exclude='.dockerignore' \
        .

    # BSD tar's --exclude='./docs' matches 'docs' at any level, not just root
    # Explicitly add back the azure-uksouth docs that were incorrectly excluded
    echo -e "${YELLOW}📋 Adding azure-uksouth docs to tarball...${NC}" >&2
    gunzip "$tarball_path"
    local uncompressed="${tarball_path%.gz}"
    tar -rf "$uncompressed" Docker/azure-uksouth/docs/
    gzip "$uncompressed"

    local size=$(ls -lh "$tarball_path" | awk '{print $5}')
    echo -e "${GREEN}✅ Context tarball created: $size${NC}" >&2
}

# Helper function to upload tarball to blob storage and get SAS URL
upload_context_to_blob() {
    local tarball_path="$1"
    local blob_name="build-context-$(date +%s).tar.gz"
    local container_name="build-contexts"

    # Get storage account from environment or auto-discover from resource group
    # Prefer non-private storage accounts (those without 'prv' in the name)
    # Private storage accounts can't be accessed by ACR build agents
    local storage_account="${AZURE_STORAGE_ACCOUNT:-}"
    if [ -z "$storage_account" ]; then
        echo -e "${YELLOW}🔍 Auto-discovering storage account in $RESOURCE_GROUP...${NC}" >&2
        # First try to find a non-private storage account (without 'prv' in name)
        storage_account=$(az storage account list \
            --resource-group "$RESOURCE_GROUP" \
            --query "[?contains(name, 'prv')==\`false\`].name | [0]" \
            --output tsv 2>/dev/null)

        # Fall back to any storage account
        if [ -z "$storage_account" ]; then
            storage_account=$(az storage account list \
                --resource-group "$RESOURCE_GROUP" \
                --query "[0].name" \
                --output tsv 2>/dev/null)
        fi

        if [ -z "$storage_account" ]; then
            echo -e "${RED}❌ No storage account found in $RESOURCE_GROUP${NC}" >&2
            return 1
        fi
        echo -e "${GREEN}✅ Using storage account: $storage_account${NC}" >&2
    fi

    # Get storage account key for authentication
    local storage_key=$(az storage account keys list \
        --account-name "$storage_account" \
        --query "[0].value" \
        --output tsv 2>/dev/null)

    if [ -z "$storage_key" ]; then
        echo -e "${RED}❌ Failed to get storage account key${NC}" >&2
        return 1
    fi

    # Create container if it doesn't exist (ignore errors if it exists)
    az storage container create \
        --name "$container_name" \
        --account-name "$storage_account" \
        --account-key "$storage_key" \
        >/dev/null 2>&1 || true

    # Upload the tarball
    echo -e "${YELLOW}📤 Uploading context to blob storage...${NC}" >&2
    az storage blob upload \
        --account-name "$storage_account" \
        --container-name "$container_name" \
        --name "$blob_name" \
        --file "$tarball_path" \
        --account-key "$storage_key" \
        --overwrite \
        --only-show-errors >&2

    if [ $? -ne 0 ]; then
        echo -e "${RED}❌ Failed to upload context to blob storage${NC}" >&2
        return 1
    fi

    # Generate SAS URL valid for 1 hour
    local expiry=$(date -u -v+1H +%Y-%m-%dT%H:%MZ 2>/dev/null || date -u -d '+1 hour' +%Y-%m-%dT%H:%MZ)
    local sas_url=$(az storage blob generate-sas \
        --account-name "$storage_account" \
        --container-name "$container_name" \
        --name "$blob_name" \
        --permissions r \
        --expiry "$expiry" \
        --account-key "$storage_key" \
        --full-uri \
        --output tsv)

    echo "$sas_url"
}

# Build and push images
if [ "$BUILD_SERVER" = true ]; then
    echo -e "${YELLOW}🔨 Building server image (with bundled CCP4)...${NC}"

    # Use Dockerfile which has CCP4 baked into the base image
    # The base image should already exist in ACR (built by build-base-image.sh)
    # The arp-warp layer should also exist (built by build-arpwarp-image.sh)
    #
    # By default, uses the same ACR as the server image (from .env.deployment)
    # To use a different ACR (e.g., cross-region), set BASE_IMAGE_ACR explicitly
    BASE_IMAGE_ACR="${BASE_IMAGE_ACR:-$ACR_LOGIN_SERVER}"
    CCP4_VERSION="${CCP4_VERSION:-ccp4-20251105}"
    BASE_IMAGE_NAME="${BASE_IMAGE_NAME:-ccp4i2/base-arpwarp}"

    echo -e "${YELLOW}📦 Using base image from: ${BASE_IMAGE_ACR}/${BASE_IMAGE_NAME}:${CCP4_VERSION}${NC}"

    # Create filtered context tarball and upload to blob storage
    CONTEXT_TARBALL="/tmp/build-context-$$.tar.gz"
    create_filtered_context "$CONTEXT_TARBALL"

    CONTEXT_URL=$(upload_context_to_blob "$CONTEXT_TARBALL")
    if [ -z "$CONTEXT_URL" ]; then
        echo -e "${RED}❌ Failed to get context URL${NC}"
        rm -f "$CONTEXT_TARBALL"
        exit 1
    fi

    echo -e "${YELLOW}🔨 Starting ACR build with filtered context...${NC}"
    az acr build \
      --registry $ACR_NAME \
      --image ccp4i2/server:$IMAGE_TAG_SERVER \
      --image ccp4i2/server:latest \
      --file Docker/server/Dockerfile \
      --platform linux/amd64 \
      --build-arg ACR_LOGIN_SERVER=$BASE_IMAGE_ACR \
      --build-arg CCP4_VERSION=$CCP4_VERSION \
      --build-arg BASE_IMAGE_NAME=$BASE_IMAGE_NAME \
      "$CONTEXT_URL"

    BUILD_RESULT=$?
    rm -f "$CONTEXT_TARBALL"

    if [ $BUILD_RESULT -ne 0 ]; then
        echo -e "${RED}❌ Failed to build server image${NC}"
        exit 1
    fi
fi

if [ "$BUILD_WEB" = true ]; then
    echo -e "${YELLOW}🔨 Building web image...${NC}"
    # Build from repo root to include icons from server/ccp4i2/{qticons,svgicons}
    if [ -f "Docker/client/Dockerfile" ]; then
        # Create filtered context tarball and upload to blob storage
        CONTEXT_TARBALL="/tmp/build-context-$$.tar.gz"
        create_filtered_context "$CONTEXT_TARBALL"

        CONTEXT_URL=$(upload_context_to_blob "$CONTEXT_TARBALL")
        if [ -z "$CONTEXT_URL" ]; then
            echo -e "${RED}❌ Failed to get context URL${NC}"
            rm -f "$CONTEXT_TARBALL"
            exit 1
        fi

        echo -e "${YELLOW}🔨 Starting ACR build with filtered context...${NC}"
        az acr build \
          --registry $ACR_NAME \
          --image ccp4i2/web:$IMAGE_TAG_WEB \
          --image ccp4i2/web:latest \
          --file Docker/client/Dockerfile \
          --platform linux/amd64 \
          --build-arg NEXT_PUBLIC_AAD_CLIENT_ID=${NEXT_PUBLIC_AAD_CLIENT_ID:-""} \
          --build-arg NEXT_PUBLIC_AAD_TENANT_ID=${NEXT_PUBLIC_AAD_TENANT_ID:-""} \
          --build-arg NEXT_PUBLIC_REQUIRE_AUTH=${NEXT_PUBLIC_REQUIRE_AUTH:-"false"} \
          "$CONTEXT_URL"

        BUILD_RESULT=$?
        rm -f "$CONTEXT_TARBALL"

        if [ $BUILD_RESULT -ne 0 ]; then
            echo -e "${RED}❌ Failed to build web image${NC}"
            exit 1
        fi
    else
        echo -e "${YELLOW}⚠️  Docker/client/Dockerfile not found, skipping web build${NC}"
        BUILD_WEB=false
    fi
fi

# Note: nginx is not used in Container Apps architecture (built-in ingress)

# Update environment file with image tags
if [ "$BUILD_WEB" = true ]; then
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' "s/^IMAGE_TAG_WEB=.*/IMAGE_TAG_WEB=$IMAGE_TAG_WEB/" "$BICEP_DIR/.env.deployment"
    else
        sed -i "s/^IMAGE_TAG_WEB=.*/IMAGE_TAG_WEB=$IMAGE_TAG_WEB/" "$BICEP_DIR/.env.deployment"
    fi
    echo -e "${YELLOW}📝 Updated .env.deployment with IMAGE_TAG_WEB: $IMAGE_TAG_WEB${NC}"
fi

if [ "$BUILD_SERVER" = true ]; then
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' "s/^IMAGE_TAG_SERVER=.*/IMAGE_TAG_SERVER=$IMAGE_TAG_SERVER/" "$BICEP_DIR/.env.deployment"
    else
        sed -i "s/^IMAGE_TAG_SERVER=.*/IMAGE_TAG_SERVER=$IMAGE_TAG_SERVER/" "$BICEP_DIR/.env.deployment"
    fi
    echo -e "${YELLOW}� Updated .env.deployment with IMAGE_TAG_SERVER: $IMAGE_TAG_SERVER${NC}"
fi

# Restore original ACR network default action
echo -e "${YELLOW}🔧 Restoring ACR network security settings...${NC}"
echo -e "${YELLOW}📝 Restoring default action to: $ORIGINAL_DEFAULT_ACTION${NC}"
az acr update --name $ACR_NAME --default-action $ORIGINAL_DEFAULT_ACTION

# Security note: Keep IP firewall rules but can optionally disable public access
echo -e "${YELLOW}🔒 Managing ACR security configuration...${NC}"
if [ "$CURRENT_IP" != "unknown" ]; then
    echo -e "${GREEN}✅ Current IP ($CURRENT_IP) will remain in firewall allow list${NC}"
else
    echo -e "${YELLOW}⚠️  Current IP could not be determined, check firewall manually if needed${NC}"
fi

# Uncomment the next line if you want to disable public access after builds
# (Note: This will block access for other IPs, keeping only the firewall allow list)
# az acr update --name $ACR_NAME --public-network-enabled false

echo -e "${GREEN}✅ Images built and pushed successfully${NC}"

if [ "$BUILD_WEB" = true ]; then
    echo -e "${YELLOW}📝 Web image tag: $IMAGE_TAG_WEB${NC}"
fi

if [ "$BUILD_SERVER" = true ]; then
    echo -e "${YELLOW}📝 Server image tag: $IMAGE_TAG_SERVER${NC}"
fi

if [ "$BUILD_SERVER" = true ] && [ "$BUILD_WEB" = true ]; then
    echo -e "${YELLOW}🔐 Security: ACR public access enabled with IP firewall restrictions${NC}"
elif [ "$BUILD_SERVER" = true ]; then
    echo -e "${YELLOW}🔐 Built: Server image${NC}"
elif [ "$BUILD_WEB" = true ]; then
    echo -e "${YELLOW}🔐 Built: Web image${NC}"
fi

# Display current ACR network configuration
echo -e "${YELLOW}📋 Current ACR network configuration:${NC}"
az acr show --name $ACR_NAME --query "{publicNetworkAccess: publicNetworkAccess, defaultAction: networkRuleSet.defaultAction, allowedIPs: networkRuleSet.ipRules[].ipAddressOrRange}" -o table