#!/bin/bash

# Docker Build and Push Script
# Uses local Docker with docker-compose for parallel builds
# Cross-platform optimized: builds natively on Apple Silicon, targets linux/amd64
#
# Usage: ./build-and-push.sh [server|web]
#   - No argument: Build both server and web images in parallel (default)
#   - server: Build only server image
#   - web: Build only web image
#
# See also: build-and-push-acr.sh for ACR-based builds (alternative)

# Ensure Homebrew paths are available
export PATH="/opt/homebrew/bin:$PATH"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

set -e

# Change to the correct working directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BICEP_DIR="$(dirname "$SCRIPT_DIR")"
WORKING_DIR="$(cd "$BICEP_DIR/../.." && pwd)"

# Parse command line arguments
BUILD_SERVER=true
BUILD_WEB=true
BUILD_TARGET=""

if [ $# -gt 0 ]; then
    case "$1" in
        "server")
            BUILD_WEB=false
            BUILD_TARGET="server"
            echo -e "${YELLOW}🔧 Building server image only${NC}"
            ;;
        "web")
            BUILD_SERVER=false
            BUILD_TARGET="web"
            echo -e "${YELLOW}🔧 Building web image only${NC}"
            ;;
        *)
            echo -e "${RED}❌ Invalid argument. Use 'server', 'web', or no argument for both${NC}"
            exit 1
            ;;
    esac
else
    echo -e "${YELLOW}🔧 Building both server and web images (parallel)${NC}"
fi

echo -e "${GREEN}🐳 Building and pushing Docker images (local build)${NC}"
echo -e "${YELLOW}📁 Working directory: $WORKING_DIR${NC}"

cd "$WORKING_DIR"

# Load environment variables
if [ -f "$BICEP_DIR/.env.deployment" ]; then
    source "$BICEP_DIR/.env.deployment"
else
    echo -e "${RED}❌ .env.deployment not found${NC}"
    exit 1
fi

# Security check
echo -e "${YELLOW}🔒 Checking authentication configuration...${NC}"
if [ "$CCP4I2_REQUIRE_AUTH" != "true" ]; then
    echo -e "${RED}╔═══════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${RED}║                   ⚠️  SECURITY WARNING ⚠️                          ║${NC}"
    echo -e "${RED}╠═══════════════════════════════════════════════════════════════════╣${NC}"
    echo -e "${RED}║  CCP4I2_REQUIRE_AUTH is not set to 'true' in .env.deployment     ║${NC}"
    echo -e "${RED}║  Current value: CCP4I2_REQUIRE_AUTH=\"$CCP4I2_REQUIRE_AUTH\"${NC}"
    echo -e "${RED}╚═══════════════════════════════════════════════════════════════════╝${NC}"

    if [ "$ALLOW_NO_AUTH_BUILD" = "true" ]; then
        echo -e "${YELLOW}⚠️  ALLOW_NO_AUTH_BUILD=true - proceeding with NO AUTH${NC}"
        sleep 5
    else
        exit 1
    fi
else
    echo -e "${GREEN}✅ Authentication required: CCP4I2_REQUIRE_AUTH=true${NC}"
fi

# Login to ACR (needed for pulling base image and pushing)
echo -e "${YELLOW}🔑 Logging into Azure Container Registry...${NC}"
az acr login --name $ACR_NAME
if [ $? -ne 0 ]; then
    echo -e "${RED}❌ ACR login failed. Try 'az login' first.${NC}"
    exit 1
fi
echo -e "${GREEN}✅ ACR login successful${NC}"

# Generate image tag
TIMESTAMP=$(date +%Y%m%d-%H%M%S)

# Export variables for docker compose
export TIMESTAMP
export ACR_LOGIN_SERVER
export ACR_NAME
export CCP4_VERSION="${CCP4_VERSION:-ccp4-20251105}"
export BASE_IMAGE_NAME="${BASE_IMAGE_NAME:-ccp4i2/base-arpwarp}"
export NEXT_PUBLIC_AAD_CLIENT_ID="${NEXT_PUBLIC_AAD_CLIENT_ID:-}"
export NEXT_PUBLIC_AAD_TENANT_ID="${NEXT_PUBLIC_AAD_TENANT_ID:-}"
export NEXT_PUBLIC_REQUIRE_AUTH="${NEXT_PUBLIC_REQUIRE_AUTH:-false}"

echo -e "${YELLOW}📦 Build configuration:${NC}"
echo -e "   ACR: $ACR_LOGIN_SERVER"
echo -e "   Tag: $TIMESTAMP"
echo -e "   Platform: linux/amd64"
if [ "$BUILD_SERVER" = true ]; then
    echo -e "   Server base: $ACR_LOGIN_SERVER/$BASE_IMAGE_NAME:$CCP4_VERSION"
fi

# Check if docker compose file exists
COMPOSE_FILE="$BICEP_DIR/docker-compose.yml"
if [ ! -f "$COMPOSE_FILE" ]; then
    echo -e "${RED}❌ docker-compose.yml not found at $COMPOSE_FILE${NC}"
    exit 1
fi

# Build images using docker compose (parallel by default)
echo -e "${YELLOW}🔨 Building images...${NC}"

if [ -z "$BUILD_TARGET" ]; then
    # Build both in parallel
    docker compose -f "$COMPOSE_FILE" build --parallel
else
    # Build specific target
    docker compose -f "$COMPOSE_FILE" build "$BUILD_TARGET"
fi

if [ $? -ne 0 ]; then
    echo -e "${RED}❌ Build failed${NC}"
    exit 1
fi

echo -e "${GREEN}✅ Build completed${NC}"

# Tag and push images
echo -e "${YELLOW}📤 Tagging and pushing images to ACR...${NC}"

if [ "$BUILD_WEB" = true ]; then
    echo -e "${YELLOW}   Tagging web image as latest...${NC}"
    docker tag $ACR_LOGIN_SERVER/ccp4i2/web:$TIMESTAMP $ACR_LOGIN_SERVER/ccp4i2/web:latest
    echo -e "${YELLOW}   Pushing web image...${NC}"
    docker push $ACR_LOGIN_SERVER/ccp4i2/web:$TIMESTAMP
    docker push $ACR_LOGIN_SERVER/ccp4i2/web:latest
fi

if [ "$BUILD_SERVER" = true ]; then
    echo -e "${YELLOW}   Tagging server image as latest...${NC}"
    docker tag $ACR_LOGIN_SERVER/ccp4i2/server:$TIMESTAMP $ACR_LOGIN_SERVER/ccp4i2/server:latest
    echo -e "${YELLOW}   Pushing server image...${NC}"
    docker push $ACR_LOGIN_SERVER/ccp4i2/server:$TIMESTAMP
    docker push $ACR_LOGIN_SERVER/ccp4i2/server:latest
fi

echo -e "${GREEN}✅ Push completed${NC}"

# Update .env.deployment with new tags
if [ "$BUILD_WEB" = true ]; then
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' "s/^IMAGE_TAG_WEB=.*/IMAGE_TAG_WEB=$TIMESTAMP/" "$BICEP_DIR/.env.deployment"
    else
        sed -i "s/^IMAGE_TAG_WEB=.*/IMAGE_TAG_WEB=$TIMESTAMP/" "$BICEP_DIR/.env.deployment"
    fi
    echo -e "${YELLOW}📝 Updated IMAGE_TAG_WEB: $TIMESTAMP${NC}"
fi

if [ "$BUILD_SERVER" = true ]; then
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' "s/^IMAGE_TAG_SERVER=.*/IMAGE_TAG_SERVER=$TIMESTAMP/" "$BICEP_DIR/.env.deployment"
    else
        sed -i "s/^IMAGE_TAG_SERVER=.*/IMAGE_TAG_SERVER=$TIMESTAMP/" "$BICEP_DIR/.env.deployment"
    fi
    echo -e "${YELLOW}📝 Updated IMAGE_TAG_SERVER: $TIMESTAMP${NC}"
fi

echo -e "${GREEN}✅ Images built and pushed successfully${NC}"
echo -e "${YELLOW}📝 Image tag: $TIMESTAMP${NC}"
echo ""
echo -e "${YELLOW}Next steps:${NC}"
echo -e "  Deploy with: ./scripts/deploy-applications.sh web"
echo -e "            or ./scripts/deploy-applications.sh server"
