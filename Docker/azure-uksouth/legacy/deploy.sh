#!/usr/bin/env bash
set -euo pipefail

#############################################
# CCP4i2 Azure Deployment Script
#
# This script builds Docker images, pushes them to Azure Container Registry,
# and deploys the infrastructure using Bicep templates.
#
# Prerequisites:
#   - Azure CLI (az) installed and logged in
#   - Docker installed and running
#   - azcopy installed (for uploading CCP4 data)
#
# Usage:
#   ./deploy.sh [command] [options]
#
# Commands:
#   build       Build Docker images locally
#   push        Push images to Azure Container Registry
#   deploy      Deploy infrastructure using Bicep
#   upload-ccp4 Upload CCP4 distribution to Azure Files
#   all         Run all steps (build, push, deploy)
#   status      Show deployment status
#
# Examples:
#   ./deploy.sh build
#   ./deploy.sh push --registry myregistry
#   ./deploy.sh deploy --resource-group ccp4i2-rg
#   ./deploy.sh all --resource-group ccp4i2-rg --registry myregistry
#############################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
ENV_FILE="${SCRIPT_DIR}/../.env"

# Default values
RESOURCE_GROUP="${AZURE_RESOURCE_GROUP:-ccp4i2-rg}"
LOCATION="${AZURE_LOCATION:-uksouth}"
BASE_NAME="${AZURE_BASE_NAME:-ccp4i2}"
REGISTRY_NAME="${AZURE_REGISTRY:-}"
SERVER_TAG="${SERVER_TAG:-latest}"
CLIENT_TAG="${CLIENT_TAG:-latest}"
CCP4_VERSION="${CCP4_VERSION:-ccp4-9}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Load .env file if it exists
load_env() {
    if [[ -f "$ENV_FILE" ]]; then
        echo -e "${BLUE}>>> Loading configuration from $ENV_FILE${NC}"
        set -a
        source <(grep -v '^#' "$ENV_FILE" | grep -v '^$' | grep '=')
        set +a
    fi
}

# Print usage
usage() {
    cat << EOF
Usage: $0 [command] [options]

Commands:
    build           Build Docker images locally
    push            Push images to Azure Container Registry
    deploy          Deploy infrastructure using Bicep
    upload-ccp4     Upload CCP4 distribution to Azure Files
    all             Run all steps (build, push, deploy)
    status          Show deployment status
    help            Show this help message

Options:
    --resource-group, -g    Azure resource group (default: ccp4i2-rg)
    --location, -l          Azure location (default: uksouth)
    --registry, -r          Azure Container Registry name
    --base-name, -n         Base name for resources (default: ccp4i2)
    --ccp4-version          CCP4 version directory name (default: ccp4-9)
    --server-tag            Server image tag (default: latest)
    --client-tag            Client image tag (default: latest)

Environment Variables (can be set in Docker/.env):
    AZURE_RESOURCE_GROUP    Resource group name
    AZURE_LOCATION          Azure region
    AZURE_REGISTRY          Container registry name
    AZURE_BASE_NAME         Base name for resources
    CCP4_VERSION            CCP4 version (e.g., ccp4-20251105)

Examples:
    $0 build
    $0 push --registry ccp4i2acr
    $0 deploy --resource-group ccp4i2-rg --location uksouth
    $0 all --resource-group ccp4i2-rg --registry ccp4i2acr
    $0 upload-ccp4 /path/to/ccp4-9
EOF
}

# Check prerequisites
check_prerequisites() {
    echo -e "${BLUE}>>> Checking prerequisites...${NC}"

    local missing=()

    command -v az >/dev/null 2>&1 || missing+=("az (Azure CLI)")
    command -v docker >/dev/null 2>&1 || missing+=("docker")

    if [[ ${#missing[@]} -gt 0 ]]; then
        echo -e "${RED}ERROR: Missing required tools:${NC}"
        for tool in "${missing[@]}"; do
            echo "  - $tool"
        done
        exit 1
    fi

    # Check if logged in to Azure
    if ! az account show >/dev/null 2>&1; then
        echo -e "${YELLOW}>>> Not logged in to Azure. Logging in...${NC}"
        az login
    fi

    echo -e "${GREEN}>>> Prerequisites OK${NC}"
}

# Build Docker images
build_images() {
    echo -e "${BLUE}>>> Building Docker images...${NC}"

    # Build server image
    echo -e "${BLUE}>>> Building server image...${NC}"
    docker build \
        --platform linux/amd64 \
        -t "ccp4i2-server:${SERVER_TAG}" \
        -f "$PROJECT_ROOT/Docker/server/Dockerfile" \
        "$PROJECT_ROOT"

    # Build client image (if Dockerfile exists)
    if [[ -f "$PROJECT_ROOT/Docker/client/Dockerfile" ]]; then
        echo -e "${BLUE}>>> Building client image...${NC}"
        docker build \
            --platform linux/amd64 \
            -t "ccp4i2-client:${CLIENT_TAG}" \
            -f "$PROJECT_ROOT/Docker/client/Dockerfile" \
            "$PROJECT_ROOT"
    else
        echo -e "${YELLOW}>>> No client Dockerfile found, skipping client build${NC}"
    fi

    echo -e "${GREEN}>>> Build complete${NC}"
}

# Create Azure Container Registry if it doesn't exist
ensure_registry() {
    if [[ -z "$REGISTRY_NAME" ]]; then
        # Generate registry name from base name (remove hyphens, must be lowercase alphanumeric)
        REGISTRY_NAME="${BASE_NAME//[^a-zA-Z0-9]/}acr"
    fi

    echo -e "${BLUE}>>> Ensuring container registry exists: $REGISTRY_NAME${NC}"

    # Create resource group if needed
    if ! az group show --name "$RESOURCE_GROUP" >/dev/null 2>&1; then
        echo -e "${BLUE}>>> Creating resource group: $RESOURCE_GROUP${NC}"
        az group create --name "$RESOURCE_GROUP" --location "$LOCATION"
    fi

    # Create registry if needed
    if ! az acr show --name "$REGISTRY_NAME" --resource-group "$RESOURCE_GROUP" >/dev/null 2>&1; then
        echo -e "${BLUE}>>> Creating container registry: $REGISTRY_NAME${NC}"
        az acr create \
            --resource-group "$RESOURCE_GROUP" \
            --name "$REGISTRY_NAME" \
            --sku Basic \
            --admin-enabled true
    fi

    echo -e "${GREEN}>>> Registry ready: $REGISTRY_NAME${NC}"
}

# Push images to Azure Container Registry
push_images() {
    ensure_registry

    echo -e "${BLUE}>>> Logging in to registry: $REGISTRY_NAME${NC}"
    az acr login --name "$REGISTRY_NAME"

    local registry_url="${REGISTRY_NAME}.azurecr.io"

    # Tag and push server image
    echo -e "${BLUE}>>> Pushing server image...${NC}"
    docker tag "ccp4i2-server:${SERVER_TAG}" "${registry_url}/ccp4i2-server:${SERVER_TAG}"
    docker push "${registry_url}/ccp4i2-server:${SERVER_TAG}"

    # Tag and push client image if it exists
    if docker image inspect "ccp4i2-client:${CLIENT_TAG}" >/dev/null 2>&1; then
        echo -e "${BLUE}>>> Pushing client image...${NC}"
        docker tag "ccp4i2-client:${CLIENT_TAG}" "${registry_url}/ccp4i2-client:${CLIENT_TAG}"
        docker push "${registry_url}/ccp4i2-client:${CLIENT_TAG}"
    fi

    echo -e "${GREEN}>>> Push complete${NC}"
    echo ""
    echo "Images available at:"
    echo "  ${registry_url}/ccp4i2-server:${SERVER_TAG}"
    if docker image inspect "ccp4i2-client:${CLIENT_TAG}" >/dev/null 2>&1; then
        echo "  ${registry_url}/ccp4i2-client:${CLIENT_TAG}"
    fi
}

# Deploy infrastructure
deploy_infrastructure() {
    ensure_registry

    echo -e "${BLUE}>>> Deploying infrastructure with Bicep...${NC}"

    local registry_url="${REGISTRY_NAME}.azurecr.io"

    # Generate secure passwords if not provided
    if [[ -z "${DB_PASSWORD:-}" ]]; then
        DB_PASSWORD=$(openssl rand -base64 24 | tr -dc 'a-zA-Z0-9' | head -c 24)
        echo -e "${YELLOW}>>> Generated database password (save this!): $DB_PASSWORD${NC}"
    fi

    if [[ -z "${SECRET_KEY:-}" ]]; then
        SECRET_KEY=$(openssl rand -base64 48 | tr -dc 'a-zA-Z0-9' | head -c 50)
        echo -e "${YELLOW}>>> Generated Django secret key${NC}"
    fi

    # Deploy using Bicep
    az deployment group create \
        --resource-group "$RESOURCE_GROUP" \
        --template-file "$SCRIPT_DIR/container-app.bicep" \
        --parameters \
            baseName="$BASE_NAME" \
            location="$LOCATION" \
            serverImage="${registry_url}/ccp4i2-server:${SERVER_TAG}" \
            clientImage="${registry_url}/ccp4i2-client:${CLIENT_TAG}" \
            dbPassword="$DB_PASSWORD" \
            secretKey="$SECRET_KEY"

    echo -e "${GREEN}>>> Deployment complete!${NC}"
    show_status
}

# Upload CCP4 distribution
upload_ccp4() {
    local ccp4_path="${1:-}"

    if [[ -z "$ccp4_path" ]]; then
        echo -e "${RED}ERROR: CCP4 path required${NC}"
        echo "Usage: $0 upload-ccp4 /path/to/ccp4-distribution"
        exit 1
    fi

    if [[ ! -d "$ccp4_path" ]]; then
        echo -e "${RED}ERROR: CCP4 path does not exist: $ccp4_path${NC}"
        exit 1
    fi

    # Use azcopy-files.sh if available
    if [[ -x "$SCRIPT_DIR/azcopy-files.sh" ]]; then
        "$SCRIPT_DIR/azcopy-files.sh" "$ccp4_path" "${CCP4_VERSION}"
    else
        echo -e "${RED}ERROR: azcopy-files.sh not found or not executable${NC}"
        exit 1
    fi
}

# Show deployment status
show_status() {
    echo ""
    echo -e "${BLUE}=== Deployment Status ===${NC}"

    # Get container app URLs
    if az containerapp show --name "${BASE_NAME}-server" --resource-group "$RESOURCE_GROUP" >/dev/null 2>&1; then
        local server_url=$(az containerapp show \
            --name "${BASE_NAME}-server" \
            --resource-group "$RESOURCE_GROUP" \
            --query "properties.configuration.ingress.fqdn" \
            --output tsv)
        echo -e "${GREEN}Server URL:${NC} https://$server_url"
    fi

    if az containerapp show --name "${BASE_NAME}-client" --resource-group "$RESOURCE_GROUP" >/dev/null 2>&1; then
        local client_url=$(az containerapp show \
            --name "${BASE_NAME}-client" \
            --resource-group "$RESOURCE_GROUP" \
            --query "properties.configuration.ingress.fqdn" \
            --output tsv)
        echo -e "${GREEN}Client URL:${NC} https://$client_url"
    fi

    # Get storage account
    if az storage account show --name "${BASE_NAME}storage" --resource-group "$RESOURCE_GROUP" >/dev/null 2>&1; then
        echo -e "${GREEN}Storage Account:${NC} ${BASE_NAME}storage"
    fi

    # Get database server
    if az postgres flexible-server show --name "${BASE_NAME}-db" --resource-group "$RESOURCE_GROUP" >/dev/null 2>&1; then
        local db_fqdn=$(az postgres flexible-server show \
            --name "${BASE_NAME}-db" \
            --resource-group "$RESOURCE_GROUP" \
            --query "fullyQualifiedDomainName" \
            --output tsv)
        echo -e "${GREEN}Database:${NC} $db_fqdn"
    fi

    echo ""
    echo -e "${BLUE}>>> Next Steps:${NC}"
    echo "1. Upload CCP4 distribution: $0 upload-ccp4 /path/to/ccp4-distribution"
    echo "2. Access the application at the Client URL above"
    echo "3. Monitor logs: az containerapp logs show --name ${BASE_NAME}-server --resource-group $RESOURCE_GROUP"
}

# Parse arguments
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --resource-group|-g)
                RESOURCE_GROUP="$2"
                shift 2
                ;;
            --location|-l)
                LOCATION="$2"
                shift 2
                ;;
            --registry|-r)
                REGISTRY_NAME="$2"
                shift 2
                ;;
            --base-name|-n)
                BASE_NAME="$2"
                shift 2
                ;;
            --ccp4-version)
                CCP4_VERSION="$2"
                shift 2
                ;;
            --server-tag)
                SERVER_TAG="$2"
                shift 2
                ;;
            --client-tag)
                CLIENT_TAG="$2"
                shift 2
                ;;
            *)
                break
                ;;
        esac
    done
    echo "$@"
}

# Main
main() {
    load_env

    local remaining_args=$(parse_args "$@")
    set -- $remaining_args

    local command="${1:-help}"
    shift || true

    case $command in
        build)
            check_prerequisites
            build_images
            ;;
        push)
            check_prerequisites
            push_images
            ;;
        deploy)
            check_prerequisites
            deploy_infrastructure
            ;;
        upload-ccp4)
            check_prerequisites
            upload_ccp4 "$@"
            ;;
        all)
            check_prerequisites
            build_images
            push_images
            deploy_infrastructure
            ;;
        status)
            check_prerequisites
            show_status
            ;;
        help|--help|-h)
            usage
            ;;
        *)
            echo -e "${RED}Unknown command: $command${NC}"
            usage
            exit 1
            ;;
    esac
}

main "$@"
