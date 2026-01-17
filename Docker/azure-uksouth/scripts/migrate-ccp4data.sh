#!/bin/bash

# Migrate CCP4 Data to Public Storage Account
#
# After deploying the updated infrastructure with separate storage accounts:
# - Private: stornekmayz3n2 (projects, media, static files)
# - Public: storpub<unique> (ccp4data for CCP4 software)
#
# This script copies ccp4data from the old private storage to the new public storage.
#
# Usage: ./migrate-ccp4data.sh [check|copy]

# Ensure Homebrew paths are available
export PATH="/opt/homebrew/bin:$PATH"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Get the directory where the script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="$SCRIPT_DIR/../.env.deployment"

# Load environment variables
if [ -f "$ENV_FILE" ]; then
    source "$ENV_FILE"
else
    echo -e "${RED}Error: .env.deployment not found${NC}"
    exit 1
fi

SHARE_NAME="ccp4data"

# Check if azcopy is installed
check_azcopy() {
    if ! command -v azcopy &> /dev/null; then
        echo -e "${RED}Error: azcopy is not installed${NC}"
        echo "Install with: brew install azcopy"
        exit 1
    fi
}

# Get storage account info by querying Azure directly
get_storage_accounts() {
    echo -e "${BLUE}Getting storage account info...${NC}"

    # Query storage accounts directly - more reliable than deployment outputs
    # Private storage: starts with 'stor' but not 'storpub'
    PRIVATE_STORAGE=$(az storage account list \
        --resource-group $RESOURCE_GROUP \
        --query "[?starts_with(name, 'stor') && !starts_with(name, 'storpub')].name" \
        -o tsv 2>/dev/null | head -1)

    # Public storage: starts with 'storpub'
    PUBLIC_STORAGE=$(az storage account list \
        --resource-group $RESOURCE_GROUP \
        --query "[?starts_with(name, 'storpub')].name" \
        -o tsv 2>/dev/null | head -1)

    if [ -z "$PRIVATE_STORAGE" ]; then
        echo -e "${RED}Error: Private storage account not found${NC}"
        exit 1
    fi

    if [ -z "$PUBLIC_STORAGE" ]; then
        echo -e "${RED}Error: Public storage account not found${NC}"
        echo -e "${YELLOW}Have you deployed the updated infrastructure.bicep?${NC}"
        echo ""
        echo "Run: ./deploy-infrastructure.sh"
        exit 1
    fi

    echo -e "${GREEN}Private storage: $PRIVATE_STORAGE${NC}"
    echo -e "${GREEN}Public storage:  $PUBLIC_STORAGE${NC}"
}

# Enable network access
enable_network_access() {
    local account=$1
    echo -e "${YELLOW}Enabling network access for $account...${NC}"

    CURRENT_IP=$(curl -s https://ipinfo.io/ip 2>/dev/null || curl -s https://api.ipify.org 2>/dev/null)

    az storage account network-rule add \
        --account-name "$account" \
        --resource-group "$RESOURCE_GROUP" \
        --ip-address "$CURRENT_IP" \
        --output none 2>/dev/null || true

    az storage account update \
        --name "$account" \
        --resource-group "$RESOURCE_GROUP" \
        --default-action Allow \
        --output none 2>/dev/null || true
}

# Restore network access
restore_network_access() {
    local account=$1
    local default_action=$2

    if [ -n "$account" ] && [ -n "$default_action" ]; then
        az storage account update \
            --name "$account" \
            --resource-group "$RESOURCE_GROUP" \
            --default-action "$default_action" \
            --output none 2>/dev/null || true
    fi
}

# Generate SAS token
generate_sas() {
    local account=$1
    local permissions=$2

    local key=$(az storage account keys list \
        --account-name "$account" \
        --resource-group "$RESOURCE_GROUP" \
        --query "[0].value" -o tsv)

    local expiry=$(date -u -v+24H "+%Y-%m-%dT%H:%MZ" 2>/dev/null || date -u -d "+24 hours" "+%Y-%m-%dT%H:%MZ")

    az storage account generate-sas \
        --account-name "$account" \
        --account-key "$key" \
        --permissions "$permissions" \
        --resource-types sco \
        --services f \
        --expiry "$expiry" \
        -o tsv
}

# Check what needs to be migrated
check_migration() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  CCP4 Data Migration Check${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    check_azcopy
    get_storage_accounts

    # Enable access to both accounts
    enable_network_access "$PRIVATE_STORAGE"
    sleep 2

    # Get keys
    PRIVATE_KEY=$(az storage account keys list \
        --account-name "$PRIVATE_STORAGE" \
        --resource-group "$RESOURCE_GROUP" \
        --query "[0].value" -o tsv)

    PUBLIC_KEY=$(az storage account keys list \
        --account-name "$PUBLIC_STORAGE" \
        --resource-group "$RESOURCE_GROUP" \
        --query "[0].value" -o tsv 2>/dev/null)

    # Check source (private storage)
    echo -e "${YELLOW}Checking source (private storage: $PRIVATE_STORAGE)...${NC}"
    SOURCE_EXISTS=$(az storage share exists \
        --name $SHARE_NAME \
        --account-name "$PRIVATE_STORAGE" \
        --account-key "$PRIVATE_KEY" \
        --query "exists" -o tsv 2>/dev/null || echo "false")

    if [ "$SOURCE_EXISTS" = "true" ]; then
        SOURCE_ITEMS=$(az storage file list \
            --share-name $SHARE_NAME \
            --account-name "$PRIVATE_STORAGE" \
            --account-key "$PRIVATE_KEY" \
            --query "length([])" -o tsv 2>/dev/null || echo "0")
        echo -e "${GREEN}Source share exists with $SOURCE_ITEMS root items${NC}"
    else
        echo -e "${YELLOW}Source share does not exist in private storage${NC}"
        SOURCE_ITEMS=0
    fi

    # Check destination (public storage)
    echo -e "${YELLOW}Checking destination (public storage: $PUBLIC_STORAGE)...${NC}"
    if [ -n "$PUBLIC_KEY" ]; then
        DEST_EXISTS=$(az storage share exists \
            --name $SHARE_NAME \
            --account-name "$PUBLIC_STORAGE" \
            --account-key "$PUBLIC_KEY" \
            --query "exists" -o tsv 2>/dev/null || echo "false")

        if [ "$DEST_EXISTS" = "true" ]; then
            DEST_ITEMS=$(az storage file list \
                --share-name $SHARE_NAME \
                --account-name "$PUBLIC_STORAGE" \
                --account-key "$PUBLIC_KEY" \
                --query "length([])" -o tsv 2>/dev/null || echo "0")
            echo -e "${GREEN}Destination share exists with $DEST_ITEMS root items${NC}"
        else
            echo -e "${YELLOW}Destination share does not exist yet${NC}"
            echo -e "${YELLOW}It will be created when you deploy infrastructure${NC}"
            DEST_ITEMS=0
        fi
    else
        echo -e "${YELLOW}Public storage account not yet deployed${NC}"
        DEST_EXISTS="false"
        DEST_ITEMS=0
    fi

    # Restore private storage to Deny
    restore_network_access "$PRIVATE_STORAGE" "Deny"

    echo ""
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Migration Summary${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""
    echo "Source:      $PRIVATE_STORAGE/$SHARE_NAME ($SOURCE_ITEMS items)"
    echo "Destination: $PUBLIC_STORAGE/$SHARE_NAME ($DEST_ITEMS items)"
    echo ""

    if [ "$SOURCE_EXISTS" = "true" ] && [ "$SOURCE_ITEMS" -gt 0 ]; then
        echo -e "${YELLOW}Migration needed: Copy ccp4data to public storage${NC}"
        echo ""
        echo "Run: $0 copy"
    elif [ "$SOURCE_EXISTS" = "false" ] || [ "$SOURCE_ITEMS" = "0" ]; then
        echo -e "${GREEN}No migration needed: Source is empty${NC}"
        echo -e "${YELLOW}Upload CCP4 data directly to: $PUBLIC_STORAGE/$SHARE_NAME${NC}"
    fi
}

# Perform the copy
do_copy() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Migrating CCP4 Data${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    check_azcopy
    get_storage_accounts

    # Enable access
    enable_network_access "$PRIVATE_STORAGE"
    enable_network_access "$PUBLIC_STORAGE"
    sleep 3

    # Generate SAS tokens
    echo -e "${YELLOW}Generating SAS tokens...${NC}"
    SOURCE_SAS=$(generate_sas "$PRIVATE_STORAGE" "rl")
    DEST_SAS=$(generate_sas "$PUBLIC_STORAGE" "rwdl")

    if [ -z "$SOURCE_SAS" ] || [ -z "$DEST_SAS" ]; then
        echo -e "${RED}Failed to generate SAS tokens${NC}"
        exit 1
    fi

    # Build URLs
    SOURCE_URL="https://${PRIVATE_STORAGE}.file.core.windows.net/${SHARE_NAME}/*?${SOURCE_SAS}"
    DEST_URL="https://${PUBLIC_STORAGE}.file.core.windows.net/${SHARE_NAME}/?${DEST_SAS}"

    echo ""
    echo -e "${YELLOW}Copying ccp4data...${NC}"
    echo -e "${BLUE}From: $PRIVATE_STORAGE/$SHARE_NAME${NC}"
    echo -e "${BLUE}To:   $PUBLIC_STORAGE/$SHARE_NAME${NC}"
    echo ""

    azcopy copy \
        "$SOURCE_URL" \
        "$DEST_URL" \
        --recursive \
        --preserve-smb-info=false \
        --skip-version-check \
        --log-level=WARNING

    if [ $? -eq 0 ]; then
        echo ""
        echo -e "${GREEN}========================================${NC}"
        echo -e "${GREEN}  Migration Complete!${NC}"
        echo -e "${GREEN}========================================${NC}"
        echo ""
        echo -e "${YELLOW}Next steps:${NC}"
        echo "1. Deploy applications to use the new storage mount"
        echo "2. Optionally delete ccp4data from private storage:"
        echo "   az storage share delete --name $SHARE_NAME --account-name $PRIVATE_STORAGE"
    else
        echo -e "${RED}Migration failed${NC}"
    fi

    # Restore private storage to Deny
    restore_network_access "$PRIVATE_STORAGE" "Deny"
}

# Show usage
show_usage() {
    echo "CCP4 Data Migration Script"
    echo ""
    echo "Migrates ccp4data from private storage to public storage account."
    echo ""
    echo "Usage: $0 <command>"
    echo ""
    echo "Commands:"
    echo "  check  Check migration status"
    echo "  copy   Perform the migration"
    echo ""
    echo "Prerequisites:"
    echo "  - Deploy updated infrastructure.bicep first"
    echo "  - This creates the new public storage account"
    echo ""
}

# Main
case "${1:-}" in
    check)
        check_migration
        ;;
    copy)
        do_copy
        ;;
    help|--help|-h)
        show_usage
        ;;
    *)
        if [ -n "$1" ]; then
            echo -e "${RED}Unknown command: $1${NC}"
            echo ""
        fi
        show_usage
        exit 1
        ;;
esac
