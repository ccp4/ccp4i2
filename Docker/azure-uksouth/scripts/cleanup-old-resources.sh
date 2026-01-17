#!/bin/bash

# Cleanup Old/Orphaned Azure Resources
#
# After migrating to the new storage architecture, this script removes:
# - Old storage mounts that pointed to stornekmayz3n2
# - Old private endpoints that can't be updated
#
# Usage: ./cleanup-old-resources.sh [check|cleanup]

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

# Container Apps Environment name
CAE_NAME="ccp4i2-bicep-env-ne"

# Old resources to clean up
OLD_MOUNTS=("ccp4data-mount" "staticfiles-mount" "mediafiles-mount" "projects-mount" "ccp4data-public" "staticfiles-private" "mediafiles-private" "projects-private")
OLD_PRIVATE_ENDPOINTS=("ccp4i2-bicep-ne-storage-pe" "ccp4i2-bicep-ne-blob-pe")

# Check what resources exist
check_resources() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Checking for Old Resources${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    # Check old storage mounts
    echo -e "${YELLOW}Checking old storage mounts...${NC}"
    EXISTING_MOUNTS=$(az containerapp env storage list \
        --name "$CAE_NAME" \
        --resource-group "$RESOURCE_GROUP" \
        --query "[].name" -o tsv 2>/dev/null)

    MOUNTS_TO_DELETE=()
    for mount in "${OLD_MOUNTS[@]}"; do
        if echo "$EXISTING_MOUNTS" | grep -q "^${mount}$"; then
            echo -e "${RED}  Found: $mount${NC}"
            MOUNTS_TO_DELETE+=("$mount")
        else
            echo -e "${GREEN}  Not found: $mount (already cleaned)${NC}"
        fi
    done

    echo ""
    echo -e "${YELLOW}Checking old private endpoints...${NC}"
    EXISTING_PES=$(az network private-endpoint list \
        --resource-group "$RESOURCE_GROUP" \
        --query "[].name" -o tsv 2>/dev/null)

    PES_TO_DELETE=()
    for pe in "${OLD_PRIVATE_ENDPOINTS[@]}"; do
        if echo "$EXISTING_PES" | grep -q "^${pe}$"; then
            echo -e "${RED}  Found: $pe${NC}"
            PES_TO_DELETE+=("$pe")
        else
            echo -e "${GREEN}  Not found: $pe (already cleaned)${NC}"
        fi
    done

    echo ""
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Summary${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    TOTAL_TO_DELETE=$((${#MOUNTS_TO_DELETE[@]} + ${#PES_TO_DELETE[@]}))

    if [ $TOTAL_TO_DELETE -eq 0 ]; then
        echo -e "${GREEN}No old resources to clean up!${NC}"
    else
        echo -e "${YELLOW}Resources to delete:${NC}"
        echo "  Storage mounts: ${#MOUNTS_TO_DELETE[@]}"
        echo "  Private endpoints: ${#PES_TO_DELETE[@]}"
        echo ""
        echo -e "${YELLOW}Run '$0 cleanup' to delete these resources${NC}"
    fi
}

# Delete old resources
do_cleanup() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Cleaning Up Old Resources${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    # Get existing mounts
    EXISTING_MOUNTS=$(az containerapp env storage list \
        --name "$CAE_NAME" \
        --resource-group "$RESOURCE_GROUP" \
        --query "[].name" -o tsv 2>/dev/null)

    # Delete old storage mounts
    echo -e "${YELLOW}Deleting old storage mounts...${NC}"
    for mount in "${OLD_MOUNTS[@]}"; do
        if echo "$EXISTING_MOUNTS" | grep -q "^${mount}$"; then
            echo -e "${BLUE}  Deleting: $mount${NC}"
            az containerapp env storage remove \
                --name "$CAE_NAME" \
                --resource-group "$RESOURCE_GROUP" \
                --storage-name "$mount" \
                --yes \
                --output none 2>/dev/null

            if [ $? -eq 0 ]; then
                echo -e "${GREEN}  Deleted: $mount${NC}"
            else
                echo -e "${RED}  Failed to delete: $mount${NC}"
            fi
        else
            echo -e "${GREEN}  Skipped: $mount (not found)${NC}"
        fi
    done

    # Get existing private endpoints
    EXISTING_PES=$(az network private-endpoint list \
        --resource-group "$RESOURCE_GROUP" \
        --query "[].name" -o tsv 2>/dev/null)

    # Delete old private endpoints
    echo ""
    echo -e "${YELLOW}Deleting old private endpoints...${NC}"
    for pe in "${OLD_PRIVATE_ENDPOINTS[@]}"; do
        if echo "$EXISTING_PES" | grep -q "^${pe}$"; then
            echo -e "${BLUE}  Deleting: $pe${NC}"
            az network private-endpoint delete \
                --name "$pe" \
                --resource-group "$RESOURCE_GROUP" \
                --output none 2>/dev/null

            if [ $? -eq 0 ]; then
                echo -e "${GREEN}  Deleted: $pe${NC}"
            else
                echo -e "${RED}  Failed to delete: $pe${NC}"
            fi
        else
            echo -e "${GREEN}  Skipped: $pe (not found)${NC}"
        fi
    done

    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}  Cleanup Complete!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
    echo -e "${YELLOW}Note: The new resources are:${NC}"
    echo "  Storage mounts: ccp4-software, staticfiles-private, mediafiles-private, projects-private"
    echo "  Private endpoints: *-storage-private-pe, *-blob-private-pe"
    echo ""
}

# Show usage
show_usage() {
    echo "Cleanup Old/Orphaned Azure Resources"
    echo ""
    echo "After migrating to the new storage architecture, this script removes:"
    echo "  - Old storage mounts (ccp4data-mount, staticfiles-mount, etc.)"
    echo "  - Old private endpoints that can't be updated"
    echo ""
    echo "Usage: $0 <command>"
    echo ""
    echo "Commands:"
    echo "  check    Check what old resources exist"
    echo "  cleanup  Delete old resources"
    echo ""
}

# Main command handler
case "${1:-}" in
    check)
        check_resources
        ;;
    cleanup)
        do_cleanup
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
