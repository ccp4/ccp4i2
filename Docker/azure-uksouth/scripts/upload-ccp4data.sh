#!/bin/bash

# Upload CCP4 Data to Azure Files (Public Storage)
#
# This script uploads a local CCP4 installation to the public Azure Files share.
# It uses tar to preserve symlinks, uploads the archive, then extracts in place.
#
# Usage: ./upload-ccp4data.sh <local-ccp4-path> [folder-name]
#
# Examples:
#   ./upload-ccp4data.sh /path/to/ccp4-20251105
#   ./upload-ccp4data.sh /path/to/ccp4-20251105 ccp4-20251105
#   ./upload-ccp4data.sh /path/to/arp-warp arp-warp

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

# Validate arguments
LOCAL_PATH="${1:-}"
FOLDER_NAME="${2:-}"

if [ -z "$LOCAL_PATH" ]; then
    echo "Upload CCP4 Data to Azure Files"
    echo ""
    echo "Usage: $0 <local-ccp4-path> [folder-name]"
    echo ""
    echo "Arguments:"
    echo "  local-ccp4-path  Path to local CCP4 installation directory"
    echo "  folder-name      Name for the folder in Azure (defaults to source folder name)"
    echo ""
    echo "Examples:"
    echo "  $0 ../ccp4-20251105"
    echo "  $0 /opt/ccp4-20251105 ccp4-20251105"
    echo "  $0 /path/to/arp-warp arp-warp"
    echo ""
    exit 1
fi

if [ ! -d "$LOCAL_PATH" ]; then
    echo -e "${RED}Error: Directory not found: $LOCAL_PATH${NC}"
    exit 1
fi

# Get absolute path and folder name
LOCAL_PATH=$(cd "$LOCAL_PATH" && pwd)
if [ -z "$FOLDER_NAME" ]; then
    FOLDER_NAME=$(basename "$LOCAL_PATH")
fi

echo -e "${CYAN}========================================${NC}"
echo -e "${CYAN}  Upload CCP4 Data to Azure Files${NC}"
echo -e "${CYAN}========================================${NC}"
echo ""
echo -e "${BLUE}Source: $LOCAL_PATH${NC}"
echo -e "${BLUE}Target: $FOLDER_NAME${NC}"
echo ""

# Get public storage account
get_storage_account() {
    echo -e "${YELLOW}Finding public storage account...${NC}"

    PUBLIC_STORAGE=$(az storage account list \
        --resource-group $RESOURCE_GROUP \
        --query "[?starts_with(name, 'storpub')].name" \
        -o tsv 2>/dev/null | head -1)

    if [ -z "$PUBLIC_STORAGE" ]; then
        echo -e "${RED}Error: Public storage account not found${NC}"
        echo -e "${YELLOW}Have you deployed the updated infrastructure.bicep?${NC}"
        exit 1
    fi

    echo -e "${GREEN}Public storage: $PUBLIC_STORAGE${NC}"
}

# Get storage key
get_storage_key() {
    STORAGE_KEY=$(az storage account keys list \
        --account-name "$PUBLIC_STORAGE" \
        --resource-group "$RESOURCE_GROUP" \
        --query "[0].value" -o tsv)

    if [ -z "$STORAGE_KEY" ]; then
        echo -e "${RED}Error: Could not get storage key${NC}"
        exit 1
    fi
}

# Check symlinks in source
check_symlinks() {
    echo -e "${YELLOW}Checking for symlinks...${NC}"

    SYMLINK_COUNT=$(find "$LOCAL_PATH" -type l 2>/dev/null | wc -l | tr -d ' ')

    if [ "$SYMLINK_COUNT" -gt 0 ]; then
        echo -e "${GREEN}Found $SYMLINK_COUNT symlinks - these will be preserved in the archive${NC}"
    else
        echo -e "${GREEN}No symlinks found${NC}"
    fi
}

# Create tar archive
create_archive() {
    ARCHIVE_NAME="${FOLDER_NAME}.tar.gz"
    ARCHIVE_PATH="/tmp/$ARCHIVE_NAME"

    echo -e "${YELLOW}Creating archive: $ARCHIVE_NAME${NC}"
    echo -e "${BLUE}This preserves symlinks...${NC}"

    # Create tar from parent directory to get correct folder structure
    PARENT_DIR=$(dirname "$LOCAL_PATH")
    SOURCE_DIR=$(basename "$LOCAL_PATH")

    # If folder name differs from source, we need to rename during archive
    if [ "$SOURCE_DIR" = "$FOLDER_NAME" ]; then
        tar -czhf "$ARCHIVE_PATH" -C "$PARENT_DIR" "$SOURCE_DIR"
    else
        # Use transform to rename the folder
        tar -czhf "$ARCHIVE_PATH" -C "$PARENT_DIR" --transform "s|^$SOURCE_DIR|$FOLDER_NAME|" "$SOURCE_DIR"
    fi

    if [ $? -ne 0 ]; then
        echo -e "${RED}Error creating archive${NC}"
        exit 1
    fi

    ARCHIVE_SIZE=$(ls -lh "$ARCHIVE_PATH" | awk '{print $5}')
    echo -e "${GREEN}Archive created: $ARCHIVE_SIZE${NC}"
}

# Upload archive to Azure Files
upload_archive() {
    echo -e "${YELLOW}Uploading archive to Azure Files...${NC}"
    echo -e "${BLUE}This may take a while for large archives...${NC}"

    az storage file upload \
        --share-name "$SHARE_NAME" \
        --source "$ARCHIVE_PATH" \
        --path "$ARCHIVE_NAME" \
        --account-name "$PUBLIC_STORAGE" \
        --account-key "$STORAGE_KEY" \
        --output none

    if [ $? -ne 0 ]; then
        echo -e "${RED}Error uploading archive${NC}"
        exit 1
    fi

    echo -e "${GREEN}Archive uploaded successfully${NC}"
}

# Extract archive on Azure (using Azure Container Instance)
extract_archive() {
    echo -e "${YELLOW}Extracting archive on Azure...${NC}"

    # Create a temporary container to extract the archive
    ACI_NAME="ccp4-extract-$(date +%s)"

    echo -e "${BLUE}Creating extraction container: $ACI_NAME${NC}"

    # Build the extraction command
    EXTRACT_CMD="cd /mnt/ccp4data && tar -xzf $ARCHIVE_NAME && rm $ARCHIVE_NAME && echo 'Extraction complete'"

    az container create \
        --resource-group "$RESOURCE_GROUP" \
        --name "$ACI_NAME" \
        --image mcr.microsoft.com/azure-cli:latest \
        --restart-policy Never \
        --azure-file-volume-share-name "$SHARE_NAME" \
        --azure-file-volume-account-name "$PUBLIC_STORAGE" \
        --azure-file-volume-account-key "$STORAGE_KEY" \
        --azure-file-volume-mount-path /mnt/ccp4data \
        --command-line "/bin/bash -c '$EXTRACT_CMD'" \
        --output none

    if [ $? -ne 0 ]; then
        echo -e "${RED}Error creating extraction container${NC}"
        echo -e "${YELLOW}You may need to extract manually${NC}"
        return 1
    fi

    echo -e "${YELLOW}Waiting for extraction to complete...${NC}"

    # Wait for container to complete
    for i in {1..60}; do
        STATE=$(az container show \
            --resource-group "$RESOURCE_GROUP" \
            --name "$ACI_NAME" \
            --query "containers[0].instanceView.currentState.state" \
            -o tsv 2>/dev/null)

        if [ "$STATE" = "Terminated" ]; then
            EXIT_CODE=$(az container show \
                --resource-group "$RESOURCE_GROUP" \
                --name "$ACI_NAME" \
                --query "containers[0].instanceView.currentState.exitCode" \
                -o tsv 2>/dev/null)

            if [ "$EXIT_CODE" = "0" ]; then
                echo -e "${GREEN}Extraction completed successfully${NC}"
            else
                echo -e "${RED}Extraction failed with exit code: $EXIT_CODE${NC}"
                echo -e "${YELLOW}Container logs:${NC}"
                az container logs --resource-group "$RESOURCE_GROUP" --name "$ACI_NAME"
            fi
            break
        fi

        echo -n "."
        sleep 5
    done
    echo ""

    # Cleanup container
    echo -e "${YELLOW}Cleaning up extraction container...${NC}"
    az container delete \
        --resource-group "$RESOURCE_GROUP" \
        --name "$ACI_NAME" \
        --yes \
        --output none 2>/dev/null || true
}

# Cleanup local archive
cleanup() {
    if [ -f "$ARCHIVE_PATH" ]; then
        echo -e "${YELLOW}Cleaning up local archive...${NC}"
        rm "$ARCHIVE_PATH"
    fi
}

# Verify upload
verify_upload() {
    echo -e "${YELLOW}Verifying upload...${NC}"

    EXISTS=$(az storage file exists \
        --share-name "$SHARE_NAME" \
        --path "$FOLDER_NAME" \
        --account-name "$PUBLIC_STORAGE" \
        --account-key "$STORAGE_KEY" \
        --query "exists" -o tsv 2>/dev/null)

    if [ "$EXISTS" = "true" ]; then
        echo -e "${GREEN}Folder $FOLDER_NAME exists on Azure Files${NC}"

        # List top-level contents
        echo -e "${BLUE}Contents:${NC}"
        az storage file list \
            --share-name "$SHARE_NAME" \
            --path "$FOLDER_NAME" \
            --account-name "$PUBLIC_STORAGE" \
            --account-key "$STORAGE_KEY" \
            --query "[].name" -o tsv 2>/dev/null | head -10
    else
        echo -e "${YELLOW}Folder not yet visible (may still be extracting)${NC}"
    fi
}

# Main
get_storage_account
get_storage_key
check_symlinks
create_archive
upload_archive
extract_archive
cleanup
verify_upload

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Upload Complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo -e "${YELLOW}Next steps:${NC}"
echo "1. Update CCP4_VERSION in .env.deployment if needed:"
echo "   CCP4_VERSION=$FOLDER_NAME"
echo ""
echo "2. Deploy applications to use the new data:"
echo "   ./deploy-applications.sh"
echo ""
