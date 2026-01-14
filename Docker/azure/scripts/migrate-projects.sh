#!/bin/bash

# Migrate CCP4i2 Projects from Legacy Storage to New Storage
#
# Source: ddudatabasestorageac / ddudatebasefileshare / CompoundDatabaseData / {CCP4I2_PROJECTS, CCP4I2_NEW_PROJECTS}
# Destination: stornekmayz3n2 / ccp4i2-projects
#
# Usage: ./migrate-projects.sh [check|copy|snapshot]
#   check    - Check for filename clashes and show what would be copied
#   copy     - Perform the actual copy
#   snapshot - Create a snapshot after copy is complete

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

# Source storage account (legacy)
SOURCE_STORAGE_ACCOUNT="ddudatabasestorageac"
SOURCE_SHARE="ddudatabasefileshare"
SOURCE_BASE_PATH="CompoundDatabaseData"
SOURCE_FOLDERS=("CCP4I2_PROJECTS" "CCP4I2_NEW_PROJECTS")

# Destination storage account (new private storage)
# Will be determined dynamically - looks for storage account starting with 'storprv'
DEST_SHARE="ccp4i2-projects"

# Get destination storage account dynamically
get_dest_storage_account() {
    DEST_STORAGE_ACCOUNT=$(az storage account list \
        --resource-group $RESOURCE_GROUP \
        --query "[?starts_with(name, 'storprv')].name" \
        -o tsv 2>/dev/null | head -1)

    if [ -z "$DEST_STORAGE_ACCOUNT" ]; then
        echo -e "${RED}Error: Private storage account (storprv*) not found${NC}"
        echo -e "${YELLOW}Have you deployed the updated infrastructure.bicep?${NC}"
        exit 1
    fi

    echo -e "${GREEN}Destination storage: $DEST_STORAGE_ACCOUNT${NC}"
}

# Check if azcopy is installed
check_azcopy() {
    if ! command -v azcopy &> /dev/null; then
        echo -e "${RED}Error: azcopy is not installed${NC}"
        echo ""
        echo "Install with:"
        echo "  brew install azcopy"
        echo ""
        echo "Or download from:"
        echo "  https://docs.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-v10"
        exit 1
    fi
    echo -e "${GREEN}azcopy found: $(which azcopy)${NC}"
}

# Get storage account key
get_storage_key() {
    local account=$1
    local rg=$2

    if [ -z "$rg" ]; then
        # Try to find resource group
        rg=$(az storage account list --query "[?name=='$account'].resourceGroup" -o tsv 2>/dev/null)
    fi

    if [ -z "$rg" ]; then
        echo -e "${RED}Could not find resource group for storage account: $account${NC}"
        return 1
    fi

    az storage account keys list \
        --account-name "$account" \
        --resource-group "$rg" \
        --query "[0].value" -o tsv
}

# Enable network access for a storage account
enable_network_access() {
    local account=$1
    local rg=$2

    echo -e "${YELLOW}Enabling network access for $account...${NC}"

    CURRENT_IP=$(curl -s https://ipinfo.io/ip 2>/dev/null || curl -s https://api.ipify.org 2>/dev/null)

    if [ -z "$CURRENT_IP" ]; then
        echo -e "${RED}Could not determine public IP${NC}"
        return 1
    fi

    echo -e "${GREEN}Current IP: $CURRENT_IP${NC}"

    # Try to add IP rule (may fail if already exists or no permissions)
    az storage account network-rule add \
        --account-name "$account" \
        --resource-group "$rg" \
        --ip-address "$CURRENT_IP" \
        --output none 2>/dev/null || true

    # Temporarily allow access
    az storage account update \
        --name "$account" \
        --resource-group "$rg" \
        --default-action Allow \
        --output none 2>/dev/null || true

    sleep 3
}

# Generate SAS token for a storage account
generate_sas() {
    local account=$1
    local key=$2
    local permissions=$3  # r=read, w=write, l=list, etc.

    # Generate SAS valid for 24 hours
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

# Check for filename clashes
check_clashes() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Checking for Filename Clashes${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    check_azcopy
    get_dest_storage_account

    # Get keys
    echo -e "${YELLOW}Getting storage account keys...${NC}"

    SOURCE_RG=$(az storage account list --query "[?name=='$SOURCE_STORAGE_ACCOUNT'].resourceGroup" -o tsv 2>/dev/null)
    DEST_RG=$RESOURCE_GROUP

    enable_network_access "$SOURCE_STORAGE_ACCOUNT" "$SOURCE_RG"
    enable_network_access "$DEST_STORAGE_ACCOUNT" "$DEST_RG"

    SOURCE_KEY=$(get_storage_key "$SOURCE_STORAGE_ACCOUNT" "$SOURCE_RG")
    DEST_KEY=$(get_storage_key "$DEST_STORAGE_ACCOUNT" "$DEST_RG")

    if [ -z "$SOURCE_KEY" ] || [ -z "$DEST_KEY" ]; then
        echo -e "${RED}Failed to get storage keys${NC}"
        exit 1
    fi

    # Create temp files for directory listings
    TEMP_DIR=$(mktemp -d)
    SOURCE_DIRS="$TEMP_DIR/source_dirs.txt"
    DEST_DIRS="$TEMP_DIR/dest_dirs.txt"

    # List source directories from both folders
    echo -e "${YELLOW}Listing source directories...${NC}"
    > "$SOURCE_DIRS"

    for folder in "${SOURCE_FOLDERS[@]}"; do
        echo -e "${BLUE}  Checking $SOURCE_BASE_PATH/$folder...${NC}"
        az storage file list \
            --share-name "$SOURCE_SHARE" \
            --account-name "$SOURCE_STORAGE_ACCOUNT" \
            --account-key "$SOURCE_KEY" \
            --path "$SOURCE_BASE_PATH/$folder" \
            --query "[?type=='dir'].name" -o tsv 2>/dev/null >> "$SOURCE_DIRS" || true
    done

    SOURCE_COUNT=$(wc -l < "$SOURCE_DIRS" | tr -d ' ')
    echo -e "${GREEN}Found $SOURCE_COUNT project directories in source${NC}"

    # List destination directories
    echo -e "${YELLOW}Listing destination directories...${NC}"
    az storage file list \
        --share-name "$DEST_SHARE" \
        --account-name "$DEST_STORAGE_ACCOUNT" \
        --account-key "$DEST_KEY" \
        --query "[?type=='dir'].name" -o tsv 2>/dev/null > "$DEST_DIRS" || true

    DEST_COUNT=$(wc -l < "$DEST_DIRS" | tr -d ' ')
    echo -e "${GREEN}Found $DEST_COUNT project directories in destination${NC}"

    # Check for clashes
    echo ""
    echo -e "${YELLOW}Checking for clashes...${NC}"
    CLASHES=$(comm -12 <(sort "$SOURCE_DIRS" | uniq) <(sort "$DEST_DIRS"))

    if [ -n "$CLASHES" ]; then
        echo -e "${RED}========================================${NC}"
        echo -e "${RED}  WARNING: Filename Clashes Found!${NC}"
        echo -e "${RED}========================================${NC}"
        echo "$CLASHES" | while read clash; do
            echo -e "${RED}  - $clash${NC}"
        done
        CLASH_COUNT=$(echo "$CLASHES" | wc -l | tr -d ' ')
        echo ""
        echo -e "${RED}$CLASH_COUNT directories would clash${NC}"
        echo -e "${YELLOW}These will be SKIPPED during copy (azcopy default behavior)${NC}"
        echo -e "${YELLOW}Use --overwrite=true to overwrite existing files${NC}"
    else
        echo -e "${GREEN}No filename clashes detected!${NC}"
    fi

    # Summary
    echo ""
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Migration Summary${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""
    echo "Source: $SOURCE_STORAGE_ACCOUNT/$SOURCE_SHARE"
    for folder in "${SOURCE_FOLDERS[@]}"; do
        count=$(az storage file list \
            --share-name "$SOURCE_SHARE" \
            --account-name "$SOURCE_STORAGE_ACCOUNT" \
            --account-key "$SOURCE_KEY" \
            --path "$SOURCE_BASE_PATH/$folder" \
            --query "[?type=='dir'] | length(@)" -o tsv 2>/dev/null || echo "0")
        echo "  - $SOURCE_BASE_PATH/$folder: $count directories"
    done
    echo ""
    echo "Destination: $DEST_STORAGE_ACCOUNT/$DEST_SHARE"
    echo "  - Current directories: $DEST_COUNT"
    echo ""
    echo -e "${YELLOW}Run '$0 copy' to perform the migration${NC}"

    # Cleanup
    rm -rf "$TEMP_DIR"
}

# Perform the copy
do_copy() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Migrating Project Data${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    check_azcopy
    get_dest_storage_account

    # Get resource groups
    SOURCE_RG=$(az storage account list --query "[?name=='$SOURCE_STORAGE_ACCOUNT'].resourceGroup" -o tsv 2>/dev/null)
    DEST_RG=$RESOURCE_GROUP

    if [ -z "$SOURCE_RG" ]; then
        echo -e "${RED}Could not find resource group for source storage account${NC}"
        exit 1
    fi

    echo -e "${BLUE}Source: $SOURCE_STORAGE_ACCOUNT (RG: $SOURCE_RG)${NC}"
    echo -e "${BLUE}Destination: $DEST_STORAGE_ACCOUNT (RG: $DEST_RG)${NC}"
    echo ""

    # Enable network access
    enable_network_access "$SOURCE_STORAGE_ACCOUNT" "$SOURCE_RG"
    enable_network_access "$DEST_STORAGE_ACCOUNT" "$DEST_RG"

    # Get storage keys
    echo -e "${YELLOW}Getting storage account keys...${NC}"
    SOURCE_KEY=$(get_storage_key "$SOURCE_STORAGE_ACCOUNT" "$SOURCE_RG")
    DEST_KEY=$(get_storage_key "$DEST_STORAGE_ACCOUNT" "$DEST_RG")

    if [ -z "$SOURCE_KEY" ] || [ -z "$DEST_KEY" ]; then
        echo -e "${RED}Failed to get storage keys${NC}"
        exit 1
    fi

    # Generate SAS tokens
    echo -e "${YELLOW}Generating SAS tokens...${NC}"
    SOURCE_SAS=$(generate_sas "$SOURCE_STORAGE_ACCOUNT" "$SOURCE_KEY" "rl")
    DEST_SAS=$(generate_sas "$DEST_STORAGE_ACCOUNT" "$DEST_KEY" "rwdl")

    if [ -z "$SOURCE_SAS" ] || [ -z "$DEST_SAS" ]; then
        echo -e "${RED}Failed to generate SAS tokens${NC}"
        exit 1
    fi

    echo -e "${GREEN}SAS tokens generated (valid for 24 hours)${NC}"
    echo ""

    # Copy each source folder
    for folder in "${SOURCE_FOLDERS[@]}"; do
        echo -e "${CYAN}----------------------------------------${NC}"
        echo -e "${CYAN}Copying: $folder${NC}"
        echo -e "${CYAN}----------------------------------------${NC}"

        SOURCE_URL="https://${SOURCE_STORAGE_ACCOUNT}.file.core.windows.net/${SOURCE_SHARE}/${SOURCE_BASE_PATH}/${folder}/*?${SOURCE_SAS}"
        DEST_URL="https://${DEST_STORAGE_ACCOUNT}.file.core.windows.net/${DEST_SHARE}/?${DEST_SAS}"

        echo -e "${YELLOW}Source: $SOURCE_STORAGE_ACCOUNT/$SOURCE_SHARE/$SOURCE_BASE_PATH/$folder${NC}"
        echo -e "${YELLOW}Dest:   $DEST_STORAGE_ACCOUNT/$DEST_SHARE${NC}"
        echo ""

        # Run azcopy with progress
        # --recursive: copy subdirectories
        # --preserve-smb-info: preserve file metadata
        # --skip-version-check: faster startup
        azcopy copy \
            "$SOURCE_URL" \
            "$DEST_URL" \
            --recursive \
            --preserve-smb-info=false \
            --skip-version-check \
            --log-level=WARNING

        if [ $? -eq 0 ]; then
            echo -e "${GREEN}Successfully copied $folder${NC}"
        else
            echo -e "${RED}Error copying $folder${NC}"
        fi
        echo ""
    done

    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}  Migration Complete!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
    echo -e "${YELLOW}Next steps:${NC}"
    echo "1. Verify the data: az storage file list --share-name $DEST_SHARE --account-name $DEST_STORAGE_ACCOUNT"
    echo "2. Create a backup snapshot: $0 snapshot"
    echo "3. Update database references if needed"
    echo ""
}

# Create snapshot after migration
create_snapshot() {
    echo -e "${CYAN}Creating backup snapshot...${NC}"
    "$SCRIPT_DIR/backup-projects.sh" create
}

# Show usage
show_usage() {
    echo "CCP4i2 Projects Migration Script"
    echo ""
    echo "Migrates project data from legacy storage to new storage account."
    echo ""
    echo "Source:"
    echo "  Storage Account: $SOURCE_STORAGE_ACCOUNT"
    echo "  File Share: $SOURCE_SHARE"
    echo "  Paths: $SOURCE_BASE_PATH/{CCP4I2_PROJECTS, CCP4I2_NEW_PROJECTS}"
    echo ""
    echo "Destination:"
    echo "  Storage Account: $DEST_STORAGE_ACCOUNT"
    echo "  File Share: $DEST_SHARE"
    echo ""
    echo "Usage: $0 <command>"
    echo ""
    echo "Commands:"
    echo "  check     Check for filename clashes (run this first)"
    echo "  copy      Perform the migration"
    echo "  snapshot  Create a backup snapshot after migration"
    echo ""
    echo "Recommended workflow:"
    echo "  1. $0 check     # Review what will be copied"
    echo "  2. $0 copy      # Perform the copy"
    echo "  3. $0 snapshot  # Create backup snapshot"
    echo ""
}

# Main command handler
case "${1:-}" in
    check)
        check_clashes
        ;;
    copy)
        do_copy
        ;;
    snapshot)
        create_snapshot
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
