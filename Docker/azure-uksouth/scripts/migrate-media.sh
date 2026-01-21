#!/bin/bash

# Migrate Media Files from Legacy Storage to New Storage (UK South)
#
# Source: ddudatabasestorageac / ddudatabasefileshare / CompoundDatabaseData/media
# Destination: storprv* (dynamically discovered in ccp4i2-bicep-rg-uksouth) / django-uploads (blob container)
#
# These media files back *File objects in the registry, assay, and construct database.
# The new infrastructure uses Blob storage (cheaper) instead of Azure Files for Django uploads.
#
# Path transformations during migration:
# Note: Django's AzureStorage uses the blob container root (no media/ prefix)
#   - Most files: media/* -> /* (copied to container root)
#   - Batch QC files: media/RegBatchQCFile_NCL-* -> RegisterCompounds/BatchQCFiles/NCL-*
#     (relocates per-compound QC directories into a dedicated subdirectory)
#   - Construct database: media/ConstructDatabase/ConstructDatabase/* -> ConstructDatabase/*
#     (flattens legacy double-nested directory structure)
#
# Environment: Docker/azure-uksouth/.env.deployment
#
# Usage: ./migrate-media.sh [check|copy]
#   check    - Check for filename clashes and show what would be copied
#   copy     - Perform the actual copy

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
SOURCE_PATH="CompoundDatabaseData/media"

# Destination storage account (new private storage)
# Will be determined dynamically - looks for storage account starting with 'storprv'
# Uses Blob storage container instead of File share (per infrastructure.bicep)
# Note: Django's AzureStorage uses the container root, not a media/ subdirectory
DEST_CONTAINER="django-uploads"
DEST_PATH=""

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

    # Wait for network rules to propagate (Azure can take 30+ seconds)
    echo -e "${YELLOW}Waiting 30 seconds for network rules to propagate...${NC}"
    sleep 30
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

# Check source and destination
check_clashes() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Checking Media Migration${NC}"
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

    # Check source media directory (File Share)
    echo -e "${YELLOW}Checking source media directory...${NC}"
    SOURCE_FILES=$(az storage file list \
        --share-name "$SOURCE_SHARE" \
        --account-name "$SOURCE_STORAGE_ACCOUNT" \
        --account-key "$SOURCE_KEY" \
        --path "$SOURCE_PATH" \
        --query "length(@)" -o tsv 2>/dev/null || echo "0")

    echo -e "${GREEN}Found $SOURCE_FILES items in source ($SOURCE_PATH)${NC}"

    # List top-level directories in source media
    echo ""
    echo -e "${YELLOW}Source media subdirectories:${NC}"
    BATCH_QC_COUNT=0
    HAS_CONSTRUCT_DB=false
    HAS_ASSAY_COMPOUNDS=false
    HAS_REGISTER_COMPOUNDS=false
    OTHER_DIRS=""
    while read dir; do
        if [[ "$dir" == RegBatchQCFile_* ]]; then
            ((BATCH_QC_COUNT++))
        elif [[ "$dir" == "ConstructDatabase" ]]; then
            HAS_CONSTRUCT_DB=true
        elif [[ "$dir" == "AssayCompounds" ]]; then
            HAS_ASSAY_COMPOUNDS=true
        elif [[ "$dir" == "RegisterCompounds" ]]; then
            HAS_REGISTER_COMPOUNDS=true
        else
            OTHER_DIRS="$OTHER_DIRS$dir\n"
        fi
    done < <(az storage file list \
        --share-name "$SOURCE_SHARE" \
        --account-name "$SOURCE_STORAGE_ACCOUNT" \
        --account-key "$SOURCE_KEY" \
        --path "$SOURCE_PATH" \
        --query "[?type=='dir'].name" -o tsv 2>/dev/null)

    # Show non-QC directories
    echo -e "$OTHER_DIRS" | while read dir; do
        if [ -n "$dir" ]; then
            echo -e "  ${BLUE}$dir${NC}"
        fi
    done

    # Summarize batch QC directories
    if [ $BATCH_QC_COUNT -gt 0 ]; then
        echo -e "  ${CYAN}RegBatchQCFile_NCL-* (${BATCH_QC_COUNT} directories)${NC}"
        echo -e "    ${YELLOW}-> Will be relocated to: RegisterCompounds/BatchQCFiles/NCL-*${NC}"
    fi

    # Note ConstructDatabase special handling
    if [ "$HAS_CONSTRUCT_DB" = true ]; then
        echo -e "  ${CYAN}ConstructDatabase${NC}"
        echo -e "    ${YELLOW}-> Will flatten and copy to container root: ConstructDatabase/ConstructDatabase/* -> ConstructDatabase/*${NC}"
    fi

    # Check destination (Blob Container)
    echo ""
    echo -e "${YELLOW}Checking destination blob container...${NC}"

    # Check if destination container exists
    CONTAINER_EXISTS=$(az storage container exists \
        --name "$DEST_CONTAINER" \
        --account-name "$DEST_STORAGE_ACCOUNT" \
        --account-key "$DEST_KEY" \
        --query "exists" -o tsv 2>/dev/null || echo "false")

    if [ "$CONTAINER_EXISTS" != "true" ]; then
        echo -e "${YELLOW}Destination container '$DEST_CONTAINER' does not exist${NC}"
        echo -e "${YELLOW}It should be created by infrastructure.bicep deployment${NC}"
        DEST_BLOBS=0
    else
        # Count blobs in the media prefix
        DEST_BLOBS=$(az storage blob list \
            --container-name "$DEST_CONTAINER" \
            --account-name "$DEST_STORAGE_ACCOUNT" \
            --account-key "$DEST_KEY" \
            --prefix "$DEST_PATH/" \
            --query "length(@)" -o tsv 2>/dev/null || echo "0")
        echo -e "${GREEN}Found $DEST_BLOBS blobs in destination ($DEST_CONTAINER/$DEST_PATH)${NC}"
    fi

    # Summary
    echo ""
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Migration Summary${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""
    echo "Source (Azure File Share):"
    echo "  Storage Account: $SOURCE_STORAGE_ACCOUNT"
    echo "  File Share: $SOURCE_SHARE"
    echo "  Path: $SOURCE_PATH"
    echo "  Items: $SOURCE_FILES"
    echo ""
    echo "Destination (Azure Blob Storage):"
    echo "  Storage Account: $DEST_STORAGE_ACCOUNT"
    echo "  Container: $DEST_CONTAINER"
    echo "  Path: $DEST_PATH/"
    echo "  Current Blobs: $DEST_BLOBS"
    echo ""
    echo -e "${YELLOW}Run '$0 copy' to perform the migration${NC}"
}

# Perform the copy
do_copy() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Migrating Media Files${NC}"
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

    # Verify destination container exists (should be created by infrastructure.bicep)
    echo -e "${YELLOW}Checking destination container exists...${NC}"
    CONTAINER_EXISTS=$(az storage container exists \
        --name "$DEST_CONTAINER" \
        --account-name "$DEST_STORAGE_ACCOUNT" \
        --account-key "$DEST_KEY" \
        --query "exists" -o tsv 2>/dev/null || echo "false")

    if [ "$CONTAINER_EXISTS" != "true" ]; then
        echo -e "${RED}Destination container '$DEST_CONTAINER' does not exist${NC}"
        echo -e "${YELLOW}Please deploy infrastructure.bicep first to create the container${NC}"
        exit 1
    fi
    echo -e "${GREEN}Container '$DEST_CONTAINER' exists${NC}"

    # Generate SAS tokens
    # Source: File share needs 'f' service type
    # Destination: Blob container needs 'b' service type
    echo -e "${YELLOW}Generating SAS tokens...${NC}"

    # Source SAS for File share
    local source_expiry=$(date -u -v+24H "+%Y-%m-%dT%H:%MZ" 2>/dev/null || date -u -d "+24 hours" "+%Y-%m-%dT%H:%MZ")
    SOURCE_SAS=$(az storage account generate-sas \
        --account-name "$SOURCE_STORAGE_ACCOUNT" \
        --account-key "$SOURCE_KEY" \
        --permissions "rl" \
        --resource-types sco \
        --services f \
        --expiry "$source_expiry" \
        -o tsv)

    # Destination SAS for Blob storage
    local dest_expiry=$(date -u -v+24H "+%Y-%m-%dT%H:%MZ" 2>/dev/null || date -u -d "+24 hours" "+%Y-%m-%dT%H:%MZ")
    DEST_SAS=$(az storage account generate-sas \
        --account-name "$DEST_STORAGE_ACCOUNT" \
        --account-key "$DEST_KEY" \
        --permissions "rwdlac" \
        --resource-types sco \
        --services b \
        --expiry "$dest_expiry" \
        -o tsv)

    if [ -z "$SOURCE_SAS" ] || [ -z "$DEST_SAS" ]; then
        echo -e "${RED}Failed to generate SAS tokens${NC}"
        exit 1
    fi

    echo -e "${GREEN}SAS tokens generated (valid for 24 hours)${NC}"
    echo ""

    # Get list of top-level directories to determine what needs relocation
    echo -e "${YELLOW}Analyzing source directories...${NC}"
    BATCH_QC_DIRS=()
    OTHER_ITEMS=()
    HAS_CONSTRUCT_DB=false

    # Track directories that have double-nested structure in legacy storage
    HAS_ASSAY_COMPOUNDS=false
    HAS_REGISTER_COMPOUNDS=false

    while read item; do
        if [[ "$item" == RegBatchQCFile_* ]]; then
            BATCH_QC_DIRS+=("$item")
        elif [[ "$item" == "ConstructDatabase" ]]; then
            HAS_CONSTRUCT_DB=true
        elif [[ "$item" == "AssayCompounds" ]]; then
            HAS_ASSAY_COMPOUNDS=true
        elif [[ "$item" == "RegisterCompounds" ]]; then
            HAS_REGISTER_COMPOUNDS=true
        else
            OTHER_ITEMS+=("$item")
        fi
    done < <(az storage file list \
        --share-name "$SOURCE_SHARE" \
        --account-name "$SOURCE_STORAGE_ACCOUNT" \
        --account-key "$SOURCE_KEY" \
        --path "$SOURCE_PATH" \
        --query "[].name" -o tsv 2>/dev/null)

    echo -e "${GREEN}Found ${#BATCH_QC_DIRS[@]} batch QC directories to relocate${NC}"
    echo -e "${GREEN}Found ${#OTHER_ITEMS[@]} other items to copy${NC}"
    if [ "$HAS_CONSTRUCT_DB" = true ]; then
        echo -e "${GREEN}Found ConstructDatabase directory (will flatten nested structure)${NC}"
    fi
    if [ "$HAS_ASSAY_COMPOUNDS" = true ]; then
        echo -e "${GREEN}Found AssayCompounds directory (will flatten nested structure)${NC}"
    fi
    if [ "$HAS_REGISTER_COMPOUNDS" = true ]; then
        echo -e "${GREEN}Found RegisterCompounds directory (will flatten nested structure)${NC}"
    fi
    echo ""

    # Step 1: Copy non-QC items (preserving paths)
    if [ ${#OTHER_ITEMS[@]} -gt 0 ]; then
        echo -e "${CYAN}----------------------------------------${NC}"
        echo -e "${CYAN}Step 1: Copying non-QC media files${NC}"
        echo -e "${CYAN}----------------------------------------${NC}"

        for item in "${OTHER_ITEMS[@]}"; do
            echo -e "${YELLOW}  Copying: $item${NC}"
            SOURCE_URL="https://${SOURCE_STORAGE_ACCOUNT}.file.core.windows.net/${SOURCE_SHARE}/${SOURCE_PATH}/${item}?${SOURCE_SAS}"
            DEST_URL="https://${DEST_STORAGE_ACCOUNT}.blob.core.windows.net/${DEST_CONTAINER}/${item}?${DEST_SAS}"

            azcopy copy \
                "$SOURCE_URL" \
                "$DEST_URL" \
                --recursive \
                --skip-version-check \
                --log-level=ERROR 2>/dev/null

            if [ $? -ne 0 ]; then
                echo -e "${RED}    Error copying $item${NC}"
            fi
        done
        echo -e "${GREEN}Non-QC files copied${NC}"
        echo ""
    fi

    # Step 2: Copy batch QC directories with path transformation
    # RegBatchQCFile_NCL-XXXXX -> RegisterCompounds/BatchQCFiles/NCL-XXXXX
    if [ ${#BATCH_QC_DIRS[@]} -gt 0 ]; then
        echo -e "${CYAN}----------------------------------------${NC}"
        echo -e "${CYAN}Step 2: Relocating batch QC directories${NC}"
        echo -e "${CYAN}  From: media/RegBatchQCFile_NCL-*${NC}"
        echo -e "${CYAN}  To:   RegisterCompounds/BatchQCFiles/NCL-*${NC}"
        echo -e "${CYAN}----------------------------------------${NC}"

        for qc_dir in "${BATCH_QC_DIRS[@]}"; do
            # Extract compound ID from directory name (e.g., RegBatchQCFile_NCL-00029551 -> NCL-00029551)
            COMPOUND_ID="${qc_dir#RegBatchQCFile_}"
            NEW_PATH="RegisterCompounds/BatchQCFiles/${COMPOUND_ID}"

            echo -e "${YELLOW}  $qc_dir -> $NEW_PATH${NC}"

            SOURCE_URL="https://${SOURCE_STORAGE_ACCOUNT}.file.core.windows.net/${SOURCE_SHARE}/${SOURCE_PATH}/${qc_dir}/*?${SOURCE_SAS}"
            DEST_URL="https://${DEST_STORAGE_ACCOUNT}.blob.core.windows.net/${DEST_CONTAINER}/${NEW_PATH}/?${DEST_SAS}"

            azcopy copy \
                "$SOURCE_URL" \
                "$DEST_URL" \
                --recursive \
                --skip-version-check \
                --log-level=ERROR 2>/dev/null

            if [ $? -ne 0 ]; then
                echo -e "${RED}    Error copying $qc_dir${NC}"
            fi
        done
        echo -e "${GREEN}Batch QC directories relocated${NC}"
        echo ""
    fi

    # Step 3: Copy ConstructDatabase with path flattening
    # Legacy data has double-nested structure: ConstructDatabase/ConstructDatabase/NCLCON-*
    # Django expects files at container root: ConstructDatabase/NCLCON-* (no media/ prefix)
    # This is because django-storages uses the container as the root, not a media/ subdirectory
    if [ "$HAS_CONSTRUCT_DB" = true ]; then
        echo -e "${CYAN}----------------------------------------${NC}"
        echo -e "${CYAN}Step 3: Flattening ConstructDatabase${NC}"
        echo -e "${CYAN}  From: media/ConstructDatabase/ConstructDatabase/*${NC}"
        echo -e "${CYAN}  To:   ConstructDatabase/* (container root, no media/ prefix)${NC}"
        echo -e "${CYAN}----------------------------------------${NC}"

        # Copy from the nested ConstructDatabase/ConstructDatabase/ to container root ConstructDatabase/
        # Note: We copy to container root (no DEST_PATH) because django-storages doesn't use media/ prefix
        SOURCE_URL="https://${SOURCE_STORAGE_ACCOUNT}.file.core.windows.net/${SOURCE_SHARE}/${SOURCE_PATH}/ConstructDatabase/ConstructDatabase/*?${SOURCE_SAS}"
        DEST_URL="https://${DEST_STORAGE_ACCOUNT}.blob.core.windows.net/${DEST_CONTAINER}/ConstructDatabase/?${DEST_SAS}"

        echo -e "${YELLOW}  Copying nested ConstructDatabase contents...${NC}"
        azcopy copy \
            "$SOURCE_URL" \
            "$DEST_URL" \
            --recursive \
            --skip-version-check \
            --log-level=ERROR 2>/dev/null

        if [ $? -ne 0 ]; then
            echo -e "${RED}    Error copying ConstructDatabase${NC}"
            echo -e "${YELLOW}    Attempting direct copy (source may not be nested)...${NC}"
            # Fallback: try direct copy in case the source isn't double-nested
            SOURCE_URL="https://${SOURCE_STORAGE_ACCOUNT}.file.core.windows.net/${SOURCE_SHARE}/${SOURCE_PATH}/ConstructDatabase/*?${SOURCE_SAS}"
            azcopy copy \
                "$SOURCE_URL" \
                "$DEST_URL" \
                --recursive \
                --skip-version-check \
                --log-level=ERROR 2>/dev/null
        fi

        echo -e "${GREEN}ConstructDatabase copied${NC}"
        echo ""
    fi

    # Step 4: Copy AssayCompounds with path flattening
    # Legacy data has double-nested structure: AssayCompounds/AssayCompounds/*
    # Django expects: AssayCompounds/*
    if [ "$HAS_ASSAY_COMPOUNDS" = true ]; then
        echo -e "${CYAN}----------------------------------------${NC}"
        echo -e "${CYAN}Step 4: Flattening AssayCompounds${NC}"
        echo -e "${CYAN}  From: media/AssayCompounds/AssayCompounds/*${NC}"
        echo -e "${CYAN}  To:   AssayCompounds/* (container root)${NC}"
        echo -e "${CYAN}----------------------------------------${NC}"

        SOURCE_URL="https://${SOURCE_STORAGE_ACCOUNT}.file.core.windows.net/${SOURCE_SHARE}/${SOURCE_PATH}/AssayCompounds/AssayCompounds/*?${SOURCE_SAS}"
        DEST_URL="https://${DEST_STORAGE_ACCOUNT}.blob.core.windows.net/${DEST_CONTAINER}/AssayCompounds/?${DEST_SAS}"

        echo -e "${YELLOW}  Copying nested AssayCompounds contents...${NC}"
        azcopy copy \
            "$SOURCE_URL" \
            "$DEST_URL" \
            --recursive \
            --skip-version-check \
            --log-level=ERROR 2>/dev/null

        if [ $? -ne 0 ]; then
            echo -e "${RED}    Error copying AssayCompounds${NC}"
            echo -e "${YELLOW}    Attempting direct copy (source may not be nested)...${NC}"
            SOURCE_URL="https://${SOURCE_STORAGE_ACCOUNT}.file.core.windows.net/${SOURCE_SHARE}/${SOURCE_PATH}/AssayCompounds/*?${SOURCE_SAS}"
            azcopy copy \
                "$SOURCE_URL" \
                "$DEST_URL" \
                --recursive \
                --skip-version-check \
                --log-level=ERROR 2>/dev/null
        fi

        echo -e "${GREEN}AssayCompounds copied${NC}"
        echo ""
    fi

    # Step 5: Copy RegisterCompounds with path flattening
    # Legacy data has double-nested structure: RegisterCompounds/RegisterCompounds/*
    # Django expects: RegisterCompounds/*
    if [ "$HAS_REGISTER_COMPOUNDS" = true ]; then
        echo -e "${CYAN}----------------------------------------${NC}"
        echo -e "${CYAN}Step 5: Flattening RegisterCompounds${NC}"
        echo -e "${CYAN}  From: media/RegisterCompounds/RegisterCompounds/*${NC}"
        echo -e "${CYAN}  To:   RegisterCompounds/* (container root)${NC}"
        echo -e "${CYAN}----------------------------------------${NC}"

        SOURCE_URL="https://${SOURCE_STORAGE_ACCOUNT}.file.core.windows.net/${SOURCE_SHARE}/${SOURCE_PATH}/RegisterCompounds/RegisterCompounds/*?${SOURCE_SAS}"
        DEST_URL="https://${DEST_STORAGE_ACCOUNT}.blob.core.windows.net/${DEST_CONTAINER}/RegisterCompounds/?${DEST_SAS}"

        echo -e "${YELLOW}  Copying nested RegisterCompounds contents...${NC}"
        azcopy copy \
            "$SOURCE_URL" \
            "$DEST_URL" \
            --recursive \
            --skip-version-check \
            --log-level=ERROR 2>/dev/null

        if [ $? -ne 0 ]; then
            echo -e "${RED}    Error copying RegisterCompounds${NC}"
            echo -e "${YELLOW}    Attempting direct copy (source may not be nested)...${NC}"
            SOURCE_URL="https://${SOURCE_STORAGE_ACCOUNT}.file.core.windows.net/${SOURCE_SHARE}/${SOURCE_PATH}/RegisterCompounds/*?${SOURCE_SAS}"
            azcopy copy \
                "$SOURCE_URL" \
                "$DEST_URL" \
                --recursive \
                --skip-version-check \
                --log-level=ERROR 2>/dev/null
        fi

        echo -e "${GREEN}RegisterCompounds copied${NC}"
        echo ""
    fi

    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}  Media Migration Complete!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
    echo -e "${YELLOW}Next steps:${NC}"
    echo "1. Verify the data: az storage blob list --container-name $DEST_CONTAINER --account-name $DEST_STORAGE_ACCOUNT"
    echo "2. Ensure Django is configured to use Azure Blob Storage (django-storages)"
    echo "3. Update AZURE_CONTAINER and DEFAULT_FILE_STORAGE settings if needed"
    echo ""
}

# List blobs after migration (for verification)
list_blobs() {
    echo -e "${CYAN}Listing migrated media blobs...${NC}"

    get_dest_storage_account
    DEST_KEY=$(get_storage_key "$DEST_STORAGE_ACCOUNT" "$RESOURCE_GROUP")

    if [ -z "$DEST_KEY" ]; then
        echo -e "${RED}Failed to get storage key${NC}"
        exit 1
    fi

    echo -e "${YELLOW}Blobs in $DEST_CONTAINER (container root):${NC}"
    az storage blob list \
        --container-name "$DEST_CONTAINER" \
        --account-name "$DEST_STORAGE_ACCOUNT" \
        --account-key "$DEST_KEY" \
        --query "[].{Name:name, Size:properties.contentLength}" \
        -o table
}

# Show usage
show_usage() {
    echo "Media Files Migration Script"
    echo ""
    echo "Migrates media files (registry, assay, construct database attachments)"
    echo "from legacy File Share storage to new Blob storage."
    echo ""
    echo "Source (Azure File Share):"
    echo "  Storage Account: $SOURCE_STORAGE_ACCOUNT"
    echo "  File Share: $SOURCE_SHARE"
    echo "  Path: $SOURCE_PATH"
    echo ""
    echo "Destination (Azure Blob Storage):"
    echo "  Storage Account: storprv* (auto-discovered)"
    echo "  Container: $DEST_CONTAINER"
    echo "  Path: (container root - Django uses root, no media/ prefix)"
    echo ""
    echo "Usage: $0 <command>"
    echo ""
    echo "Commands:"
    echo "  check     Check source and destination (run this first)"
    echo "  copy      Perform the migration"
    echo "  list      List blobs after migration (verification)"
    echo ""
    echo "Recommended workflow:"
    echo "  1. $0 check     # Review what will be copied"
    echo "  2. $0 copy      # Perform the copy"
    echo "  3. $0 list      # Verify the copied blobs"
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
    list)
        list_blobs
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
