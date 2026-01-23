#!/bin/bash

# Migrate Media Files from Legacy Storage to New Storage (UK South)
#
# Source: ddudatabasestorageac / ddudatabasefileshare / CompoundDatabaseData/media
# Destination: storprv* (dynamically discovered in ccp4i2-bicep-rg-uksouth) / django-uploads (blob container)
#
# These media files back *File objects in the registry, assay, and construct database.
# The new infrastructure uses Blob storage (cheaper) instead of Azure Files for Django uploads.
#
# Path transformations during migration (legacy -> new):
#   AssayCompounds/Experiments/*  -> compounds/assays/data/*
#   AssayCompounds/DataSeries/*   -> compounds/assays/plots/*  (not used - plots are in Experiments)
#   AssayCompounds/Protocols/*    -> compounds/assays/protocols/*
#   AssayCompounds/svg/*          -> compounds/assays/svg/*
#   RegisterCompounds/svg/*       -> compounds/registry/svg/*
#   RegBatchQCFile_NCL-*          -> compounds/registry/qc/NCL-*
#   ConstructDatabase/*           -> compounds/constructs/*
#
# Environment: Docker/azure-uksouth/.env.deployment
#
# Usage: ./migrate-media.sh [check|copy|rename]
#   check    - Check source directories and show what would be copied
#   copy     - Perform fresh migration with path transformations
#   rename   - Rename already-migrated blobs from legacy to new paths

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
    echo -e "${YELLOW}Source media subdirectories and path transformations:${NC}"
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

    # Show path transformations
    if [ "$HAS_ASSAY_COMPOUNDS" = true ]; then
        echo -e "  ${CYAN}AssayCompounds/${NC}"
        echo -e "    ${YELLOW}Experiments/*  -> compounds/assays/data/*${NC}"
        echo -e "    ${YELLOW}Protocols/*    -> compounds/assays/protocols/*${NC}"
        echo -e "    ${YELLOW}svg/*          -> compounds/assays/svg/*${NC}"
        echo -e "    ${YELLOW}DataSeries/*   -> compounds/assays/plots/*${NC}"
    fi

    if [ "$HAS_REGISTER_COMPOUNDS" = true ]; then
        echo -e "  ${CYAN}RegisterCompounds/${NC}"
        echo -e "    ${YELLOW}svg/*          -> compounds/registry/svg/*${NC}"
    fi

    if [ $BATCH_QC_COUNT -gt 0 ]; then
        echo -e "  ${CYAN}RegBatchQCFile_NCL-* (${BATCH_QC_COUNT} directories)${NC}"
        echo -e "    ${YELLOW}-> compounds/registry/qc/NCL-*${NC}"
    fi

    if [ "$HAS_CONSTRUCT_DB" = true ]; then
        echo -e "  ${CYAN}ConstructDatabase/${NC}"
        echo -e "    ${YELLOW}-> compounds/constructs/*${NC}"
    fi

    # Show other directories
    echo -e "$OTHER_DIRS" | while read dir; do
        if [ -n "$dir" ]; then
            echo -e "  ${BLUE}$dir (will be skipped)${NC}"
        fi
    done

    # Check destination (Blob Container)
    echo ""
    echo -e "${YELLOW}Checking destination blob container...${NC}"

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
        # Count blobs in compounds/ prefix
        DEST_BLOBS=$(az storage blob list \
            --container-name "$DEST_CONTAINER" \
            --account-name "$DEST_STORAGE_ACCOUNT" \
            --account-key "$DEST_KEY" \
            --prefix "compounds/" \
            --query "length(@)" -o tsv 2>/dev/null || echo "0")
        echo -e "${GREEN}Found $DEST_BLOBS blobs in destination (compounds/*)${NC}"
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
    echo "  Target prefix: compounds/"
    echo "  Current Blobs: $DEST_BLOBS"
    echo ""
    echo -e "${YELLOW}Run '$0 copy' to perform the migration${NC}"
}

# Perform the copy with path transformations to new schema
do_copy() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Migrating Media Files (Fresh Copy)${NC}"
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

    # Verify destination container exists
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
    echo -e "${YELLOW}Generating SAS tokens...${NC}"

    local source_expiry=$(date -u -v+24H "+%Y-%m-%dT%H:%MZ" 2>/dev/null || date -u -d "+24 hours" "+%Y-%m-%dT%H:%MZ")
    SOURCE_SAS=$(az storage account generate-sas \
        --account-name "$SOURCE_STORAGE_ACCOUNT" \
        --account-key "$SOURCE_KEY" \
        --permissions "rl" \
        --resource-types sco \
        --services f \
        --expiry "$source_expiry" \
        -o tsv)

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

    # Helper function to copy with path transformation
    copy_with_transform() {
        local src_subpath=$1
        local dest_subpath=$2
        local description=$3

        echo -e "${YELLOW}  $description${NC}"
        echo -e "${BLUE}    From: $src_subpath${NC}"
        echo -e "${BLUE}    To:   $dest_subpath${NC}"

        SOURCE_URL="https://${SOURCE_STORAGE_ACCOUNT}.file.core.windows.net/${SOURCE_SHARE}/${SOURCE_PATH}/${src_subpath}?${SOURCE_SAS}"
        DEST_URL="https://${DEST_STORAGE_ACCOUNT}.blob.core.windows.net/${DEST_CONTAINER}/${dest_subpath}?${DEST_SAS}"

        azcopy copy \
            "$SOURCE_URL" \
            "$DEST_URL" \
            --recursive \
            --skip-version-check \
            --log-level=WARNING 2>&1 | grep -v "^$" | head -5

        if [ ${PIPESTATUS[0]} -ne 0 ]; then
            echo -e "${RED}    Warning: Some errors during copy${NC}"
        else
            echo -e "${GREEN}    Done${NC}"
        fi
    }

    # =========================================================================
    # Step 1: AssayCompounds subdirectories -> compounds/assays/*
    # =========================================================================
    echo -e "${CYAN}----------------------------------------${NC}"
    echo -e "${CYAN}Step 1: Migrating AssayCompounds${NC}"
    echo -e "${CYAN}----------------------------------------${NC}"

    # AssayCompounds/Experiments/* -> compounds/assays/data/*
    copy_with_transform \
        "AssayCompounds/Experiments/*" \
        "compounds/assays/data/" \
        "Experiments -> compounds/assays/data"

    # AssayCompounds/Protocols/* -> compounds/assays/protocols/*
    copy_with_transform \
        "AssayCompounds/Protocols/*" \
        "compounds/assays/protocols/" \
        "Protocols -> compounds/assays/protocols"

    # AssayCompounds/svg/* -> compounds/assays/svg/*
    copy_with_transform \
        "AssayCompounds/svg/*" \
        "compounds/assays/svg/" \
        "svg -> compounds/assays/svg"

    # AssayCompounds/DataSeries/* -> compounds/assays/plots/* (if needed)
    copy_with_transform \
        "AssayCompounds/DataSeries/*" \
        "compounds/assays/plots/" \
        "DataSeries -> compounds/assays/plots"

    echo ""

    # =========================================================================
    # Step 2: RegisterCompounds subdirectories -> compounds/registry/*
    # =========================================================================
    echo -e "${CYAN}----------------------------------------${NC}"
    echo -e "${CYAN}Step 2: Migrating RegisterCompounds${NC}"
    echo -e "${CYAN}----------------------------------------${NC}"

    # RegisterCompounds/svg/* -> compounds/registry/svg/*
    copy_with_transform \
        "RegisterCompounds/svg/*" \
        "compounds/registry/svg/" \
        "svg -> compounds/registry/svg"

    echo ""

    # =========================================================================
    # Step 3: Batch QC files -> compounds/registry/qc/*
    # =========================================================================
    echo -e "${CYAN}----------------------------------------${NC}"
    echo -e "${CYAN}Step 3: Migrating Batch QC Files${NC}"
    echo -e "${CYAN}----------------------------------------${NC}"

    # Get list of RegBatchQCFile_* directories
    BATCH_QC_DIRS=()
    while read item; do
        if [[ "$item" == RegBatchQCFile_* ]]; then
            BATCH_QC_DIRS+=("$item")
        fi
    done < <(az storage file list \
        --share-name "$SOURCE_SHARE" \
        --account-name "$SOURCE_STORAGE_ACCOUNT" \
        --account-key "$SOURCE_KEY" \
        --path "$SOURCE_PATH" \
        --query "[].name" -o tsv 2>/dev/null)

    echo -e "${GREEN}Found ${#BATCH_QC_DIRS[@]} batch QC directories${NC}"

    for qc_dir in "${BATCH_QC_DIRS[@]}"; do
        # Extract compound ID (e.g., RegBatchQCFile_NCL-00029551 -> NCL-00029551)
        COMPOUND_ID="${qc_dir#RegBatchQCFile_}"
        copy_with_transform \
            "${qc_dir}/*" \
            "compounds/registry/qc/${COMPOUND_ID}/" \
            "${qc_dir} -> compounds/registry/qc/${COMPOUND_ID}"
    done

    echo ""

    # =========================================================================
    # Step 4: ConstructDatabase -> compounds/constructs/*
    # =========================================================================
    echo -e "${CYAN}----------------------------------------${NC}"
    echo -e "${CYAN}Step 4: Migrating ConstructDatabase${NC}"
    echo -e "${CYAN}----------------------------------------${NC}"

    # Check if source has nested structure (ConstructDatabase/ConstructDatabase/*)
    # or flat structure (ConstructDatabase/NCLCON-*)
    NESTED_EXISTS=$(az storage file list \
        --share-name "$SOURCE_SHARE" \
        --account-name "$SOURCE_STORAGE_ACCOUNT" \
        --account-key "$SOURCE_KEY" \
        --path "$SOURCE_PATH/ConstructDatabase/ConstructDatabase" \
        --query "length(@)" -o tsv 2>/dev/null || echo "0")

    if [ "$NESTED_EXISTS" -gt 0 ]; then
        # Nested structure - copy from ConstructDatabase/ConstructDatabase/*
        copy_with_transform \
            "ConstructDatabase/ConstructDatabase/*" \
            "compounds/constructs/" \
            "ConstructDatabase (nested) -> compounds/constructs"
    else
        # Flat structure - copy from ConstructDatabase/*
        copy_with_transform \
            "ConstructDatabase/*" \
            "compounds/constructs/" \
            "ConstructDatabase -> compounds/constructs"
    fi

    echo ""

    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}  Media Migration Complete!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
    echo -e "${YELLOW}Next steps:${NC}"
    echo "1. Verify the data: az storage blob list --container-name $DEST_CONTAINER --account-name $DEST_STORAGE_ACCOUNT --prefix compounds/"
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

    echo -e "${YELLOW}Blobs in $DEST_CONTAINER/compounds/:${NC}"
    echo ""

    # Count by subdirectory
    echo -e "${CYAN}Blob counts by directory:${NC}"
    az storage blob list \
        --container-name "$DEST_CONTAINER" \
        --account-name "$DEST_STORAGE_ACCOUNT" \
        --account-key "$DEST_KEY" \
        --prefix "compounds/" \
        --query "[].name" -o tsv 2>/dev/null | \
        sed 's|/[^/]*$||' | sort | uniq -c | sort -rn | head -20

    echo ""
    echo -e "${CYAN}Sample blobs:${NC}"
    az storage blob list \
        --container-name "$DEST_CONTAINER" \
        --account-name "$DEST_STORAGE_ACCOUNT" \
        --account-key "$DEST_KEY" \
        --prefix "compounds/" \
        --num-results 20 \
        --query "[].{Name:name, Size:properties.contentLength}" \
        -o table
}

# Show usage
show_usage() {
    echo "Media Files Migration Script"
    echo ""
    echo "Migrates media files (registry, assay, construct database attachments)"
    echo "from legacy File Share storage to new Blob storage with path transformations."
    echo ""
    echo "Source (Azure File Share):"
    echo "  Storage Account: $SOURCE_STORAGE_ACCOUNT"
    echo "  File Share: $SOURCE_SHARE"
    echo "  Path: $SOURCE_PATH"
    echo ""
    echo "Destination (Azure Blob Storage):"
    echo "  Storage Account: storprv* (auto-discovered)"
    echo "  Container: $DEST_CONTAINER"
    echo "  Target prefix: compounds/"
    echo ""
    echo "Path transformations:"
    echo "  AssayCompounds/Experiments/*  -> compounds/assays/data/*"
    echo "  AssayCompounds/Protocols/*    -> compounds/assays/protocols/*"
    echo "  AssayCompounds/svg/*          -> compounds/assays/svg/*"
    echo "  AssayCompounds/DataSeries/*   -> compounds/assays/plots/*"
    echo "  RegisterCompounds/svg/*       -> compounds/registry/svg/*"
    echo "  RegBatchQCFile_NCL-*          -> compounds/registry/qc/NCL-*"
    echo "  ConstructDatabase/*           -> compounds/constructs/*"
    echo ""
    echo "Usage: $0 <command>"
    echo ""
    echo "Commands:"
    echo "  check     Check source and destination, show path transformations"
    echo "  copy      Perform the migration with path transformations"
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
