#!/bin/bash

# Migrate Fixture Files from Legacy Storage to New Blob Storage (UK South)
#
# Source: ddudatabasestorageac / ddudatabasefileshare / CompoundDatabaseData
#         Fixture files match pattern: YYYYMMDD-HH-MM-*.json
# Destination: storprv* blob container "fixtures" folder in django-uploads
#
# These fixtures are periodic snapshots of the legacy databases:
#   - *-CCP4i2.json (CCP4i2 Django dumpdata)
#   - *-RegisterCompounds.json (Compound registry)
#   - *-AssayCompounds.json (Assay data)
#   - *-ConstructDatabase.json (Construct/plasmid data)
#   - *-auth.json (Django auth users/groups)
#   - *-reversion.json (Django reversion history)
#
# Environment: Docker/azure-uksouth/.env.deployment
#
# Usage: ./migrate-fixtures.sh [list|copy|latest]
#   list   - List available fixtures in source storage
#   copy   - Copy all fixtures to new storage
#   latest - Copy only the latest fixture of each type

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
SOURCE_PATH="CompoundDatabaseData"

# Destination storage (new blob storage)
DEST_CONTAINER="django-uploads"
DEST_PATH="fixtures"

# Fixture types we're looking for (case-sensitive, must match actual filenames)
FIXTURE_TYPES=("CCP4i2" "RegisterCompounds" "AssayCompounds" "ConstructDatabase" "auth" "reversion")

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
        exit 1
    fi
    echo -e "${GREEN}azcopy found: $(which azcopy)${NC}"
}

# Get storage account key
get_storage_key() {
    local account=$1
    local rg=$2

    if [ -z "$rg" ]; then
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

# Generate SAS token for storage account
generate_sas_token() {
    local account=$1
    local key=$2
    local service=$3  # 'blob' or 'file'
    local permissions=$4  # e.g., 'rl' for read+list or 'rwdl' for full

    local expiry=$(date -u -v+24H +"%Y-%m-%dT%H:%MZ" 2>/dev/null || date -u -d "+24 hours" +"%Y-%m-%dT%H:%MZ")

    if [ "$service" = "file" ]; then
        az storage share generate-sas \
            --name "$SOURCE_SHARE" \
            --account-name "$account" \
            --account-key "$key" \
            --permissions "$permissions" \
            --expiry "$expiry" \
            -o tsv 2>/dev/null
    else
        az storage container generate-sas \
            --name "$DEST_CONTAINER" \
            --account-name "$account" \
            --account-key "$key" \
            --permissions "$permissions" \
            --expiry "$expiry" \
            -o tsv 2>/dev/null
    fi
}

# List fixture files in source storage
list_fixtures() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Listing Fixture Files in Legacy Storage${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    # Get source storage key
    SOURCE_KEY=$(get_storage_key "$SOURCE_STORAGE_ACCOUNT")
    if [ -z "$SOURCE_KEY" ]; then
        echo -e "${RED}Failed to get source storage key${NC}"
        exit 1
    fi

    echo -e "${YELLOW}Source: ${SOURCE_STORAGE_ACCOUNT}/${SOURCE_SHARE}/${SOURCE_PATH}${NC}"
    echo ""

    # List files matching fixture pattern
    echo -e "${BLUE}Available fixture files:${NC}"
    echo ""

    for fixture_type in "${FIXTURE_TYPES[@]}"; do
        echo -e "${GREEN}=== ${fixture_type} ===${NC}"

        # List files matching pattern YYYYMMDD-HH-MM-{type}.json
        # Filter: starts with '20' (2020s dates) AND ends with -{type}.json
        az storage file list \
            --share-name "$SOURCE_SHARE" \
            --account-name "$SOURCE_STORAGE_ACCOUNT" \
            --account-key "$SOURCE_KEY" \
            --path "$SOURCE_PATH" \
            --query "[?starts_with(name, '20') && ends_with(name, '-${fixture_type}.json')].{name:name, size:properties.contentLength}" \
            -o table 2>/dev/null | tail -n +3 | sort -r | head -10

        echo ""
    done

    # Show total count
    echo -e "${YELLOW}Summary:${NC}"
    TOTAL=$(az storage file list \
        --share-name "$SOURCE_SHARE" \
        --account-name "$SOURCE_STORAGE_ACCOUNT" \
        --account-key "$SOURCE_KEY" \
        --path "$SOURCE_PATH" \
        --query "length([?ends_with(name, '.json')])" \
        -o tsv 2>/dev/null)
    echo "  Total JSON files in ${SOURCE_PATH}: ${TOTAL}"
}

# Copy fixtures to new storage
copy_fixtures() {
    local mode=$1  # 'all' or 'latest'

    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Copying Fixture Files to New Storage${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    check_azcopy
    get_dest_storage_account

    # Get storage keys
    SOURCE_KEY=$(get_storage_key "$SOURCE_STORAGE_ACCOUNT")
    DEST_KEY=$(get_storage_key "$DEST_STORAGE_ACCOUNT" "$RESOURCE_GROUP")

    if [ -z "$SOURCE_KEY" ] || [ -z "$DEST_KEY" ]; then
        echo -e "${RED}Failed to get storage keys${NC}"
        exit 1
    fi

    # Generate SAS tokens
    echo -e "${YELLOW}Generating SAS tokens...${NC}"
    SOURCE_SAS=$(generate_sas_token "$SOURCE_STORAGE_ACCOUNT" "$SOURCE_KEY" "file" "rl")
    DEST_SAS=$(generate_sas_token "$DEST_STORAGE_ACCOUNT" "$DEST_KEY" "blob" "rwdlac")

    if [ -z "$SOURCE_SAS" ] || [ -z "$DEST_SAS" ]; then
        echo -e "${RED}Failed to generate SAS tokens${NC}"
        exit 1
    fi

    echo -e "${GREEN}SAS tokens generated (valid for 24 hours)${NC}"
    echo ""

    # Ensure destination folder exists (azcopy will create it)
    echo -e "${YELLOW}Destination: ${DEST_STORAGE_ACCOUNT}/${DEST_CONTAINER}/${DEST_PATH}/${NC}"
    echo ""

    if [ "$mode" = "latest" ]; then
        # Copy only the latest fixture of each type
        echo -e "${CYAN}Copying latest fixtures only...${NC}"
        echo ""

        for fixture_type in "${FIXTURE_TYPES[@]}"; do
            echo -e "${GREEN}Finding latest ${fixture_type} fixture...${NC}"

            # Get the latest file for this type (must start with '20' for 2020s dates)
            LATEST_FILE=$(az storage file list \
                --share-name "$SOURCE_SHARE" \
                --account-name "$SOURCE_STORAGE_ACCOUNT" \
                --account-key "$SOURCE_KEY" \
                --path "$SOURCE_PATH" \
                --query "[?starts_with(name, '20') && ends_with(name, '-${fixture_type}.json')].name" \
                -o tsv 2>/dev/null | sort -r | head -1)

            if [ -n "$LATEST_FILE" ]; then
                echo -e "  ${YELLOW}Copying: ${LATEST_FILE}${NC}"

                SOURCE_URL="https://${SOURCE_STORAGE_ACCOUNT}.file.core.windows.net/${SOURCE_SHARE}/${SOURCE_PATH}/${LATEST_FILE}?${SOURCE_SAS}"
                DEST_URL="https://${DEST_STORAGE_ACCOUNT}.blob.core.windows.net/${DEST_CONTAINER}/${DEST_PATH}/${LATEST_FILE}?${DEST_SAS}"

                azcopy copy "$SOURCE_URL" "$DEST_URL" \
                    --skip-version-check \
                    --log-level=ERROR 2>/dev/null

                if [ $? -eq 0 ]; then
                    echo -e "  ${GREEN}✓ Copied successfully${NC}"
                else
                    echo -e "  ${RED}✗ Failed to copy${NC}"
                fi
            else
                echo -e "  ${YELLOW}No fixtures found for ${fixture_type}${NC}"
            fi
            echo ""
        done
    else
        # Copy all fixtures
        echo -e "${CYAN}Copying all fixture files...${NC}"
        echo ""

        for fixture_type in "${FIXTURE_TYPES[@]}"; do
            echo -e "${GREEN}Copying ${fixture_type} fixtures...${NC}"

            # Get all files for this type (must start with '20' for 2020s dates)
            FILES=$(az storage file list \
                --share-name "$SOURCE_SHARE" \
                --account-name "$SOURCE_STORAGE_ACCOUNT" \
                --account-key "$SOURCE_KEY" \
                --path "$SOURCE_PATH" \
                --query "[?starts_with(name, '20') && ends_with(name, '-${fixture_type}.json')].name" \
                -o tsv 2>/dev/null)

            COUNT=0
            while read -r file; do
                if [ -n "$file" ]; then
                    echo -e "  ${YELLOW}Copying: ${file}${NC}"

                    SOURCE_URL="https://${SOURCE_STORAGE_ACCOUNT}.file.core.windows.net/${SOURCE_SHARE}/${SOURCE_PATH}/${file}?${SOURCE_SAS}"
                    DEST_URL="https://${DEST_STORAGE_ACCOUNT}.blob.core.windows.net/${DEST_CONTAINER}/${DEST_PATH}/${file}?${DEST_SAS}"

                    azcopy copy "$SOURCE_URL" "$DEST_URL" \
                        --skip-version-check \
                        --log-level=ERROR \
                        --output-level=quiet 2>/dev/null

                    ((COUNT++))
                fi
            done <<< "$FILES"

            echo -e "  ${GREEN}Copied ${COUNT} ${fixture_type} fixtures${NC}"
            echo ""
        done
    fi

    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}  Fixture Migration Complete!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
    echo -e "${YELLOW}Fixtures are now available at:${NC}"
    echo "  https://${DEST_STORAGE_ACCOUNT}.blob.core.windows.net/${DEST_CONTAINER}/${DEST_PATH}/"
    echo ""
    echo -e "${YELLOW}To list copied fixtures:${NC}"
    echo "  az storage blob list --container-name $DEST_CONTAINER --account-name $DEST_STORAGE_ACCOUNT --prefix ${DEST_PATH}/ --query '[].name' -o table"
}

# Show usage
show_usage() {
    echo "Usage: $0 [list|copy|latest]"
    echo ""
    echo "Commands:"
    echo "  list    - List available fixture files in legacy storage"
    echo "  copy    - Copy ALL fixture files to new blob storage"
    echo "  latest  - Copy only the LATEST fixture of each type"
    echo ""
    echo "Fixture types managed:"
    echo "  - CCP4i2.json            (CCP4i2 projects, jobs, files)"
    echo "  - RegisterCompounds.json (Compound registry)"
    echo "  - AssayCompounds.json    (Assay data)"
    echo "  - ConstructDatabase.json (Construct/plasmid data)"
    echo "  - auth.json              (Django auth users/groups)"
    echo "  - reversion.json         (Django reversion history)"
    echo ""
    echo "Source: ${SOURCE_STORAGE_ACCOUNT}/${SOURCE_SHARE}/${SOURCE_PATH}"
    echo "Destination: storprv*/${DEST_CONTAINER}/${DEST_PATH}/"
}

# Main command dispatcher
case "$1" in
    list)
        list_fixtures
        ;;
    copy)
        copy_fixtures "all"
        ;;
    latest)
        copy_fixtures "latest"
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
