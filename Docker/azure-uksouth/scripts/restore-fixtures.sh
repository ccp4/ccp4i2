#!/bin/bash

# Restore Fixture Files from Blob Storage to Azure File Share (UK South)
#
# Source: storprv* blob container "django-uploads/fixtures/"
# Destination: storprv* file share "ccp4i2-projects/fixtures/"
#
# This script downloads fixture files from blob storage to the mounted
# file share, making them accessible to Django management commands
# running inside containers.
#
# Use cases:
#   - Initial population of a new deployment from legacy data
#   - Disaster recovery to restore from backups
#
# Environment: Docker/azure-uksouth/.env.deployment
#
# Usage: ./restore-fixtures.sh [list|latest|copy|download <filename>]
#   list              - List available fixtures in blob storage
#   latest            - Download only the latest fixture of each type
#   copy              - Download all fixtures to file share
#   download <file>   - Download a specific fixture file

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

# Source: Blob storage
SOURCE_CONTAINER="django-uploads"
SOURCE_PATH="fixtures"

# Destination: File share (mounted in containers at /mnt/azure-files)
DEST_SHARE="ccp4i2-projects"
DEST_PATH="fixtures"

# Fixture types we're looking for
FIXTURE_TYPES=("CCP4i2" "RegisterCompounds" "AssayCompounds" "ConstructDatabase" "auth" "reversion" "ccp4i2" "registry" "assays" "constructs" "users")

# Get storage account dynamically
get_storage_account() {
    STORAGE_ACCOUNT=$(az storage account list \
        --resource-group $RESOURCE_GROUP \
        --query "[?starts_with(name, 'storprv')].name" \
        -o tsv 2>/dev/null | head -1)

    if [ -z "$STORAGE_ACCOUNT" ]; then
        echo -e "${RED}Error: Private storage account (storprv*) not found${NC}"
        echo -e "${YELLOW}Have you deployed the infrastructure?${NC}"
        exit 1
    fi

    echo -e "${GREEN}Storage account: $STORAGE_ACCOUNT${NC}"
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
    az storage account keys list \
        --account-name "$STORAGE_ACCOUNT" \
        --resource-group "$RESOURCE_GROUP" \
        --query "[0].value" -o tsv
}

# Generate SAS token for blob container
generate_blob_sas() {
    local key=$1
    local permissions=$2  # e.g., 'rl' for read+list

    local expiry=$(date -u -v+24H +"%Y-%m-%dT%H:%MZ" 2>/dev/null || date -u -d "+24 hours" +"%Y-%m-%dT%H:%MZ")

    az storage container generate-sas \
        --name "$SOURCE_CONTAINER" \
        --account-name "$STORAGE_ACCOUNT" \
        --account-key "$key" \
        --permissions "$permissions" \
        --expiry "$expiry" \
        -o tsv 2>/dev/null
}

# Generate SAS token for file share
generate_file_sas() {
    local key=$1
    local permissions=$2  # e.g., 'rwdl' for read+write+delete+list

    local expiry=$(date -u -v+24H +"%Y-%m-%dT%H:%MZ" 2>/dev/null || date -u -d "+24 hours" +"%Y-%m-%dT%H:%MZ")

    az storage share generate-sas \
        --name "$DEST_SHARE" \
        --account-name "$STORAGE_ACCOUNT" \
        --account-key "$key" \
        --permissions "$permissions" \
        --expiry "$expiry" \
        -o tsv 2>/dev/null
}

# List fixture files in blob storage
list_fixtures() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Fixtures Available in Blob Storage${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    get_storage_account

    # Get storage key
    STORAGE_KEY=$(get_storage_key)
    if [ -z "$STORAGE_KEY" ]; then
        echo -e "${RED}Failed to get storage key${NC}"
        exit 1
    fi

    echo -e "${YELLOW}Source: ${STORAGE_ACCOUNT}/${SOURCE_CONTAINER}/${SOURCE_PATH}/${NC}"
    echo ""

    # List all fixture files (filter for dated fixtures starting with '20')
    echo -e "${BLUE}Available fixtures:${NC}"
    echo ""

    # Group by date (most recent first) - only show properly dated fixtures
    az storage blob list \
        --container-name "$SOURCE_CONTAINER" \
        --account-name "$STORAGE_ACCOUNT" \
        --account-key "$STORAGE_KEY" \
        --prefix "${SOURCE_PATH}/20" \
        --query "[].{name:name, size:properties.contentLength, modified:properties.lastModified}" \
        -o table 2>/dev/null | tail -n +3 | sort -r | head -30

    echo ""

    # Show count by type (filter for proper dated fixtures starting with '20')
    echo -e "${YELLOW}Summary by type:${NC}"
    for fixture_type in "${FIXTURE_TYPES[@]}"; do
        COUNT=$(az storage blob list \
            --container-name "$SOURCE_CONTAINER" \
            --account-name "$STORAGE_ACCOUNT" \
            --account-key "$STORAGE_KEY" \
            --prefix "${SOURCE_PATH}/20" \
            --query "length([?contains(name, '-${fixture_type}.json')])" \
            -o tsv 2>/dev/null)
        if [ "$COUNT" != "0" ] && [ -n "$COUNT" ]; then
            echo "  ${fixture_type}: ${COUNT} files"
        fi
    done
}

# Download fixtures to file share
download_fixtures() {
    local mode=$1  # 'all' or 'latest'

    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Downloading Fixtures to File Share${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    check_azcopy
    get_storage_account

    # Get storage key
    STORAGE_KEY=$(get_storage_key)
    if [ -z "$STORAGE_KEY" ]; then
        echo -e "${RED}Failed to get storage key${NC}"
        exit 1
    fi

    # Generate SAS tokens
    echo -e "${YELLOW}Generating SAS tokens...${NC}"
    SOURCE_SAS=$(generate_blob_sas "$STORAGE_KEY" "rl")
    DEST_SAS=$(generate_file_sas "$STORAGE_KEY" "rwdlc")

    if [ -z "$SOURCE_SAS" ] || [ -z "$DEST_SAS" ]; then
        echo -e "${RED}Failed to generate SAS tokens${NC}"
        exit 1
    fi

    echo -e "${GREEN}SAS tokens generated (valid for 24 hours)${NC}"
    echo ""

    # Ensure destination folder exists
    echo -e "${YELLOW}Creating destination folder if needed...${NC}"
    az storage directory create \
        --share-name "$DEST_SHARE" \
        --account-name "$STORAGE_ACCOUNT" \
        --account-key "$STORAGE_KEY" \
        --name "$DEST_PATH" \
        --output none 2>/dev/null || true

    echo -e "${YELLOW}Destination: ${STORAGE_ACCOUNT}/${DEST_SHARE}/${DEST_PATH}/${NC}"
    echo ""

    # Track downloaded filenames for each type (using regular variables for portability)
    DOWNLOADED_CCP4i2=""
    DOWNLOADED_RegisterCompounds=""
    DOWNLOADED_AssayCompounds=""
    DOWNLOADED_ConstructDatabase=""
    DOWNLOADED_auth=""
    DOWNLOADED_reversion=""

    if [ "$mode" = "latest" ]; then
        # Download only the latest fixture of each type
        echo -e "${CYAN}Downloading latest fixtures only...${NC}"
        echo ""

        for fixture_type in "${FIXTURE_TYPES[@]}"; do
            # Get the latest file for this type (filter for dated fixtures starting with '20')
            LATEST_FILE=$(az storage blob list \
                --container-name "$SOURCE_CONTAINER" \
                --account-name "$STORAGE_ACCOUNT" \
                --account-key "$STORAGE_KEY" \
                --prefix "${SOURCE_PATH}/20" \
                --query "[?contains(name, '-${fixture_type}.json')].name" \
                -o tsv 2>/dev/null | sort -r | head -1)

            if [ -n "$LATEST_FILE" ]; then
                FILENAME=$(basename "$LATEST_FILE")
                echo -e "${GREEN}Downloading: ${FILENAME}${NC}"

                SOURCE_URL="https://${STORAGE_ACCOUNT}.blob.core.windows.net/${SOURCE_CONTAINER}/${LATEST_FILE}?${SOURCE_SAS}"
                DEST_URL="https://${STORAGE_ACCOUNT}.file.core.windows.net/${DEST_SHARE}/${DEST_PATH}/${FILENAME}?${DEST_SAS}"

                azcopy copy "$SOURCE_URL" "$DEST_URL" \
                    --skip-version-check \
                    --log-level=ERROR \
                    --output-level=quiet 2>/dev/null

                if [ $? -eq 0 ]; then
                    echo -e "  ${GREEN}Done${NC}"
                    # Store filename in corresponding variable
                    case "$fixture_type" in
                        CCP4i2) DOWNLOADED_CCP4i2="$FILENAME" ;;
                        RegisterCompounds) DOWNLOADED_RegisterCompounds="$FILENAME" ;;
                        AssayCompounds) DOWNLOADED_AssayCompounds="$FILENAME" ;;
                        ConstructDatabase) DOWNLOADED_ConstructDatabase="$FILENAME" ;;
                        auth) DOWNLOADED_auth="$FILENAME" ;;
                        reversion) DOWNLOADED_reversion="$FILENAME" ;;
                    esac
                else
                    echo -e "  ${RED}Failed${NC}"
                fi
            fi
        done
    else
        # Download all fixtures
        echo -e "${CYAN}Downloading all fixture files...${NC}"
        echo ""

        SOURCE_URL="https://${STORAGE_ACCOUNT}.blob.core.windows.net/${SOURCE_CONTAINER}/${SOURCE_PATH}/*?${SOURCE_SAS}"
        DEST_URL="https://${STORAGE_ACCOUNT}.file.core.windows.net/${DEST_SHARE}/${DEST_PATH}/?${DEST_SAS}"

        azcopy copy "$SOURCE_URL" "$DEST_URL" \
            --skip-version-check \
            --log-level=ERROR \
            --recursive 2>/dev/null

        if [ $? -eq 0 ]; then
            echo -e "${GREEN}All fixtures copied successfully${NC}"
        else
            echo -e "${RED}Some files may have failed to copy${NC}"
        fi
    fi

    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}  Download Complete!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
    echo -e "${YELLOW}Fixtures are now available at:${NC}"
    echo "  File share: ${DEST_SHARE}/${DEST_PATH}/"
    echo "  Container path: /mnt/azure-files/${DEST_PATH}/"
    echo ""

    # Use actual filenames if available, otherwise use placeholders
    AUTH_FILE="${DOWNLOADED_auth:-YYYYMMDD-HH-MM-auth.json}"
    REGISTRY_FILE="${DOWNLOADED_RegisterCompounds:-YYYYMMDD-HH-MM-RegisterCompounds.json}"
    ASSAYS_FILE="${DOWNLOADED_AssayCompounds:-YYYYMMDD-HH-MM-AssayCompounds.json}"
    CONSTRUCTS_FILE="${DOWNLOADED_ConstructDatabase:-YYYYMMDD-HH-MM-ConstructDatabase.json}"
    CCP4I2_FILE="${DOWNLOADED_CCP4i2:-YYYYMMDD-HH-MM-CCP4i2.json}"

    echo -e "${YELLOW}To import fixtures, connect to container and run:${NC}"
    echo ""
    echo -e "${CYAN}# 1. Import compounds registry AND assays (includes auth users):${NC}"
    echo "python manage.py import_legacy_compounds \\"
    echo "    --auth-fixture /mnt/azure-files/${DEST_PATH}/${AUTH_FILE} \\"
    echo "    --registry-fixture /mnt/azure-files/${DEST_PATH}/${REGISTRY_FILE} \\"
    echo "    --assays-fixture /mnt/azure-files/${DEST_PATH}/${ASSAYS_FILE}"
    echo ""
    echo -e "${CYAN}# 2. Import construct/plasmid database:${NC}"
    echo "python manage.py import_legacy_constructs \\"
    echo "    --auth-fixture /mnt/azure-files/${DEST_PATH}/${AUTH_FILE} \\"
    echo "    --constructs-fixture /mnt/azure-files/${DEST_PATH}/${CONSTRUCTS_FILE}"
    echo ""
    echo -e "${CYAN}# 3. Import CCP4i2 projects/jobs (optional):${NC}"
    echo "python manage.py import_legacy_ccp4i2 \\"
    echo "    /mnt/azure-files/${DEST_PATH}/${CCP4I2_FILE}"
}

# Download a specific file
download_file() {
    local filename=$1

    echo -e "${CYAN}Downloading: ${filename}${NC}"
    echo ""

    check_azcopy
    get_storage_account

    STORAGE_KEY=$(get_storage_key)
    SOURCE_SAS=$(generate_blob_sas "$STORAGE_KEY" "rl")
    DEST_SAS=$(generate_file_sas "$STORAGE_KEY" "rwdlc")

    # Ensure destination folder exists
    az storage directory create \
        --share-name "$DEST_SHARE" \
        --account-name "$STORAGE_ACCOUNT" \
        --account-key "$STORAGE_KEY" \
        --name "$DEST_PATH" \
        --output none 2>/dev/null || true

    SOURCE_URL="https://${STORAGE_ACCOUNT}.blob.core.windows.net/${SOURCE_CONTAINER}/${SOURCE_PATH}/${filename}?${SOURCE_SAS}"
    DEST_URL="https://${STORAGE_ACCOUNT}.file.core.windows.net/${DEST_SHARE}/${DEST_PATH}/${filename}?${DEST_SAS}"

    azcopy copy "$SOURCE_URL" "$DEST_URL" \
        --skip-version-check \
        --log-level=ERROR 2>/dev/null

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}Downloaded to: ${DEST_SHARE}/${DEST_PATH}/${filename}${NC}"
        echo -e "${GREEN}Container path: /mnt/azure-files/${DEST_PATH}/${filename}${NC}"
    else
        echo -e "${RED}Failed to download ${filename}${NC}"
        exit 1
    fi
}

# Show usage
show_usage() {
    echo "Usage: $0 [list|latest|copy|download <filename>]"
    echo ""
    echo "Commands:"
    echo "  list              - List available fixture files in blob storage"
    echo "  latest            - Download only the LATEST fixture of each type"
    echo "  copy              - Download ALL fixture files to file share"
    echo "  download <file>   - Download a specific fixture file"
    echo ""
    echo "Source: storprv*/${SOURCE_CONTAINER}/${SOURCE_PATH}/"
    echo "Destination: storprv*/${DEST_SHARE}/${DEST_PATH}/"
    echo "Container mount: /mnt/azure-files/${DEST_PATH}/"
    echo ""
    echo "Fixture types:"
    echo "  Legacy: CCP4i2, RegisterCompounds, AssayCompounds, ConstructDatabase, auth, reversion"
    echo "  New:    ccp4i2, registry, assays, constructs, users"
    echo ""
    echo "Example workflow:"
    echo "  1. List available fixtures:"
    echo "     ./restore-fixtures.sh list"
    echo ""
    echo "  2. Download latest fixtures:"
    echo "     ./restore-fixtures.sh latest"
    echo ""
    echo "  3. Connect to container:"
    echo "     az containerapp exec --name ccp4i2-bicep-server --resource-group \$RESOURCE_GROUP"
    echo ""
    echo "  4. Import fixtures (run these in order):"
    echo ""
    echo "     # Compounds registry + assays (includes auth users):"
    echo "     python manage.py import_legacy_compounds \\"
    echo "         --auth-fixture /mnt/azure-files/fixtures/YYYYMMDD-HH-MM-auth.json \\"
    echo "         --registry-fixture /mnt/azure-files/fixtures/YYYYMMDD-HH-MM-RegisterCompounds.json \\"
    echo "         --assays-fixture /mnt/azure-files/fixtures/YYYYMMDD-HH-MM-AssayCompounds.json"
    echo ""
    echo "     # Construct/plasmid database:"
    echo "     python manage.py import_legacy_constructs \\"
    echo "         --auth-fixture /mnt/azure-files/fixtures/YYYYMMDD-HH-MM-auth.json \\"
    echo "         --constructs-fixture /mnt/azure-files/fixtures/YYYYMMDD-HH-MM-ConstructDatabase.json"
    echo ""
    echo "     # CCP4i2 projects/jobs (optional, large file):"
    echo "     python manage.py import_legacy_ccp4i2 \\"
    echo "         /mnt/azure-files/fixtures/YYYYMMDD-HH-MM-CCP4i2.json"
}

# Main command dispatcher
case "$1" in
    list)
        list_fixtures
        ;;
    copy)
        download_fixtures "all"
        ;;
    latest)
        download_fixtures "latest"
        ;;
    download)
        if [ -z "$2" ]; then
            echo -e "${RED}Error: filename required${NC}"
            echo "Usage: $0 download <filename>"
            exit 1
        fi
        download_file "$2"
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
