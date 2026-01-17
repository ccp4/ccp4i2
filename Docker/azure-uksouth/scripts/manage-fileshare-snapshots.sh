#!/bin/bash

# Azure File Share Snapshot Management Script
# Usage: ./manage-fileshare-snapshots.sh <command> [options]
#
# Commands:
#   create [share_name]       Create a snapshot of the file share
#   list [share_name]         List all snapshots for a file share
#   restore <snapshot> [share_name]  Restore files from a snapshot
#   delete <snapshot> [share_name]   Delete a snapshot
#   diff <snapshot> [share_name]     Show files added/modified since snapshot
#
# Default share: ccp4i2-projects

# Ensure Homebrew paths are available
export PATH="/opt/homebrew/bin:$PATH"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Get the directory where the script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="$SCRIPT_DIR/../.env.deployment"

# Load environment variables
if [ -f "$ENV_FILE" ]; then
    source "$ENV_FILE"
else
    echo -e "${RED}Error: .env.deployment not found at $ENV_FILE${NC}"
    echo -e "${YELLOW}Run deploy-infrastructure.sh first.${NC}"
    exit 1
fi

# Default file share name
DEFAULT_SHARE="ccp4i2-projects"

# Error handling
cleanup_on_error() {
    echo -e "${RED}Operation failed. Cleaning up...${NC}"
    restore_network_access
    exit 1
}

trap cleanup_on_error ERR

# Function to get private storage account (for projects/media)
get_storage_account() {
    echo -e "${BLUE}Getting private storage account...${NC}"

    # Query for private storage account (storprv*)
    STORAGE_ACCOUNT_NAME=$(az storage account list \
        --resource-group $RESOURCE_GROUP \
        --query "[?starts_with(name, 'storprv')].name" \
        -o tsv 2>/dev/null | head -1)

    if [ -z "$STORAGE_ACCOUNT_NAME" ]; then
        echo -e "${RED}Error: Private storage account (storprv*) not found${NC}"
        echo -e "${YELLOW}Have you deployed the updated infrastructure.bicep?${NC}"
        exit 1
    fi

    echo -e "${GREEN}Storage account: $STORAGE_ACCOUNT_NAME${NC}"
}

# Function to temporarily enable public network access
enable_network_access() {
    echo -e "${YELLOW}Configuring network access for storage account...${NC}"

    # Get current public IP
    CURRENT_IP=$(curl -s https://ipinfo.io/ip 2>/dev/null || curl -s https://api.ipify.org 2>/dev/null || echo "unknown")

    if [ "$CURRENT_IP" = "unknown" ]; then
        echo -e "${RED}Error: Could not determine public IP address${NC}"
        echo -e "${YELLOW}You may need to manually configure network access${NC}"
        exit 1
    fi

    echo -e "${GREEN}Current public IP: $CURRENT_IP${NC}"

    # Store original network settings
    ORIGINAL_DEFAULT_ACTION=$(az storage account show \
        --name $STORAGE_ACCOUNT_NAME \
        --resource-group $RESOURCE_GROUP \
        --query "networkRuleSet.defaultAction" -o tsv 2>/dev/null || echo "Allow")

    # Add current IP to firewall rules
    echo -e "${YELLOW}Adding IP to storage account firewall...${NC}"

    # Check if IP already exists
    EXISTING_IP=$(az storage account network-rule list \
        --account-name $STORAGE_ACCOUNT_NAME \
        --resource-group $RESOURCE_GROUP \
        --query "ipRules[?contains(ipAddressOrRange, '$CURRENT_IP')].ipAddressOrRange" -o tsv 2>/dev/null)

    if [ -z "$EXISTING_IP" ]; then
        az storage account network-rule add \
            --account-name $STORAGE_ACCOUNT_NAME \
            --resource-group $RESOURCE_GROUP \
            --ip-address $CURRENT_IP \
            --output none 2>/dev/null || true
        IP_ADDED=true
        echo -e "${GREEN}IP address added to firewall${NC}"
    else
        echo -e "${GREEN}IP already in firewall allow list${NC}"
        IP_ADDED=false
    fi

    # Temporarily set default action to Allow for this operation
    echo -e "${YELLOW}Temporarily allowing network access...${NC}"
    az storage account update \
        --name $STORAGE_ACCOUNT_NAME \
        --resource-group $RESOURCE_GROUP \
        --default-action Allow \
        --output none

    # Wait for network rules to propagate
    echo -e "${YELLOW}Waiting for network rules to propagate...${NC}"
    sleep 5

    # Get storage account key
    STORAGE_KEY=$(az storage account keys list \
        --account-name $STORAGE_ACCOUNT_NAME \
        --resource-group $RESOURCE_GROUP \
        --query "[0].value" -o tsv)

    if [ -z "$STORAGE_KEY" ]; then
        echo -e "${RED}Error: Could not retrieve storage account key${NC}"
        restore_network_access
        exit 1
    fi
}

# Function to restore network access settings
restore_network_access() {
    if [ -n "$STORAGE_ACCOUNT_NAME" ] && [ -n "$ORIGINAL_DEFAULT_ACTION" ]; then
        echo -e "${YELLOW}Restoring network security settings...${NC}"
        az storage account update \
            --name $STORAGE_ACCOUNT_NAME \
            --resource-group $RESOURCE_GROUP \
            --default-action $ORIGINAL_DEFAULT_ACTION \
            --output none 2>/dev/null || true
        echo -e "${GREEN}Network settings restored${NC}"
    fi
}

# Function to create a snapshot
create_snapshot() {
    local share_name="${1:-$DEFAULT_SHARE}"

    echo -e "${BLUE}Creating snapshot of file share: $share_name${NC}"

    get_storage_account
    enable_network_access

    # Create snapshot
    SNAPSHOT_TIME=$(az storage share snapshot \
        --name $share_name \
        --account-name $STORAGE_ACCOUNT_NAME \
        --account-key "$STORAGE_KEY" \
        --query "snapshot" -o tsv)

    if [ -n "$SNAPSHOT_TIME" ]; then
        echo -e "${GREEN}Snapshot created successfully!${NC}"
        echo -e "${YELLOW}Snapshot timestamp: $SNAPSHOT_TIME${NC}"
        echo ""
        echo -e "${BLUE}To restore from this snapshot, run:${NC}"
        echo -e "  $0 restore \"$SNAPSHOT_TIME\" $share_name"
        echo ""
        echo -e "${BLUE}To see what changed since this snapshot:${NC}"
        echo -e "  $0 diff \"$SNAPSHOT_TIME\" $share_name"
    else
        echo -e "${RED}Failed to create snapshot${NC}"
    fi

    restore_network_access
}

# Function to list snapshots
list_snapshots() {
    local share_name="${1:-$DEFAULT_SHARE}"

    echo -e "${BLUE}Listing snapshots for file share: $share_name${NC}"

    get_storage_account
    enable_network_access

    # List all shares including snapshots
    echo ""
    echo -e "${YELLOW}Available snapshots:${NC}"
    echo "----------------------------------------"

    az storage share list \
        --account-name $STORAGE_ACCOUNT_NAME \
        --account-key "$STORAGE_KEY" \
        --include-snapshots \
        --query "[?name=='$share_name' && snapshot!=null].{Snapshot:snapshot, Quota:properties.quota}" \
        -o table

    # Count snapshots
    SNAPSHOT_COUNT=$(az storage share list \
        --account-name $STORAGE_ACCOUNT_NAME \
        --account-key "$STORAGE_KEY" \
        --include-snapshots \
        --query "[?name=='$share_name' && snapshot!=null] | length(@)" -o tsv)

    echo ""
    echo -e "${GREEN}Total snapshots: $SNAPSHOT_COUNT${NC}"
    echo -e "${YELLOW}Note: Azure allows up to 200 snapshots per share${NC}"

    restore_network_access
}

# Function to show diff between snapshot and current state
show_diff() {
    local snapshot="$1"
    local share_name="${2:-$DEFAULT_SHARE}"

    if [ -z "$snapshot" ]; then
        echo -e "${RED}Error: Snapshot timestamp required${NC}"
        echo "Usage: $0 diff <snapshot_timestamp> [share_name]"
        exit 1
    fi

    echo -e "${BLUE}Comparing snapshot to current state...${NC}"
    echo -e "${YELLOW}Snapshot: $snapshot${NC}"
    echo -e "${YELLOW}Share: $share_name${NC}"

    get_storage_account
    enable_network_access

    # Create temporary files for comparison
    TEMP_DIR=$(mktemp -d)
    SNAPSHOT_FILES="$TEMP_DIR/snapshot_files.txt"
    CURRENT_FILES="$TEMP_DIR/current_files.txt"

    echo -e "${YELLOW}Listing files in snapshot...${NC}"
    az storage file list \
        --share-name $share_name \
        --account-name $STORAGE_ACCOUNT_NAME \
        --account-key "$STORAGE_KEY" \
        --snapshot "$snapshot" \
        --query "[].name" -o tsv > "$SNAPSHOT_FILES" 2>/dev/null || true

    echo -e "${YELLOW}Listing current files...${NC}"
    az storage file list \
        --share-name $share_name \
        --account-name $STORAGE_ACCOUNT_NAME \
        --account-key "$STORAGE_KEY" \
        --query "[].name" -o tsv > "$CURRENT_FILES" 2>/dev/null || true

    echo ""
    echo -e "${GREEN}=== Files added since snapshot ===${NC}"
    comm -13 <(sort "$SNAPSHOT_FILES") <(sort "$CURRENT_FILES") | head -50

    echo ""
    echo -e "${RED}=== Files removed since snapshot ===${NC}"
    comm -23 <(sort "$SNAPSHOT_FILES") <(sort "$CURRENT_FILES") | head -50

    # Count changes
    ADDED=$(comm -13 <(sort "$SNAPSHOT_FILES") <(sort "$CURRENT_FILES") | wc -l | tr -d ' ')
    REMOVED=$(comm -23 <(sort "$SNAPSHOT_FILES") <(sort "$CURRENT_FILES") | wc -l | tr -d ' ')

    echo ""
    echo -e "${YELLOW}Summary: $ADDED files added, $REMOVED files removed${NC}"
    echo -e "${YELLOW}Note: This shows root-level items only. Use --recursive for deep comparison.${NC}"

    # Cleanup
    rm -rf "$TEMP_DIR"

    restore_network_access
}

# Function to restore from snapshot
restore_snapshot() {
    local snapshot="$1"
    local share_name="${2:-$DEFAULT_SHARE}"

    if [ -z "$snapshot" ]; then
        echo -e "${RED}Error: Snapshot timestamp required${NC}"
        echo "Usage: $0 restore <snapshot_timestamp> [share_name]"
        exit 1
    fi

    echo -e "${BLUE}Restore options from snapshot${NC}"
    echo -e "${YELLOW}Snapshot: $snapshot${NC}"
    echo -e "${YELLOW}Share: $share_name${NC}"
    echo ""

    get_storage_account
    enable_network_access

    echo -e "${YELLOW}Choose restore method:${NC}"
    echo "1) Restore entire share (creates new share from snapshot)"
    echo "2) Restore specific file or directory"
    echo "3) Cancel"
    echo ""
    read -p "Enter choice (1-3): " choice

    case $choice in
        1)
            # Create a new share from the snapshot
            local new_share_name="${share_name}-restored-$(date +%Y%m%d-%H%M%S)"
            echo -e "${YELLOW}Creating new share: $new_share_name${NC}"

            az storage share-rm create \
                --storage-account $STORAGE_ACCOUNT_NAME \
                --resource-group $RESOURCE_GROUP \
                --name $new_share_name \
                --quota 100 \
                --output none

            echo -e "${GREEN}New share created: $new_share_name${NC}"
            echo -e "${YELLOW}Now copying files from snapshot...${NC}"
            echo -e "${RED}Note: For large restores, consider using AzCopy or Azure Portal${NC}"
            echo ""
            echo -e "${BLUE}To copy files manually using AzCopy:${NC}"
            echo "azcopy copy 'https://$STORAGE_ACCOUNT_NAME.file.core.windows.net/$share_name?snapshot=$snapshot&<SAS_TOKEN>' 'https://$STORAGE_ACCOUNT_NAME.file.core.windows.net/$new_share_name?<SAS_TOKEN>' --recursive"
            ;;
        2)
            read -p "Enter file or directory path to restore: " restore_path
            if [ -z "$restore_path" ]; then
                echo -e "${RED}No path specified${NC}"
            else
                echo -e "${YELLOW}Copying $restore_path from snapshot...${NC}"

                # Copy file from snapshot to current share
                az storage file copy start \
                    --account-name $STORAGE_ACCOUNT_NAME \
                    --account-key "$STORAGE_KEY" \
                    --destination-share $share_name \
                    --destination-path "$restore_path.restored" \
                    --source-share $share_name \
                    --source-path "$restore_path" \
                    --source-snapshot "$snapshot"

                echo -e "${GREEN}File copied to: $restore_path.restored${NC}"
                echo -e "${YELLOW}Rename manually if needed${NC}"
            fi
            ;;
        3)
            echo "Cancelled"
            ;;
        *)
            echo -e "${RED}Invalid choice${NC}"
            ;;
    esac

    restore_network_access
}

# Function to delete a snapshot
delete_snapshot() {
    local snapshot="$1"
    local share_name="${2:-$DEFAULT_SHARE}"

    if [ -z "$snapshot" ]; then
        echo -e "${RED}Error: Snapshot timestamp required${NC}"
        echo "Usage: $0 delete <snapshot_timestamp> [share_name]"
        exit 1
    fi

    echo -e "${YELLOW}Deleting snapshot: $snapshot${NC}"
    echo -e "${YELLOW}Share: $share_name${NC}"

    read -p "Are you sure you want to delete this snapshot? (y/N): " confirm
    if [ "$confirm" != "y" ] && [ "$confirm" != "Y" ]; then
        echo "Cancelled"
        exit 0
    fi

    get_storage_account
    enable_network_access

    az storage share delete \
        --name $share_name \
        --account-name $STORAGE_ACCOUNT_NAME \
        --account-key "$STORAGE_KEY" \
        --snapshot "$snapshot"

    echo -e "${GREEN}Snapshot deleted${NC}"

    restore_network_access
}

# Function to show usage
show_usage() {
    echo "Azure File Share Snapshot Management"
    echo ""
    echo "Usage: $0 <command> [options]"
    echo ""
    echo "Commands:"
    echo "  create [share_name]              Create a snapshot (default: $DEFAULT_SHARE)"
    echo "  list [share_name]                List all snapshots"
    echo "  diff <snapshot> [share_name]     Show changes since snapshot"
    echo "  restore <snapshot> [share_name]  Restore from snapshot"
    echo "  delete <snapshot> [share_name]   Delete a snapshot"
    echo ""
    echo "Available file shares:"
    echo "  ccp4i2-projects  - CCP4i2 project data (default)"
    echo "  ccp4data         - CCP4 installation data"
    echo "  staticfiles      - Django static files"
    echo "  mediafiles       - Django media files"
    echo ""
    echo "Examples:"
    echo "  $0 create                        # Snapshot ccp4i2-projects"
    echo "  $0 create ccp4data               # Snapshot ccp4data share"
    echo "  $0 list                          # List ccp4i2-projects snapshots"
    echo "  $0 diff \"2024-01-14T10:30:00.0000000Z\""
    echo ""
}

# Main command handler
case "${1:-}" in
    create)
        create_snapshot "$2"
        ;;
    list)
        list_snapshots "$2"
        ;;
    diff)
        show_diff "$2" "$3"
        ;;
    restore)
        restore_snapshot "$2" "$3"
        ;;
    delete)
        delete_snapshot "$2" "$3"
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
