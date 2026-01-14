#!/bin/bash

# Project Data Backup Script
# Creates a snapshot of the ccp4i2-projects file share for backup/migration purposes
#
# Workflow:
#   1. Copy data to the file share (via your migration process)
#   2. Run this script to create a snapshot
#   3. Continue using the system (new files accrue in DB and filesystem)
#   4. Later: Restore snapshot and merge any new files
#
# Usage: ./backup-projects.sh [create|status|incremental-diff]

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
SNAPSHOT_FILE="$SCRIPT_DIR/../.last-backup-snapshot"

# Load environment variables
if [ -f "$ENV_FILE" ]; then
    source "$ENV_FILE"
else
    echo -e "${RED}Error: .env.deployment not found${NC}"
    exit 1
fi

# File share to backup
SHARE_NAME="ccp4i2-projects"

# Function to get private storage account (for projects/media)
get_storage_info() {
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
}

# Function to enable network access
enable_network_access() {
    CURRENT_IP=$(curl -s https://ipinfo.io/ip 2>/dev/null || curl -s https://api.ipify.org 2>/dev/null)

    ORIGINAL_DEFAULT_ACTION=$(az storage account show \
        --name $STORAGE_ACCOUNT_NAME \
        --resource-group $RESOURCE_GROUP \
        --query "networkRuleSet.defaultAction" -o tsv 2>/dev/null || echo "Allow")

    az storage account network-rule add \
        --account-name $STORAGE_ACCOUNT_NAME \
        --resource-group $RESOURCE_GROUP \
        --ip-address $CURRENT_IP \
        --output none 2>/dev/null || true

    az storage account update \
        --name $STORAGE_ACCOUNT_NAME \
        --resource-group $RESOURCE_GROUP \
        --default-action Allow \
        --output none

    sleep 3

    STORAGE_KEY=$(az storage account keys list \
        --account-name $STORAGE_ACCOUNT_NAME \
        --resource-group $RESOURCE_GROUP \
        --query "[0].value" -o tsv)
}

# Function to restore network settings
restore_network_access() {
    if [ -n "$STORAGE_ACCOUNT_NAME" ] && [ -n "$ORIGINAL_DEFAULT_ACTION" ]; then
        az storage account update \
            --name $STORAGE_ACCOUNT_NAME \
            --resource-group $RESOURCE_GROUP \
            --default-action $ORIGINAL_DEFAULT_ACTION \
            --output none 2>/dev/null || true
    fi
}

trap restore_network_access EXIT

# Function to create backup snapshot
create_backup() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  CCP4i2 Projects Backup${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    get_storage_info
    echo -e "${BLUE}Storage Account: $STORAGE_ACCOUNT_NAME${NC}"
    echo -e "${BLUE}File Share: $SHARE_NAME${NC}"
    echo ""

    enable_network_access

    # Get current share stats
    echo -e "${YELLOW}Getting share statistics...${NC}"
    SHARE_QUOTA=$(az storage share show \
        --name $SHARE_NAME \
        --account-name $STORAGE_ACCOUNT_NAME \
        --account-key "$STORAGE_KEY" \
        --query "properties.quota" -o tsv 2>/dev/null || echo "unknown")

    echo -e "${GREEN}Share quota: ${SHARE_QUOTA}GB${NC}"

    # Count existing snapshots
    SNAPSHOT_COUNT=$(az storage share list \
        --account-name $STORAGE_ACCOUNT_NAME \
        --account-key "$STORAGE_KEY" \
        --include-snapshots \
        --query "[?name=='$SHARE_NAME' && snapshot!=null] | length(@)" -o tsv)

    echo -e "${GREEN}Existing snapshots: $SNAPSHOT_COUNT/200${NC}"
    echo ""

    # Create snapshot
    echo -e "${YELLOW}Creating backup snapshot...${NC}"
    SNAPSHOT_TIME=$(az storage share snapshot \
        --name $SHARE_NAME \
        --account-name $STORAGE_ACCOUNT_NAME \
        --account-key "$STORAGE_KEY" \
        --query "snapshot" -o tsv)

    if [ -n "$SNAPSHOT_TIME" ]; then
        # Save snapshot timestamp for later reference
        echo "$SNAPSHOT_TIME" > "$SNAPSHOT_FILE"

        echo ""
        echo -e "${GREEN}========================================${NC}"
        echo -e "${GREEN}  Backup Created Successfully!${NC}"
        echo -e "${GREEN}========================================${NC}"
        echo ""
        echo -e "${YELLOW}Snapshot: $SNAPSHOT_TIME${NC}"
        echo -e "${YELLOW}Saved to: $SNAPSHOT_FILE${NC}"
        echo ""
        echo -e "${BLUE}Next steps:${NC}"
        echo "1. Continue using the system normally"
        echo "2. To see new files since backup:"
        echo -e "   ${CYAN}$0 incremental-diff${NC}"
        echo "3. To restore this snapshot:"
        echo -e "   ${CYAN}$SCRIPT_DIR/manage-fileshare-snapshots.sh restore \"$SNAPSHOT_TIME\"${NC}"
        echo ""
    else
        echo -e "${RED}Failed to create backup snapshot${NC}"
        exit 1
    fi
}

# Function to show backup status
show_status() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Backup Status${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    get_storage_info
    enable_network_access

    # Check for last backup
    if [ -f "$SNAPSHOT_FILE" ]; then
        LAST_SNAPSHOT=$(cat "$SNAPSHOT_FILE")
        echo -e "${GREEN}Last backup snapshot: $LAST_SNAPSHOT${NC}"

        # Parse timestamp for human-readable format
        # Azure timestamps are like: 2024-01-14T10:30:00.0000000Z
        BACKUP_DATE=$(echo "$LAST_SNAPSHOT" | cut -d'T' -f1)
        BACKUP_TIME=$(echo "$LAST_SNAPSHOT" | cut -d'T' -f2 | cut -d'.' -f1)
        echo -e "${GREEN}Backup date: $BACKUP_DATE $BACKUP_TIME UTC${NC}"
    else
        echo -e "${YELLOW}No backup snapshot recorded${NC}"
        echo -e "${YELLOW}Run '$0 create' to create a backup${NC}"
    fi

    echo ""
    echo -e "${BLUE}All snapshots for $SHARE_NAME:${NC}"
    echo "----------------------------------------"

    az storage share list \
        --account-name $STORAGE_ACCOUNT_NAME \
        --account-key "$STORAGE_KEY" \
        --include-snapshots \
        --query "[?name=='$SHARE_NAME' && snapshot!=null].snapshot" \
        -o tsv | while read snap; do
            if [ "$snap" = "$LAST_SNAPSHOT" ]; then
                echo -e "${GREEN}$snap (current backup)${NC}"
            else
                echo "$snap"
            fi
        done

    SNAPSHOT_COUNT=$(az storage share list \
        --account-name $STORAGE_ACCOUNT_NAME \
        --account-key "$STORAGE_KEY" \
        --include-snapshots \
        --query "[?name=='$SHARE_NAME' && snapshot!=null] | length(@)" -o tsv)

    echo ""
    echo -e "${YELLOW}Total snapshots: $SNAPSHOT_COUNT/200${NC}"
}

# Function to show incremental diff since last backup
show_incremental_diff() {
    if [ ! -f "$SNAPSHOT_FILE" ]; then
        echo -e "${RED}No backup snapshot found${NC}"
        echo -e "${YELLOW}Run '$0 create' first to create a backup${NC}"
        exit 1
    fi

    LAST_SNAPSHOT=$(cat "$SNAPSHOT_FILE")

    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Changes Since Last Backup${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""
    echo -e "${YELLOW}Backup snapshot: $LAST_SNAPSHOT${NC}"
    echo ""

    # Delegate to the main snapshot management script
    "$SCRIPT_DIR/manage-fileshare-snapshots.sh" diff "$LAST_SNAPSHOT" "$SHARE_NAME"
}

# Main command handler
case "${1:-create}" in
    create)
        create_backup
        ;;
    status)
        show_status
        ;;
    incremental-diff|diff)
        show_incremental_diff
        ;;
    help|--help|-h)
        echo "CCP4i2 Projects Backup Script"
        echo ""
        echo "Usage: $0 [command]"
        echo ""
        echo "Commands:"
        echo "  create           Create a backup snapshot (default)"
        echo "  status           Show backup status and list snapshots"
        echo "  incremental-diff Show files changed since last backup"
        echo ""
        echo "Workflow:"
        echo "  1. Migrate/copy data to the file share"
        echo "  2. Run '$0 create' to snapshot current state"
        echo "  3. System continues operating (new files added)"
        echo "  4. Run '$0 incremental-diff' to see new files"
        echo "  5. Restore snapshot if needed and merge new files"
        echo ""
        ;;
    *)
        echo -e "${RED}Unknown command: $1${NC}"
        echo "Run '$0 help' for usage"
        exit 1
        ;;
esac
