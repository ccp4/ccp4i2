#!/bin/bash

# Scheduled Backup Script for Azure Container Apps Job
#
# Creates fixture backups and writes them to the mounted Azure Files share.
# Designed to run as an Azure Container Apps Job on a schedule.
#
# The server container mounts Azure Files at /mnt/azure-files, so backups
# are written directly to: /mnt/azure-files/fixtures/YYYYMMDD-HH-MM-{app}.json
#
# Setup as Container Apps Job:
#   See deploy-backup-job.sh for deployment instructions
#
# Environment variables (inherited from server container):
#   DATABASE_URL - PostgreSQL connection string
#   DJANGO_SETTINGS_MODULE - azure_extensions.settings

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

echo -e "${CYAN}========================================${NC}"
echo -e "${CYAN}  CCP4i2 Scheduled Database Backup${NC}"
echo -e "${CYAN}========================================${NC}"
echo ""
echo "Timestamp: $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
echo ""

# Backup destination (mounted Azure Files)
BACKUP_DIR="/mnt/azure-files/fixtures"

# Verify mount exists
if [ ! -d "/mnt/azure-files" ]; then
    echo -e "${RED}Error: Azure Files not mounted at /mnt/azure-files${NC}"
    echo "This script must run in the Azure Container Apps environment"
    exit 1
fi

# Create fixtures directory if needed
mkdir -p "$BACKUP_DIR"

# Verify database connection
if [ -z "$DATABASE_URL" ]; then
    echo -e "${RED}Error: DATABASE_URL not set${NC}"
    exit 1
fi

echo -e "${GREEN}Backup directory: $BACKUP_DIR${NC}"
echo ""

# Set Django settings
export DJANGO_SETTINGS_MODULE="${DJANGO_SETTINGS_MODULE:-azure_extensions.settings}"

# Navigate to server directory
cd /app/server

# Create timestamp
TIMESTAMP=$(date -u '+%Y%m%d-%H-%M')

echo -e "${YELLOW}Creating backups with timestamp: $TIMESTAMP${NC}"
echo ""

# Apps to backup with their Django app labels
declare -A APPS
APPS["ccp4i2"]="ccp4i2.db"
APPS["registry"]="registry"
APPS["assays"]="assays"
APPS["constructs"]="constructs"
APPS["users"]="users auth"

BACKUP_COUNT=0
TOTAL_SIZE=0

for app_name in "${!APPS[@]}"; do
    app_labels="${APPS[$app_name]}"
    filename="${TIMESTAMP}-${app_name}.json"
    filepath="${BACKUP_DIR}/${filename}"

    echo -e "${CYAN}Backing up: $app_name${NC}"

    # Run dumpdata
    if python manage.py dumpdata $app_labels --indent 2 > "$filepath" 2>/dev/null; then
        size=$(stat -f%z "$filepath" 2>/dev/null || stat -c%s "$filepath" 2>/dev/null || echo "0")
        size_mb=$(echo "scale=2; $size / 1048576" | bc 2>/dev/null || echo "?")
        echo -e "  ${GREEN}✓ Created: $filename (${size_mb} MB)${NC}"
        ((BACKUP_COUNT++))
        TOTAL_SIZE=$((TOTAL_SIZE + size))
    else
        echo -e "  ${YELLOW}⚠ Skipped: $app_name (app not installed or empty)${NC}"
        rm -f "$filepath"  # Remove empty file
    fi
done

echo ""

# Cleanup old backups (keep last 30 days)
echo -e "${YELLOW}Cleaning up old backups (keeping last 30 days)...${NC}"
find "$BACKUP_DIR" -name "*.json" -mtime +30 -delete 2>/dev/null || true

# Count remaining backups
REMAINING=$(ls -1 "$BACKUP_DIR"/*.json 2>/dev/null | wc -l || echo "0")

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Backup Complete${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Backups created: $BACKUP_COUNT"
echo "Total backups on disk: $REMAINING"
echo "Location: $BACKUP_DIR"
echo ""
