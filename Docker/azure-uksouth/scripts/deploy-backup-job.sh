#!/bin/bash

# Deploy Scheduled Backup Job to Azure Container Apps
#
# Creates a Container Apps Job that runs daily to backup the database
# to the mounted Azure Files share.
#
# Usage: ./deploy-backup-job.sh [create|update|delete|trigger|logs]
#
# Environment: Docker/azure-uksouth/.env.deployment

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

# Get script directory and load environment
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="$SCRIPT_DIR/../.env.deployment"

if [ -f "$ENV_FILE" ]; then
    source "$ENV_FILE"
else
    echo -e "${RED}Error: .env.deployment not found${NC}"
    exit 1
fi

# Job configuration
JOB_NAME="ccp4i2-backup-job"
SCHEDULE="0 2 * * *"  # Daily at 2am UTC

# Get the Container Apps environment name
get_environment() {
    CONTAINER_ENV=$(az containerapp env list \
        --resource-group "$RESOURCE_GROUP" \
        --query "[0].name" -o tsv 2>/dev/null)

    if [ -z "$CONTAINER_ENV" ]; then
        echo -e "${RED}Error: No Container Apps environment found${NC}"
        exit 1
    fi
    echo -e "${GREEN}Environment: $CONTAINER_ENV${NC}"
}

# Get the latest server image
get_server_image() {
    SERVER_IMAGE=$(az containerapp show \
        --name ccp4i2-bicep-server \
        --resource-group "$RESOURCE_GROUP" \
        --query "properties.template.containers[0].image" -o tsv 2>/dev/null)

    if [ -z "$SERVER_IMAGE" ]; then
        echo -e "${RED}Error: Could not get server image${NC}"
        exit 1
    fi
    echo -e "${GREEN}Server image: $SERVER_IMAGE${NC}"
}

# Get storage account for file share mount
get_storage_account() {
    STORAGE_ACCOUNT=$(az storage account list \
        --resource-group "$RESOURCE_GROUP" \
        --query "[?starts_with(name, 'storprv')].name" -o tsv 2>/dev/null | head -1)

    if [ -z "$STORAGE_ACCOUNT" ]; then
        echo -e "${RED}Error: Storage account not found${NC}"
        exit 1
    fi
    echo -e "${GREEN}Storage account: $STORAGE_ACCOUNT${NC}"
}

# Create the backup job
create_job() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Creating Backup Job${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    get_environment
    get_server_image
    get_storage_account
    echo ""

    echo -e "${YELLOW}Creating Container Apps Job: $JOB_NAME${NC}"
    echo "Schedule: $SCHEDULE (daily at 2am UTC)"
    echo ""

    # Get environment variables from server container
    SERVER_ENVS=$(az containerapp show \
        --name ccp4i2-bicep-server \
        --resource-group "$RESOURCE_GROUP" \
        --query "properties.template.containers[0].env" -o json)

    # Create the job
    az containerapp job create \
        --name "$JOB_NAME" \
        --resource-group "$RESOURCE_GROUP" \
        --environment "$CONTAINER_ENV" \
        --trigger-type "Schedule" \
        --cron-expression "$SCHEDULE" \
        --replica-timeout 1800 \
        --replica-retry-limit 1 \
        --replica-completion-count 1 \
        --parallelism 1 \
        --image "$SERVER_IMAGE" \
        --cpu 0.5 \
        --memory 1Gi \
        --registry-server "${ACR_NAME}.azurecr.io" \
        --env-vars "DJANGO_SETTINGS_MODULE=azure_extensions.settings" \
        --command "/bin/bash" "-c" "/app/Docker/azure-uksouth/scripts/run-scheduled-backup.sh"

    echo ""
    echo -e "${GREEN}Job created successfully!${NC}"
    echo ""
    echo -e "${YELLOW}Note: You need to configure environment variables and storage mounts${NC}"
    echo "manually in the Azure Portal or via additional az commands."
    echo ""
    echo "Required environment variables:"
    echo "  - DATABASE_URL"
    echo "  - DJANGO_SETTINGS_MODULE=azure_extensions.settings"
    echo ""
    echo "Required storage mount:"
    echo "  - Mount Azure Files share to /mnt/azure-files"
}

# Update the job with latest server image
update_job() {
    echo -e "${CYAN}Updating backup job with latest server image...${NC}"

    get_server_image

    az containerapp job update \
        --name "$JOB_NAME" \
        --resource-group "$RESOURCE_GROUP" \
        --image "$SERVER_IMAGE"

    echo -e "${GREEN}Job updated!${NC}"
}

# Delete the job
delete_job() {
    echo -e "${YELLOW}Deleting backup job: $JOB_NAME${NC}"

    az containerapp job delete \
        --name "$JOB_NAME" \
        --resource-group "$RESOURCE_GROUP" \
        --yes

    echo -e "${GREEN}Job deleted!${NC}"
}

# Trigger the job manually
trigger_job() {
    echo -e "${CYAN}Triggering backup job manually...${NC}"

    az containerapp job start \
        --name "$JOB_NAME" \
        --resource-group "$RESOURCE_GROUP"

    echo -e "${GREEN}Job triggered! Use '$0 logs' to view progress.${NC}"
}

# View job execution logs
view_logs() {
    echo -e "${CYAN}Fetching recent job executions...${NC}"

    # List recent executions
    az containerapp job execution list \
        --name "$JOB_NAME" \
        --resource-group "$RESOURCE_GROUP" \
        --query "[].{Name:name, Status:properties.status, StartTime:properties.startTime}" \
        -o table

    echo ""
    echo -e "${YELLOW}To view logs for a specific execution:${NC}"
    echo "az containerapp job logs show --name $JOB_NAME --resource-group $RESOURCE_GROUP --execution <execution-name>"
}

# Show job status
show_status() {
    echo -e "${CYAN}Backup Job Status${NC}"
    echo "=" "60"

    az containerapp job show \
        --name "$JOB_NAME" \
        --resource-group "$RESOURCE_GROUP" \
        --query "{Name:name, Status:properties.provisioningState, Schedule:properties.configuration.scheduleTriggerConfig.cronExpression, Image:properties.template.containers[0].image}" \
        -o table 2>/dev/null || echo -e "${YELLOW}Job not found${NC}"
}

# Show usage
show_usage() {
    echo "Deploy Scheduled Backup Job"
    echo ""
    echo "Usage: $0 <command>"
    echo ""
    echo "Commands:"
    echo "  create    Create the backup job"
    echo "  update    Update job with latest server image"
    echo "  delete    Delete the backup job"
    echo "  trigger   Trigger the job manually"
    echo "  logs      View recent job executions"
    echo "  status    Show job status"
    echo ""
    echo "The job runs daily at 2am UTC and creates backups to:"
    echo "  /mnt/azure-files/fixtures/YYYYMMDD-HH-MM-{app}.json"
}

# Main
case "${1:-}" in
    create)
        create_job
        ;;
    update)
        update_job
        ;;
    delete)
        delete_job
        ;;
    trigger)
        trigger_job
        ;;
    logs)
        view_logs
        ;;
    status)
        show_status
        ;;
    *)
        show_usage
        exit 1
        ;;
esac
