#!/bin/bash

# View clean, chronological logs from Container Apps Job
# This queries Log Analytics and shows only the log messages without metadata

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Load environment
if [ -f .env.deployment ]; then
    source .env.deployment
else
    echo "❌ .env.deployment not found. Run deploy-infrastructure.sh first."
    exit 1
fi

JOB_NAME="${1:-ccp4i2-bicep-maintenance-job}"
TAIL_LINES="${2:-100}"

echo -e "${BLUE}📋 Fetching clean logs for: $JOB_NAME${NC}"
echo -e "${YELLOW}Last $TAIL_LINES lines...${NC}"
echo ""

# Get the Log Analytics workspace ID
WORKSPACE_ID=$(az containerapp env show \
    --name ccp4i2-bicep-env-ne \
    --resource-group $RESOURCE_GROUP \
    --query "properties.appLogsConfiguration.logAnalyticsConfiguration.customerId" \
    -o tsv)

if [ -z "$WORKSPACE_ID" ]; then
    echo "❌ Could not find Log Analytics workspace"
    exit 1
fi

# Query Log Analytics for clean logs
az monitor log-analytics query \
    --workspace $WORKSPACE_ID \
    --analytics-query "
        ContainerAppConsoleLogs_CL
        | where ContainerAppName_s == '$JOB_NAME'
        | project TimeGenerated, Log_s
        | order by TimeGenerated asc
        | take $TAIL_LINES
    " \
    --query "tables[0].rows[]" \
    -o tsv | while IFS=$'\t' read -r timestamp log; do
        # Format timestamp to be more readable
        formatted_time=$(date -j -f "%Y-%m-%dT%H:%M:%S" "$(echo $timestamp | cut -d'.' -f1)" "+%H:%M:%S" 2>/dev/null || echo "$timestamp")
        echo "[$formatted_time] $log"
    done

echo ""
echo -e "${GREEN}✅ End of logs${NC}"
echo -e "${YELLOW}💡 To follow logs in real-time:${NC}"
echo -e "   az containerapp job logs show -n $JOB_NAME -g $RESOURCE_GROUP --container server --follow --tail 50"
