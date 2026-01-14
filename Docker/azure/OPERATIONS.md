# Azure Deployment Operations Guide

This guide covers day-to-day operations for the CCP4i2 Azure deployment.

## Quick Reference

```bash
# Always source environment first
source ./Docker/azure/.env.deployment

# Check deployment status
az containerapp show --name ccp4i2-bicep-server --resource-group "$RESOURCE_GROUP" \
  --query "{image: properties.template.containers[0].image, state: properties.runningState}" -o json

# View logs
az containerapp logs show --name ccp4i2-bicep-server --resource-group "$RESOURCE_GROUP" --type console --tail 100

# Build and deploy
./Docker/azure/scripts/build-and-push.sh server
./Docker/azure/scripts/deploy-applications.sh server
```

## Container Apps Overview

| App Name | Purpose | Scaling |
|----------|---------|---------|
| `ccp4i2-bicep-server` | Django REST API | 1-10 replicas |
| `ccp4i2-bicep-web` | Next.js frontend | 1-5 replicas |
| `ccp4i2-bicep-worker` | Background job processor | 0-5 replicas (queue-based) |

## Deployment Scripts

All scripts are in `Docker/azure/scripts/`. Source `.env.deployment` before running.

### Build and Push Images

```bash
# Build and push a single image
./Docker/azure/scripts/build-and-push.sh server   # Django backend
./Docker/azure/scripts/build-and-push.sh web      # Next.js frontend
./Docker/azure/scripts/build-and-push.sh worker   # Same as server (different entrypoint)
```

### Deploy Applications

```bash
# Deploy single application
./Docker/azure/scripts/deploy-applications.sh server
./Docker/azure/scripts/deploy-applications.sh web

# Full deployment (infrastructure + applications)
./Docker/azure/scripts/deploy-complete.sh
```

### Other Scripts

```bash
# Infrastructure management
./Docker/azure/scripts/deploy-infrastructure.sh   # Core Azure resources
./Docker/azure/scripts/deploy-servicebus.sh       # Service Bus queue
./Docker/azure/scripts/deploy-management.sh       # Management VM

# Management access
./Docker/azure/scripts/connect-management.sh      # SSH to management VM

# Maintenance jobs
./Docker/azure/scripts/run-maintenance-job.sh     # Run maintenance
./Docker/azure/scripts/view-job-logs.sh           # View job output
```

## Environment Variables

### Viewing Current Environment

```bash
# List all env vars for server app
az containerapp show --name ccp4i2-bicep-server --resource-group "$RESOURCE_GROUP" \
  --query "properties.template.containers[0].env[].{name:name, value:value}" -o table
```

### Adding/Updating Environment Variables

```bash
# Add or update a single env var
az containerapp update \
  --name ccp4i2-bicep-server \
  --resource-group "$RESOURCE_GROUP" \
  --set-env-vars MY_VAR=my_value

# Add multiple env vars
az containerapp update \
  --name ccp4i2-bicep-server \
  --resource-group "$RESOURCE_GROUP" \
  --set-env-vars VAR1=value1 VAR2=value2
```

### Removing Environment Variables

```bash
az containerapp update \
  --name ccp4i2-bicep-server \
  --resource-group "$RESOURCE_GROUP" \
  --remove-env-vars MY_VAR
```

### Permanent Environment Variables

For permanent env vars, add them to `infrastructure/applications.bicep` in the `env` array:

```bicep
env: [
  // ... existing vars ...
  {
    name: 'MY_NEW_VAR'
    value: 'my_value'
  }
]
```

Then redeploy: `./Docker/azure/scripts/deploy-applications.sh server`

## Admin API Endpoints

### Import Status

Check current database counts:

```bash
curl -X GET "https://$WEB_URL/api/compounds/admin/import-status/" \
  -H "Authorization: Bearer $TOKEN"
```

### Import Legacy Fixtures

Import legacy CCP4i2 fixtures:

```bash
curl -X POST "https://$WEB_URL/api/compounds/admin/import-legacy/" \
  -H "Authorization: Bearer $TOKEN" \
  -F "users_fixture=@auth.json" \
  -F "registry_fixture=@RegisterCompounds.json" \
  -F "assays_fixture=@AssayCompounds.json"

# Dry run (validate without importing)
curl -X POST "https://$WEB_URL/api/compounds/admin/import-legacy/" \
  -H "Authorization: Bearer $TOKEN" \
  -F "registry_fixture=@RegisterCompounds.json" \
  -F "dry_run=true"
```

### Reset Compounds Data (Migration Phase Only)

**WARNING: This permanently deletes all compounds data!**

This endpoint is gated by an environment variable and requires explicit confirmation.

#### Enable Reset (Temporary)

```bash
# Enable the reset endpoint
az containerapp update \
  --name ccp4i2-bicep-server \
  --resource-group "$RESOURCE_GROUP" \
  --set-env-vars CCP4I2_ALLOW_DB_RESET=true
```

#### Perform Reset

```bash
# Reset compounds data only
curl -X POST "https://$WEB_URL/api/compounds/admin/reset-data/" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"confirmation": "DELETE ALL COMPOUNDS DATA"}'

# Reset compounds data AND clear legacy user mappings
curl -X POST "https://$WEB_URL/api/compounds/admin/reset-data/" \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"confirmation": "DELETE ALL COMPOUNDS DATA", "include_users": "true"}'
```

#### Disable Reset (After Migration)

```bash
# Remove the env var to disable reset endpoint
az containerapp update \
  --name ccp4i2-bicep-server \
  --resource-group "$RESOURCE_GROUP" \
  --remove-env-vars CCP4I2_ALLOW_DB_RESET
```

## Monitoring

### Check Revision Health

```bash
# Get latest revision name
REVISION=$(az containerapp show --name ccp4i2-bicep-server --resource-group "$RESOURCE_GROUP" \
  --query "properties.latestRevisionName" -o tsv)

# Check revision status
az containerapp revision show \
  --name ccp4i2-bicep-server \
  --resource-group "$RESOURCE_GROUP" \
  --revision "$REVISION" \
  --query "{active: properties.active, replicas: properties.replicas, healthState: properties.healthState, runningState: properties.runningState}"
```

### View Logs

```bash
# Console logs (stdout/stderr)
az containerapp logs show \
  --name ccp4i2-bicep-server \
  --resource-group "$RESOURCE_GROUP" \
  --type console \
  --tail 100

# System logs (container events)
az containerapp logs show \
  --name ccp4i2-bicep-server \
  --resource-group "$RESOURCE_GROUP" \
  --type system \
  --tail 50

# Follow logs in real-time
az containerapp logs show \
  --name ccp4i2-bicep-server \
  --resource-group "$RESOURCE_GROUP" \
  --type console \
  --follow
```

### Check Image Tags

```bash
# List recent image tags
az acr repository show-tags \
  --name "$ACR_NAME" \
  --repository ccp4i2/server \
  --orderby time_desc \
  --top 10

# Check currently deployed image
az containerapp show --name ccp4i2-bicep-server --resource-group "$RESOURCE_GROUP" \
  --query "properties.template.containers[0].image" -o tsv
```

## Troubleshooting

### Revision Failed to Start

1. Check logs for errors:
   ```bash
   az containerapp logs show --name ccp4i2-bicep-server --resource-group "$RESOURCE_GROUP" --type console --tail 200
   ```

2. Common issues:
   - **Migration error**: Check for Django migration issues in logs
   - **Import error**: Python import failures (check module paths)
   - **Database connection**: Check DB_HOST and credentials
   - **Volume mount**: Check Azure Files connectivity

### Rollback to Previous Revision

```bash
# List revisions
az containerapp revision list \
  --name ccp4i2-bicep-server \
  --resource-group "$RESOURCE_GROUP" \
  --query "[].{name:name, active:properties.active, created:properties.createdTime}" \
  -o table

# Activate a previous revision (replace with actual revision name)
az containerapp revision activate \
  --name ccp4i2-bicep-server \
  --resource-group "$RESOURCE_GROUP" \
  --revision ccp4i2-bicep-server--0000081

# Route traffic to it
az containerapp ingress traffic set \
  --name ccp4i2-bicep-server \
  --resource-group "$RESOURCE_GROUP" \
  --revision-weight ccp4i2-bicep-server--0000081=100
```

### Database Access

Connect to PostgreSQL via management VM:

```bash
# SSH to management VM
./Docker/azure/scripts/connect-management.sh

# On the VM, connect to PostgreSQL
psql "host=$DB_HOST dbname=postgres user=ccp4i2 password=$DB_PASSWORD sslmode=require"
```

### Force Restart

```bash
# Restart by updating a dummy env var
az containerapp update \
  --name ccp4i2-bicep-server \
  --resource-group "$RESOURCE_GROUP" \
  --set-env-vars RESTART_TIMESTAMP=$(date +%s)
```

## Data Migration Workflow

Typical workflow for importing legacy data:

1. **Enable reset endpoint** (if needed):
   ```bash
   az containerapp update --name ccp4i2-bicep-server --resource-group "$RESOURCE_GROUP" \
     --set-env-vars CCP4I2_ALLOW_DB_RESET=true
   ```

2. **Check current state**:
   ```bash
   curl -X GET "https://$WEB_URL/api/compounds/admin/import-status/" -H "Authorization: Bearer $TOKEN"
   ```

3. **Reset if needed**:
   ```bash
   curl -X POST "https://$WEB_URL/api/compounds/admin/reset-data/" \
     -H "Authorization: Bearer $TOKEN" \
     -H "Content-Type: application/json" \
     -d '{"confirmation": "DELETE ALL COMPOUNDS DATA", "include_users": "true"}'
   ```

4. **Import fixtures**:
   ```bash
   # Import users first
   curl -X POST "https://$WEB_URL/api/compounds/admin/import-legacy/" \
     -H "Authorization: Bearer $TOKEN" \
     -F "users_fixture=@auth.json"

   # Then import registry and assays
   curl -X POST "https://$WEB_URL/api/compounds/admin/import-legacy/" \
     -H "Authorization: Bearer $TOKEN" \
     -F "registry_fixture=@RegisterCompounds.json" \
     -F "assays_fixture=@AssayCompounds.json"
   ```

5. **Verify**:
   ```bash
   curl -X GET "https://$WEB_URL/api/compounds/admin/import-status/" -H "Authorization: Bearer $TOKEN"
   ```

6. **Disable reset when done**:
   ```bash
   az containerapp update --name ccp4i2-bicep-server --resource-group "$RESOURCE_GROUP" \
     --remove-env-vars CCP4I2_ALLOW_DB_RESET
   ```

## Security Checklist

Before going to production:

- [ ] Remove `CCP4I2_ALLOW_DB_RESET` environment variable
- [ ] Verify `DEBUG=false` is set
- [ ] Confirm `PLATFORM_ADMIN_EMAILS` is correctly configured
- [ ] Check CORS origins are restricted to actual frontend domain
- [ ] Ensure Key Vault access policies are correct
- [ ] Review Azure AD app registration permissions
