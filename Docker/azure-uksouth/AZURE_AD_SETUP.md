# Azure AD Configuration for CCP4i2

## Current Status: ✅ Configured

Groups claim is configured and Teams authorization is working.

## Configuration Summary

**App Registration:** `386da83f-1bf4-4ad8-b742-79b600e2208b`
**Tenant:** `9c5012c9-b616-44c2-a917-66814fbe3e87`
**Authorized Group:** Newcastle Drug Discovery Unit (`6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d`)

### Token Configuration

| Claim | Type | Tokens |
|-------|------|--------|
| groups | Security groups, Group ID | ID, Access |
| idtyp | App-only token signal | Access |

## Adding New Authorized Groups

To grant access to additional Teams/groups:

1. **Find Group ID:**
   - Azure Portal → Microsoft Entra ID → Groups
   - Search for the group name
   - Copy the Object ID

2. **Update Environment:**
   - Edit `Docker/azure-uksouth/.env.deployment`
   - Add to `ALLOWED_AZURE_AD_GROUPS` (comma-separated):
     ```bash
     ALLOWED_AZURE_AD_GROUPS=6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d,new-group-id
     ```

3. **Redeploy:**
   ```bash
   ./scripts/deploy-applications.sh server
   ./scripts/deploy-applications.sh worker
   ```

## Troubleshooting

### Users Not Authorized Despite Being in Team

1. **User must log out and log back in** to get a new token with groups claim
2. Verify group ID is correct in `.env.deployment`
3. Check Azure AD Token Configuration has groups claim

### Groups Claim Not in Token

If `groups` array is missing from JWT:

1. Azure Portal → App Registrations → Your App → Token Configuration
2. Verify "groups" claim is listed
3. If missing, add: Groups claim → Security groups → Group ID

### Group Overage (>200 Groups)

If a user has >200 group memberships:
- Azure AD uses overage indicator instead of full list
- Middleware returns 403 with specific message
- Solution: Create dedicated smaller access group

## Verifying Configuration

Check deployed environment:
```bash
az containerapp show \
  --name ccp4i2-bicep-server \
  --resource-group ccp4i2-bicep-rg-uksouth \
  --query "properties.template.containers[0].env[?name=='ALLOWED_AZURE_AD_GROUPS'].value" \
  -o tsv
```

Check logs:
```bash
az containerapp logs show \
  --name ccp4i2-bicep-server \
  --resource-group ccp4i2-bicep-rg-uksouth \
  --follow | grep -E "(authorized|denied)"
```
