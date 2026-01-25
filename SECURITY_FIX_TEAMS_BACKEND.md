# Teams Authorization - Backend Implementation

## Status: ✅ IMPLEMENTED

This fix has been implemented and deployed. Teams/Groups authorization is now enforced on the backend.

## What Was Implemented

### Backend Middleware (`server/ccp4i2/middleware/azure_auth.py`)

The Django middleware now validates the `groups` claim in JWT tokens:

1. Extracts `groups` array from validated JWT token
2. Compares against `ALLOWED_AZURE_AD_GROUPS` environment variable
3. Returns 403 Forbidden if user is not in an authorized group
4. Handles group overage case (>200 groups)

### Configuration

**Azure AD App Registration:**
- Token Configuration → Groups claim → Security groups → Group ID

**Environment Variable:**
```bash
ALLOWED_AZURE_AD_GROUPS=6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d
```

## Security Posture

### Before Fix (Vulnerable)

| Attack | Protected? |
|--------|-----------|
| External attacker | ✓ (Wrong tenant) |
| Token from another app | ✓ (Wrong audience) |
| Tenant user, not in team, using UI | ✓ (Frontend blocks) |
| Tenant user, not in team, console attack | ✗ VULNERABLE |

### After Fix (Current)

| Attack | Protected? |
|--------|-----------|
| External attacker | ✓ (Wrong tenant) |
| Token from another app | ✓ (Wrong audience) |
| Tenant user, not in team, using UI | ✓ (Frontend blocks) |
| Tenant user, not in team, console attack | ✓ (Backend blocks) |

## Testing

To verify the fix is working:

1. **Authorized users** (in Newcastle Drug Discovery Unit team) can use the app normally
2. **Unauthorized users** (not in team) get 403 Forbidden on API calls
3. **Console bypass attempts** return 403 instead of 200

Check logs for:
- ✅ `✅ User aad_xxx authorized via Teams/Groups membership`
- ❌ `Access denied for user aad_xxx - not in authorized groups`

## Adding Additional Groups

To allow additional teams/groups access:

1. Find the Group Object ID in Azure Portal → Entra ID → Groups
2. Add to `ALLOWED_AZURE_AD_GROUPS` as comma-separated list:
   ```bash
   ALLOWED_AZURE_AD_GROUPS=6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d,another-group-id
   ```
3. Redeploy server and worker

## Remaining Considerations

While Teams authorization is now enforced, note that:

- All authorized team members can see all projects (no per-project ownership)
- For true multi-tenancy, consider implementing project ownership in the future

## Files Changed

- `server/ccp4i2/middleware/azure_auth.py` - Teams authorization logic
- `Docker/azure-uksouth/.env.deployment` - `ALLOWED_AZURE_AD_GROUPS` config
- `Docker/azure-uksouth/scripts/build-and-push.sh` - Auth requirement check
- `Docker/azure-uksouth/scripts/deploy-applications.sh` - Auth requirement check
