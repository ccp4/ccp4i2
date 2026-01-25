# Security Fix Deployment Summary

## What Was Fixed

**Critical Vulnerability:** Frontend-only Teams authorization could be bypassed via browser console

**Attack Vector:**
```javascript
// Any tenant user could bypass the UI and access all data directly:
fetch('/api/ccp4i2/projects/', {
  headers: { 'Authorization': 'Bearer ' + validToken }
})
// Before fix: 200 OK with all projects ✗
// After fix: 403 Forbidden ✓
```

## Changes Made

### 1. Backend Middleware ([azure_auth.py](server/ccp4i2/middleware/azure_auth.py))

**Added:** Teams/Groups authorization check after JWT validation (line ~341)

**Logic:**
- Extracts `groups` claim from validated JWT token
- Compares against `ALLOWED_AZURE_AD_GROUPS` environment variable
- Returns 403 Forbidden if user not in allowed groups
- Handles group overage case (>200 groups)

**Code added:** ~70 lines of authorization logic

### 2. Deployment Configuration ([.env.deployment](Docker/azure-uksouth/.env.deployment))

**Added:**
```bash
ALLOWED_AZURE_AD_GROUPS=6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d
```

This is the Object ID of the "Newcastle Drug Discovery Unit" Microsoft Teams group (matches frontend configuration).

## Deployment Steps

### Step 1: Configure Azure AD (One-Time Setup)

**What:** Enable groups claim in JWT tokens

**How:**
1. Go to [Azure Portal](https://portal.azure.com) → Azure Active Directory → App registrations
2. Find your app: Client ID `386da83f-1bf4-4ad8-b742-79b600e2208b`
3. Click **Token configuration** → **+ Add groups claim**
4. Select **Security groups**, check **Group ID** for ID and Access tokens
5. Click **Add**

**Why:** Puts the user's group memberships into the JWT token so backend can validate

**Details:** See [AZURE_AD_SETUP.md](Docker/azure-uksouth/AZURE_AD_SETUP.md)

### Step 2: Build Server Image

```bash
cd Docker/azure-uksouth
./scripts/build-and-push.sh server
```

**What this does:**
- Creates filtered build context (respecting .dockerignore)
- Uploads to blob storage
- Builds server image in ACR with new middleware code
- Tags as both `latest` and timestamped version
- Updates `IMAGE_TAG_SERVER` in `.env.deployment`

**Time:** ~10-15 minutes

### Step 3: Deploy Server and Worker

```bash
# Deploy server (handles API requests)
./scripts/deploy-applications.sh server

# Deploy worker (handles background jobs)
./scripts/deploy-applications.sh worker
```

**What this does:**
- Deploys new container images to Azure Container Apps
- Sets `ALLOWED_AZURE_AD_GROUPS` environment variable
- Restarts containers with new configuration

**Time:** ~5 minutes total

### Step 4: Verify Deployment

```bash
# Check environment variable is set
az containerapp show \
  --name ccp4i2-bicep-server \
  --resource-group ccp4i2-bicep-rg-uksouth \
  --query "properties.template.containers[0].env[?name=='ALLOWED_AZURE_AD_GROUPS'].value" \
  -o tsv

# Expected output: 6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d
```

### Step 5: Test the Fix

See detailed testing guide: [TEST_SECURITY_FIX.md](Docker/azure-uksouth/TEST_SECURITY_FIX.md)

**Quick test:**
1. Log in as a user NOT in the "Newcastle Drug Discovery Unit" team
2. Should see "Access Denied" on frontend ✓
3. Open browser console, try: `fetch('/api/ccp4i2/projects/')`
4. Should get **403 Forbidden** (not 200 OK) ✓

## What This Fixes ✅

- ✓ **Console bypass attacks** - Non-team members can't call APIs directly
- ✓ **Token reuse attacks** - Valid tenant tokens from non-team users rejected
- ✓ **Frontend/backend parity** - Both layers enforce same authorization
- ✓ **Security in depth** - Backend validates even if frontend is compromised

## What This DOESN'T Fix ❌

- ✗ **Cross-user data access** - Team members can still see each other's projects
- ✗ **Project ownership** - No per-project access control
- ✗ **API authorization** - ViewSets still use `permission_classes=[]`

**For full multi-tenancy:** Implement project ownership (see [SECURITY_VULNERABILITIES.md](Docker/SECURITY_VULNERABILITIES.md#option-2-project-ownership-high-effort-best-security))

## Rollback Plan

If something goes wrong:

**Quick disable (no redeploy):**
```bash
az containerapp update \
  --name ccp4i2-bicep-server \
  --resource-group ccp4i2-bicep-rg-uksouth \
  --remove-env-vars ALLOWED_AZURE_AD_GROUPS
```

**Full rollback:**
```bash
# Revert code
git checkout HEAD -- server/ccp4i2/middleware/azure_auth.py
git checkout HEAD -- Docker/azure-uksouth/.env.deployment

# Rebuild and redeploy
cd Docker/azure-uksouth
./scripts/build-and-push.sh server
./scripts/deploy-applications.sh server worker
```

## Files Changed

| File | Change | Lines |
|------|--------|-------|
| `server/ccp4i2/middleware/azure_auth.py` | Added Teams/Groups authorization | +70 |
| `Docker/azure-uksouth/.env.deployment` | Added `ALLOWED_AZURE_AD_GROUPS` | +6 |

## Documentation Added

| File | Purpose |
|------|---------|
| [SECURITY.md](Docker/SECURITY.md) | Overall security architecture |
| [SECURITY_VULNERABILITIES.md](Docker/SECURITY_VULNERABILITIES.md) | Detailed vulnerability analysis |
| [SECURITY_FIX_TEAMS_BACKEND.md](SECURITY_FIX_TEAMS_BACKEND.md) | Implementation guide (3 options) |
| [SECURITY_QUICK_REFERENCE.md](SECURITY_QUICK_REFERENCE.md) | Quick copy/paste fix guide |
| [AZURE_AD_SETUP.md](Docker/azure-uksouth/AZURE_AD_SETUP.md) | Azure AD configuration steps |
| [TEST_SECURITY_FIX.md](Docker/azure-uksouth/TEST_SECURITY_FIX.md) | Testing procedures |

## Build Script Security Enhancement

Also added security checks to prevent accidental no-auth deployments:

- [build-and-push.sh](Docker/azure-uksouth/scripts/build-and-push.sh) - Refuses to build if `CCP4I2_REQUIRE_AUTH != true`
- [deploy-applications.sh](Docker/azure-uksouth/scripts/deploy-applications.sh) - Refuses to deploy if `CCP4I2_REQUIRE_AUTH != true`

Override for local dev only: `ALLOW_NO_AUTH_BUILD=true`

## Success Metrics

After deployment, verify:

- [ ] Authorized users can use app normally
- [ ] Unauthorized users blocked at frontend (existing)
- [ ] Unauthorized users blocked at backend (NEW - this fix)
- [ ] Console attacks return 403 Forbidden (NEW - this fix)
- [ ] Logs show authorization messages
- [ ] No performance degradation

## Security Posture

### Before Fix

| Attack | Protected? |
|--------|-----------|
| External attacker | ✓ (Wrong tenant) |
| Token from another app | ✓ (Wrong audience) |
| Tenant user, not in team, using UI | ✓ (Frontend blocks) |
| Tenant user, not in team, console attack | ✗ VULNERABLE |

### After Fix

| Attack | Protected? |
|--------|-----------|
| External attacker | ✓ (Wrong tenant) |
| Token from another app | ✓ (Wrong audience) |
| Tenant user, not in team, using UI | ✓ (Frontend blocks) |
| Tenant user, not in team, console attack | ✓ (Backend blocks) |

## Monitoring

After deployment, monitor for:

```bash
# Watch authorization denials
az containerapp logs show \
  --name ccp4i2-bicep-server \
  --resource-group ccp4i2-bicep-rg-uksouth \
  --follow | grep "Access denied"

# Watch successful authorizations
az containerapp logs show \
  --name ccp4i2-bicep-server \
  --resource-group ccp4i2-bicep-rg-uksouth \
  --follow | grep "authorized via Teams"
```

**Red flags:**
- Many "Access denied" messages → Users trying to access who shouldn't
- No authorization messages at all → Groups claim not working

## Next Steps

1. **Deploy this fix** (critical - closes security hole)
2. **Test thoroughly** (see TEST_SECURITY_FIX.md)
3. **Monitor logs** for authorization patterns
4. **Consider project ownership** for true multi-tenancy (future work)
5. **Add automated tests** for authorization logic

## Questions?

- **Azure AD setup:** See [AZURE_AD_SETUP.md](Docker/azure-uksouth/AZURE_AD_SETUP.md)
- **Testing:** See [TEST_SECURITY_FIX.md](Docker/azure-uksouth/TEST_SECURITY_FIX.md)
- **Detailed vulnerability info:** See [SECURITY_VULNERABILITIES.md](Docker/SECURITY_VULNERABILITIES.md)
- **Quick reference:** See [SECURITY_QUICK_REFERENCE.md](SECURITY_QUICK_REFERENCE.md)

## Timeline

**Total time to deploy:** ~30-40 minutes

- Azure AD configuration: 5 minutes (one-time)
- Build server image: 10-15 minutes
- Deploy server + worker: 5-10 minutes
- Verify and test: 10 minutes

## Support

If you encounter issues:

1. Check logs: `az containerapp logs show --name ccp4i2-bicep-server ...`
2. Verify env var: `az containerapp show --name ccp4i2-bicep-server ... --query env`
3. Test token has groups: Browser console → decode JWT token
4. Check Azure AD config: Token configuration has groups claim
5. Rollback if needed: Remove `ALLOWED_AZURE_AD_GROUPS` env var
