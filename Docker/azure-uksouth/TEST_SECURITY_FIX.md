# Testing the Teams Authorization Security Fix

## What Was Fixed

**Vulnerability:** Frontend-only Teams authorization could be bypassed via browser console
**Fix:** Added backend validation of `groups` claim in JWT tokens

## Pre-Deployment Testing (Local)

If you want to test locally before deploying to Azure:

### 1. Test Middleware Logic

Create a test file to verify the middleware logic:

```python
# tests/test_teams_auth.py
import os
import json
import pytest
from django.http import HttpRequest, JsonResponse
from server.ccp4i2.middleware.azure_auth import AzureADAuthMiddleware

def test_groups_authorization_blocks_non_members():
    """Test that users not in allowed groups are blocked"""
    os.environ['ALLOWED_AZURE_AD_GROUPS'] = '6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d'
    os.environ['CCP4I2_REQUIRE_AUTH'] = 'true'

    # Simulate a JWT with user not in the allowed group
    test_claims = {
        'sub': 'test-user-id-12345',
        'email': 'test@example.com',
        'groups': ['other-group-id-99999']  # Wrong group
    }

    # The middleware should return 403 Forbidden
    # (actual test implementation would mock JWT validation)

def test_groups_authorization_allows_members():
    """Test that users in allowed groups are granted access"""
    os.environ['ALLOWED_AZURE_AD_GROUPS'] = '6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d'

    test_claims = {
        'sub': 'test-user-id-12345',
        'email': 'test@example.com',
        'groups': ['6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d']  # Correct group
    }

    # Should create user and allow access
```

Run with:
```bash
cd server
ccp4-python -m pytest tests/test_teams_auth.py -v
```

## Post-Deployment Testing (Production)

After deploying to Azure, perform these tests:

### Test 1: Verify Environment Variable

Check that the group ID is set correctly:

```bash
az containerapp show \
  --name ccp4i2-bicep-server \
  --resource-group ccp4i2-bicep-rg-uksouth \
  --query "properties.template.containers[0].env[?name=='ALLOWED_AZURE_AD_GROUPS']" \
  -o json
```

**Expected output:**
```json
[
  {
    "name": "ALLOWED_AZURE_AD_GROUPS",
    "value": "6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d"
  }
]
```

### Test 2: Authorized User Access (Should Work)

**Prerequisites:**
- Test with a user who IS in the "Newcastle Drug Discovery Unit" team

**Steps:**
1. Visit your application: `https://ccp4i2-bicep-web.uksouth.azurecontainerapps.io`
2. Log in with an authorized user account
3. Should see the application UI normally ✓
4. Open browser console and run:

```javascript
// This should work (200 OK) because user is in the team
fetch('/api/ccp4i2/health/')
  .then(r => r.json())
  .then(console.log);

// This should also work
fetch('/api/ccp4i2/projects/')
  .then(r => r.json())
  .then(console.log);
```

**Expected result:**
- ✓ Both requests return 200 OK
- ✓ Data is returned normally
- ✓ Application works as expected

### Test 3: Unauthorized User Blocked (Critical Test)

**Prerequisites:**
- Test with a user who is in your Azure AD tenant BUT NOT in the team
- Or create a test user specifically for this

**Steps:**
1. Visit your application with the unauthorized user account
2. Try to log in
3. Frontend should show "Access Denied" page ✓
4. Open browser console and try the console bypass attack:

```javascript
// Attempt to bypass frontend authorization
fetch('/api/ccp4i2/projects/')
  .then(r => {
    console.log('Status:', r.status);
    return r.json();
  })
  .then(console.log)
  .catch(console.error);
```

**Expected result (AFTER fix):**
```javascript
// Status: 403
// {
//   "success": false,
//   "error": "Access denied: You are not a member of an authorized group.
//             This application requires membership in the Newcastle Drug Discovery Unit team.
//             Please contact your administrator to request access."
// }
```

**Before fix (vulnerable):**
```javascript
// Status: 200  ✗ BAD
// [ { id: 1, name: "Project A", ... }, ... ]  ✗ LEAKED DATA
```

### Test 4: Direct API Call with cURL

Test the API directly (simulating an attacker with a valid token):

```bash
# Get an access token for a user NOT in the team
# (Use Azure CLI or MSAL to obtain token)
TOKEN="eyJ0eXAiOiJKV1QiLCJhbGc..."  # Unauthorized user's token

# Try to access the API
curl -H "Authorization: Bearer $TOKEN" \
     https://ccp4i2-bicep-web.uksouth.azurecontainerapps.io/api/ccp4i2/projects/

# Expected: 403 Forbidden
# {
#   "success": false,
#   "error": "Access denied: You are not a member of an authorized group..."
# }
```

### Test 5: Check Application Logs

View real-time logs to see authorization checks:

```bash
# Server logs
az containerapp logs show \
  --name ccp4i2-bicep-server \
  --resource-group ccp4i2-bicep-rg-uksouth \
  --follow

# Look for these log messages:

# ✓ Authorized user:
# [INFO] ✅ User aad_12345678 authorized via Teams/Groups membership

# ✗ Unauthorized user:
# [WARNING] Access denied for user aad_87654321 - not in authorized groups.
#           User groups: ['other-group-id'], Required: ['6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d']
```

### Test 6: Token Inspection

Verify that authorized users have the correct group claim in their tokens:

**In browser console (after login as authorized user):**

```javascript
// Get the current MSAL account
const accounts = window.msal?.getAllAccounts();
if (!accounts || accounts.length === 0) {
  console.error('No MSAL accounts found');
} else {
  const account = accounts[0];
  console.log('Account:', account.username);

  // Get ID token
  const response = await window.msal.acquireTokenSilent({
    scopes: ['openid', 'profile'],
    account: account
  });

  // Decode JWT (middle part is the payload)
  const parts = response.idToken.split('.');
  const payload = JSON.parse(atob(parts[1]));

  console.log('Full JWT claims:', payload);
  console.log('Groups claim:', payload.groups);

  // Verify Newcastle Drug Discovery Unit group is present
  const expectedGroup = '6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d';
  if (payload.groups && payload.groups.includes(expectedGroup)) {
    console.log('✅ User has required group membership in token');
  } else {
    console.error('❌ Required group NOT in token!');
    console.log('This means Azure AD is not emitting groups claim correctly');
  }
}
```

**Expected for authorized users:**
```javascript
Groups claim: [
  "6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d",  // Newcastle Drug Discovery Unit
  "other-group-id-1",
  "other-group-id-2"
]
✅ User has required group membership in token
```

**Expected for unauthorized users:**
```javascript
Groups claim: [
  "other-group-id-1",
  "other-group-id-2"
]
❌ Required group NOT in token!
```

## Test Results Checklist

- [ ] Environment variable `ALLOWED_AZURE_AD_GROUPS` is set in container app
- [ ] Azure AD app registration has "groups claim" configured
- [ ] Authorized users can log in and use the application normally
- [ ] Authorized users' tokens contain the correct group ID
- [ ] Unauthorized users see "Access Denied" on frontend
- [ ] Unauthorized users get 403 Forbidden when trying console bypass
- [ ] Logs show authorization messages for both authorized/unauthorized users
- [ ] Direct API calls (cURL) with unauthorized tokens return 403

## Common Issues

### Issue 1: Authorized users being blocked

**Symptoms:**
- User is in the team but getting 403 Forbidden
- Logs show: `User groups: [], Required: [6f35cbeb...]`

**Cause:** Azure AD not emitting groups claim

**Fix:**
1. Check Azure AD app registration Token Configuration
2. Ensure "groups claim" is added
3. User must log out and log back in (to get new token with groups)

### Issue 2: No groups in token

**Symptoms:**
- Token decoded but `groups` field is empty or missing

**Possible causes:**
1. Azure AD Free (not Premium P1/P2) - groups claim requires Premium
2. Token Configuration not saved properly
3. User needs to re-authenticate

**Fix:**
1. Verify Azure AD license level
2. Re-configure Token Configuration
3. Clear browser cache and re-login

### Issue 3: All users blocked

**Symptoms:**
- Even authorized users get 403 Forbidden
- Logs show: `Groups authorization not configured`

**Cause:** `ALLOWED_AZURE_AD_GROUPS` env var not set or not deployed

**Fix:**
```bash
# Check if env var is set
az containerapp show \
  --name ccp4i2-bicep-server \
  --resource-group ccp4i2-bicep-rg-uksouth \
  --query "properties.template.containers[0].env[?name=='ALLOWED_AZURE_AD_GROUPS']"

# If missing, redeploy
cd Docker/azure-uksouth
./scripts/deploy-applications.sh server worker
```

## Rollback Procedure

If the fix causes issues and you need to rollback:

### Option 1: Quick Disable (No Redeploy)

Remove the environment variable to disable the check:

```bash
az containerapp update \
  --name ccp4i2-bicep-server \
  --resource-group ccp4i2-bicep-rg-uksouth \
  --remove-env-vars ALLOWED_AZURE_AD_GROUPS

az containerapp update \
  --name ccp4i2-bicep-worker \
  --resource-group ccp4i2-bicep-rg-uksouth \
  --remove-env-vars ALLOWED_AZURE_AD_GROUPS
```

**Effect:** Backend Teams check disabled, all authenticated users allowed (frontend still checks)

### Option 2: Revert Code Changes

1. Revert the middleware change:
   ```bash
   git diff server/ccp4i2/middleware/azure_auth.py
   git checkout HEAD -- server/ccp4i2/middleware/azure_auth.py
   ```

2. Rebuild and redeploy:
   ```bash
   cd Docker/azure-uksouth
   ./scripts/build-and-push.sh server
   ./scripts/deploy-applications.sh server worker
   ```

## Success Criteria

The fix is working correctly when:

1. ✅ Authorized users (in team) can use the app normally
2. ✅ Unauthorized users (not in team) are blocked at frontend
3. ✅ Unauthorized users CANNOT bypass via console (403 Forbidden)
4. ✅ Logs show authorization checks happening
5. ✅ Direct API calls with unauthorized tokens fail (403)
6. ✅ Performance is not noticeably affected (groups check is fast)

## Performance Impact

**Expected impact:** Negligible
- Groups validation adds ~1ms per request (reading from already-validated JWT)
- No external API calls (groups are in the token)
- No database queries

**Monitoring:**
```bash
# Check response times before/after
time curl -H "Authorization: Bearer $TOKEN" \
     https://ccp4i2-bicep-web.uksouth.azurecontainerapps.io/api/ccp4i2/health/
```

## Next Steps After Testing

Once all tests pass:

1. ✅ Security fix is confirmed working
2. Consider adding automated tests for authorization logic
3. Consider adding project ownership for true multi-tenancy (see SECURITY_VULNERABILITIES.md)
4. Document the authorization process for new team members
5. Set up monitoring/alerts for authorization failures

## Questions?

- See: [AZURE_AD_SETUP.md](AZURE_AD_SETUP.md) for Azure AD configuration
- See: [SECURITY_VULNERABILITIES.md](../SECURITY_VULNERABILITIES.md) for detailed vulnerability analysis
- See: [SECURITY_FIX_TEAMS_BACKEND.md](../../SECURITY_FIX_TEAMS_BACKEND.md) for implementation details
