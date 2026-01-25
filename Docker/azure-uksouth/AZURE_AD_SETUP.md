# Azure AD Configuration for Groups Claim

## Overview

To enable backend Teams authorization, you need to configure your Azure AD app registration to include the `groups` claim in JWT tokens.

## Prerequisites

- Azure AD Premium P1 or P2 license (required for groups claim)
- Access to Azure Portal with permissions to modify app registrations
- App registration already exists: Client ID `386da83f-1bf4-4ad8-b742-79b600e2208b`

## Step-by-Step Instructions

### 1. Access App Registration

1. Go to [Azure Portal](https://portal.azure.com)
2. Navigate to: **Azure Active Directory** → **App registrations**
3. Search for your app using client ID: `386da83f-1bf4-4ad8-b742-79b600e2208b`
4. Click on the app to open its configuration

### 2. Configure Token Configuration

1. In the left sidebar, click **Token configuration**
2. Click **+ Add groups claim** button
3. In the "Edit groups claim" dialog:
   - Select **Security groups** checkbox
   - Under "Customize token properties by type":
     - For **ID**: Check **Group ID**
     - For **Access**: Check **Group ID**
     - For **SAML**: Check **Group ID** (if using SAML)
4. Click **Add** to save

**What this does:**
- Azure AD will now include a `groups` claim in JWT tokens
- The claim contains an array of Group Object IDs (UUIDs)
- These are the security groups the user is a member of

### 3. Verify Configuration

After adding the groups claim, you should see:

```
Token configuration
├─ Groups claim
│  ├─ Security groups
│  ├─ Emit as Group ID
│  └─ ID tokens, Access tokens
```

### 4. Test Token Contains Groups

To verify the configuration works:

1. Log in to your application as a user in the "Newcastle Drug Discovery Unit" team
2. Open browser console
3. Run this to decode the JWT token (only works if you have the token):

```javascript
// Get the token from MSAL
const accounts = window.msal?.getAllAccounts();
if (accounts?.[0]) {
  const response = await window.msal.acquireTokenSilent({
    scopes: ["openid", "profile"],
    account: accounts[0]
  });

  // Decode JWT (middle part between dots)
  const parts = response.idToken.split('.');
  const payload = JSON.parse(atob(parts[1]));

  console.log('Groups in token:', payload.groups);
  // Should show: ["6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d", ...]
}
```

**Expected result:**
- `groups` array should be present in the token
- Should contain: `"6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d"` (Newcastle Drug Discovery Unit)

### 5. Important Notes

#### Group Overage (>200 Groups)

If a user is a member of more than 200 groups:
- Azure AD won't include the full `groups` array
- Instead, it adds `_claim_names` and `_claim_sources` with an overage indicator
- The middleware will detect this and return a 403 error
- **Solution:** Create a dedicated smaller security group just for app access

#### Group ID vs Display Name

- We use **Group Object IDs** (UUIDs), not display names
- Group IDs are stable and don't change when the team is renamed
- Group ID for "Newcastle Drug Discovery Unit": `6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d`

#### Finding Group IDs

To find the Object ID of any group:

1. Azure Portal → **Azure Active Directory** → **Groups**
2. Search for the group name (e.g., "Newcastle Drug Discovery Unit")
3. Click on the group
4. Copy the **Object ID** (UUID format)

## Deployment

After configuring Azure AD:

1. **Build the server image:**
   ```bash
   cd Docker/azure-uksouth
   ./scripts/build-and-push.sh server
   ```

2. **Deploy the server and worker:**
   ```bash
   ./scripts/deploy-applications.sh server
   ./scripts/deploy-applications.sh worker
   ```

3. **Verify deployment:**
   ```bash
   # Check that ALLOWED_AZURE_AD_GROUPS is set
   az containerapp show \
     --name ccp4i2-bicep-server \
     --resource-group ccp4i2-bicep-rg-uksouth \
     --query "properties.template.containers[0].env[?name=='ALLOWED_AZURE_AD_GROUPS'].value" \
     -o tsv

   # Should output: 6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d
   ```

## Testing

### Test 1: Authorized User (In Team)

1. Log in as a user who IS in the "Newcastle Drug Discovery Unit" team
2. Should see the app UI normally ✓
3. API calls should work ✓
4. Check logs for: `✅ User aad_xxx authorized via Teams/Groups membership`

### Test 2: Unauthorized User (Not in Team)

1. Create a test user in your tenant who is NOT in the team
2. Have them try to log in to the app
3. Frontend should block them with "Access Denied" ✓
4. If they try console bypass:
   ```javascript
   fetch('/api/ccp4i2/projects/')
   ```
5. Should get **403 Forbidden** with error message ✓
6. Check logs for: `Access denied for user aad_xxx - not in authorized groups`

### Test 3: Console Bypass Protection

**Before the fix:**
```javascript
// User not in team opens console
fetch('/api/ccp4i2/projects/', {
  headers: { 'Authorization': 'Bearer ' + token }
})
// Result: 200 OK with all projects ✗
```

**After the fix:**
```javascript
// Same attack attempt
fetch('/api/ccp4i2/projects/', {
  headers: { 'Authorization': 'Bearer ' + token }
})
// Result: 403 Forbidden ✓
// Response: {"success": false, "error": "Access denied: You are not a member of an authorized group..."}
```

## Troubleshooting

### Problem: "Groups claim not in token"

**Symptoms:**
- Users who should have access are being denied
- Logs show: `User aad_xxx has groups: []`

**Solutions:**
1. Verify Token Configuration in Azure AD has groups claim enabled
2. Check you're using ID token (not access token) - groups are in ID token
3. Wait 5 minutes after configuration for token cache to refresh
4. Have user log out and log back in to get fresh token

### Problem: "Group overage detected"

**Symptoms:**
- Error: "Your account has too many group memberships to verify automatically"

**Solutions:**
1. Create a dedicated "CCP4i2 Access" security group with <200 members
2. Add that group's Object ID to `ALLOWED_AZURE_AD_GROUPS`
3. Add users to this smaller group instead of relying on Teams membership directly

### Problem: "User in team but still denied"

**Symptoms:**
- User is definitely in the Teams group
- Still getting 403 Forbidden

**Debugging:**
1. Check the token has the correct group ID:
   ```javascript
   // In browser console after login
   const accounts = window.msal?.getAllAccounts();
   const token = await window.msal.acquireTokenSilent({
     scopes: ["openid"], account: accounts[0]
   });
   const payload = JSON.parse(atob(token.idToken.split('.')[1]));
   console.log('Groups:', payload.groups);
   ```

2. Verify the group ID matches:
   - Token contains: `6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d`
   - `.env.deployment` has: `ALLOWED_AZURE_AD_GROUPS=6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d`

3. Check deployment has the env var:
   ```bash
   az containerapp show \
     --name ccp4i2-bicep-server \
     --resource-group ccp4i2-bicep-rg-uksouth \
     --query "properties.template.containers[0].env" \
     -o json | grep -A 2 ALLOWED_AZURE_AD_GROUPS
   ```

## Security Implications

### What This Protects Against ✅

- ✓ **Console bypass attacks** - Non-team members can't call APIs directly
- ✓ **Token reuse** - Valid tokens from non-team users are rejected
- ✓ **Unauthorized access** - Only team members can use the app

### What This DOESN'T Protect Against ❌

- ✗ **Cross-user data access** - Team members can still see each other's projects
- ✗ **Project ownership** - No per-project access control
- ✗ **Token theft** - If an authorized user's token is stolen, it works

**For full multi-tenancy, you'll also need project ownership (see SECURITY_VULNERABILITIES.md).**

## Documentation References

- [Microsoft Docs: Configure groups claim](https://learn.microsoft.com/en-us/azure/active-directory/develop/optional-claims#configuring-groups-optional-claims)
- [Microsoft Docs: Groups overage claim](https://learn.microsoft.com/en-us/azure/active-directory/develop/id-tokens#groups-overage-claim)
- [Azure AD Premium licensing](https://www.microsoft.com/en-us/security/business/identity-access-management/azure-ad-pricing)

## Rollback Plan

If you need to disable Teams authorization (emergency):

1. **Quick rollback (no redeploy):**
   ```bash
   az containerapp update \
     --name ccp4i2-bicep-server \
     --resource-group ccp4i2-bicep-rg-uksouth \
     --remove-env-vars ALLOWED_AZURE_AD_GROUPS
   ```

2. **Permanent rollback:**
   - Edit `.env.deployment` and remove `ALLOWED_AZURE_AD_GROUPS=...` line
   - Redeploy server and worker

**Note:** Removing the env var disables the check - all authenticated users will be allowed (less secure).
