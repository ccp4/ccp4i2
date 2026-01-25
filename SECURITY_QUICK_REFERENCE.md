# Security Quick Reference

## Current Status: VULNERABLE ⚠️

### What Works ✅
- **Frontend authorization** - Teams membership checked via Microsoft Graph API
- **JWT validation** - Signature, expiration, audience, issuer all verified
- **Token audience check** - Rejects tokens from other Azure apps

### What's Broken 🚨
- **Backend has NO Teams validation** - Console attacks bypass frontend
- **No project ownership** - All users see all projects
- **No authorization on API endpoints** - All use `permission_classes=[]`

## Attack Vectors

### ✅ BLOCKED: Token from Another App
```javascript
// Attacker has token for different app → REJECTED by Django
// Reason: aud claim doesn't match AZURE_AD_CLIENT_ID
```

### 🚨 WIDE OPEN: Console Bypass
```javascript
// Attacker in your tenant, not in Teams → Frontend blocks them
// BUT they can bypass via console:
fetch('/api/ccp4i2/projects/', {
  headers: { 'Authorization': 'Bearer ' + validToken }
})
// Result: 200 OK with all projects (SHOULD BE 403!)
```

### 🚨 WIDE OPEN: Cross-User Access
```javascript
// User A creates project ID=123
// User B (different tenant user) can access it:
fetch('/api/ccp4i2/projects/123/')
// Result: 200 OK with User A's data (SHOULD BE 403!)
```

## Immediate Fix (15 Minutes)

Add Teams validation to Django backend:

### File: `server/ccp4i2/middleware/azure_auth.py`

**Add after line 340** (after `request.azure_ad_user_id = azure_ad_sub`):

```python
# --- TEAMS AUTHORIZATION CHECK ---
ALLOWED_AZURE_AD_GROUPS = os.environ.get("ALLOWED_AZURE_AD_GROUPS", "").split(",")
user_groups = claims.get("groups", [])

if ALLOWED_AZURE_AD_GROUPS and ALLOWED_AZURE_AD_GROUPS[0]:
    # Check for overage (>200 groups)
    if "_claim_names" in claims or "_claim_sources" in claims:
        logger.warning(f"Group claims overage for {azure_ad_sub[:8]}")
        return JsonResponse(
            {"success": False, "error": "Too many groups - cannot verify membership"},
            status=403,
        )

    # Validate user is in allowed groups
    has_access = any(group_id in ALLOWED_AZURE_AD_GROUPS for group_id in user_groups)

    if not has_access:
        logger.warning(
            f"Access denied for {azure_ad_sub[:8]}. "
            f"User groups: {user_groups}, Required: {ALLOWED_AZURE_AD_GROUPS}"
        )
        return JsonResponse(
            {
                "success": False,
                "error": "Access denied: You are not a member of an authorized group. "
                         "Contact your administrator to request access.",
            },
            status=403,
        )

    logger.info(f"✅ User {azure_ad_sub[:8]} authorized via Teams membership")
# --- END TEAMS AUTHORIZATION CHECK ---
```

### File: `Docker/azure-uksouth/.env.deployment`

Add this line (your Teams group ID from frontend code):

```bash
ALLOWED_AZURE_AD_GROUPS=6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d
```

### Azure AD Configuration

1. Go to: Azure Portal → App Registrations → Your App → Token Configuration
2. Click "Add groups claim"
3. Select "Security groups"
4. Check "Group ID"
5. Apply to both ID tokens and Access tokens

### Deploy

```bash
cd Docker/azure-uksouth
./scripts/build-and-push.sh server
./scripts/deploy-applications.sh server worker
```

## Testing the Fix

### Test 1: Console Attack Should Fail
```javascript
// User NOT in "Newcastle Drug Discovery Unit" team
// Opens console on your app and runs:
fetch('/api/ccp4i2/projects/')

// BEFORE fix: 200 OK with all projects ✗
// AFTER fix: 403 Forbidden ✓
```

### Test 2: Frontend Still Works
```javascript
// User IN the team logs in normally
// Should see app UI and be able to use it
// API calls via UI should work: 200 OK
```

### Test 3: Check Logs
```bash
# Should see in Django logs:
# ✅ User aad_xxxx authorized via Teams membership
# OR for non-members:
# Access denied for aad_xxxx. User groups: [...], Required: [...]
```

## Future Work (Not Urgent)

### Project Ownership (For Multi-Tenancy)
- Add `owner` field to Project model
- Filter querysets by `owner=request.user`
- Implement sharing mechanism

See: [SECURITY_VULNERABILITIES.md](Docker/SECURITY_VULNERABILITIES.md#option-2-project-ownership-high-effort-best-security)

## Questions?

- **Do I have Azure AD Premium?**
  - Check: Azure Portal → Azure AD → Licenses
  - Need P1 or P2 for `groups` claim in tokens

- **What if I don't have Premium?**
  - Use App Roles instead (simpler) or Microsoft Graph API (more complex)
  - See: [SECURITY_FIX_TEAMS_BACKEND.md](SECURITY_FIX_TEAMS_BACKEND.md#option-3-app-roles-alternative-to-groups)

- **How do I find my Team's Group ID?**
  - Already in your code: `6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d`
  - Or: Azure Portal → Azure AD → Groups → Search "Newcastle Drug Discovery Unit" → Object ID

- **Will this break existing users?**
  - No - users in the Teams group will continue working
  - Users NOT in the group will be blocked (as intended)
  - Frontend already does this, backend will now enforce it too

## Build Script Protection

Your build/deploy scripts now refuse to proceed without `CCP4I2_REQUIRE_AUTH=true`:

```bash
# This will FAIL with security warning:
CCP4I2_REQUIRE_AUTH=false ./scripts/build-and-push.sh

# Override for local dev only:
ALLOW_NO_AUTH_BUILD=true ./scripts/build-and-push.sh
```

## Documentation

- [SECURITY.md](Docker/SECURITY.md) - Overall security architecture
- [SECURITY_VULNERABILITIES.md](Docker/SECURITY_VULNERABILITIES.md) - Detailed vulnerability analysis
- [SECURITY_FIX_TEAMS_BACKEND.md](SECURITY_FIX_TEAMS_BACKEND.md) - Implementation guide for backend Teams auth

## Summary

**Priority 1 (CRITICAL):** Add Teams validation to backend (15 min fix above)

**Priority 2 (HIGH):** Add project ownership for true multi-tenancy (future work)

**Priority 3 (MEDIUM):** Remove token-in-query-parameter support (security hygiene)
