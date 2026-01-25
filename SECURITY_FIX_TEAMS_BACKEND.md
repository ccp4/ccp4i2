# Adding Teams Authorization to Django Backend

## Problem

Currently, Teams membership is only checked in the **frontend** (`require-auth.tsx`). This can be bypassed via direct API calls from the browser console.

**Attack:** Any tenant user can bypass the frontend UI and call APIs directly:
```javascript
fetch('/api/ccp4i2/projects/', {
  headers: { 'Authorization': 'Bearer ' + validToken }
})
```

## Solution: Backend Teams Validation

Add Teams membership validation to Django middleware using the `groups` claim in the JWT token.

### Option 1: Using JWT `groups` Claim (Recommended)

**Prerequisites:**
1. Azure AD Premium P1 or P2 license (required for group claims in tokens)
2. Configure Azure AD app to emit `groups` claim

**Steps:**

#### 1. Configure Azure AD App Registration

In Azure Portal → App Registrations → Your App → Token Configuration:

1. Click "Add groups claim"
2. Select "Security groups"
3. Check "Group ID" (not display name - IDs are stable)
4. Apply to both ID tokens and Access tokens

**Note:** If you have >200 groups, Azure AD won't include them in the token. Instead, it adds an `_claim_names` overage indicator. For most orgs, this isn't an issue.

#### 2. Update Django Middleware

Edit `server/ccp4i2/middleware/azure_auth.py`:

```python
# At the top, add environment variable for allowed groups
ALLOWED_AZURE_AD_GROUPS = os.environ.get("ALLOWED_AZURE_AD_GROUPS", "").split(",")
# Format: comma-separated Azure AD Group Object IDs
# Example: "6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d,other-group-id"

# In process_token() method, after line 339 (after getting azure_ad_sub):
# Add this right after: request.azure_ad_user_id = azure_ad_sub

# --- Add Teams/Groups authorization check ---
# Extract groups claim from validated JWT
user_groups = claims.get("groups", [])
logger.info(f"User {azure_ad_sub[:8]}... has groups: {user_groups}")

# Check if groups validation is enabled and required
if ALLOWED_AZURE_AD_GROUPS and ALLOWED_AZURE_AD_GROUPS[0]:  # If configured
    # Check for group claims overage (happens with >200 groups)
    if "_claim_names" in claims or "_claim_sources" in claims:
        logger.warning(
            f"Group claims overage detected for {azure_ad_sub[:8]}. "
            "User has >200 groups. Cannot validate Teams membership from token. "
            "Consider using smaller security groups or Microsoft Graph API."
        )
        return JsonResponse(
            {
                "success": False,
                "error": "Group membership cannot be verified - user has too many groups. "
                         "Contact administrator to be added to a dedicated app access group.",
            },
            status=403,
        )

    # Validate user is in at least one allowed group
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
                         "Please contact your administrator to request access to the "
                         "Newcastle Drug Discovery Unit team.",
                "user_groups": user_groups,  # For debugging (remove in production)
                "required_groups": ALLOWED_AZURE_AD_GROUPS,  # For debugging (remove in production)
            },
            status=403,
        )

    logger.info(f"✅ User {azure_ad_sub[:8]} authorized via group membership")
# --- End Teams/Groups authorization check ---
```

#### 3. Update `.env.deployment`

```bash
# Add to Docker/azure-uksouth/.env.deployment:

# Teams/Groups Authorization
# Comma-separated list of Azure AD Group Object IDs that can access the app
# Get this from: Azure Portal → Azure AD → Groups → [Your Team] → Object ID
ALLOWED_AZURE_AD_GROUPS=6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d
```

**Finding Group IDs:**
- Azure Portal → Azure Active Directory → Groups
- Search for "Newcastle Drug Discovery Unit"
- Click on it → Copy "Object ID" (UUID format)

#### 4. Deploy

```bash
cd Docker/azure-uksouth
./scripts/build-and-push.sh server
./scripts/deploy-applications.sh server worker
```

### Option 2: Using Microsoft Graph API (No Premium License Required)

If you don't have Azure AD Premium P1/P2, you can validate Teams membership by calling Microsoft Graph API from Django.

**Pros:**
- Works with Azure AD Free
- Can check Teams membership dynamically
- No token size concerns

**Cons:**
- Adds latency (external API call on every request)
- Requires server-side Graph API permissions
- More complex implementation

**Implementation:**

```python
# server/ccp4i2/middleware/azure_auth.py

import requests
from functools import lru_cache
import time

# Cache Graph API tokens for 50 minutes (Azure tokens valid for 60 minutes)
@lru_cache(maxsize=1)
def get_graph_api_token(cache_key: str) -> str:
    """
    Get an app-only access token for Microsoft Graph API
    Uses client credentials flow (server-to-server auth)

    cache_key is just for cache invalidation (current hour timestamp)
    """
    tenant_id = os.environ.get("AZURE_AD_TENANT_ID")
    client_id = os.environ.get("AZURE_AD_CLIENT_ID")
    client_secret = os.environ.get("AZURE_AD_CLIENT_SECRET")  # NEW: Need this

    if not all([tenant_id, client_id, client_secret]):
        raise ValueError("Missing Graph API credentials in environment")

    token_url = f"https://login.microsoftonline.com/{tenant_id}/oauth2/v2.0/token"

    response = requests.post(
        token_url,
        data={
            "grant_type": "client_credentials",
            "client_id": client_id,
            "client_secret": client_secret,
            "scope": "https://graph.microsoft.com/.default",
        },
    )

    if response.status_code != 200:
        raise Exception(f"Failed to get Graph API token: {response.text}")

    return response.json()["access_token"]


def check_user_teams_membership(user_id: str, allowed_team_ids: list[str]) -> bool:
    """
    Check if user is a member of any allowed Teams using Microsoft Graph API

    Args:
        user_id: Azure AD user object ID (from 'sub' claim)
        allowed_team_ids: List of Team IDs to check

    Returns:
        True if user is member of at least one allowed team
    """
    # Get app-only token (cache key invalidates every hour)
    cache_key = str(int(time.time() // 3000))  # Changes every ~50 minutes
    token = get_graph_api_token(cache_key)

    # Call Graph API to get user's joined teams
    headers = {"Authorization": f"Bearer {token}"}

    # Use /users/{id}/joinedTeams endpoint
    graph_url = f"https://graph.microsoft.com/v1.0/users/{user_id}/joinedTeams"

    response = requests.get(graph_url, headers=headers)

    if response.status_code != 200:
        logger.error(f"Graph API call failed: {response.text}")
        # Fail closed - deny access if we can't verify
        return False

    user_teams = response.json().get("value", [])
    user_team_ids = [team["id"] for team in user_teams]

    # Check if user is in any allowed team
    return any(team_id in allowed_team_ids for team_id in user_team_ids)


# In process_token() method, add after getting azure_ad_sub:

allowed_team_ids = os.environ.get("ALLOWED_TEAMS_IDS", "").split(",")
if allowed_team_ids and allowed_team_ids[0]:
    logger.info(f"Checking Teams membership for {azure_ad_sub[:8]}...")

    if not check_user_teams_membership(azure_ad_sub, allowed_team_ids):
        logger.warning(f"Access denied - not in allowed Teams")
        return JsonResponse(
            {
                "success": False,
                "error": "Access denied: You are not a member of an authorized Microsoft Team. "
                         "Please contact your administrator to request access.",
            },
            status=403,
        )

    logger.info(f"✅ User authorized via Teams membership")
```

**Additional Requirements for Option 2:**

1. **Create Client Secret:**
   - Azure Portal → App Registrations → Your App → Certificates & secrets
   - New client secret → Copy the value
   - Add to `.env.deployment`: `AZURE_AD_CLIENT_SECRET=your-secret-value`

2. **Grant API Permissions:**
   - Azure Portal → App Registrations → Your App → API permissions
   - Add permission → Microsoft Graph → Application permissions
   - Select `Team.ReadBasic.All` (application permission, not delegated)
   - Click "Grant admin consent"

3. **Store Team IDs:**
   ```bash
   # In .env.deployment:
   ALLOWED_TEAMS_IDS=6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d
   ```

### Option 3: App Roles (Alternative to Groups)

If you prefer role-based access instead of group-based:

**Azure AD Configuration:**
1. Azure Portal → App Registrations → Your App → App roles
2. Create new app role:
   - Display name: "CCP4i2 User"
   - Allowed member types: Users/Groups
   - Value: `User`
   - Description: "Standard user access to CCP4i2"
3. Assign users/groups to this role via Enterprise Applications

**Django Middleware:**
```python
# In process_token(), check for app roles:
roles = claims.get("roles", [])

if not roles or "User" not in roles:
    return JsonResponse(
        {"success": False, "error": "Access denied: Missing required app role"},
        status=403,
    )
```

**Pros:**
- Simpler than groups
- Works with Azure AD Free
- No Graph API calls needed
- Token already contains roles

**Cons:**
- Requires manual role assignment for each user
- Less flexible than Teams-based access

## Recommended Approach

**For your use case (Newcastle Drug Discovery Unit):**

Use **Option 1: JWT `groups` Claim** if you have Azure AD Premium.

**Why:**
- ✓ Simplest implementation
- ✓ No external API calls (lower latency)
- ✓ Matches your existing frontend logic
- ✓ Team membership managed automatically via Teams UI
- ✓ Secure (group IDs in cryptographically signed JWT)

**If no Azure AD Premium:**
Use **Option 3: App Roles** (simpler than Option 2).

## Testing

After implementing, test the bypass is fixed:

```javascript
// 1. Create a test user NOT in the Newcastle Drug Discovery Unit team
// 2. Have them log in normally - frontend should block them ✓
// 3. Have them open console and run:

fetch('/api/ccp4i2/projects/', {
  headers: { 'Authorization': 'Bearer ' + await getToken() }
})

// Expected result (AFTER fix): 403 Forbidden
// Current result (BEFORE fix): 200 OK with all projects
```

## Migration Plan

1. **Phase 1:** Implement backend Groups validation (Option 1)
2. **Phase 2:** Test with a few users
3. **Phase 3:** Deploy to production
4. **Phase 4:** Remove debugging fields (`user_groups`, `required_groups`) from error response

## Environment Variables Summary

Add to `.env.deployment`:

```bash
# Backend Teams/Groups Authorization
ALLOWED_AZURE_AD_GROUPS=6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d

# OR for App Roles approach (alternative):
# ALLOWED_APP_ROLES=User

# OR for Graph API approach (if no Premium license):
# AZURE_AD_CLIENT_SECRET=your-client-secret
# ALLOWED_TEAMS_IDS=6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d
```
