# CCP4i2 Docker Security Guide

## Authentication Architecture

CCP4i2 uses a **two-layer authentication system**:

### Layer 1: Next.js Proxy (Frontend)
- **Location:** `client/renderer/app/api/proxy/*/route.ts`
- **Purpose:** Token presence enforcement
- **Validation:** Checks if token exists (does NOT validate signature/claims)
- **Forwards:** All requests with auth headers to Django backend

### Layer 2: Django Middleware (Backend)
- **Location:** `server/ccp4i2/middleware/azure_auth.py`
- **Purpose:** Full JWT cryptographic validation
- **Validation:**
  - ✓ RS256 signature verification against Azure AD public keys
  - ✓ Token expiration (`exp` claim)
  - ✓ Audience validation (`aud` must match client ID)
  - ✓ Issuer validation (`iss` must match Azure AD tenant)
  - ✓ User identity extraction from verified `sub` claim

## Why No Validation at Next.js Proxy?

**Short answer:** The backend server is not publicly exposed, so it's redundant.

**Architecture justification:**
1. **Backend is isolated** - Only accessible via the Next.js proxy container
2. **Network security** - Django's `ALLOWED_HOSTS` restricts access to the frontend
3. **Single source of truth** - All JWT validation happens in Django (easier to maintain)
4. **Performance** - Avoids duplicate cryptographic operations
5. **Key rotation** - Only Django needs to update JWKS cache

The proxy's role is to:
- Enforce whether auth is required (yes/no decision)
- Extract and forward tokens
- Let Django handle all security-critical validation

## Deployment Modes

### 1. Local Development (No Auth) - **DEFAULT**

**Configuration:**
```bash
CCP4I2_REQUIRE_AUTH=false
NEXT_PUBLIC_REQUIRE_AUTH=false
```

**Behavior:**
- Django middleware creates a `dev_admin` superuser account
- **ALL requests get full superuser access** without any credentials
- Suitable ONLY for local development on your machine

**Security implications:**
- ⚠️ **No authentication** - Anyone can access the application
- ⚠️ **No authorization** - All users have admin privileges
- ⚠️ **No multi-tenancy** - All users see all projects

**When to use:**
- Local development on your laptop/desktop
- Behind a firewall with no external access
- Testing features that don't require auth flow

### 2. Local Development (With Auth)

**Configuration:**
```bash
CCP4I2_REQUIRE_AUTH=true
AZURE_AD_CLIENT_ID=your-app-registration-client-id
AZURE_AD_TENANT_ID=your-tenant-id
NEXT_PUBLIC_REQUIRE_AUTH=true
NEXT_PUBLIC_AAD_CLIENT_ID=your-app-registration-client-id
NEXT_PUBLIC_AAD_TENANT_ID=your-tenant-id
```

**Behavior:**
- Full Azure AD authentication flow
- JWT tokens validated against Azure AD
- Users created/authenticated via verified claims

**When to use:**
- Testing authentication workflows
- Developing auth-related features
- Simulating production environment locally

### 3. Production Deployment (Azure Container Apps)

**Configuration:**
- **Enforced by build scripts** - Will not proceed without `CCP4I2_REQUIRE_AUTH=true`
- Uses `.env.deployment` in `Docker/azure-uksouth/`

**Behavior:**
- Full Azure AD authentication required
- Django validates all JWT tokens
- Azure Easy Auth provides additional security layer (optional)

## Security Safeguards

### Build-Time Protection

Both `build-and-push.sh` and `deploy-applications.sh` now include security checks:

```bash
# Security check: Ensure authentication is required for production deployments
if [ "$CCP4I2_REQUIRE_AUTH" != "true" ]; then
    echo "⚠️  SECURITY WARNING ⚠️"
    echo "CCP4I2_REQUIRE_AUTH is not set to 'true'"
    echo "This would deploy WITHOUT authentication!"
    exit 1
fi
```

**To override for local testing only:**
```bash
export ALLOW_NO_AUTH_BUILD=true  # For build-and-push.sh
export ALLOW_NO_AUTH_DEPLOY=true # For deploy-applications.sh
./Docker/azure-uksouth/scripts/build-and-push.sh
```

### Runtime Protection

Django middleware checks `CCP4I2_REQUIRE_AUTH` on every request:
- If `false`: Creates `dev_admin` superuser (development mode)
- If `true`: Validates JWT and rejects invalid/missing tokens with 401

## Current Security Vulnerabilities

### CRITICAL: No Authorization on Projects/Jobs

**Issue:** All viewsets use `permission_classes=[]` (AllowAny)
- Any authenticated user can access ANY project/job
- No multi-tenancy isolation
- No ownership checks

**Impact:** In production, all users can see all other users' data

**Mitigation needed:**
```python
# In ProjectViewSet, JobViewSet, etc.
permission_classes = [IsAuthenticated, IsProjectOwner]
```

### HIGH: Token Exposure in URLs

**Issue:** Tokens can be passed as query parameters: `?access_token=xyz`
- Intentional for file downloads (anchor clicks can't set headers)
- Tokens logged in access logs, browser history, proxy logs

**Mitigation:**
- Use short-lived tokens (Azure AD default: 1 hour)
- Consider alternative download mechanisms (signed URLs, download tokens)
- Rate-limit download endpoints

### MEDIUM: Email Fallback from Untrusted Header

**Issue:** If JWT lacks email claim, falls back to `X-User-Email` header
- Headers can be spoofed if not behind trusted proxy
- Only affects display name, not primary identity (`sub` is still verified)

**Mitigation:**
```python
# Don't accept X-User-Email if JWT has 'sub' claim
if azure_ad_sub and not email:
    # Don't fall back to header - leave email empty
    email = f"{azure_ad_sub}@no-email-provided.local"
```

## Checklist for Production Deployment

Before deploying to production, ensure:

- [ ] `CCP4I2_REQUIRE_AUTH=true` in `.env.deployment`
- [ ] `NEXT_PUBLIC_REQUIRE_AUTH=true` in `.env.deployment`
- [ ] Azure AD app registration configured with correct redirect URIs
- [ ] `PLATFORM_ADMIN_EMAILS` set to admin users
- [ ] `SECRET_KEY` changed from default value
- [ ] `DEBUG=false` in production settings
- [ ] HTTPS enabled (handled by Azure Container Apps ingress)
- [ ] Database backups configured
- [ ] Monitor logs for auth failures

## Testing Authentication Locally

### 1. Test No-Auth Mode (Default)
```bash
cd Docker
cp .env.example .env
# Leave CCP4I2_REQUIRE_AUTH=false
docker compose up
# Access http://localhost:3000 - no login required
```

### 2. Test Azure AD Auth
```bash
cd Docker
cp .env.example .env
# Edit .env and set:
#   CCP4I2_REQUIRE_AUTH=true
#   AZURE_AD_CLIENT_ID=<your-app-id>
#   AZURE_AD_TENANT_ID=<your-tenant-id>
#   NEXT_PUBLIC_REQUIRE_AUTH=true
#   NEXT_PUBLIC_AAD_CLIENT_ID=<your-app-id>
#   NEXT_PUBLIC_AAD_TENANT_ID=<your-tenant-id>
docker compose up --build
# Access http://localhost:3000 - redirects to Microsoft login
```

### 3. Test Production Build Rejection
```bash
cd Docker/azure-uksouth
# Edit .env.deployment and set CCP4I2_REQUIRE_AUTH=false
./scripts/build-and-push.sh
# Should fail with security warning
```

## File Reference

### Authentication Configuration
- **Azure production config:** `Docker/azure-uksouth/.env.deployment`
- **Local docker-compose config:** `Docker/.env` (copy from `.env.example`)
- **Example with docs:** `Docker/.env.example`

### Authentication Implementation
- **Django middleware:** `server/ccp4i2/middleware/azure_auth.py`
- **Next.js proxy (CCP4i2):** `client/renderer/app/api/proxy/ccp4i2/[...path]/route.ts`
- **Next.js proxy (Compounds):** `apps/compounds/frontend/app/api/proxy/compounds/[...path]/route.ts`

### Security Checks
- **Build script:** `Docker/azure-uksouth/scripts/build-and-push.sh`
- **Deploy script:** `Docker/azure-uksouth/scripts/deploy-applications.sh`

## Common Questions

### Q: Why is auth disabled by default for local development?
**A:** To simplify the developer experience. Most local development doesn't require realistic auth flows. Developers can enable auth when testing auth-specific features.

### Q: Can I deploy to my own server without Azure AD?
**A:** Yes, but you'll need to implement alternative authentication:
1. Set `CCP4I2_REQUIRE_AUTH=false` (creates dev_admin)
2. Put the application behind a reverse proxy with HTTP basic auth
3. Or implement a different auth backend in Django

### Q: How do I add a new admin user?
**A:** Add their email to `PLATFORM_ADMIN_EMAILS` in `.env.deployment` and redeploy. The middleware automatically grants admin privileges to matching emails.

### Q: What happens if Azure AD keys rotate?
**A:** The Django middleware caches Azure AD's JWKS (public keys) for 1 hour. Key rotation is handled automatically - the middleware will fetch new keys when the cache expires.

### Q: Can I use this with a different identity provider?
**A:** Yes, but you'd need to modify `azure_auth.py` to:
1. Change the JWKS URL to your provider's
2. Adjust claim names (`sub`, `email`, etc.)
3. Update issuer/audience validation logic
