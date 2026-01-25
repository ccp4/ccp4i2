# CCP4i2 Security Documentation

## Current Security Status: ✅ Teams Authorization Implemented

The critical console bypass vulnerability has been fixed. Backend now validates Teams/Groups membership.

## Authentication & Authorization Layers

### Layer 1: Next.js Proxy (Frontend)
- **Location:** `client/renderer/app/api/proxy/*/route.ts`
- **Validates:** Token presence (not validity)
- **Forwards:** Requests with auth headers to Django

### Layer 2: Django Middleware (Backend)
- **Location:** `server/ccp4i2/middleware/azure_auth.py`
- **Validates:**
  - ✓ JWT signature (RS256 against Azure AD keys)
  - ✓ Token expiration (`exp` claim)
  - ✓ Audience (`aud` must match client ID)
  - ✓ Issuer (`iss` must match Azure AD tenant)
  - ✓ **Teams/Groups membership** (`groups` claim)

### Layer 3: Frontend Authorization
- **Location:** `client/renderer/components/require-auth.tsx`
- **Validates:** Teams membership via Microsoft Graph API
- **Shows:** "Access Denied" UI for non-members

## Protection Summary

| Attack Vector | Protected? | How |
|---------------|-----------|-----|
| External attacker (wrong tenant) | ✓ | Issuer validation |
| Token from another app | ✓ | Audience validation |
| Expired token | ✓ | Expiration validation |
| Invalid signature | ✓ | RS256 signature check |
| Tenant user not in team (UI) | ✓ | Frontend blocks |
| Tenant user not in team (console) | ✓ | Backend groups check |

## Configuration

### Required Environment Variables

```bash
# Backend auth (Django)
CCP4I2_REQUIRE_AUTH=true
AZURE_AD_TENANT_ID=9c5012c9-b616-44c2-a917-66814fbe3e87
AZURE_AD_CLIENT_ID=386da83f-1bf4-4ad8-b742-79b600e2208b

# Teams/Groups authorization
ALLOWED_AZURE_AD_GROUPS=6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d

# Frontend auth (Next.js)
NEXT_PUBLIC_REQUIRE_AUTH=true
NEXT_PUBLIC_AAD_CLIENT_ID=386da83f-1bf4-4ad8-b742-79b600e2208b
NEXT_PUBLIC_AAD_TENANT_ID=9c5012c9-b616-44c2-a917-66814fbe3e87
```

### Azure AD Configuration

Token Configuration must include:
- Groups claim → Security groups → Group ID
- Applied to: ID tokens, Access tokens

## Build Script Protections

Both `build-and-push.sh` and `deploy-applications.sh` will refuse to proceed if `CCP4I2_REQUIRE_AUTH != true`.

Override for local development only:
```bash
ALLOW_NO_AUTH_BUILD=true ./scripts/build-and-push.sh
```

## Future Considerations

### Project Ownership (Not Yet Implemented)

Currently, all authorized team members can see all projects. For true multi-tenancy:

1. Add `owner` field to Project model
2. Filter querysets by owner
3. Implement sharing mechanism

This is a future enhancement if per-user data isolation is needed.

## Rollback Procedure

If needed, disable Teams authorization:

```bash
# Quick disable (no redeploy)
az containerapp update \
  --name ccp4i2-bicep-server \
  --resource-group ccp4i2-bicep-rg-uksouth \
  --remove-env-vars ALLOWED_AZURE_AD_GROUPS
```

**Note:** Removing `ALLOWED_AZURE_AD_GROUPS` allows all authenticated tenant users (less secure).

## Log Messages

**Authorized user:**
```
✅ User aad_12345678 authorized via Teams/Groups membership
```

**Unauthorized user:**
```
Access denied for user aad_87654321 - not in authorized groups.
User groups: ['other-group'], Required: ['6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d']
```
