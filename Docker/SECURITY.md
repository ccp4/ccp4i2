# CCP4i2 Security Architecture

## Overview

CCP4i2 uses a multi-layer authentication and authorization system with Azure AD (Microsoft Entra ID).

## Authentication Flow

```
User → Next.js Frontend → Next.js Proxy → Django Backend
         │                    │                │
         │                    │                ├─ Validate JWT signature
         │                    │                ├─ Check expiration
         │                    │                ├─ Verify audience (app ID)
         │                    │                ├─ Verify issuer (tenant)
         │                    │                └─ Check groups membership ✓
         │                    │
         │                    └─ Forward token to backend
         │
         └─ MSAL authentication, Teams membership check (UI)
```

## Security Layers

### 1. JWT Token Validation (Backend)
- **RS256 signature** verification against Azure AD public keys
- **Expiration** check (`exp` claim)
- **Audience** validation (`aud` must match app client ID)
- **Issuer** validation (`iss` must match Azure AD tenant)

### 2. Teams/Groups Authorization (Backend)
- **Groups claim** extracted from validated JWT
- **Membership check** against `ALLOWED_AZURE_AD_GROUPS`
- **403 Forbidden** returned for non-members

### 3. Frontend Authorization (UI)
- **Microsoft Graph API** check for Teams membership
- **Access Denied** page shown to non-members
- **Defense in depth** - backend also validates

## Configuration

### Environment Variables

| Variable | Purpose |
|----------|---------|
| `CCP4I2_REQUIRE_AUTH` | Enable authentication (must be `true` for production) |
| `AZURE_AD_TENANT_ID` | Your Azure AD tenant ID |
| `AZURE_AD_CLIENT_ID` | App registration client ID |
| `ALLOWED_AZURE_AD_GROUPS` | Comma-separated group IDs for authorization |

### Azure AD App Registration

Required configuration in Token Configuration:
- Groups claim → Security groups → Group ID
- Applied to ID and Access tokens

## Deployment Protections

Build and deploy scripts refuse to proceed unless `CCP4I2_REQUIRE_AUTH=true`.

Override for local development:
```bash
ALLOW_NO_AUTH_BUILD=true ./scripts/build-and-push.sh
```

## Local Development

For local development without authentication:
```bash
CCP4I2_REQUIRE_AUTH=false
```

This creates a `dev_admin` superuser automatically. **Never use in production.**

## Adding New Authorized Groups

1. Find Group Object ID in Azure Portal → Entra ID → Groups
2. Add to `ALLOWED_AZURE_AD_GROUPS` (comma-separated)
3. Redeploy server and worker

## Monitoring

Check logs for authorization events:
- ✅ `User aad_xxx authorized via Teams/Groups membership`
- ❌ `Access denied for user aad_xxx - not in authorized groups`
