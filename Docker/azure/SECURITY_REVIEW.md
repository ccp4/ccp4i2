# Security Review - CCP4i2 Azure Deployment

Last reviewed: 2026-01-01

## Summary

The application has good security fundamentals around network isolation and secrets management, but there are configuration issues that should be addressed.

---

## Issues to Address

### HIGH Priority

- [x] **DEBUG=true in production Bicep** ✅ Fixed 2026-01-01
  - Location: `infrastructure/applications.bicep:228`
  - Risk: Exposes stack traces and sensitive debugging info to users
  - Fix: Changed `value: 'true'` to `value: 'false'`

- [x] **Hardcoded insecure SECRET_KEY fallback** ✅ Fixed 2026-01-01
  - Location: `server/ccp4i2/config/settings.py:17-28`
  - Risk: If SECRET_KEY env var is not set, an insecure default is used
  - Fix: Now raises ValueError in production (DEBUG=false) if SECRET_KEY not set

### MEDIUM Priority

- [ ] **Server container runs as root**
  - Location: `Docker/server/Dockerfile`
  - Risk: Container escape could give root access
  - Fix: Add non-root user similar to client Dockerfile

- [ ] **No authorization/RBAC implementation**
  - Location: Middleware layer
  - Risk: Any authenticated user can access any project/data
  - Fix: Implement per-project authorization checks

- [ ] **ACR uses admin credentials instead of managed identity**
  - Location: `infrastructure/applications.bicep:100-109`
  - Risk: Admin credentials are more powerful than needed
  - Fix: Use managed identity for ACR authentication

### LOW Priority

- [ ] **Storage has public network access enabled**
  - Location: `infrastructure/infrastructure.bicep:244`
  - Risk: Storage is publicly accessible (though SAS-protected)
  - Note: Required for external SAS URL uploads - acceptable with current SAS controls

- [ ] **DB SSL certificate verification disabled**
  - Location: `infrastructure/applications.bicep:224-225`
  - Risk: MITM attacks (mitigated by private endpoint)
  - Fix: Enable SSL cert verification

---

## What's Done Well

### Network Security (Strong)
- All backend services use private endpoints (PostgreSQL, Key Vault, Service Bus, ACR)
- Container Apps run in dedicated VNet with proper subnet isolation
- Server API is internal-only (`external: false`)
- Private DNS zones properly configured for all services
- NSGs configured for Container Apps and private endpoints

### Secrets Management (Strong)
- Azure Key Vault with RBAC and private access only
- Database passwords generated randomly and stored in Key Vault
- Container Apps use managed identity for Key Vault access
- User delegation SAS tokens (short-lived, limited permissions)
- No secrets in source code or Docker images

### Authentication (Good)
- Azure AD JWT validation with proper checks:
  - Signature verification (RS256)
  - Token expiration
  - Audience (client ID)
  - Issuer (tenant ID)
- JWKS key caching with 1-hour TTL
- SSL certificate handling with certifi
- Supports both Bearer token and Azure Easy Auth headers

### Container Security (Partial)
- Client container runs as non-root user
- Minimal base images (Alpine, slim)
- Health checks configured
- Liveness, readiness, and startup probes

---

## Architecture Notes

```
External Users
     |
     v
[Azure AD] <-- Authentication
     |
     v
[Web Container App] -- external, HTTPS only
     |
     v (internal VNet)
[Server Container App] -- internal only
     |
     +---> [PostgreSQL] -- private endpoint
     +---> [Key Vault] -- private endpoint
     +---> [Storage] -- private endpoint (files), public (blobs for SAS)
     +---> [Service Bus] -- private endpoint
     |
     v
[Worker Container App] -- no ingress, background jobs
```

---

## Recommendations for Future

1. **Rate Limiting** - Add API rate limiting to prevent abuse
2. **WAF** - Consider Azure Front Door or Application Gateway for WAF protection
3. **Image Scanning** - Enable Azure Container Registry vulnerability scanning
4. **Secret Rotation** - Implement automatic secret rotation policies
5. **Audit Logging** - Add security-specific alerting and SIEM integration
6. **File Validation** - Add file type validation for uploads
