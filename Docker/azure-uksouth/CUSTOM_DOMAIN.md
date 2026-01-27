# Custom Domain & Certificate Management

This guide covers setting up, maintaining, and recovering custom domains and SSL certificates for the CCP4i2 Azure Container Apps deployment.

## Table of Contents

1. [Current Configuration](#current-configuration)
2. [Initial Custom Domain Setup](#initial-custom-domain-setup)
3. [Certificate Management](#certificate-management)
4. [Routine Maintenance](#routine-maintenance)
5. [Troubleshooting](#troubleshooting)
6. [Disaster Recovery](#disaster-recovery)
7. [Quick Reference Commands](#quick-reference-commands)

---

## Current Configuration

| Setting | Value |
|---------|-------|
| Custom Domain | `ddudatabase.ncl.ac.uk` |
| Container App | `ccp4i2-bicep-web` |
| Environment | `ccp4i2-bicep-env-uk` |
| Resource Group | `ccp4i2-bicep-rg-uksouth` |
| Default FQDN | `ccp4i2-bicep-web.ambitiousrock-87fff394.uksouth.azurecontainerapps.io` |

---

## Initial Custom Domain Setup

### Prerequisites

Before setting up a custom domain, you need:

1. **DNS access** - Ability to create CNAME and TXT records for your domain
2. **Azure CLI** - Installed and authenticated (`az login`)
3. **Environment variables** - Source the deployment config

```bash
# Always start by sourcing environment variables
. ./Docker/azure-uksouth/.env.deployment
```

### Step 1: Get Domain Verification Token

Azure requires a TXT record to prove domain ownership. Get the verification token:

```bash
# Get the Container Apps environment's custom domain verification ID
az containerapp env show \
  --name ccp4i2-bicep-env-uk \
  --resource-group "$RESOURCE_GROUP" \
  --query "properties.customDomainConfiguration.customDomainVerificationId" \
  -o tsv
```

This returns a value like: `E055531894B68972D34C9700D6042863796A62AF940D092099E30F2EF022CF73`

### Step 2: Configure DNS Records

Request your IT team to add these DNS records:

| Record Type | Host/Name | Value |
|-------------|-----------|-------|
| CNAME | `ddudatabase.ncl.ac.uk` | `ccp4i2-bicep-web.ambitiousrock-87fff394.uksouth.azurecontainerapps.io` |
| TXT | `asuid.ddudatabase.ncl.ac.uk` | `<verification-id-from-step-1>` |

**Why CNAME instead of A record?**
- Azure Container Apps IP addresses can change during scaling or maintenance
- CNAME automatically follows IP changes
- A records require manual updates when IPs change

### Step 3: Verify DNS Propagation

Wait for DNS propagation (can take 5-60 minutes), then verify:

```bash
# Check CNAME record
nslookup ddudatabase.ncl.ac.uk

# Check TXT verification record
nslookup -type=TXT asuid.ddudatabase.ncl.ac.uk

# Alternative using dig (more detailed)
dig ddudatabase.ncl.ac.uk CNAME +short
dig asuid.ddudatabase.ncl.ac.uk TXT +short
```

Expected output for TXT record:
```
"E055531894B68972D34C9700D6042863796A62AF940D092099E30F2EF022CF73"
```

### Step 4: Add Hostname to Container App

```bash
az containerapp hostname add \
  --name ccp4i2-bicep-web \
  --resource-group "$RESOURCE_GROUP" \
  --hostname ddudatabase.ncl.ac.uk
```

Expected output:
```json
[
  {
    "bindingType": "Disabled",
    "name": "ddudatabase.ncl.ac.uk"
  }
]
```

### Step 5: Bind Managed Certificate

```bash
az containerapp hostname bind \
  --name ccp4i2-bicep-web \
  --resource-group "$RESOURCE_GROUP" \
  --hostname ddudatabase.ncl.ac.uk \
  --environment ccp4i2-bicep-env-uk \
  --validation-method CNAME
```

This command:
1. Validates domain ownership via the TXT record
2. Provisions a free managed SSL certificate
3. Binds the certificate to your hostname
4. Enables HTTPS traffic

Certificate provisioning typically takes 1-5 minutes.

### Step 6: Verify Configuration

```bash
# List all hostnames and their binding status
az containerapp hostname list \
  --name ccp4i2-bicep-web \
  --resource-group "$RESOURCE_GROUP" \
  -o table

# Test HTTPS access
curl -I https://ddudatabase.ncl.ac.uk
```

Expected hostname list output:
```
Name                      BindingType
------------------------  -------------
ddudatabase.ncl.ac.uk     SniEnabled
```

### Step 7: Update CORS Configuration (if needed)

If your application requires CORS, update the environment variables to include the custom domain:

```bash
# Add to .env.deployment
CUSTOM_DOMAIN=ddudatabase.ncl.ac.uk
```

Then redeploy the application to pick up the new CORS origins.

---

## Certificate Management

### Azure Managed Certificates

Azure Container Apps provides **free managed certificates** that:
- Auto-renew before expiration (typically 60 days before)
- Require no manual intervention
- Support TLS 1.2 and 1.3

### Checking Certificate Status

```bash
# View certificate details via Azure Portal
# Navigate to: Container Apps > ccp4i2-bicep-web > Custom domains

# Or via CLI - list managed certificates in the environment
az containerapp env certificate list \
  --name ccp4i2-bicep-env-uk \
  --resource-group "$RESOURCE_GROUP" \
  -o table
```

### Certificate Expiration Monitoring

Set up Azure Monitor alerts for certificate expiration:

```bash
# Create an action group for notifications (one-time setup)
az monitor action-group create \
  --name ccp4i2-cert-alerts \
  --resource-group "$RESOURCE_GROUP" \
  --short-name CertAlert \
  --email-receiver name=Admin email=nmemn@newcastle.ac.uk
```

### Using Your Own Certificate (Optional)

If you need to use a custom certificate (e.g., EV certificate for branding):

```bash
# Upload PFX certificate to the environment
az containerapp env certificate upload \
  --name ccp4i2-bicep-env-uk \
  --resource-group "$RESOURCE_GROUP" \
  --certificate-file /path/to/certificate.pfx \
  --certificate-name my-custom-cert \
  --password "<pfx-password>"

# Bind the custom certificate to hostname
az containerapp hostname bind \
  --name ccp4i2-bicep-web \
  --resource-group "$RESOURCE_GROUP" \
  --hostname ddudatabase.ncl.ac.uk \
  --environment ccp4i2-bicep-env-uk \
  --certificate my-custom-cert
```

**Important**: Custom certificates require manual renewal before expiration.

---

## Routine Maintenance

### Monthly Health Check

Run these checks monthly to ensure domain/certificate health:

```bash
#!/bin/bash
# Monthly domain health check script

. ./Docker/azure-uksouth/.env.deployment

echo "=== Domain Health Check ==="
echo "Date: $(date)"
echo ""

echo "1. Checking DNS resolution..."
nslookup ddudatabase.ncl.ac.uk

echo ""
echo "2. Checking hostname binding..."
az containerapp hostname list \
  --name ccp4i2-bicep-web \
  --resource-group "$RESOURCE_GROUP" \
  -o table

echo ""
echo "3. Checking certificate status..."
az containerapp env certificate list \
  --name ccp4i2-bicep-env-uk \
  --resource-group "$RESOURCE_GROUP" \
  -o table

echo ""
echo "4. Testing HTTPS connectivity..."
curl -sI https://ddudatabase.ncl.ac.uk | head -5

echo ""
echo "5. Checking certificate expiration..."
echo | openssl s_client -servername ddudatabase.ncl.ac.uk -connect ddudatabase.ncl.ac.uk:443 2>/dev/null | openssl x509 -noout -dates
```

### Before Major Changes

Before making infrastructure changes that might affect the domain:

1. **Document current state**
   ```bash
   az containerapp hostname list --name ccp4i2-bicep-web --resource-group "$RESOURCE_GROUP" -o json > hostname-backup.json
   ```

2. **Note the verification ID** (needed if re-adding domain)
   ```bash
   az containerapp env show --name ccp4i2-bicep-env-uk --resource-group "$RESOURCE_GROUP" \
     --query "properties.customDomainConfiguration.customDomainVerificationId" -o tsv
   ```

3. **Keep DNS records** - Don't remove CNAME or TXT records during maintenance

---

## Troubleshooting

### Common Issues

#### 1. "Domain ownership validation failed"

**Cause**: TXT record not found or incorrect

**Solution**:
```bash
# Verify TXT record exists
dig asuid.ddudatabase.ncl.ac.uk TXT +short

# Compare with expected verification ID
az containerapp env show --name ccp4i2-bicep-env-uk --resource-group "$RESOURCE_GROUP" \
  --query "properties.customDomainConfiguration.customDomainVerificationId" -o tsv
```

If mismatched, request IT to update the TXT record.

#### 2. "Certificate provisioning failed"

**Cause**: DNS not pointing to Container App or CAA records blocking issuance

**Solution**:
```bash
# Verify CNAME points to correct target
dig ddudatabase.ncl.ac.uk CNAME +short
# Should return: ccp4i2-bicep-web.ambitiousrock-87fff394.uksouth.azurecontainerapps.io

# Check for CAA records that might block Azure
dig ddudatabase.ncl.ac.uk CAA +short
# If CAA exists, ensure it allows: digicert.com, microsoft.com
```

#### 3. "Hostname already exists"

**Cause**: Hostname bound to another app or leftover from previous config

**Solution**:
```bash
# Remove existing hostname binding
az containerapp hostname delete \
  --name ccp4i2-bicep-web \
  --resource-group "$RESOURCE_GROUP" \
  --hostname ddudatabase.ncl.ac.uk \
  --yes

# Re-add hostname
az containerapp hostname add ...
```

#### 4. SSL Certificate Errors in Browser

**Cause**: Certificate not yet provisioned or propagated

**Solution**:
```bash
# Check certificate binding status
az containerapp hostname list --name ccp4i2-bicep-web --resource-group "$RESOURCE_GROUP" -o table

# BindingType should be "SniEnabled", not "Disabled"
# If "Disabled", re-run the bind command
```

#### 5. Mixed Content Warnings

**Cause**: Application making HTTP requests instead of HTTPS

**Solution**: Update application configuration to use HTTPS URLs, especially for:
- API endpoints
- Asset URLs
- WebSocket connections

---

## Disaster Recovery

### Scenario 1: Container App Deleted/Recreated

If the Container App is deleted and recreated:

1. **DNS records remain valid** - No changes needed to CNAME or TXT records

2. **Re-add hostname**:
   ```bash
   az containerapp hostname add \
     --name ccp4i2-bicep-web \
     --resource-group "$RESOURCE_GROUP" \
     --hostname ddudatabase.ncl.ac.uk
   ```

3. **Re-bind certificate**:
   ```bash
   az containerapp hostname bind \
     --name ccp4i2-bicep-web \
     --resource-group "$RESOURCE_GROUP" \
     --hostname ddudatabase.ncl.ac.uk \
     --environment ccp4i2-bicep-env-uk \
     --validation-method CNAME
   ```

### Scenario 2: Environment Deleted/Recreated

If the Container Apps Environment is recreated:

1. **Get new verification ID**:
   ```bash
   az containerapp env show --name ccp4i2-bicep-env-uk --resource-group "$RESOURCE_GROUP" \
     --query "properties.customDomainConfiguration.customDomainVerificationId" -o tsv
   ```

2. **Request IT to update TXT record** with new verification ID

3. **Wait for DNS propagation** (5-60 minutes)

4. **Add and bind hostname** (same as initial setup, Steps 4-5)

### Scenario 3: Complete Region Failure

If UK South region is unavailable:

1. **Deploy to new region** using existing Bicep templates
   ```bash
   # Update .env.deployment with new region settings
   # Run infrastructure deployment
   ./Docker/azure-uksouth/scripts/deploy-infrastructure.sh
   ./Docker/azure-uksouth/scripts/deploy-applications.sh
   ```

2. **Get new Container App FQDN**:
   ```bash
   az containerapp show --name ccp4i2-bicep-web --resource-group "$NEW_RESOURCE_GROUP" \
     --query "properties.configuration.ingress.fqdn" -o tsv
   ```

3. **Update DNS CNAME** to point to new FQDN

4. **Get new verification ID and update TXT record**

5. **Add and bind hostname** to new Container App

### Scenario 4: Domain Transfer/Change

If changing from `ddudatabase.ncl.ac.uk` to a new domain:

1. **Remove old hostname**:
   ```bash
   az containerapp hostname delete \
     --name ccp4i2-bicep-web \
     --resource-group "$RESOURCE_GROUP" \
     --hostname ddudatabase.ncl.ac.uk \
     --yes
   ```

2. **Request IT to create new DNS records** for new domain

3. **Add new hostname** following initial setup process

4. **Update application configuration** (CORS, redirects, etc.)

5. **Optionally keep old domain** as redirect:
   - Keep old CNAME pointing to app
   - Add both hostnames to Container App
   - Configure application-level redirect from old to new domain

---

## Quick Reference Commands

```bash
# Source environment (always run first)
. ./Docker/azure-uksouth/.env.deployment

# Get verification ID for DNS TXT record
az containerapp env show --name ccp4i2-bicep-env-uk --resource-group "$RESOURCE_GROUP" \
  --query "properties.customDomainConfiguration.customDomainVerificationId" -o tsv

# List current hostnames
az containerapp hostname list --name ccp4i2-bicep-web --resource-group "$RESOURCE_GROUP" -o table

# Add hostname
az containerapp hostname add --name ccp4i2-bicep-web --resource-group "$RESOURCE_GROUP" \
  --hostname ddudatabase.ncl.ac.uk

# Bind managed certificate
az containerapp hostname bind --name ccp4i2-bicep-web --resource-group "$RESOURCE_GROUP" \
  --hostname ddudatabase.ncl.ac.uk --environment ccp4i2-bicep-env-uk --validation-method CNAME

# Remove hostname
az containerapp hostname delete --name ccp4i2-bicep-web --resource-group "$RESOURCE_GROUP" \
  --hostname ddudatabase.ncl.ac.uk --yes

# List certificates
az containerapp env certificate list --name ccp4i2-bicep-env-uk --resource-group "$RESOURCE_GROUP" -o table

# Check certificate expiration
echo | openssl s_client -servername ddudatabase.ncl.ac.uk -connect ddudatabase.ncl.ac.uk:443 2>/dev/null | openssl x509 -noout -dates

# Get Container App default FQDN
az containerapp show --name ccp4i2-bicep-web --resource-group "$RESOURCE_GROUP" \
  --query "properties.configuration.ingress.fqdn" -o tsv
```

---

## Related Documentation

- [OPERATIONS.md](./OPERATIONS.md) - General operations guide
- [AZURE_AD_SETUP.md](./AZURE_AD_SETUP.md) - Authentication configuration
- [SECURITY_REVIEW.md](./SECURITY_REVIEW.md) - Security assessment

## Support

For issues with:
- **DNS records**: Contact Newcastle University IT
- **Azure infrastructure**: Contact Azure administrators or raise Azure support ticket
- **Application issues**: Check application logs via `az containerapp logs show`
