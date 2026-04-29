# Multi-Instance Deployment Guide

This guide covers deploying additional CCP4i2 instances that are fully isolated from the main production deployment. Each instance gets its own resource group, database, storage, and (optionally) its own Azure AD app registration.

## Security Hardening (baked into Bicep)

The infrastructure template applies these settings by default:

| Setting | Value | Rationale |
|---------|-------|-----------|
| PostgreSQL `publicNetworkAccess` | `Disabled` | Private endpoint only; no public IPs can reach the DB |
| Private storage `defaultAction` | `Deny` | Only VNet and trusted Azure services |
| Public storage (CCP4 software) | Default Allow | Read-only CCP4 binaries, no user data |
| Key Vault `publicNetworkAccess` | `Disabled` | Private endpoint only |
| TLS minimum | `1.2` | On all storage, DB, Key Vault |
| HTTPS enforcement | Enabled | Container Apps `allowInsecure: false` |
| Log Analytics retention | 90 days | Supports security investigations |
| Storage blob public access | `false` | No anonymous blob access |
| Server ingress | Internal only | Not reachable from internet; behind web proxy |

If an existing instance was deployed before these were added, apply them manually via `az` commands or re-run the infrastructure deployment (incremental).

## Image-lineage discipline

The shared ACR (`ccp4acrshareduk14fb`) holds two parallel image lineages during the materia transition. The rule below is load-bearing for DDUDatabase safety: a build aimed at materia development must be physically incapable of replacing DDU's image.

| Lineage | ACR repos | Source branch | Consumed by |
|---------|-----------|--------------|-------------|
| `ccp4i2/*` (stable) | `ccp4i2/web`, `ccp4i2/server` | `django` | DDUDatabase, kawamura, ccp4i2-demo (current `ccp4i2-bicep-*` apps) |
| `materia/*` (development) | `materia/web`, `materia/server` | `django-sliced` | demo only (`materia-demo-*` apps) |

How the separation is enforced:

1. **Env files** select the lineage via `IMAGE_REPO_WEB` / `IMAGE_REPO_SERVER`. The stable env files (`.env.deployment`, `.env.demo`, `.env.kawamura`) leave these unset, so the default `ccp4i2/{web,server}` applies. Only `.env.demo-materia` sets the materia values.
2. **`--env` is mandatory** in [build-and-push.sh](scripts/build-and-push.sh) and [deploy-applications.sh](scripts/deploy-applications.sh). There is no default, so a script invocation without `--env` errors out instead of silently aiming at DDU.
3. **No `:latest` tag is pushed.** Every consumer references an explicit timestamp tag in its env file, so a stray re-tag in ACR cannot move any live container app.
4. **DDU and kawamura must never be deployed with `.env.demo-materia`.** The container-app prefix in that env file (`materia-demo`) and the resource group (`ccp4i2-demo-rg-uksouth`) make a misfire visible, but the rule is the discipline, not the safeguard.

DDU cutover to materia is an explicit, deferred event: rebuild from the `django-sliced` branch, validate on demo, run data-migration parity tests, then deliberately point `.env.deployment` at the materia lineage inside an agreed maintenance window with rollback ready. (`django-sliced` is the CCP4i2-side branch where the slice-prep work happens during the transition; the materia platform itself lives in a separate repo.)

## Prerequisites

- Azure CLI (`az`) installed and logged in
- Owner of an Azure AD app registration (or ability to create one)
- Permission to create resource groups in the Azure subscription

## 1. App Registration Setup

You can reuse an existing app registration or create a new one. Each instance needs its own app registration if you want independent access control.

### Option A: Create a new app registration

```bash
az ad app create --display-name "ccp4i2-<instance-name>" \
  --sign-in-audience AzureADMultipleOrgs
```

### Option B: Reuse an existing registration

If reusing a registration that was created for another purpose, ensure it is correctly configured (see below).

### Required app registration settings

These settings are on the Azure AD app registration object (not in Bicep/ARM).
Run all commands using the app's **application (client) ID**.

```bash
APP_CLIENT_ID="<your-client-id>"
APP_OBJECT_ID=$(az ad app show --id $APP_CLIENT_ID --query id -o tsv)
```

#### a. Access token version (required)

The backend validates v2.0 JWT tokens. Without this, audience and issuer checks will fail with 401.

```bash
az rest --method PATCH \
  --url "https://graph.microsoft.com/v1.0/applications/$APP_OBJECT_ID" \
  --headers "Content-Type=application/json" \
  --body '{"api": {"requestedAccessTokenVersion": 2}}'
```

#### b. Sign-in audience

- `AzureADMyOrg` — single tenant (only your organisation)
- `AzureADMultipleOrgs` — multi-tenant (any Entra ID organisation)

For demos with external collaborators, use multi-tenant:

```bash
az ad app update --id $APP_CLIENT_ID --sign-in-audience AzureADMultipleOrgs
```

#### c. Groups claim (required for group-based access control)

Enables the `groups` claim in JWT tokens so the backend can enforce group membership.

```bash
az rest --method PATCH \
  --url "https://graph.microsoft.com/v1.0/applications/$APP_OBJECT_ID" \
  --headers "Content-Type=application/json" \
  --body '{
    "groupMembershipClaims": "SecurityGroup",
    "optionalClaims": {
      "accessToken": [{"name": "groups", "essential": false}],
      "idToken": [{"name": "groups", "essential": false}]
    }
  }'
```

#### d. Redirect URIs

Add after deployment, once you know the Container App URL:

```bash
INSTANCE_URL="https://<your-container-app-fqdn>"
az rest --method PATCH \
  --url "https://graph.microsoft.com/v1.0/applications/$APP_OBJECT_ID" \
  --headers "Content-Type=application/json" \
  --body "{
    \"spa\": {
      \"redirectUris\": [
        \"${INSTANCE_URL}/\",
        \"${INSTANCE_URL}/auth/callback\",
        \"http://localhost:3000/\",
        \"http://localhost:3000/auth/callback\"
      ]
    }
  }"
```

### App roles (optional)

The app registration should have `User` and `Task.ReadWrite` roles. If creating a new registration, add them:

```bash
az rest --method PATCH \
  --url "https://graph.microsoft.com/v1.0/applications/$APP_OBJECT_ID" \
  --headers "Content-Type=application/json" \
  --body '{
    "appRoles": [
      {
        "allowedMemberTypes": ["User"],
        "description": "A User with regular privilege",
        "displayName": "User",
        "isEnabled": true,
        "value": "User"
      },
      {
        "allowedMemberTypes": ["User"],
        "description": "ReadWrite access allows to view and amend projects",
        "displayName": "Read and Write",
        "isEnabled": true,
        "value": "Task.ReadWrite"
      }
    ]
  }'
```

## 2. Security Group for Access Control

Create a security group to whitelist users. The backend checks JWT group claims against `ALLOWED_AZURE_AD_GROUPS`.

```bash
az ad group create \
  --display-name "ccp4i2-<instance-name>-users" \
  --mail-nickname "ccp4i2-<instance-name>-users" \
  --description "Users authorised for CCP4i2 <instance-name> instance"
```

Note the group ID from the output — you'll need it for `.env` configuration.

Add yourself:

```bash
MY_OID=$(az ad signed-in-user show --query id -o tsv)
az ad group member add --group "ccp4i2-<instance-name>-users" --member-id "$MY_OID"
```

## 3. Environment Configuration

Create an `.env.<instance>` file based on `.env.deployment`. Key fields to change:

```bash
# Resource Group
RESOURCE_GROUP=ccp4i2-<instance>-rg-uksouth

# AAD Configuration (use the new app registration)
NEXT_PUBLIC_AAD_CLIENT_ID=<app-client-id>
NEXT_PUBLIC_AAD_TENANT_ID=<tenant-id>
AZURE_AD_CLIENT_ID=<app-client-id>
AZURE_AD_TENANT_ID=<tenant-id>

# Access control group
ALLOWED_AZURE_AD_GROUPS=<group-id>
```

The remaining fields (ACR, identity IDs, storage account names) are populated after infrastructure deployment.

## 4. Deploy Infrastructure

```bash
# Create resource group
az group create --name ccp4i2-<instance>-rg-uksouth --location uksouth

# Deploy infrastructure (Standard tier Service Bus)
DEMO_DB_PASSWORD=$(openssl rand -base64 24 | tr -d '/+=' | head -c 24)
az deployment group create \
  --resource-group ccp4i2-<instance>-rg-uksouth \
  --template-file Docker/azure-uksouth/infrastructure/infrastructure.bicep \
  --parameters \
    prefix=ccp4i2-<instance> \
    environment=uk \
    postgresAdminPassword="$DEMO_DB_PASSWORD" \
    skipCcp4Storage=true
```

Note the outputs — you'll need ACR name, Key Vault name, etc. for the env file and app deployment.

## 5. Build and Push Images

The frontend image is instance-agnostic (auth config is loaded at runtime). You can share images across instances. The **preferred pattern** for new instances is to consume an existing shared ACR (e.g. `ccp4acrukbwmx`) rather than provisioning and maintaining a per-instance ACR — see [docs/SHARED_ACR.md](docs/SHARED_ACR.md) for the cross-RG wiring (AcrPull role, private endpoint, DNS zone group). The instructions below describe the self-contained per-instance-ACR path, which is still useful for fully isolated instances.

```bash
# Login to the instance's ACR
az acr login --name <acr-name>

# Build web image (no AAD baked in)
docker build --platform linux/amd64 \
  -t <acr-login-server>/ccp4i2/web:latest \
  -f Docker/client/Dockerfile \
  --build-arg NEXT_PUBLIC_REQUIRE_AUTH=true \
  .

# Build server image
docker build --platform linux/amd64 \
  -t <acr-login-server>/ccp4i2/server:latest \
  -f Docker/server/Dockerfile \
  --build-arg ACR_LOGIN_SERVER=<acr-login-server> \
  --build-arg CCP4_VERSION=ccp4-20251105 \
  --build-arg BASE_IMAGE_NAME=ccp4i2/base-arpwarp \
  .

# Push both
docker push <acr-login-server>/ccp4i2/web:latest
docker push <acr-login-server>/ccp4i2/server:latest
```

## 6. Deploy Container Apps

```bash
az deployment group create \
  --resource-group ccp4i2-<instance>-rg-uksouth \
  --template-file Docker/azure-uksouth/infrastructure/applications.bicep \
  --parameters \
    prefix=ccp4i2-<instance> \
    environment=uk \
    aadClientId=<app-client-id> \
    aadTenantId=<tenant-id> \
    allowedAzureAdGroups=<group-id> \
    imageTagWeb=latest \
    imageTagServer=latest \
    containerAppsEnvironmentId=<from-infra-output> \
    acrLoginServer=<from-infra-output> \
    acrName=<from-infra-output> \
    postgresServerFqdn=<from-infra-output> \
    keyVaultName=<from-infra-output> \
    storageAccountName=<from-infra-output> \
    containerAppsIdentityId=<from-infra-output> \
    containerAppsIdentityPrincipalId=<from-infra-output> \
    containerAppsIdentityClientId=<from-infra-output> \
    platformAdminEmails="you@example.com" \
    skipCcp4Storage=true
```

## 7. Add Redirect URIs

After deployment, get the web URL and add redirect URIs (see step 1d above):

```bash
az containerapp show \
  --name ccp4i2-<instance>-web \
  --resource-group ccp4i2-<instance>-rg-uksouth \
  --query "properties.configuration.ingress.fqdn" -o tsv
```

## 8. Inviting External Users

Use the provided script to invite collaborators:

```bash
./Docker/azure-uksouth/scripts/invite-user.sh colleague@york.ac.uk
```

This sends a B2B guest invitation and adds the user to the access control group. Users from any organisation (including non-Entra ID orgs like Google Workspace) can be invited — they'll authenticate via email one-time passcode.

## 9. Managing Instances with Scripts

The deployment scripts **require** an `--env` flag — there is no implicit default. Forgetting it is the most likely way to deploy to the wrong instance, so it is a hard error.

```bash
# DDUDatabase (production)
./scripts/build-and-push.sh --env .env.deployment web
./scripts/deploy-applications.sh --env .env.deployment server

# Demo (ccp4i2 stable)
./scripts/build-and-push.sh --env .env.demo web
./scripts/deploy-applications.sh --env .env.demo server

# Kawamura
./scripts/build-and-push.sh --env .env.kawamura web

# Materia development → demo only (separate image lineage)
./scripts/build-and-push.sh --env .env.demo-materia web
./scripts/deploy-applications.sh --env .env.demo-materia server
```

See **Image-lineage discipline** below for the rule that protects DDUDatabase from materia development work.

## 10. Future: In-App User Management

The `invite-user.sh` script manages B2B guest invitations from the command line. For a more accessible workflow, the CCP4i2 admin interface could be extended to manage guest users directly through the app.

The Microsoft Graph REST API supports all the required operations:
- `POST /invitations` — invite a user by email
- `POST /groups/{id}/members/$ref` — add to the access control group
- `GET /groups/{id}/members` — list current members
- `DELETE /groups/{id}/members/{userId}/$ref` — remove access

Implementation considerations:
- The Django backend would need a **confidential client credential** (client secret or certificate) or **managed identity with Graph API permissions** to call the Graph API server-side
- Access to user management should be restricted to platform admins (`PLATFORM_ADMIN_EMAILS`)
- The admin interface at `/admin/` already exists — a "Manage Users" section could wrap these Graph API calls
- This would allow instance administrators to invite and manage collaborators without CLI access

## Instance Isolation

Each instance is fully isolated:

| Resource | Isolated per instance |
|----------|----------------------|
| Resource Group | Yes |
| Database (PostgreSQL) | Yes |
| Storage Account | Yes |
| Key Vault | Yes |
| Service Bus | Yes |
| Container Apps | Yes |
| Managed Identity | Yes |
| App Registration | Yes (if using separate) |
| ACR | Yes (created by Bicep) |

A token issued for one instance's app registration will fail audience validation on another instance's server. Managed identities have no cross-instance RBAC roles.

## Tearing Down an Instance

```bash
# Delete the resource group (removes all Azure resources)
az group delete --name ccp4i2-<instance>-rg-uksouth --yes

# Optionally clean up Azure AD objects
az ad group delete --group "ccp4i2-<instance>-users"
az ad app delete --id <app-client-id>
```
