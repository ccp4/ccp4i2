# Azure AD Configuration for CCP4i2

## Current Status: ✅ Configured

Groups claim is configured and Teams authorization is working.

## Configuration Summary

**App Registration:** `386da83f-1bf4-4ad8-b742-79b600e2208b`
**Tenant:** `9c5012c9-b616-44c2-a917-66814fbe3e87`
**Authorized Group:** Newcastle Drug Discovery Unit (`6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d`)

### Token Configuration

| Claim | Type | Tokens |
|-------|------|--------|
| groups | Security groups, Group ID | ID, Access |
| idtyp | App-only token signal | Access |

## Adding New Authorized Groups

To grant access to additional Teams/groups:

1. **Find Group ID:**
   - Azure Portal → Microsoft Entra ID → Groups
   - Search for the group name
   - Copy the Object ID

2. **Update Environment:**
   - Edit `Docker/azure-uksouth/.env.deployment`
   - Add to `ALLOWED_AZURE_AD_GROUPS` (comma-separated):
     ```bash
     ALLOWED_AZURE_AD_GROUPS=6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d,new-group-id
     ```

3. **Redeploy:**
   ```bash
   ./scripts/deploy-applications.sh server
   ./scripts/deploy-applications.sh worker
   ```

## Troubleshooting

### Users Not Authorized Despite Being in Team

1. **User must log out and log back in** to get a new token with groups claim
2. Verify group ID is correct in `.env.deployment`
3. Check Azure AD Token Configuration has groups claim

### Groups Claim Not in Token

If `groups` array is missing from JWT:

1. Azure Portal → App Registrations → Your App → Token Configuration
2. Verify "groups" claim is listed
3. If missing, add: Groups claim → Security groups → Group ID

### Group Overage (>200 Groups)

If a user has >200 group memberships:
- Azure AD uses overage indicator instead of full list
- Middleware returns 403 with specific message
- Solution: Create dedicated smaller access group

## Verifying Configuration

Check deployed environment:
```bash
az containerapp show \
  --name ccp4i2-bicep-server \
  --resource-group ccp4i2-bicep-rg-uksouth \
  --query "properties.template.containers[0].env[?name=='ALLOWED_AZURE_AD_GROUPS'].value" \
  -o tsv
```

Check logs:
```bash
az containerapp logs show \
  --name ccp4i2-bicep-server \
  --resource-group ccp4i2-bicep-rg-uksouth \
  --follow | grep -E "(authorized|denied)"
```

---

## Microsoft Teams SSO Configuration

For seamless authentication when the app is embedded in Microsoft Teams, additional Azure AD configuration is required.

### Why This Is Needed

When the app runs inside Teams (as a tab or personal app), normal browser redirects are blocked by the iframe security boundary. Teams SSO allows the app to obtain an authentication token without user interaction by leveraging the user's existing Teams session.

### Step 1: Set Application ID URI

1. **Azure Portal** → App Registrations → `386da83f-1bf4-4ad8-b742-79b600e2208b`
2. Go to **Expose an API**
3. Click **Set** next to "Application ID URI"
4. Set the URI to:
   ```
   api://ccp4i2-bicep-web.ambitiousrock-87fff394.uksouth.azurecontainerapps.io/386da83f-1bf4-4ad8-b742-79b600e2208b
   ```
   Format: `api://<your-app-domain>/<client-id>`

### Step 2: Add the `access_as_user` Scope

1. In **Expose an API**, click **Add a scope**
2. Configure:
   | Field | Value |
   |-------|-------|
   | Scope name | `access_as_user` |
   | Who can consent | Admins and users |
   | Admin consent display name | Access CCP4i2 as user |
   | Admin consent description | Allow Teams to access CCP4i2 on behalf of the user |
   | User consent display name | Access CCP4i2 |
   | User consent description | Allow Teams to access CCP4i2 on your behalf |
   | State | Enabled |

### Step 3: Pre-authorize Teams Client Applications

1. In **Expose an API**, click **Add a client application**
2. Add **Teams desktop/mobile** client:
   - Client ID: `1fec8e78-bce4-4aaf-ab1b-5451cc387264`
   - Select the `access_as_user` scope
3. Add **Teams web** client:
   - Client ID: `5e3ce6c0-2b1f-4285-8d4b-75ee78787346`
   - Select the `access_as_user` scope

### Step 4: Create Teams App Manifest

Create a Teams app package with a `manifest.json` containing:

```json
{
  "$schema": "https://developer.microsoft.com/en-us/json-schemas/teams/v1.16/MicrosoftTeams.schema.json",
  "manifestVersion": "1.16",
  "version": "1.0.0",
  "id": "386da83f-1bf4-4ad8-b742-79b600e2208b",
  "packageName": "com.ccp4.ccp4i2",
  "developer": {
    "name": "CCP4",
    "websiteUrl": "https://www.ccp4.ac.uk",
    "privacyUrl": "https://www.ccp4.ac.uk/privacy",
    "termsOfUseUrl": "https://www.ccp4.ac.uk/terms"
  },
  "name": {
    "short": "CCP4i2",
    "full": "CCP4i2 Crystallographic Interface"
  },
  "description": {
    "short": "Crystallographic computing environment",
    "full": "CCP4i2 provides an environment for crystallographic computing including structure determination and refinement."
  },
  "icons": {
    "outline": "outline.png",
    "color": "color.png"
  },
  "accentColor": "#1976d2",
  "staticTabs": [
    {
      "entityId": "ccp4i2-home",
      "name": "CCP4i2",
      "contentUrl": "https://ccp4i2-bicep-web.ambitiousrock-87fff394.uksouth.azurecontainerapps.io",
      "websiteUrl": "https://ccp4i2-bicep-web.ambitiousrock-87fff394.uksouth.azurecontainerapps.io",
      "scopes": ["personal"]
    }
  ],
  "permissions": ["identity"],
  "validDomains": [
    "ccp4i2-bicep-web.ambitiousrock-87fff394.uksouth.azurecontainerapps.io",
    "login.microsoftonline.com"
  ],
  "webApplicationInfo": {
    "id": "386da83f-1bf4-4ad8-b742-79b600e2208b",
    "resource": "api://ccp4i2-bicep-web.ambitiousrock-87fff394.uksouth.azurecontainerapps.io/386da83f-1bf4-4ad8-b742-79b600e2208b"
  }
}
```

**Critical:** The `webApplicationInfo` section must match exactly:
- `id` = Your Azure AD app client ID
- `resource` = The Application ID URI from Step 1

### Step 5: Deploy Teams App

1. Create app package:
   ```bash
   # Create a folder with manifest.json + icons
   mkdir teams-app
   cp manifest.json teams-app/
   cp outline.png color.png teams-app/
   cd teams-app && zip -r ../ccp4i2-teams.zip *
   ```

2. Upload to Teams Admin Center:
   - Go to https://admin.teams.microsoft.com
   - Navigate to **Teams apps** → **Manage apps**
   - Click **Upload new app** → **Upload an app to your org's app catalog**
   - Upload `ccp4i2-teams.zip`

3. Or sideload for testing:
   - In Teams, click **Apps** → **Manage your apps**
   - Click **Upload an app** → **Upload a custom app**

### Troubleshooting Teams SSO

#### "Redirecting to sign in..." hangs indefinitely
- Verify `webApplicationInfo` in manifest matches Azure AD config
- Check browser console for `[LOGIN]` or `[Teams SSO]` log messages
- Ensure Teams client app IDs are pre-authorized in Azure AD

#### "User consent required" error
- User needs to consent to the app scopes
- The app will show a consent popup - user must approve

#### "Invalid resource" error
- Application ID URI doesn't match manifest `resource` field
- Double-check Step 1 configuration

#### Teams SSO times out, falls back to popup
- This is expected behavior if Azure AD isn't fully configured
- Popup auth will still work, just requires user interaction
