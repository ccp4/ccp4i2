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

### Custom Domain Support

**Recommended:** Use a custom domain (e.g., `ddudatabase.ncl.ac.uk`) instead of the default Azure Container Apps URL. This provides:
- Stable URLs that don't change with deployments
- Shorter, memorable URLs for users
- Professional branding

If using a custom domain:
1. Configure the custom domain in Azure Container Apps
2. Set `CUSTOM_DOMAIN` in `.env.deployment`
3. Use the custom domain in all configuration below

### Step 1: Set Application ID URI

1. **Azure Portal** → App Registrations → `386da83f-1bf4-4ad8-b742-79b600e2208b`
2. Go to **Expose an API**
3. Click **Set** next to "Application ID URI"
4. Set the URI to match your domain:

   **With custom domain (recommended):**
   ```
   api://ddudatabase.ncl.ac.uk/386da83f-1bf4-4ad8-b742-79b600e2208b
   ```

   **Without custom domain:**
   ```
   api://ccp4i2-bicep-web.ambitiousrock-87fff394.uksouth.azurecontainerapps.io/386da83f-1bf4-4ad8-b742-79b600e2208b
   ```

   Format: `api://<your-app-domain>/<client-id>`

   **Important:** The domain in the URI must exactly match the domain users access the app from.

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

Create a Teams app package with a `manifest.json`. The example below uses a custom domain:

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
      "contentUrl": "https://ddudatabase.ncl.ac.uk",
      "websiteUrl": "https://ddudatabase.ncl.ac.uk",
      "scopes": ["personal"]
    }
  ],
  "permissions": ["identity"],
  "validDomains": [
    "ddudatabase.ncl.ac.uk",
    "login.microsoftonline.com"
  ],
  "webApplicationInfo": {
    "id": "386da83f-1bf4-4ad8-b742-79b600e2208b",
    "resource": "api://ddudatabase.ncl.ac.uk/386da83f-1bf4-4ad8-b742-79b600e2208b"
  }
}
```

**Critical Configuration Requirements:**

| Field | Must Match |
|-------|------------|
| `staticTabs[].contentUrl` | Your app's actual URL (custom domain or Azure Container Apps) |
| `staticTabs[].websiteUrl` | Same as contentUrl |
| `validDomains[]` | Your app's domain (without protocol) |
| `webApplicationInfo.id` | Azure AD app client ID |
| `webApplicationInfo.resource` | Application ID URI from Step 1 (must match exactly) |

**Common Mistake:** If `contentUrl` and `webApplicationInfo.resource` use different domains, you'll get:
```
"App resource defined in manifest and iframe origin do not match"
```

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

### How Teams SSO Token Storage Works

Teams SSO returns tokens differently from standard MSAL browser authentication:
- **Normal browser:** MSAL stores tokens in its internal cache, and `acquireTokenSilent()` retrieves them
- **Teams iframe:** Teams SDK provides tokens directly, but MSAL's account cache isn't populated

The app handles this by:
1. Storing the Teams token in `sessionStorage` (keys: `ccp4i2-teams-token`, `ccp4i2-teams-token-expires`)
2. The `getAccessToken()` function checks for a Teams token first before falling back to MSAL
3. Tokens auto-refresh when they expire using the Teams SDK

**Session Cookie:** The auth-session cookie uses `sameSite: "none"` in production to work within the Teams iframe (third-party context).

**Security:** This is equivalent to MSAL's default storage (also sessionStorage). Tokens are:
- Only accessible from the same origin
- Cleared when the browser tab closes
- Protected by HTTPS in production

**Files involved:**
- `client/renderer/utils/auth-token.ts` - Token storage and retrieval
- `client/renderer/components/auth-provider.tsx` - Token refresh setup
- `client/renderer/app/auth/login/login-content.tsx` - Initial Teams SSO flow

### Troubleshooting Teams SSO

#### "App resource defined in manifest and iframe origin do not match"
- The domain in `webApplicationInfo.resource` doesn't match where the app is running
- Ensure `contentUrl` and `resource` use the same domain
- If using a custom domain, update both to match

#### "Redirecting to sign in..." hangs indefinitely
- Verify `webApplicationInfo` in manifest matches Azure AD config
- Check browser console for `[LOGIN]` or `[AUTH]` log messages
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

#### API calls fail after successful login (401/403)
- Check browser console for `[AUTH-TOKEN]` log messages
- Verify the Teams token includes the `groups` claim (see Token Configuration above)
- The token should be stored in sessionStorage - check DevTools → Application → Session Storage

#### Login works but then redirects back to login (loop)
- Likely the auth-session cookie isn't being set or sent
- In Teams context, cookies must have `sameSite: "none"` and `secure: true`
- Check the cookie exists in DevTools → Application → Cookies

### Updating Teams App After Domain Change

If you change the domain (e.g., from Azure Container Apps URL to custom domain):

1. **Update Azure AD:**
   - Change Application ID URI in "Expose an API"
   - Old URI must be removed before adding new one

2. **Update Teams Manifest:**
   - Update `staticTabs[].contentUrl` and `websiteUrl`
   - Update `validDomains[]`
   - Update `webApplicationInfo.resource`
   - Increment `version` number

3. **Update Teams App:**
   - Download updated manifest from Teams Developer Portal
   - Or re-upload the modified manifest.json

4. **Rebuild and Deploy:**
   ```bash
   ./scripts/build-and-push.sh web
   ./scripts/deploy-applications.sh web
   ```
