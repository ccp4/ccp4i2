/**
 * Microsoft Teams SSO (Single Sign-On) utilities
 *
 * This module provides seamless authentication for users accessing the app
 * through Microsoft Teams. When running in Teams, we can use the Teams SDK
 * to get an auth token without any user interaction.
 *
 * The flow:
 * 1. Detect if running in Teams context
 * 2. Initialize Teams SDK
 * 3. Get SSO token from Teams
 * 4. Exchange SSO token for access token (On-Behalf-Of flow)
 *
 * ## Azure AD App Registration Requirements for Teams SSO
 *
 * To enable Teams SSO, you must configure your Azure AD app registration:
 *
 * ### 1. Set Application ID URI
 *    - Go to Azure Portal > App registrations > Your App > Expose an API
 *    - Set Application ID URI to: `api://<your-app-domain>/<client-id>`
 *    - Example: `api://ccp4i2-bicep-web.ambitiousrock-87fff394.uksouth.azurecontainerapps.io/your-client-id`
 *
 * ### 2. Add a Scope
 *    - In "Expose an API", click "Add a scope"
 *    - Scope name: `access_as_user`
 *    - Who can consent: Admins and users
 *    - Admin consent display name: "Access CCP4i2 as user"
 *    - Admin consent description: "Allow Teams to access CCP4i2 on behalf of the user"
 *    - User consent display name: "Access CCP4i2"
 *    - User consent description: "Allow Teams to access CCP4i2 on your behalf"
 *    - State: Enabled
 *
 * ### 3. Pre-authorize Teams Client Applications
 *    - In "Expose an API", click "Add a client application"
 *    - Add these Microsoft Teams client IDs:
 *      - `1fec8e78-bce4-4aaf-ab1b-5451cc387264` (Teams desktop/mobile)
 *      - `5e3ce6c0-2b1f-4285-8d4b-75ee78787346` (Teams web)
 *    - Select the `access_as_user` scope for each
 *
 * ### 4. Update Teams App Manifest
 *    - In your Teams app manifest (manifest.json), add:
 *    ```json
 *    "webApplicationInfo": {
 *      "id": "<your-client-id>",
 *      "resource": "api://<your-app-domain>/<your-client-id>"
 *    }
 *    ```
 *
 * ### 5. (Optional) On-Behalf-Of Flow for Backend
 *    - If you need to call Microsoft Graph or other APIs from your backend
 *    - Add a client secret to your app registration
 *    - Exchange the SSO token for an access token using the OBO flow
 *    - See: https://learn.microsoft.com/en-us/azure/active-directory/develop/v2-oauth2-on-behalf-of-flow
 */

// Dynamic import to avoid bundling Teams SDK in Electron builds
let teamsModule: typeof import("@microsoft/teams-js") | null = null;

/**
 * Lazily load the Teams SDK
 */
async function getTeamsSDK(): Promise<typeof import("@microsoft/teams-js")> {
  if (!teamsModule) {
    teamsModule = await import("@microsoft/teams-js");
  }
  return teamsModule;
}

/**
 * Helper to add timeout to a promise
 */
function withTimeout<T>(promise: Promise<T>, ms: number, errorMessage: string): Promise<T> {
  return Promise.race([
    promise,
    new Promise<T>((_, reject) =>
      setTimeout(() => reject(new Error(errorMessage)), ms)
    ),
  ]);
}

/**
 * Check if we're running inside Microsoft Teams
 * Includes timeout protection to prevent indefinite hangs
 */
export async function isRunningInTeams(): Promise<boolean> {
  try {
    const teams = await getTeamsSDK();

    // Add timeout to prevent hanging if Teams context never responds
    await withTimeout(
      teams.app.initialize(),
      3000,
      "Teams SDK initialization timeout"
    );

    const context = await withTimeout(
      teams.app.getContext(),
      2000,
      "Teams context timeout"
    );

    // Check if we have a valid Teams context
    return !!(context && context.app && context.app.host);
  } catch (error) {
    // Teams SDK initialization failed or timed out - not running in Teams
    console.log("[Teams SSO] Not in Teams context:", error);
    return false;
  }
}

/**
 * Result of Teams SSO authentication attempt
 */
export interface TeamsSSOResult {
  success: boolean;
  token?: string;
  error?: string;
  needsConsent?: boolean;
}

/**
 * Attempt to get an SSO token from Teams
 *
 * This uses the Teams JavaScript SDK to request an authentication token.
 * The token is a JWT that can be exchanged for an access token using
 * the On-Behalf-Of (OBO) flow on the backend.
 *
 * @param clientId - The Azure AD application client ID
 * @returns SSO result with token or error information
 */
export async function getTeamsSSOToken(clientId: string): Promise<TeamsSSOResult> {
  try {
    const teams = await getTeamsSDK();

    // Ensure Teams SDK is initialized
    await teams.app.initialize();

    // Request SSO token from Teams
    // The token will be scoped to our Azure AD app
    const token = await teams.authentication.getAuthToken({
      resources: [`api://${clientId}`],
      silent: true,
    });

    console.log("[Teams SSO] Successfully obtained SSO token");
    return {
      success: true,
      token,
    };
  } catch (error: any) {
    console.error("[Teams SSO] Failed to get SSO token:", error);

    // Check for specific error codes
    if (error.errorCode === "resourceRequiresConsent") {
      return {
        success: false,
        error: "User consent required for SSO",
        needsConsent: true,
      };
    }

    if (error.errorCode === "invalid_resource") {
      return {
        success: false,
        error: "Invalid resource - check Azure AD app configuration",
      };
    }

    return {
      success: false,
      error: error.message || "Unknown SSO error",
    };
  }
}

/**
 * Trigger interactive consent flow in Teams
 *
 * This opens a popup within Teams to allow the user to grant consent.
 * Use this when getTeamsSSOToken returns needsConsent: true.
 *
 * @param clientId - The Azure AD application client ID
 * @param tenantId - The Azure AD tenant ID
 * @returns SSO result after consent
 */
export async function triggerTeamsConsent(
  clientId: string,
  tenantId: string
): Promise<TeamsSSOResult> {
  try {
    const teams = await getTeamsSDK();

    await teams.app.initialize();

    // Use Teams authentication popup for consent
    await teams.authentication.authenticate({
      url: `https://login.microsoftonline.com/${tenantId}/oauth2/v2.0/authorize?` +
        `client_id=${clientId}` +
        `&response_type=token` +
        `&scope=openid%20profile%20${encodeURIComponent(`api://${clientId}/access_as_user`)}` +
        `&redirect_uri=${encodeURIComponent(window.location.origin + "/auth-end")}`,
      width: 600,
      height: 535,
    });

    // After consent, try SSO again
    return getTeamsSSOToken(clientId);
  } catch (error: any) {
    console.error("[Teams SSO] Consent flow failed:", error);
    return {
      success: false,
      error: error.message || "Consent flow failed",
    };
  }
}

/**
 * Get Teams context information
 *
 * Useful for logging and debugging.
 */
export async function getTeamsContext(): Promise<{
  userPrincipalName?: string;
  tenantId?: string;
  teamName?: string;
  channelName?: string;
} | null> {
  try {
    const teams = await getTeamsSDK();
    await teams.app.initialize();
    const context = await teams.app.getContext();

    return {
      userPrincipalName: context.user?.userPrincipalName,
      tenantId: context.user?.tenant?.id,
      teamName: context.team?.displayName,
      channelName: context.channel?.displayName,
    };
  } catch {
    return null;
  }
}

/**
 * Complete SSO flow with fallback
 *
 * Attempts SSO, handles consent if needed, and returns the result.
 * The caller should fall back to MSAL popup auth if this fails.
 *
 * @param clientId - The Azure AD application client ID
 * @param tenantId - The Azure AD tenant ID
 * @returns SSO result
 */
export async function attemptTeamsSSO(
  clientId: string,
  tenantId: string
): Promise<TeamsSSOResult> {
  // First, check if we're in Teams
  const inTeams = await isRunningInTeams();
  if (!inTeams) {
    return {
      success: false,
      error: "Not running in Teams",
    };
  }

  console.log("[Teams SSO] Running in Teams context, attempting SSO...");

  // Try to get SSO token silently
  const ssoResult = await getTeamsSSOToken(clientId);

  if (ssoResult.success) {
    return ssoResult;
  }

  // If consent is needed, trigger consent flow
  if (ssoResult.needsConsent) {
    console.log("[Teams SSO] Consent required, triggering consent flow...");
    return triggerTeamsConsent(clientId, tenantId);
  }

  return ssoResult;
}
