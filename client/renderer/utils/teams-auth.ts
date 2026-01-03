/**
 * Microsoft Teams membership authorization utilities
 */

import { IPublicClientApplication } from "@azure/msal-browser";
import { Client } from "@microsoft/microsoft-graph-client";
import { AuthenticationProvider } from "@microsoft/microsoft-graph-client";

/**
 * Custom authentication provider for Microsoft Graph using MSAL
 */
export class MsalAuthProvider implements AuthenticationProvider {
  constructor(private msalInstance: IPublicClientApplication) {}

  async getAccessToken(): Promise<string> {
    const accounts = this.msalInstance.getAllAccounts();
    if (accounts.length === 0) {
      throw new Error("No accounts found");
    }

    // Since portal shows Teams permission doesn't require admin consent,
    // we can try to get it directly for Graph API calls
    try {
      console.log("Trying to acquire token with Teams permissions");
      const response = await this.msalInstance.acquireTokenSilent({
        scopes: ["User.Read", "Team.ReadBasic.All"],
        account: accounts[0],
      });
      console.log("Successfully acquired token with Teams permissions");
      return response.accessToken;
    } catch (silentError) {
      console.log(
        "Silent token acquisition failed, trying interactive consent:",
        silentError
      );

      try {
        // Try interactive consent with Teams permissions using redirect
        await this.msalInstance.acquireTokenRedirect({
          scopes: ["User.Read", "Team.ReadBasic.All"],
          prompt: "consent", // Force consent dialog
        });
        // Note: This will redirect, so we won't return from here
        // The token will be available after redirect completes
        throw new Error("Redirect initiated for Teams consent");
      } catch (interactiveError) {
        console.error(
          "Interactive Teams token acquisition failed:",
          interactiveError
        );

        // Fallback to basic User.Read if Teams fails
        console.log("Falling back to basic User.Read scope");
        try {
          const basicResponse = await this.msalInstance.acquireTokenSilent({
            scopes: ["User.Read"],
            account: accounts[0],
          });
          console.log("Fallback to basic scope successful");
          return basicResponse.accessToken;
        } catch (basicError) {
          console.error("Even basic token acquisition failed:", basicError);
          throw new Error(
            `All token acquisition attempts failed: ${basicError}`
          );
        }
      }
    }
  }

  /**
   * Attempt to get Teams-capable token
   * Since portal shows "Admin Consent Required: No", we can try user consent
   */
  async getTeamsAccessToken(): Promise<string> {
    const accounts = this.msalInstance.getAllAccounts();
    if (accounts.length === 0) {
      throw new Error("No accounts found");
    }

    try {
      // Try Teams permissions silently first
      console.log("Attempting to acquire Teams permissions silently");
      const response = await this.msalInstance.acquireTokenSilent({
        scopes: ["User.Read", "Team.ReadBasic.All"],
        account: accounts[0],
      });
      console.log(
        "Successfully acquired token with Teams permissions silently"
      );
      return response.accessToken;
    } catch (silentError) {
      console.log("Silent Teams token acquisition failed, error details:", {
        errorCode: silentError.errorCode,
        errorMessage: silentError.errorMessage,
        message: silentError.message,
      });

      // Since admin consent is not required according to portal,
      // this should be a user consent issue that interactive login can resolve
      try {
        console.log(
          "Trying interactive consent for Teams permissions via redirect..."
        );
        await this.msalInstance.acquireTokenRedirect({
          scopes: ["User.Read", "Team.ReadBasic.All"],
          prompt: "consent", // Force consent screen to ensure user grants Teams permission
        });
        // Note: This will redirect, so we won't return from here
        // The token will be available after redirect completes
        throw new Error("Redirect initiated for Teams consent");
      } catch (interactiveError) {
        console.error("Teams interactive token acquisition failed:", {
          errorCode: interactiveError.errorCode,
          errorMessage: interactiveError.errorMessage,
          message: interactiveError.message,
        });

        // More specific error handling
        if (interactiveError.message?.includes("AADSTS65001")) {
          throw new Error(
            "User declined consent for Teams access or additional permissions needed"
          );
        } else if (interactiveError.message?.includes("AADSTS70011")) {
          throw new Error(
            "Invalid scope: Team.ReadBasic.All may not be configured properly in Azure AD"
          );
        } else {
          throw new Error(
            `Teams token acquisition failed: ${interactiveError.message || interactiveError}`
          );
        }
      }
    }
  }
}

/**
 * Configuration for Teams-based authorization
 */
interface TeamsAuthConfig {
  /**
   * List of Team IDs that should have access to the application
   * You can find Team IDs in Microsoft Graph Explorer or Teams Admin Center
   */
  authorizedTeamIds?: string[];

  /**
   * List of Team display names that should have access
   * Less reliable than IDs but easier to configure
   */
  authorizedTeamNames?: string[];

  /**
   * If true, any team membership grants access
   * Useful for organizations where all Teams users should have access
   */
  anyTeamMembership?: boolean;
}

/**
 * Check if user has access based on Microsoft Teams membership
 */
export async function checkTeamsMembership(
  msalInstance: IPublicClientApplication,
  config: TeamsAuthConfig
): Promise<{
  hasAccess: boolean;
  reason: string;
  teamsFound?: Array<{ id: string; displayName: string }>;
  error?: string;
}> {
  try {
    console.log("=== TEAMS MEMBERSHIP CHECK STARTING ===");

    const authProvider = new MsalAuthProvider(msalInstance);

    // Try to get a Teams-capable token
    let accessToken: string;
    try {
      accessToken = await authProvider.getTeamsAccessToken();
      console.log("Teams-capable token acquired successfully");
    } catch (tokenError) {
      console.warn("Teams token acquisition failed:", tokenError);

      // Return helpful information about the consent requirement
      if (
        tokenError.message?.includes("consent") ||
        tokenError.message?.includes("AADSTS65001")
      ) {
        return {
          hasAccess: false,
          reason:
            "Teams access requires your permission to read basic team membership information. Click 'Retry Authorization' to grant this permission.",
          error: "admin_consent_required", // Keep same error code but different messaging
        };
      }

      // For other token errors, fall back to basic access
      return {
        hasAccess: false,
        reason: `Teams authentication failed: ${tokenError.message}`,
        error: "teams_auth_failed",
      };
    }

    const graphClient = Client.initWithMiddleware({ authProvider });

    console.log("Graph client initialized, attempting to call /me/joinedTeams");

    // Try to get user's joined teams
    let teamsResponse;
    try {
      teamsResponse = await graphClient.api("/me/joinedTeams").get();
      console.log("Raw Teams API response:", teamsResponse);
    } catch (apiError: any) {
      console.error("Teams API call failed:", apiError);
      console.error("API Error details:", {
        code: apiError.code,
        message: apiError.message,
        statusCode: apiError.statusCode,
      });

      // Provide specific guidance for common API errors
      if (apiError.code === "Forbidden" || apiError.statusCode === 403) {
        return {
          hasAccess: false,
          reason:
            "Access to Teams data is forbidden. The application may need additional permissions or admin consent.",
          error: "teams_api_forbidden",
        };
      }

      throw apiError;
    }

    const userTeams = teamsResponse.value || [];

    console.log("=== TEAMS MEMBERSHIP DEBUG ===");
    console.log("Teams API call successful!");
    console.log("User's teams count:", userTeams.length);
    console.log("User's teams:", userTeams);
    console.log("Auth config:", config);

    // Check if user has any team membership (if configured)
    if (config.anyTeamMembership && userTeams.length > 0) {
      return {
        hasAccess: true,
        reason: "User is member of at least one Microsoft Team",
        teamsFound: userTeams.map((team: any) => ({
          id: team.id,
          displayName: team.displayName,
        })),
      };
    }

    // Check specific team IDs
    if (config.authorizedTeamIds && config.authorizedTeamIds.length > 0) {
      const matchingTeams = userTeams.filter((team: any) =>
        config.authorizedTeamIds!.includes(team.id)
      );

      if (matchingTeams.length > 0) {
        return {
          hasAccess: true,
          reason: `User is member of authorized team(s): ${matchingTeams
            .map((t: any) => t.displayName)
            .join(", ")}`,
          teamsFound: matchingTeams.map((team: any) => ({
            id: team.id,
            displayName: team.displayName,
          })),
        };
      }
    }

    // Check team display names
    if (config.authorizedTeamNames && config.authorizedTeamNames.length > 0) {
      const matchingTeams = userTeams.filter((team: any) =>
        config.authorizedTeamNames!.some((name) =>
          team.displayName?.toLowerCase().includes(name.toLowerCase())
        )
      );

      if (matchingTeams.length > 0) {
        return {
          hasAccess: true,
          reason: `User is member of authorized team(s): ${matchingTeams
            .map((t: any) => t.displayName)
            .join(", ")}`,
          teamsFound: matchingTeams.map((team: any) => ({
            id: team.id,
            displayName: team.displayName,
          })),
        };
      }
    }

    // No matching teams found
    return {
      hasAccess: false,
      reason: "User is not a member of any authorized Microsoft Teams",
      teamsFound: userTeams.map((team: any) => ({
        id: team.id,
        displayName: team.displayName,
      })),
    };
  } catch (error: any) {
    console.error("Teams membership check failed:", error);
    return {
      hasAccess: false,
      reason: "Failed to check Teams membership",
      error: error.message || "Unknown error",
    };
  }
}

/**
 * Example Teams configuration
 * You can customize this based on your organization's Teams setup
 */
export const DEFAULT_TEAMS_CONFIG: TeamsAuthConfig = {
  // Option 1: Specific team IDs (most secure) - Newcastle Drug Discovery Unit only
  authorizedTeamIds: [
    "6f35cbeb-5f5d-4cf3-9b93-fe0b6eb6306d", // Newcastle Drug Discovery Unit
  ],

  // Option 2: Backup - Team name matching (less secure but more flexible)
  authorizedTeamNames: ["Newcastle Drug Discovery Unit"],

  // Option 3: Disabled - no longer allow any team membership
  anyTeamMembership: false,
};

/**
 * Generate admin consent URL for Teams permissions
 * Administrators can use this URL to grant necessary permissions
 */
export function generateAdminConsentUrl(
  clientId: string,
  tenantId: string,
  redirectUri: string = window.location.origin
): string {
  const scopes = ["User.Read", "Team.ReadBasic.All"];

  const consentUrl = new URL(
    "https://login.microsoftonline.com/common/adminconsent"
  );
  consentUrl.searchParams.set("client_id", clientId);
  consentUrl.searchParams.set("redirect_uri", redirectUri);
  consentUrl.searchParams.set("scope", scopes.join(" "));

  return consentUrl.toString();
}

/**
 * Provide helpful guidance for resolving consent issues
 */
export function getConsentGuidance(clientId: string, tenantId: string) {
  return {
    adminConsentUrl: generateAdminConsentUrl(clientId, tenantId),
    instructions: [
      "Contact your IT administrator or Azure AD administrator",
      "Ask them to grant 'Team.ReadBasic.All' permission to the CCP4i2 application",
      "Provide them with the admin consent URL above",
      "Once granted, all users in your organization will have access to Teams-based authorization",
    ],
    troubleshooting: [
      "If you are an administrator, you can click the admin consent URL to grant permissions",
      "The application needs 'Team.ReadBasic.All' to read basic team information",
      "This is a delegated permission that works on behalf of signed-in users",
      "No sensitive team data is accessed - only basic membership information",
    ],
  };
}
