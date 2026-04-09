"use client";

import { ReactNode, useEffect, useState, useRef } from "react";
import { useMsal } from "@azure/msal-react";
import { InteractionStatus } from "@azure/msal-browser";
import {
  Box,
  Typography,
  CircularProgress,
  Button,
  Alert,
} from "@mui/material";

// Teams SSO types for dynamic import
type TeamsSSOResult = {
  success: boolean;
  token?: string;
  error?: string;
  needsConsent?: boolean;
};

/**
 * Dynamically attempt Teams SSO - only loads Teams SDK when actually in Teams
 * Returns null if not in Teams or if Teams SSO module fails to load
 */
async function tryTeamsSSO(clientId: string, tenantId: string): Promise<TeamsSSOResult | null> {
  // Only attempt Teams SSO in web builds and when in an iframe
  if (!isRunningInIframe()) {
    return null;
  }

  try {
    // Dynamically import Teams SSO module - this avoids bundling in Electron
    const teamsSSO = await import("../utils/teams-sso");
    const inTeams = await teamsSSO.isRunningInTeams();

    if (!inTeams) {
      return null;
    }

    return await teamsSSO.attemptTeamsSSO(clientId, tenantId);
  } catch (error) {
    return null;
  }
}

/**
 * Detect if the app is running inside an iframe (e.g., Microsoft Teams)
 */
function isRunningInIframe(): boolean {
  try {
    return window.self !== window.top;
  } catch {
    // If we can't access window.top due to cross-origin restrictions, we're in an iframe
    return true;
  }
}

interface RequireAuthProps {
  children: ReactNode;
}

/**
 * RequireAuth - Authentication wrapper for UX purposes
 *
 * SECURITY NOTE: This component is for UX only, NOT for security enforcement.
 * All actual authorization is handled by the backend middleware:
 * - Azure AD JWT validation
 * - Group membership checks (ALLOWED_AZURE_AD_GROUPS)
 * - Token expiration validation
 *
 * This component:
 * 1. Initiates the MSAL login flow if user isn't signed in
 * 2. Shows loading states during authentication
 * 3. Sets up the token getter for API calls
 *
 * If a user bypasses this component, API calls will fail with 401/403 from
 * the backend middleware - that's the actual security boundary.
 */
export default function RequireAuth({ children }: RequireAuthProps) {
  const { instance, accounts, inProgress } = useMsal();
  const [isInitializing, setIsInitializing] = useState(true);
  const [authError, setAuthError] = useState<string | null>(null);

  // Prevent multiple login attempts
  const hasInitialized = useRef(false);

  const attemptLogin = () => {
    setAuthError(null);
    hasInitialized.current = false;
  };

  useEffect(() => {
    // Wait for MSAL to finish any in-progress operations
    if (inProgress !== InteractionStatus.None) {
      return;
    }

    // If no accounts and not yet initialized, trigger login
    if (accounts.length === 0 && !hasInitialized.current) {
      hasInitialized.current = true;
      // Use popup auth when running in an iframe (e.g., Teams) since redirects don't work
      if (isRunningInIframe()) {
        // Try Teams SSO first - this provides seamless auth for Teams users
        const clientId = process.env.NEXT_PUBLIC_AAD_CLIENT_ID || "";
        const tenantId = process.env.NEXT_PUBLIC_AAD_TENANT_ID || "";

        tryTeamsSSO(clientId, tenantId)
          .then((ssoResult) => {
            if (ssoResult?.success && ssoResult.token) {
              return instance.ssoSilent({
                scopes: ["openid", "profile"],
                loginHint: undefined,
              }).catch(() => {
                return instance.loginPopup({ scopes: ["openid", "profile"] });
              });
            } else {
              return instance.loginPopup({ scopes: ["openid", "profile"] });
            }
          })
          .catch((error) => {
            console.error("[Auth] Login failed:", error);
            setAuthError(
              error?.errorCode === "user_cancelled"
                ? "Sign-in was cancelled. Please try again to use CCP4i2."
                : "Sign-in failed. Please try again or contact your administrator."
            );
          });
      } else {
        // Standard redirect flow for normal browser usage
        instance.loginRedirect({ scopes: ["openid", "profile"], redirectUri: "/auth/callback" })
          .catch((error) => {
            console.error("[Auth] Redirect failed:", error);
            setAuthError("Unable to redirect to sign-in. Please refresh the page.");
          });
      }
      return;
    }

    // User is authenticated
    if (accounts.length > 0) {
      setAuthError(null);
      setIsInitializing(false);
    }
  }, [accounts, inProgress, instance]);

  // Show error state with retry
  if (authError) {
    return (
      <Box
        display="flex"
        justifyContent="center"
        alignItems="center"
        minHeight="100vh"
        flexDirection="column"
        gap={2}
        px={3}
      >
        <Alert severity="error" sx={{ maxWidth: 480, width: "100%" }}>
          {authError}
        </Alert>
        <Button variant="contained" onClick={attemptLogin}>
          Try again
        </Button>
      </Box>
    );
  }

  // Show loading while MSAL is working or during initialization
  if (accounts.length === 0 || isInitializing) {
    return (
      <Box
        display="flex"
        justifyContent="center"
        alignItems="center"
        minHeight="100vh"
        flexDirection="column"
        gap={2}
      >
        <CircularProgress size={40} />
        <Typography variant="body1" color="text.secondary">
          {accounts.length === 0 ? "Signing you in..." : "Loading..."}
        </Typography>
      </Box>
    );
  }

  // User is authenticated - render children
  // Backend will enforce actual authorization on API calls
  return <>{children}</>;
}
