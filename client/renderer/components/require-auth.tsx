"use client";

import { ReactNode, useEffect, useState, useRef } from "react";
import { useMsal } from "@azure/msal-react";
import { InteractionStatus } from "@azure/msal-browser";
import {
  checkTeamsMembership,
  DEFAULT_TEAMS_CONFIG,
} from "../utils/teams-auth";

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
      console.log("[Auth] In iframe but not in Teams context");
      return null;
    }

    console.log("[Auth] Teams context detected, attempting SSO...");
    return await teamsSSO.attemptTeamsSSO(clientId, tenantId);
  } catch (error) {
    console.log("[Auth] Teams SSO not available:", error);
    return null;
  }
}
import {
  Box,
  Typography,
  Paper,
  CircularProgress,
  Alert,
  Button,
} from "@mui/material";
import { Security, Warning, Info } from "@mui/icons-material";
import { useTheme } from "../theme/theme-provider";

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

interface AuthState {
  isChecking: boolean;
  hasAccess: boolean;
  reason: string;
  error?: string;
  requiresAdminConsent?: boolean;
}

export default function RequireAuth({ children }: RequireAuthProps) {
  const { customColors } = useTheme();
  const { instance, accounts, inProgress } = useMsal();
  const [authState, setAuthState] = useState<AuthState>({
    isChecking: true,
    hasAccess: false,
    reason: "Initializing...",
  });

  // Prevent infinite loops during auth checks
  const isCheckingRef = useRef(false);
  const hasInitialized = useRef(false);

  useEffect(() => {
    // Prevent running multiple times or during interactions
    if (isCheckingRef.current || inProgress !== InteractionStatus.None) {
      return;
    }

    if (accounts.length === 0) {
      console.log("No accounts found, initiating login");

      // Use popup auth when running in an iframe (e.g., Teams) since redirects don't work
      if (isRunningInIframe()) {
        // Try Teams SSO first - this provides seamless auth for Teams users
        const clientId = process.env.NEXT_PUBLIC_AAD_CLIENT_ID || "";
        const tenantId = process.env.NEXT_PUBLIC_AAD_TENANT_ID || "";

        tryTeamsSSO(clientId, tenantId)
          .then((ssoResult) => {
            if (ssoResult?.success && ssoResult.token) {
              console.log("Teams SSO successful, token obtained");
              // SSO succeeded - the token can be used directly or exchanged via OBO flow
              // For now, we still need MSAL to handle the session, so trigger a silent login
              // using the SSO token hint
              return instance.ssoSilent({
                scopes: ["openid", "profile"],
                loginHint: undefined, // Teams SSO doesn't give us the hint directly
              }).catch(() => {
                // If ssoSilent fails, fall back to popup
                console.log("SSO silent failed, falling back to popup");
                return instance.loginPopup({ scopes: ["openid", "profile"] });
              });
            } else {
              // Teams SSO not available or failed, use popup
              console.log("Teams SSO not available, using popup login");
              return instance.loginPopup({ scopes: ["openid", "profile"] });
            }
          })
          .then(() => {
            console.log("Login successful");
            hasInitialized.current = false;
          })
          .catch((error) => {
            console.error("Login failed:", error);
          });
      } else {
        // Standard redirect flow for normal browser usage
        console.log("Running in browser, using redirect login");
        instance.loginRedirect({ scopes: ["openid", "profile"], redirectUri: "/" });
      }
      return;
    }

    if (accounts.length > 0 && !hasInitialized.current) {
      console.log("Starting authorization check");
      hasInitialized.current = true;
      checkUserAccess();
    }
  }, [accounts, inProgress, instance]);

  const checkUserAccess = async () => {
    if (isCheckingRef.current) {
      console.log("Auth check already in progress, skipping");
      return;
    }

    try {
      isCheckingRef.current = true;
      setAuthState({
        isChecking: true,
        hasAccess: false,
        reason: "Checking authorization...",
      });

      console.log("=== STARTING AUTHORIZATION CHECK ===");
      const account = accounts[0];
      console.log("Account details:", {
        username: account.username,
        name: account.name,
        tenantId: account.tenantId,
        claims: account.idTokenClaims,
      });

      // Strategy 1: Check for app roles in ID token
      const roles = account.idTokenClaims?.roles as string[] | undefined;
      console.log("App roles found in token:", roles);

      if (roles && roles.includes("User")) {
        console.log("✅ Access granted via app role 'User'");
        setAuthState({
          isChecking: false,
          hasAccess: true,
          reason: "Authorized via app role",
        });
        isCheckingRef.current = false;
        return;
      }

      // Strategy 2: t check
      console.log(
        "No 'User' app role found, checking Teams membership for Newcastle Drug Discovery Unit"
      );

      try {
        const teamsResult = await checkTeamsMembership(
          instance,
          DEFAULT_TEAMS_CONFIG
        );
        console.log("Teams check result:", teamsResult);

        if (teamsResult.hasAccess) {
          console.log("✅ Access granted via Teams membership");
          setAuthState({
            isChecking: false,
            hasAccess: true,
            reason: teamsResult.reason,
          });
          isCheckingRef.current = false;
          return;
        }

        // Teams check failed
        console.log("❌ Access denied - not a member of required team");
        setAuthState({
          isChecking: false,
          hasAccess: false,
          reason: teamsResult.reason,
        });
      } catch (teamsError) {
        console.error("Teams membership check failed:", teamsError);
        setAuthState({
          isChecking: false,
          hasAccess: false,
          reason: `Teams membership check failed: ${teamsError}`,
          error: "teams_check_failed",
        });
      }
    } catch (error) {
      console.error("Authorization check failed:", error);
      setAuthState({
        isChecking: false,
        hasAccess: false,
        reason: `Authorization check failed: ${error}`,
        error: "check_failed",
      });
    } finally {
      isCheckingRef.current = false;
    }
  };

  const handleRetryAuth = () => {
    // Reset the checking flag and allow retry
    isCheckingRef.current = false;
    hasInitialized.current = false;

    setAuthState({
      isChecking: true,
      hasAccess: false,
      reason: "Retrying authorization...",
    });

    checkUserAccess();
  };

  if (accounts.length === 0) {
    return (
      <Box
        display="flex"
        justifyContent="center"
        alignItems="center"
        minHeight="100vh"
      >
        <CircularProgress size={40} />
        <Typography variant="body1" sx={{ ml: 2 }}>
          Signing you in...
        </Typography>
      </Box>
    );
  }

  if (authState.isChecking) {
    return (
      <Box
        display="flex"
        justifyContent="center"
        alignItems="center"
        minHeight="100vh"
      >
        <CircularProgress size={40} />
        <Typography variant="body1" sx={{ ml: 2 }}>
          {authState.reason}
        </Typography>
      </Box>
    );
  }

  if (!authState.hasAccess) {
    return (
      <Box
        display="flex"
        justifyContent="center"
        alignItems="center"
        minHeight="100vh"
        bgcolor={customColors.ui.veryLightGray}
        p={3}
      >
        <Paper
          elevation={3}
          sx={{
            p: 4,
            maxWidth: 600,
            textAlign: "center",
            borderRadius: 2,
          }}
        >
          <Security
            sx={{
              fontSize: 60,
              color: authState.requiresAdminConsent ? "orange" : "error.main",
              mb: 2,
            }}
          />

          <Typography
            variant="h4"
            gutterBottom
            color="text.primary"
            fontWeight="bold"
          >
            {authState.error === "admin_consent_required"
              ? "Permission Required"
              : "Access Denied"}
          </Typography>

          <Typography
            variant="body1"
            color="text.secondary"
            sx={{ mb: 3, lineHeight: 1.6 }}
          >
            {authState.reason}
          </Typography>

          {authState.error === "admin_consent_required" && (
            <Alert severity="info" sx={{ mb: 3, textAlign: "left" }}>
              <Typography variant="body2">
                <strong>What this means:</strong> The application needs your
                permission to access Microsoft Teams data to verify your
                membership. This is a one-time consent that you can grant by
                clicking "Retry Authorization" below.
              </Typography>
              <Typography variant="body2" sx={{ mt: 1 }}>
                <strong>Next steps:</strong> Click "Retry Authorization" and
                when prompted, grant permission for the application to read
                basic Teams information.
              </Typography>
            </Alert>
          )}

          {authState.error === "teams_check_failed" && (
            <Alert severity="warning" sx={{ mb: 3, textAlign: "left" }}>
              <Typography variant="body2">
                <strong>Teams Access Required:</strong> This application
                requires membership in the "Newcastle Drug Discovery Unit"
                Microsoft Team to access the system.
              </Typography>
              <Typography variant="body2" sx={{ mt: 1 }}>
                <strong>If you believe this is an error:</strong> Please contact
                your team administrator or try "Retry Authorization" to refresh
                your Teams membership status.
              </Typography>
            </Alert>
          )}

          {(!authState.error || authState.error === "check_failed") &&
            authState.reason?.includes("not a member") && (
              <Alert severity="info" sx={{ mb: 3, textAlign: "left" }}>
                <Typography variant="body2">
                  <strong>Team Membership Required:</strong> Access to this
                  application is restricted to members of the "Newcastle Drug
                  Discovery Unit" Microsoft Team.
                </Typography>
                <Typography variant="body2" sx={{ mt: 1 }}>
                  <strong>To request access:</strong> Please contact your team
                  administrator to be added to the Newcastle Drug Discovery Unit
                  team in Microsoft Teams.
                </Typography>
              </Alert>
            )}

          <Typography variant="body2" color="text.secondary" sx={{ mb: 3 }}>
            Signed in as: <strong>{accounts[0].username}</strong>
          </Typography>

          <Box
            sx={{
              display: "flex",
              gap: 2,
              justifyContent: "center",
              flexWrap: "wrap",
            }}
          >
            <Button
              variant="outlined"
              onClick={handleRetryAuth}
              startIcon={<Info />}
            >
              Retry Authorization
            </Button>

            <Button
              variant="outlined"
              onClick={() => isRunningInIframe() ? instance.logoutPopup() : instance.logoutRedirect()}
              startIcon={<Warning />}
            >
              Sign Out
            </Button>
          </Box>

          {process.env.NODE_ENV === "development" && (
            <Alert severity="info" sx={{ mt: 3, textAlign: "left" }}>
              <Typography variant="body2">
                <strong>Development Info:</strong>
              </Typography>
              <Typography variant="body2" component="div">
                • Error Code: {authState.error || "access_denied"}
                <br />• Tenant ID: {accounts[0].tenantId}
                <br />• Account Type: {accounts[0].idTokenClaims?.aud}
              </Typography>
            </Alert>
          )}
        </Paper>
      </Box>
    );
  }

  return <>{children}</>;
}
