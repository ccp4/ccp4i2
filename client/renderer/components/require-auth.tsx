"use client";

import { ReactNode, useEffect, useState, useRef } from "react";
import { useMsal } from "@azure/msal-react";
import { InteractionStatus } from "@azure/msal-browser";
import { useRouter } from "next/navigation";
import {
  checkTeamsMembership,
  DEFAULT_TEAMS_CONFIG,
} from "../utils/teams-auth";
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
  const router = useRouter();
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
      console.log("No accounts found, redirecting to login");
      // Always redirect to root "/" which is registered in Azure AD
      // After auth, user will be at root and can navigate to their destination
      instance.loginRedirect({ scopes: ["openid", "profile"], redirectUri: "/" });
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

      // Strategy 2: Teams membership check
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
              onClick={() => instance.logoutRedirect()}
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
