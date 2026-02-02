"use client";

import { useEffect, useRef, useState } from "react";
import { useMsal } from "@azure/msal-react";
import { InteractionStatus } from "@azure/msal-browser";
import { useSearchParams } from "next/navigation";

/**
 * Detect if the app is running inside an iframe (e.g., Microsoft Teams)
 */
function isRunningInIframe(): boolean {
  if (typeof window === "undefined") return false;
  try {
    return window.self !== window.top;
  } catch {
    // If we can't access window.top due to cross-origin restrictions, we're in an iframe
    return true;
  }
}

/**
 * Try Teams authentication using the Teams SDK authentication dialog.
 * This works even when popups are blocked because Teams handles it specially.
 */
async function tryTeamsAuth(clientId: string, tenantId: string): Promise<{
  success: boolean;
  token?: string;
  error?: string;
  stage?: string;
}> {
  try {
    const teamsModule = await import("@microsoft/teams-js");

    // Initialize Teams SDK with timeout
    await Promise.race([
      teamsModule.app.initialize(),
      new Promise((_, reject) => setTimeout(() => reject(new Error("Teams init timeout")), 3000))
    ]);

    // Check if we have a Teams context
    const context = await Promise.race([
      teamsModule.app.getContext(),
      new Promise<null>((_, reject) => setTimeout(() => reject(new Error("Teams context timeout")), 2000))
    ]);

    if (!context?.app?.host) {
      return { success: false, error: "Not in Teams context", stage: "context" };
    }

    // Try silent SSO first
    try {
      const token = await teamsModule.authentication.getAuthToken({
        resources: [`api://${window.location.host}/${clientId}`],
        silent: true,
      });
      return { success: true, token, stage: "sso-silent" };
    } catch (ssoError: any) {
      console.log("[LOGIN] Silent SSO failed, trying Teams auth dialog:", ssoError?.message);
    }

    // Fall back to Teams authentication dialog (not a browser popup)
    // Teams requires the auth URL to start from our domain, then redirect to Azure AD
    const startUrl = `${window.location.origin}/auth/teams-start?client_id=${clientId}&tenant_id=${tenantId}`;

    console.log("[LOGIN] Teams auth start URL:", startUrl);

    try {
      const result = await teamsModule.authentication.authenticate({
        url: startUrl,
        width: 600,
        height: 535,
      });
      return { success: true, token: result, stage: "teams-dialog" };
    } catch (dialogError: any) {
      // Return more specific error about the dialog failure
      return {
        success: false,
        error: `Dialog: ${dialogError?.message || dialogError}`,
        stage: "teams-dialog"
      };
    }
  } catch (error: any) {
    return {
      success: false,
      error: error?.message || "Teams auth failed",
      stage: "error"
    };
  }
}

/**
 * Login content that handles authentication.
 * - In normal browser: Uses redirect-based auth
 * - In iframe/Teams: Uses popup-based auth (redirects don't work in iframes)
 */
export default function LoginContent() {
  const { instance, inProgress } = useMsal();
  const searchParams = useSearchParams();
  const hasTriggeredLogin = useRef(false);
  const [statusMessage, setStatusMessage] = useState("Redirecting to sign in...");
  const [debugInfo, setDebugInfo] = useState<string | null>(null);
  const [isInIframe] = useState(() => isRunningInIframe());

  console.log("[LOGIN] Component rendered, inProgress:", inProgress, "hasTriggered:", hasTriggeredLogin.current, "inIframe:", isInIframe);

  useEffect(() => {
    console.log("[LOGIN] useEffect, inProgress:", inProgress);

    // Wait for MSAL to be ready
    if (inProgress !== InteractionStatus.None) {
      console.log("[LOGIN] MSAL busy, waiting...");
      return;
    }

    // Check if we already have accounts (shouldn't be on login page)
    const accounts = instance.getAllAccounts();
    if (accounts.length > 0) {
      console.log("[LOGIN] Already have accounts, setting cookie and redirecting");
      // User is already logged in, set cookie and redirect
      fetch("/api/auth/session", { method: "POST", credentials: "include" })
        .then(() => {
          const returnUrl = searchParams.get("returnUrl") || "/";
          window.location.replace(returnUrl);
        });
      return;
    }

    // Only trigger login once
    if (hasTriggeredLogin.current) {
      console.log("[LOGIN] Already triggered login, skipping");
      return;
    }
    hasTriggeredLogin.current = true;

    // Get the return URL from query params (set by middleware)
    const returnUrl = searchParams.get("returnUrl") || "/";
    console.log("[LOGIN] Storing returnUrl:", returnUrl);

    // Store return URL in session storage for post-login redirect
    sessionStorage.setItem("auth-return-url", returnUrl);

    // Handle authentication based on context
    if (isInIframe) {
      // Running in iframe (e.g., Teams) - redirects don't work
      console.log("[LOGIN] Running in iframe, attempting Teams auth...");
      setStatusMessage("Detecting Teams environment...");

      const clientId = process.env.NEXT_PUBLIC_AAD_CLIENT_ID || "";
      const tenantId = process.env.NEXT_PUBLIC_AAD_TENANT_ID || "";

      tryTeamsAuth(clientId, tenantId)
        .then((result) => {
          if (result.success) {
            console.log("[LOGIN] Teams auth successful via:", result.stage);
            setStatusMessage("Authenticated, loading...");
            // Set session cookie and redirect
            return fetch("/api/auth/session", { method: "POST", credentials: "include" })
              .then(() => {
                window.location.replace(returnUrl);
              });
          } else {
            // Teams auth failed - show detailed error for debugging
            console.error("[LOGIN] Teams auth failed:", result);
            const errorDetail = `${result.error || "Unknown error"} (stage: ${result.stage || "unknown"})`;
            setStatusMessage(`Sign in failed: ${errorDetail}`);
            setDebugInfo(`Origin: ${window.location.origin}`);
            hasTriggeredLogin.current = false;
          }
        })
        .catch((error: Error) => {
          console.error("[LOGIN] Teams auth exception:", error);
          setStatusMessage(`Sign in error: ${error.message}`);
          hasTriggeredLogin.current = false;
        });
    } else {
      // Normal browser - use redirect flow
      console.log("[LOGIN] Running in browser, using redirect auth...");
      instance.loginRedirect({
        scopes: ["openid", "profile"],
        redirectUri: "/auth/callback",
      }).catch((error) => {
        console.error("[LOGIN] loginRedirect failed:", error);
        hasTriggeredLogin.current = false; // Allow retry
      });
    }
  }, [instance, inProgress, searchParams, isInIframe]);

  // Minimal UI shown briefly before Azure AD redirect
  return (
    <div
      style={{
        display: "flex",
        justifyContent: "center",
        alignItems: "center",
        minHeight: "100vh",
        fontFamily: "system-ui, sans-serif",
        backgroundColor: "#f5f5f5",
      }}
    >
      <div style={{ textAlign: "center" }}>
        <div
          style={{
            width: 40,
            height: 40,
            border: "3px solid #e0e0e0",
            borderTopColor: "#1976d2",
            borderRadius: "50%",
            animation: "spin 1s linear infinite",
            margin: "0 auto 16px",
          }}
        />
        <p style={{ color: "#666", margin: 0 }}>{statusMessage}</p>
        {debugInfo && (
          <p style={{ color: "#999", margin: "8px 0 0", fontSize: "12px" }}>{debugInfo}</p>
        )}
        <style>{`
          @keyframes spin {
            to { transform: rotate(360deg); }
          }
        `}</style>
      </div>
    </div>
  );
}
