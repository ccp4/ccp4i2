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
 * Try Teams SSO with timeout protection
 * Returns null if not in Teams or if Teams SSO fails
 */
async function tryTeamsSSO(clientId: string, tenantId: string): Promise<{
  success: boolean;
  token?: string;
  error?: string;
} | null> {
  // Only attempt in iframe context
  if (!isRunningInIframe()) {
    return null;
  }

  try {
    // Add timeout to prevent indefinite hang
    const timeoutPromise = new Promise<null>((_, reject) =>
      setTimeout(() => reject(new Error("Teams SSO timeout")), 5000)
    );

    const ssoPromise = (async () => {
      const teamsSSO = await import("../../../utils/teams-sso");
      const inTeams = await teamsSSO.isRunningInTeams();

      if (!inTeams) {
        console.log("[LOGIN] In iframe but not Teams context");
        return null;
      }

      console.log("[LOGIN] Teams context detected, attempting SSO...");
      return await teamsSSO.attemptTeamsSSO(clientId, tenantId);
    })();

    return await Promise.race([ssoPromise, timeoutPromise]);
  } catch (error) {
    console.log("[LOGIN] Teams SSO not available or timed out:", error);
    return null;
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
      // Running in iframe (e.g., Teams) - redirects don't work, use popup
      console.log("[LOGIN] Running in iframe, attempting Teams SSO or popup auth...");
      setStatusMessage("Signing you in...");

      const clientId = process.env.NEXT_PUBLIC_AAD_CLIENT_ID || "";
      const tenantId = process.env.NEXT_PUBLIC_AAD_TENANT_ID || "";

      tryTeamsSSO(clientId, tenantId)
        .then((ssoResult) => {
          if (ssoResult?.success && ssoResult.token) {
            console.log("[LOGIN] Teams SSO successful");
            // Try to use the SSO token with MSAL
            return instance.ssoSilent({
              scopes: ["openid", "profile"],
            }).catch(() => {
              console.log("[LOGIN] ssoSilent failed, falling back to popup");
              return instance.loginPopup({ scopes: ["openid", "profile"] });
            });
          } else {
            console.log("[LOGIN] Teams SSO not available, using popup");
            return instance.loginPopup({ scopes: ["openid", "profile"] });
          }
        })
        .then(() => {
          console.log("[LOGIN] Login successful, setting session cookie");
          return fetch("/api/auth/session", { method: "POST", credentials: "include" });
        })
        .then(() => {
          window.location.replace(returnUrl);
        })
        .catch((error) => {
          console.error("[LOGIN] Popup login failed:", error);
          setStatusMessage("Sign in failed. Please try again.");
          hasTriggeredLogin.current = false; // Allow retry
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
        <style>{`
          @keyframes spin {
            to { transform: rotate(360deg); }
          }
        `}</style>
      </div>
    </div>
  );
}
