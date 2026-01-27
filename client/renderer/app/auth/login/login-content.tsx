"use client";

import { useEffect, useRef } from "react";
import { useMsal } from "@azure/msal-react";
import { InteractionStatus } from "@azure/msal-browser";
import { useSearchParams } from "next/navigation";

/**
 * Login content that triggers MSAL redirect.
 */
export default function LoginContent() {
  const { instance, inProgress } = useMsal();
  const searchParams = useSearchParams();
  const hasTriggeredLogin = useRef(false);

  console.log("[LOGIN] Component rendered, inProgress:", inProgress, "hasTriggered:", hasTriggeredLogin.current);

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

    // Trigger Azure AD login redirect
    console.log("[LOGIN] Triggering loginRedirect...");
    instance.loginRedirect({
      scopes: ["openid", "profile"],
      redirectUri: "/auth/callback",
    }).catch((error) => {
      console.error("[LOGIN] loginRedirect failed:", error);
      hasTriggeredLogin.current = false; // Allow retry
    });
  }, [instance, inProgress, searchParams]);

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
        <p style={{ color: "#666", margin: 0 }}>Redirecting to sign in...</p>
        <style>{`
          @keyframes spin {
            to { transform: rotate(360deg); }
          }
        `}</style>
      </div>
    </div>
  );
}
