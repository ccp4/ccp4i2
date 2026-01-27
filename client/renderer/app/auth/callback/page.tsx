"use client";

import { useEffect, useState, useRef } from "react";
import { useMsal } from "@azure/msal-react";

/**
 * Auth callback page - handles MSAL redirect response.
 *
 * This page is exempt from middleware auth to allow the MSAL flow to complete.
 *
 * IMPORTANT: We must call handleRedirectPromise() ourselves and wait for it,
 * rather than relying on auth-provider, to avoid race conditions.
 */
export default function AuthCallbackPage() {
  const { instance } = useMsal();
  const [status, setStatus] = useState("Completing sign in...");
  const hasProcessed = useRef(false);

  useEffect(() => {
    // Only process once
    if (hasProcessed.current) {
      return;
    }
    hasProcessed.current = true;

    const processAuth = async () => {
      try {
        // Wait for MSAL to process the redirect response
        // This exchanges the auth code for tokens
        const response = await instance.handleRedirectPromise();

        if (response && response.account) {
          console.log("[AUTH CALLBACK] Login successful:", response.account.username);
          setStatus("Setting up session...");

          // Set the auth session cookie
          await fetch("/api/auth/session", {
            method: "POST",
            credentials: "include",
          });

          // Get return URL from sessionStorage
          const returnUrl = sessionStorage.getItem("auth-return-url") || "/";
          sessionStorage.removeItem("auth-return-url");

          console.log("[AUTH CALLBACK] Redirecting to:", returnUrl);
          setStatus("Redirecting...");

          // Small delay to ensure cookie is set
          await new Promise(resolve => setTimeout(resolve, 100));

          // Redirect to the original destination
          window.location.replace(returnUrl);
        } else {
          // Check if we already have accounts (e.g., from a previous session)
          const accounts = instance.getAllAccounts();
          if (accounts.length > 0) {
            console.log("[AUTH CALLBACK] Existing session found");
            setStatus("Session found, setting up...");

            await fetch("/api/auth/session", {
              method: "POST",
              credentials: "include",
            });

            const returnUrl = sessionStorage.getItem("auth-return-url") || "/";
            sessionStorage.removeItem("auth-return-url");
            window.location.replace(returnUrl);
          } else {
            // No response and no accounts - user needs to log in
            console.log("[AUTH CALLBACK] No auth response, redirecting to login");
            setStatus("Please sign in...");
            setTimeout(() => {
              window.location.replace("/auth/login");
            }, 1000);
          }
        }
      } catch (error) {
        console.error("[AUTH CALLBACK] Error processing auth:", error);
        setStatus("Authentication error. Redirecting to login...");
        setTimeout(() => {
          window.location.replace("/auth/login");
        }, 2000);
      }
    };

    processAuth();
  }, [instance]);

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
        <p style={{ color: "#666", margin: 0 }}>{status}</p>
        <style>{`
          @keyframes spin {
            to { transform: rotate(360deg); }
          }
        `}</style>
      </div>
    </div>
  );
}
