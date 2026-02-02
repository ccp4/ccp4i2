"use client";

import { useEffect } from "react";

/**
 * Teams authentication callback page.
 *
 * This page is loaded inside the Teams authentication dialog after Azure AD
 * redirects back. It extracts the token from the URL fragment and notifies
 * the parent Teams app of the result.
 */
export default function TeamsCallbackPage() {
  useEffect(() => {
    async function handleCallback() {
      try {
        // Parse the hash fragment for the token
        const hash = window.location.hash.substring(1);
        const params = new URLSearchParams(hash);

        const idToken = params.get("id_token");
        const error = params.get("error");
        const errorDescription = params.get("error_description");

        if (error) {
          console.error("[Teams Callback] Auth error:", error, errorDescription);
          // Notify Teams of failure
          const teams = await import("@microsoft/teams-js");
          await teams.app.initialize();
          teams.authentication.notifyFailure(errorDescription || error);
          return;
        }

        if (idToken) {
          console.log("[Teams Callback] Received id_token, notifying success");
          // Notify Teams of success
          const teams = await import("@microsoft/teams-js");
          await teams.app.initialize();
          teams.authentication.notifySuccess(idToken);
        } else {
          console.error("[Teams Callback] No token in response");
          const teams = await import("@microsoft/teams-js");
          await teams.app.initialize();
          teams.authentication.notifyFailure("No token received");
        }
      } catch (err) {
        console.error("[Teams Callback] Error:", err);
        // Try to notify failure even if Teams SDK fails
        try {
          const teams = await import("@microsoft/teams-js");
          await teams.app.initialize();
          teams.authentication.notifyFailure("Callback processing failed");
        } catch {
          // Nothing more we can do
        }
      }
    }

    handleCallback();
  }, []);

  return (
    <div
      style={{
        display: "flex",
        justifyContent: "center",
        alignItems: "center",
        minHeight: "100vh",
        fontFamily: "system-ui, sans-serif",
      }}
    >
      <p>Completing sign in...</p>
    </div>
  );
}
