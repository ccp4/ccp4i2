"use client";

import { Suspense, useEffect } from "react";
import { useSearchParams } from "next/navigation";

/**
 * Inner component that uses useSearchParams.
 * Must be wrapped in Suspense for Next.js 15 static generation.
 */
function TeamsStartContent() {
  const searchParams = useSearchParams();

  useEffect(() => {
    const clientId = searchParams.get("client_id") || process.env.NEXT_PUBLIC_AAD_CLIENT_ID || "";
    const tenantId = searchParams.get("tenant_id") || process.env.NEXT_PUBLIC_AAD_TENANT_ID || "";
    const redirectUri = `${window.location.origin}/auth/teams-callback`;

    const authUrl = new URL(`https://login.microsoftonline.com/${tenantId}/oauth2/v2.0/authorize`);
    authUrl.searchParams.set("client_id", clientId);
    authUrl.searchParams.set("response_type", "id_token");
    authUrl.searchParams.set("scope", "openid profile email");
    authUrl.searchParams.set("redirect_uri", redirectUri);
    authUrl.searchParams.set("nonce", Math.random().toString(36).substring(2));
    authUrl.searchParams.set("response_mode", "fragment");

    // Redirect to Azure AD
    window.location.replace(authUrl.toString());
  }, [searchParams]);

  return <p>Redirecting to sign in...</p>;
}

/**
 * Teams authentication start page.
 *
 * This page runs on our domain and redirects to Azure AD.
 * Teams requires auth URLs to start from a valid domain in the manifest.
 */
export default function TeamsStartPage() {
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
      <Suspense fallback={<p>Loading...</p>}>
        <TeamsStartContent />
      </Suspense>
    </div>
  );
}
