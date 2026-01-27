"use client";
import { ReactNode, useEffect, useState } from "react";
import { MsalProvider } from "@azure/msal-react";
import { PublicClientApplication } from "@azure/msal-browser";
import { setTokenGetter, setEmailGetter, setLogoutHandler, clearTokenGetter } from "../utils/auth-token";

const clientId = process.env.NEXT_PUBLIC_AAD_CLIENT_ID || "";
const tenantId = process.env.NEXT_PUBLIC_AAD_TENANT_ID || "";

const msalConfig = {
  auth: {
    clientId,
    authority: `https://login.microsoftonline.com/${tenantId}`,
    // Redirect to /auth/callback which is exempt from middleware auth check
    // This allows the MSAL flow to complete and set the auth-session cookie
    redirectUri: "/auth/callback",
  },
};

const pca = new PublicClientApplication(msalConfig);

/**
 * Set the auth-session cookie via API route.
 * This cookie allows the middleware to gate requests server-side.
 */
async function setAuthSessionCookie(): Promise<void> {
  try {
    await fetch("/api/auth/session", {
      method: "POST",
      credentials: "include",
    });
  } catch (error) {
    console.error("[AUTH] Failed to set auth session cookie:", error);
  }
}

/**
 * Clear the auth-session cookie via API route.
 * Called during logout to ensure middleware gates future requests.
 */
async function clearAuthSessionCookie(): Promise<void> {
  try {
    await fetch("/api/auth/session", {
      method: "DELETE",
      credentials: "include",
    });
  } catch (error) {
    console.error("[AUTH] Failed to clear auth session cookie:", error);
  }
}

/**
 * Handle post-login redirect to the original URL.
 * The return URL is stored in sessionStorage by /auth/login page.
 */
function handleReturnUrlRedirect(): void {
  if (typeof window === "undefined") return;

  const returnUrl = sessionStorage.getItem("auth-return-url");
  if (returnUrl && returnUrl !== "/" && returnUrl !== window.location.pathname) {
    sessionStorage.removeItem("auth-return-url");
    // Use replace to avoid adding to history
    window.location.replace(returnUrl);
  }
}

/**
 * Get the current user's email from MSAL account.
 * Returns null if no account is available.
 */
function getAccountEmail(): string | null {
  const accounts = pca.getAllAccounts();
  if (accounts.length === 0) {
    return null;
  }
  // Account username is typically the email/UPN
  return accounts[0].username || null;
}

/**
 * Get an access token for API calls.
 * Uses the .default scope which requests all configured permissions.
 */
async function getApiAccessToken(): Promise<string | null> {
  const accounts = pca.getAllAccounts();
  if (accounts.length === 0) {
    return null;
  }

  try {
    // Use .default scope to get token for our API
    const response = await pca.acquireTokenSilent({
      scopes: [`${clientId}/.default`],
      account: accounts[0],
    });
    return response.accessToken;
  } catch (error: any) {
    // Try interactive login if silent fails
    try {
      const response = await pca.acquireTokenPopup({
        scopes: [`${clientId}/.default`],
        account: accounts[0],
      });
      return response.accessToken;
    } catch (interactiveError: any) {
      console.error("[AUTH] Token acquisition failed:", interactiveError?.message || interactiveError);
      return null;
    }
  }
}

interface AuthProviderProps {
  children: ReactNode;
}

export default function AuthProvider({ children }: AuthProviderProps) {
  const [initialized, setInitialized] = useState(false);

  useEffect(() => {
    pca
      .initialize()
      .then(() => {
        // Handle redirect promises when app loads (for redirect-based auth flows)
        return pca.handleRedirectPromise();
      })
      .then(async (response) => {
        // If we have authenticated accounts, ensure the session cookie is set
        // Note: The /auth/callback page handles redirect completion and cookie setting,
        // but we also set it here to ensure cookie exists for subsequent page loads
        if (response && response.account) {
          console.log("[AUTH] Login redirect completed");
          // Cookie is set by /auth/callback page, but ensure it's set here too
          await setAuthSessionCookie();
        } else if (pca.getAllAccounts().length > 0) {
          // User already has accounts (session exists), ensure cookie is set
          await setAuthSessionCookie();
        }

        // Set up the token and email getters for API calls
        setTokenGetter(getApiAccessToken);
        setEmailGetter(getAccountEmail);
        // Logout handler clears cookie before MSAL logout
        setLogoutHandler(async () => {
          await clearAuthSessionCookie();
          pca.logoutRedirect();
        });
        setInitialized(true);
      })
      .catch((error) => {
        console.error(
          "MSAL initialization or redirect handling failed:",
          error
        );
        setInitialized(true); // Initialize anyway to prevent blocking
      });

    // Cleanup on unmount
    return () => {
      clearTokenGetter();
    };
  }, []);

  if (!initialized) return null; // or a loading spinner

  return <MsalProvider instance={pca}>{children}</MsalProvider>;
}
