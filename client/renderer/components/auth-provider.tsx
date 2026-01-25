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
    // Default to "/" - loginRedirect() in require-auth.tsx overrides this
    // with the current path so users return to the page they were accessing
    redirectUri: "/",
  },
};

const pca = new PublicClientApplication(msalConfig);

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
      .then(() => {
        // Set up the token and email getters for API calls
        setTokenGetter(getApiAccessToken);
        setEmailGetter(getAccountEmail);
        setLogoutHandler(() => pca.logoutRedirect());
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
