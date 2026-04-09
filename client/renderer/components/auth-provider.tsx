/*
 * Copyright (C) 2025-2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
"use client";
import { ReactNode, useEffect, useState } from "react";
import { MsalProvider } from "@azure/msal-react";
import { PublicClientApplication } from "@azure/msal-browser";
import {
  setTokenGetter,
  setEmailGetter,
  setLogoutHandler,
  clearTokenGetter,
  loadTeamsToken,
  setTeamsTokenRefresher,
  setTeamsToken,
  clearTeamsToken,
} from "../utils/auth-token";
import { getAuthConfig } from "../utils/auth-config";

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
 * Check if running in an iframe (Teams context).
 */
function isRunningInIframe(): boolean {
  if (typeof window === "undefined") return false;
  try {
    return window.self !== window.top;
  } catch {
    return true;
  }
}

interface AuthProviderProps {
  children: ReactNode;
}

export default function AuthProvider({ children }: AuthProviderProps) {
  const [msalInstance, setMsalInstance] = useState<PublicClientApplication | null>(null);
  const [initialized, setInitialized] = useState(false);

  useEffect(() => {
    getAuthConfig().then((config) => {
      const pca = new PublicClientApplication({
        auth: {
          clientId: config.clientId,
          authority: `https://login.microsoftonline.com/${config.tenantId}`,
          redirectUri: "/auth/callback",
        },
      });

      pca
        .initialize()
        .then(() => {
          return pca.handleRedirectPromise();
        })
        .then(async (response) => {
          const hasStoredTeamsToken = loadTeamsToken();

          if (hasStoredTeamsToken && isRunningInIframe()) {
            // Running in Teams with stored token - set up refresher
            const refreshTeamsToken = async (): Promise<string | null> => {
              try {
                const teamsModule = await import("@microsoft/teams-js");
                await Promise.race([
                  teamsModule.app.initialize(),
                  new Promise((_, reject) => setTimeout(() => reject(new Error("Teams init timeout")), 3000))
                ]);
                const token = await teamsModule.authentication.getAuthToken({
                  resources: [`api://${window.location.host}/${config.clientId}`],
                  silent: true,
                });
                setTeamsToken(token, 3600);
                return token;
              } catch (error) {
                console.error("[AUTH] Failed to refresh Teams token:", error);
                return null;
              }
            };
            setTeamsTokenRefresher(refreshTeamsToken);
            await setAuthSessionCookie();
          } else if (response && response.account) {
            await setAuthSessionCookie();
          } else if (pca.getAllAccounts().length > 0) {
            await setAuthSessionCookie();
          }

          // Set up the token and email getters for API calls
          setTokenGetter(async () => {
            const accounts = pca.getAllAccounts();
            if (accounts.length === 0) return null;
            try {
              const resp = await pca.acquireTokenSilent({
                scopes: [`${config.clientId}/.default`],
                account: accounts[0],
              });
              return resp.accessToken;
            } catch {
              try {
                const resp = await pca.acquireTokenPopup({
                  scopes: [`${config.clientId}/.default`],
                  account: accounts[0],
                });
                return resp.accessToken;
              } catch (interactiveError: any) {
                console.error("[AUTH] Token acquisition failed:", interactiveError?.message || interactiveError);
                return null;
              }
            }
          });

          setEmailGetter(() => {
            const accounts = pca.getAllAccounts();
            if (accounts.length === 0) return null;
            return accounts[0].username || null;
          });

          setLogoutHandler(async () => {
            await clearAuthSessionCookie();
            clearTeamsToken();
            if (!isRunningInIframe()) {
              pca.logoutRedirect();
            } else {
              window.location.replace("/auth/login");
            }
          });

          setMsalInstance(pca);
          setInitialized(true);
        })
        .catch((error) => {
          console.error(
            "MSAL initialization or redirect handling failed:",
            error
          );
          setInitialized(true);
        });
    });

    return () => {
      clearTokenGetter();
    };
  }, []);

  if (!initialized || !msalInstance) return null;

  return <MsalProvider instance={msalInstance}>{children}</MsalProvider>;
}
