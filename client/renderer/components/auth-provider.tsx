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
  hasLocalSessionToken,
  createLocalSessionTokenGetter,
  createLocalSessionEmailGetter,
} from "@ccp4/ccp4i2-api";
import { getAuthConfig } from "../utils/auth-config";
import { setReauthHandler } from "../utils/reauth";

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
 * Re-stamp the auth-session cookie after a successful token acquisition.
 *
 * The cookie is a fixed 8-hour window set once at login, while MSAL's refresh
 * token outlives it. A user still happily working at hour nine had a valid
 * token and an expired cookie, and the middleware bounced their next
 * navigation to /auth/login -- a sign-in prompt caused by nothing but the
 * clock. Every silent acquisition is evidence the session is alive, so it is
 * also the moment to extend the cookie.
 *
 * Throttled: tokens come from a 4-minute cache, but Moorhen can still ask
 * often, and this is a same-origin round trip on the request path.
 */
const SESSION_COOKIE_REFRESH_MS = 5 * 60 * 1000;
let sessionCookieRefreshedAt = 0;

async function keepSessionCookieAlive(): Promise<void> {
  const now = Date.now();
  if (now - sessionCookieRefreshedAt < SESSION_COOKIE_REFRESH_MS) return;
  sessionCookieRefreshedAt = now;
  await setAuthSessionCookie();
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
    // LocalSession path (CCP4i2 desktop): the Electron preload exposed a
    // per-launch token via window.ccp4i2LocalSession. Use that getter and
    // skip MSAL entirely — desktop has no Azure AD account to bind to.
    if (hasLocalSessionToken()) {
      setTokenGetter(createLocalSessionTokenGetter());
      setEmailGetter(createLocalSessionEmailGetter());
      // The desktop session lives as long as the app process, so there is
      // nothing to renew and nothing to sign out of.
      setReauthHandler(null);
      setLogoutHandler(() => {
        // Desktop session lives until the app process dies; logout is a no-op.
        console.log("[AUTH] Local session active; logout is a no-op.");
      });
      setInitialized(true);
      return;
    }

    // MSAL path (cloud build, or Electron without preload-injected token).
    getAuthConfig().then((config) => {
      const pca = new PublicClientApplication({
        auth: {
          clientId: config.clientId,
          authority: `https://login.microsoftonline.com/${config.tenantId}`,
          redirectUri: "/auth/callback",
        },
        cache: {
          // localStorage (not the MSAL default of sessionStorage) so MSAL's
          // account + refresh-token cache survives across tabs/windows of the
          // same origin. Moorhen and Materia open job/file viewers in new
          // windows via window.open(_, "_blank"); with sessionStorage MSAL
          // the new window starts with zero accounts -> tokenGetter returns
          // null -> proxy forwards without Authorization -> Django 401 ->
          // auth-error-handler bounces the user to the AAD "choose account
          // to log out" page. localStorage scope is per-origin, so the new
          // window finds the same cached account.
          cacheLocation: "localStorage",
          storeAuthStateInCookie: false,
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

          // Set up the token and email getters for API calls.
          //
          // Silent-only path with one retry. We deliberately do NOT fall back
          // to acquireTokenPopup here: tokenGetter is invoked from inside
          // fetch wrappers, which are not a user-gesture context, so popup
          // calls would be blocked by the browser anyway. When silent fails
          // we return null; the api-fetch layer turns the resulting 401 into
          // an AUTH_ERROR_EVENT, which AuthErrorHandler surfaces as a
          // snackbar + auto-redirect to sign-in.
          //
          // The single retry catches the common transient case: parallel
          // fetches racing the underlying MSAL refresh, where the second
          // call sees the in-flight refresh and silent throws even though
          // a fresh token is about to land. Discovered when long sequences
          // (Moorhen panel mounts, R-group decomposition on Materia's
          // AggregationPage) triggered "Your session has expired"
          // snackbars: the popup fallback couldn't fire from the fetch
          // path so the request went out tokenless and got a real 401
          // anyway.
          setTokenGetter(async (options) => {
            const accounts = pca.getAllAccounts();
            if (accounts.length === 0) return null;
            const params = {
              scopes: [`${config.clientId}/.default`],
              account: accounts[0],
              // Set when a request has just been refused: go past MSAL's own
              // cache rather than re-presenting the token the server rejected.
              forceRefresh: options?.forceRefresh ?? false,
            };
            try {
              const resp = await pca.acquireTokenSilent(params);
              void keepSessionCookieAlive();
              return resp.accessToken;
            } catch (firstError: any) {
              // Brief delay lets any in-flight refresh on a sibling call
              // finish populating MSAL's account cache before we retry.
              await new Promise((resolve) => setTimeout(resolve, 250));
              try {
                const resp = await pca.acquireTokenSilent(params);
                void keepSessionCookieAlive();
                return resp.accessToken;
              } catch (secondError: any) {
                console.error(
                  "[AUTH] Silent token acquisition failed (after retry):",
                  secondError?.message || secondError
                );
                return null;
              }
            }
          });

          setEmailGetter(() => {
            const accounts = pca.getAllAccounts();
            if (accounts.length === 0) return null;
            return accounts[0].username || null;
          });

          // Signing back IN. MSAL keeps the account, so with a live AAD
          // session this round-trips without a prompt; the callback returns
          // the user to the page they were on. Contrast setLogoutHandler
          // below, which tears the session down.
          setReauthHandler(async () => {
            const accounts = pca.getAllAccounts();
            if (accounts.length === 0) return false;
            await pca.acquireTokenRedirect({
              scopes: [`${config.clientId}/.default`],
              account: accounts[0],
            });
            return true;
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
      setReauthHandler(null);
    };
  }, []);

  if (!initialized) return null;

  // LocalSession mode: render children directly. MsalProvider is not
  // needed because no MSAL hooks should be reached on the desktop auth
  // path (the renderer never sees a login flow).
  if (hasLocalSessionToken()) {
    return <>{children}</>;
  }

  if (!msalInstance) return null;
  return <MsalProvider instance={msalInstance}>{children}</MsalProvider>;
}
