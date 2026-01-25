/**
 * Authentication token provider for API requests.
 *
 * This module provides a way to inject MSAL access tokens into API requests.
 * It uses a singleton pattern to store the token getter function, which is
 * set by the AuthProvider component when the app initializes.
 */

type TokenGetter = () => Promise<string | null>;
type EmailGetter = () => string | null;
type LogoutHandler = () => void;

let tokenGetter: TokenGetter | null = null;
let emailGetter: EmailGetter | null = null;
let logoutHandler: LogoutHandler | null = null;

// Token cache to avoid calling acquireTokenSilent on every request
// This dramatically improves performance for Moorhen which makes many API calls
let cachedToken: string | null = null;
let tokenExpiresAt: number = 0;
// Cache tokens for 4 minutes (tokens typically valid for 5+ minutes)
const TOKEN_CACHE_MS = 4 * 60 * 1000;

/**
 * Set the token getter function.
 * Called by AuthProvider when MSAL is initialized.
 */
export function setTokenGetter(getter: TokenGetter): void {
  if (process.env.NODE_ENV === "development" && process.env.DEBUG_AUTH) {
    console.log("[AUTH-TOKEN] setTokenGetter called - auth provider is now configured");
  }
  tokenGetter = getter;
}

/**
 * Set the email getter function.
 * Called by AuthProvider when MSAL is initialized.
 */
export function setEmailGetter(getter: EmailGetter): void {
  emailGetter = getter;
}

/**
 * Set the logout handler function.
 * Called by AuthProvider when MSAL is initialized.
 */
export function setLogoutHandler(handler: LogoutHandler): void {
  logoutHandler = handler;
}

/**
 * Clear the token getter and cache (for logout).
 */
export function clearTokenGetter(): void {
  tokenGetter = null;
  emailGetter = null;
  logoutHandler = null;
  cachedToken = null;
  tokenExpiresAt = 0;
}

/**
 * Get the current access token.
 * Returns null if no token getter is set or if token acquisition fails.
 *
 * Uses a short-lived cache to avoid calling acquireTokenSilent on every request,
 * which dramatically improves performance for applications like Moorhen that
 * make many API calls in rapid succession.
 */
export async function getAccessToken(): Promise<string | null> {
  if (!tokenGetter) {
    // No warning - this is expected in local dev mode without Azure AD
    return null;
  }

  // Return cached token if still valid
  const now = Date.now();
  if (cachedToken && now < tokenExpiresAt) {
    return cachedToken;
  }

  try {
    const token = await tokenGetter();
    if (process.env.NODE_ENV === "development" && process.env.DEBUG_AUTH) {
      console.log("[AUTH-TOKEN] tokenGetter returned:", token ? `token(${token.length} chars)` : "null");
    }

    // Cache the token
    if (token) {
      cachedToken = token;
      tokenExpiresAt = now + TOKEN_CACHE_MS;
    }

    return token;
  } catch (error) {
    console.error("[AUTH-TOKEN] Failed to get access token:", error);
    return null;
  }
}

/**
 * Check if authentication is configured.
 */
export function isAuthConfigured(): boolean {
  return tokenGetter !== null;
}

/**
 * Get the current user's email from MSAL account info.
 * This is useful for sending to the backend as a fallback when
 * the access token doesn't include email claims.
 */
export function getUserEmail(): string | null {
  if (!emailGetter) {
    return null;
  }
  return emailGetter();
}

/**
 * Logout the current user.
 * Triggers MSAL logout redirect if configured.
 */
export function logout(): void {
  if (logoutHandler) {
    logoutHandler();
  }
}
