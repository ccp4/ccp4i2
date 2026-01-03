/**
 * Authentication token provider for API requests.
 *
 * This module provides a way to inject MSAL access tokens into API requests.
 * It uses a singleton pattern to store the token getter function, which is
 * set by the AuthProvider component when the app initializes.
 */

type TokenGetter = () => Promise<string | null>;

let tokenGetter: TokenGetter | null = null;

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
 * Clear the token getter (for logout).
 */
export function clearTokenGetter(): void {
  tokenGetter = null;
}

/**
 * Get the current access token.
 * Returns null if no token getter is set or if token acquisition fails.
 */
export async function getAccessToken(): Promise<string | null> {
  if (!tokenGetter) {
    // No warning - this is expected in local dev mode without Azure AD
    return null;
  }

  try {
    const token = await tokenGetter();
    if (process.env.NODE_ENV === "development" && process.env.DEBUG_AUTH) {
      console.log("[AUTH-TOKEN] tokenGetter returned:", token ? `token(${token.length} chars)` : "null");
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
