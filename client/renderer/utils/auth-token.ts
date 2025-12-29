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
  console.warn("[AUTH-TOKEN] setTokenGetter called - auth provider is now configured");
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
    console.warn("[AUTH-TOKEN] No tokenGetter configured - auth not initialized yet");
    return null;
  }

  try {
    console.warn("[AUTH-TOKEN] Calling tokenGetter...");
    const token = await tokenGetter();
    console.warn("[AUTH-TOKEN] tokenGetter returned:", token ? `token(${token.length} chars)` : "null");
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
