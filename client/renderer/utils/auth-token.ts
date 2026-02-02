/**
 * Authentication token provider for API requests.
 *
 * This module provides a way to inject MSAL access tokens into API requests.
 * It uses a singleton pattern to store the token getter function, which is
 * set by the AuthProvider component when the app initializes.
 *
 * For Teams iframe context, we also support storing the Teams SSO token
 * directly, since MSAL's account cache isn't populated by Teams SSO.
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

// Teams SSO token storage (used when running in Teams iframe)
// Teams SSO returns a token directly via SDK, not through MSAL
let teamsToken: string | null = null;
let teamsTokenExpiresAt: number = 0;
let isTeamsContext: boolean = false;

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
  teamsToken = null;
  teamsTokenExpiresAt = 0;
  isTeamsContext = false;
}

const TEAMS_TOKEN_KEY = "ccp4i2-teams-token";
const TEAMS_TOKEN_EXPIRES_KEY = "ccp4i2-teams-token-expires";

/**
 * Store a Teams SSO token for use in API calls.
 * Called by login-content.tsx after successful Teams SSO.
 * Persists to sessionStorage to survive page navigation.
 *
 * @param token The access token from Teams SSO
 * @param expiresInSeconds Optional token lifetime (default: 1 hour)
 */
export function setTeamsToken(token: string, expiresInSeconds: number = 3600): void {
  teamsToken = token;
  teamsTokenExpiresAt = Date.now() + (expiresInSeconds * 1000) - 60000; // Expire 1 min early
  isTeamsContext = true;

  // Persist to sessionStorage
  if (typeof window !== "undefined") {
    try {
      sessionStorage.setItem(TEAMS_TOKEN_KEY, token);
      sessionStorage.setItem(TEAMS_TOKEN_EXPIRES_KEY, teamsTokenExpiresAt.toString());
    } catch (e) {
      console.warn("[AUTH-TOKEN] Failed to persist Teams token to sessionStorage:", e);
    }
  }

  console.log("[AUTH-TOKEN] Teams token stored, expires in", expiresInSeconds, "seconds");
}

/**
 * Load Teams token from sessionStorage (called on app init).
 */
export function loadTeamsToken(): boolean {
  if (typeof window === "undefined") return false;

  try {
    const storedToken = sessionStorage.getItem(TEAMS_TOKEN_KEY);
    const storedExpires = sessionStorage.getItem(TEAMS_TOKEN_EXPIRES_KEY);

    if (storedToken && storedExpires) {
      const expiresAt = parseInt(storedExpires, 10);
      if (Date.now() < expiresAt) {
        teamsToken = storedToken;
        teamsTokenExpiresAt = expiresAt;
        isTeamsContext = true;
        console.log("[AUTH-TOKEN] Loaded Teams token from sessionStorage");
        return true;
      } else {
        // Token expired, clear it
        sessionStorage.removeItem(TEAMS_TOKEN_KEY);
        sessionStorage.removeItem(TEAMS_TOKEN_EXPIRES_KEY);
      }
    }
  } catch (e) {
    console.warn("[AUTH-TOKEN] Failed to load Teams token from sessionStorage:", e);
  }

  return false;
}

/**
 * Check if we're running in Teams context with a stored token.
 */
export function hasTeamsToken(): boolean {
  // Try loading from sessionStorage if not in memory
  if (!teamsToken && typeof window !== "undefined") {
    loadTeamsToken();
  }
  if (!teamsToken || !isTeamsContext) return false;
  return Date.now() < teamsTokenExpiresAt;
}

/**
 * Clear only the Teams token (for re-auth scenarios).
 */
export function clearTeamsToken(): void {
  teamsToken = null;
  teamsTokenExpiresAt = 0;

  if (typeof window !== "undefined") {
    try {
      sessionStorage.removeItem(TEAMS_TOKEN_KEY);
      sessionStorage.removeItem(TEAMS_TOKEN_EXPIRES_KEY);
    } catch (e) {
      // Ignore
    }
  }
}

/**
 * Set a function to refresh the Teams token.
 * This is called when the stored token expires.
 */
let teamsTokenRefresher: (() => Promise<string | null>) | null = null;

export function setTeamsTokenRefresher(refresher: () => Promise<string | null>): void {
  teamsTokenRefresher = refresher;
}

/**
 * Get the current access token.
 * Returns null if no token getter is set or if token acquisition fails.
 *
 * In Teams context, uses the stored Teams SSO token.
 * Otherwise, uses MSAL via the token getter function.
 *
 * Uses a short-lived cache to avoid calling acquireTokenSilent on every request,
 * which dramatically improves performance for applications like Moorhen that
 * make many API calls in rapid succession.
 */
export async function getAccessToken(): Promise<string | null> {
  const now = Date.now();

  // Try loading Teams token from sessionStorage if not in memory
  if (!teamsToken && typeof window !== "undefined") {
    loadTeamsToken();
  }

  // In Teams context, use the stored Teams token
  if (isTeamsContext && teamsToken) {
    // Check if Teams token is still valid
    if (now < teamsTokenExpiresAt) {
      return teamsToken;
    }

    // Token expired - try to refresh it
    if (teamsTokenRefresher) {
      try {
        const newToken = await teamsTokenRefresher();
        if (newToken) {
          teamsToken = newToken;
          teamsTokenExpiresAt = now + 3600000 - 60000; // 1 hour minus 1 min buffer
          return newToken;
        }
      } catch (error) {
        console.error("[AUTH-TOKEN] Failed to refresh Teams token:", error);
      }
    }

    // Teams token expired and couldn't refresh - clear it
    console.warn("[AUTH-TOKEN] Teams token expired, clearing");
    teamsToken = null;
    teamsTokenExpiresAt = 0;
    // Don't clear isTeamsContext - we're still in Teams, just need re-auth
  }

  // Fall back to MSAL token getter
  if (!tokenGetter) {
    // No warning - this is expected in local dev mode without Azure AD
    return null;
  }

  // Return cached token if still valid
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
