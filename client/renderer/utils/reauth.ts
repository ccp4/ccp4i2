/**
 * Re-authentication: signing back *in*, as distinct from signing out.
 *
 * When a session goes stale the app used to offer the user `logout()`, which
 * is the opposite of what they want. Signing out destroys the AAD SSO session
 * cookie -- the one thing that would have made getting back in silent -- and
 * with no account hint it first asks which account to sign out of. The user
 * then signs in again from scratch, credentials and MFA included, and lands
 * on `/` rather than the page they were reading.
 *
 * Re-authenticating instead asks AAD to re-issue a token for the account the
 * app already knows about. With a live SSO session that round-trips without a
 * prompt, and the return-url convention (see `/auth/callback`) puts the user
 * back where they were.
 *
 * Why a registry rather than `useMsal()` at the call site: on the desktop
 * LocalSession path, AuthProvider renders its children *outside*
 * `MsalProvider`, so the hook would throw. AuthProvider registers whatever
 * re-auth means for the path it took -- and on desktop that is nothing at
 * all, because the session lives as long as the app process.
 */

/** Returns true if it started a re-auth, false if it cannot. */
export type ReauthHandler = () => Promise<boolean>;

let reauthHandler: ReauthHandler | null = null;

/** Registered by AuthProvider once it knows which auth path is in use. */
export function setReauthHandler(handler: ReauthHandler | null): void {
  reauthHandler = handler;
}

/** Whether anything can re-authenticate, for UI that offers it. */
export function canReauthenticate(): boolean {
  return reauthHandler !== null;
}

/**
 * Remember where the user is, so the callback can put them back.
 *
 * Path, query and hash: a Moorhen page without its `?job=` is a different
 * page. Uses the same sessionStorage key as the login flow, which
 * `/auth/callback` already reads.
 */
export function stashReturnUrl(): void {
  if (typeof window === "undefined") return;
  const here = `${window.location.pathname}${window.location.search}${window.location.hash}`;
  try {
    sessionStorage.setItem("auth-return-url", here);
  } catch {
    // Private mode or blocked storage: the user lands on / instead. Worth
    // losing the destination rather than the sign-in.
  }
}

/**
 * Send the user through a re-auth round trip, returning them here afterwards.
 *
 * Returns false when no handler is registered, so a caller can fall back to
 * signing out rather than leaving a button that does nothing.
 */
export async function reauthenticate(): Promise<boolean> {
  if (!reauthHandler) return false;
  stashReturnUrl();
  try {
    return await reauthHandler();
  } catch (error) {
    console.error("[AUTH] Re-authentication failed to start:", error);
    return false;
  }
}
