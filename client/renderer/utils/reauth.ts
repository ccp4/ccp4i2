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

/**
 * What came of asking to re-authenticate.
 *
 * Deliberately not a boolean. "Cannot ever" and "cannot right now" are
 * different answers and only the first has anything to do with signing out:
 * a redirect already under way is the sign-in working, not failing.
 */
export type ReauthOutcome =
  /** A redirect has been started; the page is on its way to AAD. */
  | "started"
  /** One is already running. Nothing to do but let it finish. */
  | "in-progress"
  /** It could not be started. Worth telling the user; not worth signing out. */
  | "failed"
  /** Nothing can re-authenticate here at all (the desktop session). */
  | "unavailable";

let reauthHandler: ReauthHandler | null = null;
/**
 * Set once a redirect has been asked for. MSAL throws
 * ``interaction_in_progress`` if a second interaction starts while one is
 * running, and a snackbar button is easy to press twice, so the second press
 * is answered here rather than by an exception.
 */
let reauthStarted = false;

/** Registered by AuthProvider once it knows which auth path is in use. */
export function setReauthHandler(handler: ReauthHandler | null): void {
  reauthHandler = handler;
  reauthStarted = false;
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

/** MSAL's code for "an interaction is already running in this tab". */
function isInteractionInProgress(error: unknown): boolean {
  const code = (error as { errorCode?: string })?.errorCode;
  if (code === "interaction_in_progress") return true;
  const message = (error as { message?: string })?.message ?? "";
  return message.includes("interaction_in_progress");
}

/**
 * Send the user through a re-auth round trip, returning them here afterwards.
 *
 * Only "unavailable" means there is nothing to try -- and only that should
 * ever lead a caller to offer signing out. A redirect that is already running,
 * or one that failed to start, both leave the session exactly as it was, and
 * signing the user out over either would be the logout/login cycle this
 * exists to remove.
 */
export async function reauthenticate(): Promise<ReauthOutcome> {
  if (!reauthHandler) return "unavailable";
  if (reauthStarted) return "in-progress";
  stashReturnUrl();
  try {
    reauthStarted = true;
    const started = await reauthHandler();
    if (!started) {
      reauthStarted = false;
      return "failed";
    }
    return "started";
  } catch (error) {
    reauthStarted = false;
    if (isInteractionInProgress(error)) {
      // A redirect is already going, or one was abandoned with Back and left
      // its status behind. Either way this is not a reason to sign out.
      console.warn("[AUTH] A sign-in is already in progress.");
      return "in-progress";
    }
    console.error("[AUTH] Re-authentication failed to start:", error);
    return "failed";
  }
}
