/**
 * Local-session token provider.
 *
 * In CCP4i2 desktop (Electron), the main process generates a per-launch
 * secret at startup, passes it to the Django backend via an env var on the
 * spawned child process, and exposes it to the renderer via the preload
 * script (contextBridge.exposeInMainWorld). The renderer asks for the
 * token through the API surface below; the LocalSessionAuthMiddleware on
 * the server side validates the same secret on every request.
 *
 * In non-Electron contexts (web build, tests) `window.ccp4i2LocalSession`
 * is undefined and ``hasLocalSessionToken()`` returns false, letting the
 * consumer fall back to a different provider (MSAL for cloud) without any
 * code branching at the auth-token layer.
 *
 * The token never touches disk. It lives in three process memory spaces
 * (Electron main, preload, Django) and travels in one HTTP header per
 * request. Process death = secret death.
 */

import type { TokenGetter } from "../auth-token.js";

interface LocalSessionWindowSurface {
  /** The per-launch secret, exposed read-only by the Electron preload. */
  readonly token: string;
}

declare global {
  // eslint-disable-next-line @typescript-eslint/consistent-type-definitions
  interface Window {
    ccp4i2LocalSession?: LocalSessionWindowSurface;
  }
}

/** True iff the Electron preload has exposed a local-session token. */
export function hasLocalSessionToken(): boolean {
  if (typeof window === "undefined") return false;
  return typeof window.ccp4i2LocalSession?.token === "string";
}

/**
 * Build a TokenGetter that reads the per-launch secret exposed on
 * `window` by the Electron preload script.
 *
 * Returns a getter that yields null in non-Electron contexts (web build,
 * tests) so consumers can detect "no local session" and fall through to
 * another provider.
 */
export function createLocalSessionTokenGetter(): TokenGetter {
  return async () => {
    if (typeof window === "undefined") return null;
    const token = window.ccp4i2LocalSession?.token;
    return typeof token === "string" ? token : null;
  };
}
