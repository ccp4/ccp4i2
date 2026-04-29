/**
 * Per-launch local-session secret for CCP4i2 desktop authentication.
 *
 * At Electron startup we generate a cryptographically random token. The
 * same token is:
 *
 * 1. installed in ``process.env.CCP4I2_LOCAL_SESSION_TOKEN`` so that:
 *    - the Django child process inherits it via the existing
 *      ``...process.env`` spread in ccp4i2-django-server.ts, and
 *    - the Electron preload (a Node context launched per-window) reads
 *      it on load and exposes it to the renderer via contextBridge.
 * 2. consumed by ``LocalSessionAuthMiddleware`` on the Django side, which
 *    HMAC-compares the bearer token on every request.
 *
 * The token never touches disk. It lives in three process memory spaces
 * (main, preload, Django) and travels in one HTTP header per request.
 * Process death = secret death.
 *
 * Side-effect import: this module sets the env vars on first import.
 * The first import must therefore happen *before* the Django spawn and
 * *before* any BrowserWindow is created.
 */

import { randomBytes } from "node:crypto";
import os from "node:os";

/**
 * Sanitise an OS username so it forms the local-part of a valid email.
 *
 * On domain-joined Windows boxes, ``os.userInfo().username`` can return
 * ``DOMAIN\martin`` form, whose backslash is not a legal email local-
 * part character. We take the part after the last separator (handles
 * both backslash and forward-slash variants) and strip anything that
 * isn't [A-Za-z0-9._-].
 *
 * If sanitisation collapses to empty (extremely unusual), fall back to
 * the static "desktop" so the address is still well-formed.
 */
function sanitiseUsername(raw: string): string {
  const afterSeparator = raw.split(/[\\/]/).pop() ?? "";
  const sanitised = afterSeparator.replace(/[^A-Za-z0-9._-]/g, "_");
  return sanitised || "desktop";
}

const username = sanitiseUsername(os.userInfo().username);

/** The per-launch secret. 64-char hex, ASCII-only, safe in HTTP headers. */
export const LOCAL_SESSION_TOKEN: string = randomBytes(32).toString("hex");

/**
 * The Django user email the desktop session authenticates as.
 * RFC 6761 reserves ``.invalid`` so this can never collide with a real
 * deliverable address on any platform.
 */
export const LOCAL_USER_EMAIL: string = `${username}@ccp4i2.invalid`;

// Side-effect: install in process.env so the Django child process and
// the Electron preload both see the same values.
process.env.CCP4I2_LOCAL_SESSION_TOKEN = LOCAL_SESSION_TOKEN;
process.env.CCP4I2_LOCAL_USER_EMAIL = LOCAL_USER_EMAIL;
