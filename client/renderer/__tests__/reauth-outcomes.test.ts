/**
 * What asking to re-authenticate can answer, and which answer means "sign out".
 *
 * Only one of them does. The original bug was a stale session costing a full
 * logout/login cycle; a recovery path that signs the user out whenever a
 * redirect fails to *start* reintroduces it by the back door, because MSAL
 * throws interaction_in_progress for things as ordinary as pressing the
 * snackbar button twice.
 */

import { describe, it, expect, beforeEach } from "vitest";
import {
  canReauthenticate,
  reauthenticate,
  setReauthHandler,
} from "../utils/reauth";

beforeEach(() => {
  setReauthHandler(null);
});

const msalError = (errorCode: string) =>
  Object.assign(new Error(`${errorCode}: something`), { errorCode });

describe("reauthenticate", () => {
  it("reports 'unavailable' where nothing can re-authenticate", async () => {
    // The desktop session: it lives as long as the app process, so there is
    // nothing to renew. This is the only outcome a caller may answer by
    // signing out.
    expect(await reauthenticate()).toBe("unavailable");
    expect(canReauthenticate()).toBe(false);
  });

  it("reports 'started' when a redirect is under way", async () => {
    setReauthHandler(async () => true);
    expect(await reauthenticate()).toBe("started");
  });

  it("reports 'in-progress' when MSAL says one is already running", async () => {
    // BrowserCacheManager.setInteractionInProgress throws this outright.
    // A sign-in that is already happening is the thing working, not failing.
    setReauthHandler(async () => {
      throw msalError("interaction_in_progress");
    });
    expect(await reauthenticate()).toBe("in-progress");
  });

  it("recognises interaction_in_progress by message alone", async () => {
    // Not every MSAL version hands back a shaped errorCode.
    setReauthHandler(async () => {
      throw new Error("BrowserAuthError: interaction_in_progress");
    });
    expect(await reauthenticate()).toBe("in-progress");
  });

  it("answers a second press without troubling MSAL", async () => {
    // The button is a snackbar action and easy to double-click, which is the
    // commonest way to reach interaction_in_progress at all.
    let calls = 0;
    setReauthHandler(async () => {
      calls += 1;
      return true;
    });
    expect(await reauthenticate()).toBe("started");
    expect(await reauthenticate()).toBe("in-progress");
    expect(calls).toBe(1);
  });

  it("reports 'failed' for any other error, and allows a retry", async () => {
    // Failing to start leaves the session exactly as it was, so the next
    // press must be able to try again rather than being told one is running.
    let attempts = 0;
    setReauthHandler(async () => {
      attempts += 1;
      if (attempts === 1) throw msalError("network_error");
      return true;
    });
    expect(await reauthenticate()).toBe("failed");
    expect(await reauthenticate()).toBe("started");
  });

  it("registering a handler clears a previous run's state", async () => {
    setReauthHandler(async () => true);
    await reauthenticate();
    setReauthHandler(async () => true);
    expect(await reauthenticate()).toBe("started");
  });
});
