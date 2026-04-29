import { describe, it, expect, afterEach } from "vitest";

import {
  createLocalSessionEmailGetter,
  createLocalSessionTokenGetter,
  hasLocalSessionToken,
} from "../../src/providers/local-session";

afterEach(() => {
  // Ensure no test leaks the surface into the next test.
  delete (window as { ccp4i2LocalSession?: unknown }).ccp4i2LocalSession;
});

describe("hasLocalSessionToken", () => {
  it("returns false when window.ccp4i2LocalSession is unset", () => {
    expect(hasLocalSessionToken()).toBe(false);
  });

  it("returns true when window.ccp4i2LocalSession.token is a string", () => {
    (window as { ccp4i2LocalSession: { token: string } }).ccp4i2LocalSession = {
      token: "deadbeef",
    };
    expect(hasLocalSessionToken()).toBe(true);
  });

  it("returns false when ccp4i2LocalSession exists but token is not a string", () => {
    (window as { ccp4i2LocalSession: unknown }).ccp4i2LocalSession = {
      token: undefined,
    };
    expect(hasLocalSessionToken()).toBe(false);
  });
});

describe("createLocalSessionTokenGetter", () => {
  it("yields the token from window.ccp4i2LocalSession", async () => {
    (window as { ccp4i2LocalSession: { token: string } }).ccp4i2LocalSession = {
      token: "secret-from-electron",
    };
    const getter = createLocalSessionTokenGetter();
    await expect(getter()).resolves.toBe("secret-from-electron");
  });

  it("yields null when no local-session surface is exposed", async () => {
    const getter = createLocalSessionTokenGetter();
    await expect(getter()).resolves.toBeNull();
  });
});

describe("createLocalSessionEmailGetter", () => {
  it("yields the userEmail when set", () => {
    (window as { ccp4i2LocalSession: { token: string; userEmail: string } }).ccp4i2LocalSession = {
      token: "t",
      userEmail: "dev@ccp4i2.invalid",
    };
    expect(createLocalSessionEmailGetter()()).toBe("dev@ccp4i2.invalid");
  });

  it("yields null when no email is exposed", () => {
    expect(createLocalSessionEmailGetter()()).toBeNull();
  });
});
