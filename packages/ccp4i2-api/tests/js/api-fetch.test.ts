import { describe, it, expect, vi, afterEach, beforeEach } from "vitest";

import {
  AUTH_ERROR_EVENT,
  createApiFetch,
} from "../../src/api-fetch";
import {
  invalidateAccessToken,
  setTokenGetter,
  setEmailGetter,
} from "../../src/auth-token";

afterEach(() => {
  vi.restoreAllMocks();
});

describe("AUTH_ERROR_EVENT", () => {
  it("has the canonical event-name string consumers listen for", () => {
    expect(AUTH_ERROR_EVENT).toBe("ccp4i2:auth-error");
  });
});

describe("createApiFetch", () => {
  it("returns an ApiFetcher with all the expected method names", () => {
    const fetcher = createApiFetch({ baseUrl: "/api/proxy/test/" });
    const expected = [
      "apiFetch",
      "apiJson",
      "apiText",
      "apiBlob",
      "apiArrayBuffer",
      "apiPost",
      "apiPut",
      "apiPatch",
      "apiDelete",
      "apiGet",
      "apiUpload",
      "swrFetcher",
      "swrPostFetcher",
    ] as const;
    for (const name of expected) {
      expect(typeof (fetcher as Record<string, unknown>)[name]).toBe("function");
    }
  });

  it("returns independent fetchers when called twice with different baseUrls", () => {
    const ccp4i2 = createApiFetch({ baseUrl: "/api/proxy/ccp4i2/" });
    const compounds = createApiFetch({ baseUrl: "/api/proxy/compounds/" });
    // Independence: each call returns its own bound functions, not the
    // same singleton. (Smoke test — proves the factory pattern, not
    // a hardened isolation property.)
    expect(ccp4i2.apiFetch).not.toBe(compounds.apiFetch);
  });

  it("injects X-User-Email header when injectUserEmail callback returns a value", async () => {
    setTokenGetter(async () => "fake-token");
    const fetchMock = vi.spyOn(globalThis, "fetch").mockResolvedValue(
      new Response("[]", { status: 200, headers: { "content-type": "application/json" } }),
    );

    const fetcher = createApiFetch({
      baseUrl: "/api/proxy/compounds/",
      injectUserEmail: () => "user@example.com",
    });
    await fetcher.apiFetch("things");

    expect(fetchMock).toHaveBeenCalledOnce();
    const [, init] = fetchMock.mock.calls[0];
    const headers = (init as RequestInit).headers as Record<string, string>;
    expect(headers["Authorization"]).toBe("Bearer fake-token");
    expect(headers["X-User-Email"]).toBe("user@example.com");
  });

  it("omits X-User-Email when injectUserEmail returns null", async () => {
    setTokenGetter(async () => "fake-token");
    const fetchMock = vi.spyOn(globalThis, "fetch").mockResolvedValue(
      new Response("[]", { status: 200, headers: { "content-type": "application/json" } }),
    );

    const fetcher = createApiFetch({
      baseUrl: "/api/proxy/compounds/",
      injectUserEmail: () => null,
    });
    await fetcher.apiFetch("things");

    const [, init] = fetchMock.mock.calls[0];
    const headers = (init as RequestInit).headers as Record<string, string>;
    expect(headers["X-User-Email"]).toBeUndefined();
  });

  it("does not inject X-User-Email when option not provided", async () => {
    setTokenGetter(async () => "fake-token");
    setEmailGetter(() => "fallback@example.com"); // would be set by AuthProvider; should NOT leak
    const fetchMock = vi.spyOn(globalThis, "fetch").mockResolvedValue(
      new Response("[]", { status: 200, headers: { "content-type": "application/json" } }),
    );

    const fetcher = createApiFetch({ baseUrl: "/api/proxy/ccp4i2/" });
    await fetcher.apiFetch("things");

    const [, init] = fetchMock.mock.calls[0];
    const headers = (init as RequestInit).headers as Record<string, string>;
    expect(headers["X-User-Email"]).toBeUndefined();
  });
});

describe("responses with no body", () => {
  it("resolves apiDelete on a 204 rather than failing to parse it", async () => {
    // A successful DELETE is the common case, and JSON.parse of an empty
    // body throws -- so before this, every successful withdraw/delete
    // surfaced to the caller as an error over work the server had done.
    setTokenGetter(async () => "fake-token");
    vi.spyOn(globalThis, "fetch").mockResolvedValue(
      new Response(null, {
        status: 204,
        headers: { "content-type": "application/json" },
      }),
    );

    const fetcher = createApiFetch({ baseUrl: "/api/proxy/ccp4i2/" });
    await expect(fetcher.apiDelete("things/1/")).resolves.toBeUndefined();
  });

  it("still reports a DELETE the server refused", async () => {
    setTokenGetter(async () => "fake-token");
    vi.spyOn(globalThis, "fetch").mockResolvedValue(
      new Response(JSON.stringify({ error: "Site not found in this campaign" }), {
        status: 404,
        headers: { "content-type": "application/json" },
      }),
    );

    const fetcher = createApiFetch({ baseUrl: "/api/proxy/ccp4i2/" });
    await expect(fetcher.apiDelete("things/1/")).rejects.toThrow(
      "Site not found in this campaign",
    );
  });
});

describe("recovering from a 401", () => {
  // The token cache is module-global by design (one session, many fetchers),
  // so a token cached by an earlier test would be presented by this one.
  beforeEach(() => {
    invalidateAccessToken();
  });

  /** A fetch mock that answers 401 until a fresh token arrives, then 200. */
  function serverRejectingStaleTokens(freshToken: string) {
    return vi.fn(async (_url: string, init?: RequestInit) => {
      const auth = (init?.headers as Record<string, string>)?.["Authorization"];
      if (auth === `Bearer ${freshToken}`) {
        return new Response(JSON.stringify({ ok: true }), {
          status: 200,
          headers: { "content-type": "application/json" },
        });
      }
      return new Response(JSON.stringify({ error: "expired" }), {
        status: 401,
        headers: { "content-type": "application/json" },
      });
    });
  }

  it("refreshes and replays once rather than reporting the failure", async () => {
    // The whole point: a stale token should cost the user nothing. Before
    // this, the first 401 went straight to a snackbar telling them to sign in.
    let issued = 0;
    setTokenGetter(async (options) => {
      // The cache is invalidated before a forced refresh, so the getter is
      // asked again; a plain call would otherwise be served from cache.
      if (options?.forceRefresh) return "fresh";
      issued += 1;
      return "stale";
    });
    const fetchMock = serverRejectingStaleTokens("fresh");
    vi.spyOn(globalThis, "fetch").mockImplementation(fetchMock as any);

    const fetcher = createApiFetch({ baseUrl: "/api/proxy/ccp4i2/" });
    await expect(fetcher.apiGet("things/")).resolves.toEqual({ ok: true });
    expect(fetchMock).toHaveBeenCalledTimes(2);
    expect(issued).toBeGreaterThan(0);
  });

  it("reports the error when a fresh token is refused too", async () => {
    // A genuinely revoked session must still reach the user, and must not
    // spin: exactly one replay.
    setTokenGetter(async () => "no-good");
    const errors: number[] = [];
    const listener = (e: Event) =>
      errors.push((e as CustomEvent<AuthErrorDetail>).detail.status);
    window.addEventListener(AUTH_ERROR_EVENT, listener);

    const fetchMock = serverRejectingStaleTokens("never-issued");
    vi.spyOn(globalThis, "fetch").mockImplementation(fetchMock as any);

    const fetcher = createApiFetch({ baseUrl: "/api/proxy/ccp4i2/" });
    await expect(fetcher.apiGet("things/")).rejects.toThrow();
    expect(fetchMock).toHaveBeenCalledTimes(2);
    expect(errors).toEqual([401]);
    window.removeEventListener(AUTH_ERROR_EVENT, listener);
  });

  it("refreshes once however many requests fail at the same moment", async () => {
    // A stale session fails every hook on the page at once. Without
    // single-flighting, that is a dozen refreshes and a dozen chances to
    // trip an Azure AD throttle.
    let refreshes = 0;
    setTokenGetter(async (options) => {
      if (options?.forceRefresh) {
        refreshes += 1;
        await new Promise((r) => setTimeout(r, 10));
        return "fresh";
      }
      return "stale";
    });
    vi.spyOn(globalThis, "fetch").mockImplementation(
      serverRejectingStaleTokens("fresh") as any,
    );

    const fetcher = createApiFetch({ baseUrl: "/api/proxy/ccp4i2/" });
    await Promise.all(
      Array.from({ length: 8 }, (_, i) => fetcher.apiGet(`things/${i}/`)),
    );
    expect(refreshes).toBe(1);
  });

  it("does not replay a request whose body cannot be sent twice", async () => {
    // A stream is consumed by the first attempt; replaying it would send an
    // empty body, which is worse than the 401.
    setTokenGetter(async () => "stale");
    const fetchMock = serverRejectingStaleTokens("fresh");
    vi.spyOn(globalThis, "fetch").mockImplementation(fetchMock as any);

    const fetcher = createApiFetch({ baseUrl: "/api/proxy/ccp4i2/" });
    const stream = new ReadableStream({
      start(controller) {
        controller.enqueue(new TextEncoder().encode("{}"));
        controller.close();
      },
    });
    await expect(
      fetcher.apiFetch("things/", { method: "POST", body: stream as any }),
    ).rejects.toThrow();
    expect(fetchMock).toHaveBeenCalledTimes(1);
  });
});
