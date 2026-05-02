import { describe, it, expect, vi, afterEach } from "vitest";

import {
  AUTH_ERROR_EVENT,
  createApiFetch,
} from "../../src/api-fetch";
import { setTokenGetter, setEmailGetter } from "../../src/auth-token";

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
