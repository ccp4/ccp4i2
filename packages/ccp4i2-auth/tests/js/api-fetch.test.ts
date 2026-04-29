import { describe, it, expect } from "vitest";

import {
  AUTH_ERROR_EVENT,
  createApiFetch,
} from "../../src/api-fetch";

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
});
