/**
 * The prefetch is inert on the server.
 *
 * These wrappers render under SSR, where there is no window and no HTTP cache
 * to warm. Kept as a *.test.ts so it runs in the node environment with no DOM
 * -- which is exactly the condition being checked. (The behavioural tests live
 * in moorhen-wasm-prefetch.test.tsx, under jsdom, because a prefetch needs a
 * window to do anything at all.)
 */
import { describe, expect, it, vi } from "vitest";
import { prefetchMoorhenWasm } from "../lib/moorhen-wasm-prefetch";

describe("prefetchMoorhenWasm under SSR", () => {
  it("does nothing and returns a no-op teardown when there is no window", () => {
    expect(typeof window).toBe("undefined");

    const fetchMock = vi.fn();
    vi.stubGlobal("fetch", fetchMock);

    const teardown = prefetchMoorhenWasm("/api/moorhen/MoorhenAssets");

    expect(fetchMock).not.toHaveBeenCalled();
    // Still safe for an effect cleanup to call.
    expect(() => teardown()).not.toThrow();

    vi.unstubAllGlobals();
  });
});
