/**
 * The WASM download starts at mount, not after the data archives.
 *
 * MoorhenCommandCentre.init() awaits fetch(data.tar.gz),
 * rota_other-data.tar.gz and rota_arg_lys-data.tar.gz before posting
 * CootInitialize -- and only that message makes the worker probe for 64-bit
 * support, importScripts the loader and fetch the WASM. The biggest file on
 * the page is therefore last in the queue: a trace over a poor link showed
 * 50 s for the first archive with the WASM still on 360 KB.
 */
import { afterEach, describe, expect, it, vi } from "vitest";
import { moorhenWasmUrls, prefetchMoorhenWasm } from "../lib/moorhen-wasm-prefetch";

const PREFIX = "/api/moorhen/v/1.0.1-dev.g10d4c0b00/MoorhenAssets";

afterEach(() => {
  vi.unstubAllGlobals();
});

describe("moorhenWasmUrls", () => {
  it("asks for exactly what the worker will ask for", () => {
    // The worker does importScripts("./moorhen64.js") relative to
    // <prefix>/wasm/CootWorker.js. If these URLs differ by a character the
    // cache entry is not shared and the prefetch achieves nothing.
    const [js, wasm] = moorhenWasmUrls(PREFIX);
    const workerScript = `${PREFIX}/wasm/CootWorker.js`;
    const resolve = (rel: string) =>
      new URL(rel, `http://host${workerScript}`).pathname;

    const stem = js.includes("moorhen64") ? "moorhen64" : "moorhen";
    expect(js).toBe(resolve(`./${stem}.js`));
    expect(wasm).toBe(resolve(`./${stem}.wasm`));
  });

  it("warms the loader as well as the WASM", () => {
    const urls = moorhenWasmUrls(PREFIX);
    expect(urls.some((u) => u.endsWith(".js"))).toBe(true);
    expect(urls.some((u) => u.endsWith(".wasm"))).toBe(true);
  });

  it("picks one build, not both", () => {
    // Fetching both would put ~17 MB of waste on the link this is meant to
    // relieve.
    const urls = moorhenWasmUrls(PREFIX);
    expect(urls).toHaveLength(2);
    const stems = new Set(urls.map((u) => u.split("/").pop()!.split(".")[0]));
    expect(stems.size).toBe(1);
  });
});

describe("prefetchMoorhenWasm", () => {
  it("starts both requests immediately", () => {
    const fetchMock = vi.fn((_url: string, _init?: RequestInit) =>
      Promise.resolve(new Response())
    );
    vi.stubGlobal("fetch", fetchMock);

    prefetchMoorhenWasm(PREFIX);

    expect(fetchMock).toHaveBeenCalledTimes(2);
    const requested = fetchMock.mock.calls.map(([url]) => url);
    expect(requested).toEqual(moorhenWasmUrls(PREFIX));
  });

  it("returns an abort so a closed viewer stops downloading", () => {
    const signals: AbortSignal[] = [];
    vi.stubGlobal(
      "fetch",
      vi.fn((_u: string, init: RequestInit) => {
        signals.push(init.signal as AbortSignal);
        return Promise.resolve(new Response());
      })
    );

    const abort = prefetchMoorhenWasm(PREFIX);
    expect(signals.every((s) => !s.aborted)).toBe(true);
    abort();
    expect(signals.every((s) => s.aborted)).toBe(true);
  });

  it("survives fetch rejecting", async () => {
    // The worker fetches these files itself regardless, so a failed warm-up
    // must never surface as an error.
    vi.stubGlobal("fetch", vi.fn(() => Promise.reject(new Error("offline"))));
    expect(() => prefetchMoorhenWasm(PREFIX)).not.toThrow();
    await Promise.resolve();
  });

  it("survives fetch throwing synchronously", () => {
    vi.stubGlobal(
      "fetch",
      vi.fn(() => {
        throw new Error("blocked");
      })
    );
    expect(() => prefetchMoorhenWasm(PREFIX)).not.toThrow();
  });
});
