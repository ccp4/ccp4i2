/**
 * Warm the HTTP cache for Moorhen's WASM while the data archives download.
 *
 * MoorhenCommandCentre.init() creates the Coot worker, then *awaits*
 * fetch(data.tar.gz) -- and rota_other-data.tar.gz, and
 * rota_arg_lys-data.tar.gz -- and only then posts CootInitialize. That message
 * is what makes the worker run its 64-bit probe, importScripts the loader and
 * fetch the WASM. So the largest download on the page cannot begin until three
 * archives have finished: measured on a poor link, 50 s for the first archive
 * alone, with the WASM still on 360 KB when the trace was taken.
 *
 * Nothing here changes Moorhen. The wrapper simply asks for the same files, by
 * the same URLs, as soon as it mounts; the worker's later fetch then joins the
 * in-flight request or hits the cache entry it left behind. The cleaner fix is
 * in Moorhen -- start both fetches together -- and this is not a substitute for
 * it, but it recovers most of the time without waiting for a Moorhen release.
 *
 * Deliberately best-effort: every failure is swallowed. A prefetch that 404s,
 * is blocked, or is aborted must never break the viewer, because the worker
 * fetches these files again regardless and is the one that actually matters.
 */

/**
 * Does this browser take the 64-bit Moorhen build?
 *
 * The same probe the worker runs (a wasm module using a 64-bit memory index),
 * so the wrapper prefetches the file the worker will actually ask for. Getting
 * this wrong is not fatal -- it warms the wrong entry and wastes the bytes --
 * but getting it right is the whole point.
 */
export function prefersWasm64(): boolean {
  try {
    return WebAssembly.validate(
      new Uint8Array([0, 97, 115, 109, 1, 0, 0, 0, 5, 3, 1, 4, 1])
    );
  } catch {
    return false;
  }
}

/**
 * The files the Coot worker loads once CootInitialize reaches it.
 *
 * The worker does `importScripts("./moorhen64.js")` relative to itself, i.e.
 * `<urlPrefix>/wasm/`, and Emscripten's loader then fetches the matching
 * .wasm. Both are worth warming: the .js is small but sits on the critical
 * path ahead of the .wasm.
 *
 * The worker falls back to the 32-bit build if importScripts throws, and that
 * fallback is not prefetched: it is the rare path, and fetching both builds
 * would put 17 MB of waste on the link this exists to relieve.
 */
export function moorhenWasmUrls(urlPrefix: string): string[] {
  const stem = prefersWasm64() ? "moorhen64" : "moorhen";
  return [`${urlPrefix}/wasm/${stem}.js`, `${urlPrefix}/wasm/${stem}.wasm`];
}

/**
 * Start the prefetches. Returns a function that aborts them.
 *
 * Call on mount and abort on unmount, so a viewer closed early does not leave
 * 17 MB downloading in the background.
 *
 * `cache: "force-cache"` asks for the cache entry the worker will want rather
 * than a revalidation, and the response body is discarded without being read:
 * the point is to fill the HTTP cache, not to hold 17 MB in JS memory.
 */
export function prefetchMoorhenWasm(urlPrefix: string): () => void {
  if (typeof window === "undefined" || typeof fetch !== "function") {
    return () => {};
  }

  const controller = new AbortController();

  for (const url of moorhenWasmUrls(urlPrefix)) {
    try {
      fetch(url, {
        signal: controller.signal,
        cache: "force-cache",
        credentials: "same-origin",
      }).catch(() => {
        // Best-effort: the worker fetches these itself regardless.
      });
    } catch {
      // Ditto -- never let a warm-up break the viewer.
    }
  }

  return () => controller.abort();
}
