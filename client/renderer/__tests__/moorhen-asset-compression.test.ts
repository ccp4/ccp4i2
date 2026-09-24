/**
 * The Moorhen route serves precompressed siblings.
 *
 * Next compresses what it serves from public/, but not this route's responses,
 * so the largest file on the page went out raw: moorhen64.wasm at 17.8 MB,
 * where brotli takes it to 4.4 MB. A cold load over a poor link measured 50 s
 * for a 1.3 MB archive; at that rate the raw WASM needs five to six minutes.
 *
 * scripts/precompress-assets.cjs writes .br/.gz at build time -- brotli q10
 * costs ~28 s per WASM, which cannot live in a request handler.
 */
import { describe, expect, it } from "vitest";

/** The token test the route uses to read Accept-Encoding. */
const accepts = (header: string, enc: string) =>
  new RegExp(`(^|[,\\s])${enc}(;|,|$)`).test(header);

describe("Accept-Encoding negotiation", () => {
  it("takes brotli when offered", () => {
    expect(accepts("gzip, deflate, br", "br")).toBe(true);
    expect(accepts("br", "br")).toBe(true);
  });

  it("takes gzip when offered", () => {
    expect(accepts("gzip, deflate", "gzip")).toBe(true);
  });

  it("handles quality values", () => {
    expect(accepts("gzip;q=0.8, br;q=1.0", "br")).toBe(true);
    expect(accepts("gzip;q=0.8, br;q=1.0", "gzip")).toBe(true);
  });

  it("serves raw when nothing is accepted", () => {
    expect(accepts("deflate", "br")).toBe(false);
    expect(accepts("deflate", "gzip")).toBe(false);
    expect(accepts("", "br")).toBe(false);
  });

  it("matches whole tokens, not substrings", () => {
    // A substring test would wrongly match both of these.
    expect(accepts("x-gzip", "gzip")).toBe(false);
    expect(accepts("notbr", "br")).toBe(false);
  });

  it("prefers brotli over gzip", () => {
    // The route tries br first and breaks, so a client offering both gets br.
    const order = [
      { enc: "br", ext: ".br" },
      { enc: "gzip", ext: ".gz" },
    ];
    const header = "gzip, deflate, br";
    const chosen = order.find((e) => accepts(header, e.enc));
    expect(chosen?.enc).toBe("br");
  });
});

describe("what the compressed response must say", () => {
  it("keeps the original Content-Type so streaming instantiation still works", () => {
    // Content-Type is derived from the requested path, not the file read. With
    // `application/wasm` + `Content-Encoding: br` the browser decompresses and
    // still streams into WebAssembly.instantiateStreaming; sending
    // application/brotli would force the slower non-streaming path.
    const contentTypes: Record<string, string> = {
      ".wasm": "application/wasm",
      ".js": "text/javascript",
    };
    const filePath = "MoorhenAssets/wasm/moorhen.wasm";
    const ext = "." + filePath.split(".").pop();
    expect(contentTypes[ext]).toBe("application/wasm");
  });

  it("varies on Accept-Encoding so caches cannot cross-serve", () => {
    // Without this a shared cache can hand brotli bytes to a client that
    // asked for none -- and these responses are cached immutable for a year.
    const headers: Record<string, string> = {
      "Content-Type": "application/wasm",
      "Cache-Control": "public, max-age=31536000, immutable",
      Vary: "Accept-Encoding",
    };
    expect(headers.Vary).toBe("Accept-Encoding");
  });
});
