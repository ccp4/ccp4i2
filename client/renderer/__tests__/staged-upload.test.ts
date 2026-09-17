import { describe, it, expect, vi, afterEach } from "vitest";
import { setTokenGetter, clearTokenGetter } from "@ccp4/ccp4i2-api";

import { stageFile } from "../lib/staged-upload";

const cap = { chunk_bytes: 8, max_bytes: 1000, threshold_bytes: 4 };

const file = (bytes: string) => new File([bytes], "m.mrc");

describe("stageFile", () => {
  afterEach(() => clearTokenGetter());

  it("sends the bearer token on begin, every chunk and finish", async () => {
    setTokenGetter(async () => "tok-123");
    const headers: Array<Record<string, string> | undefined> = [];
    global.fetch = vi.fn(async (url: any, opts: any) => {
      const u = String(url);
      headers.push(opts?.headers);
      if (u.endsWith("staged-uploads/")) return { ok: true, json: async () => ({ upload_id: "abc", chunk_bytes: 8 }) } as any;
      if (u.includes("/finish/")) return { ok: true, json: async () => ({ state: "ready" }) } as any;
      return { ok: true } as any;
    }) as any;

    await stageFile(file("0123456789abcdefXY"), cap);

    expect(headers.length).toBe(5); // begin + 3 chunks + finish
    for (const h of headers) expect(h?.Authorization).toBe("Bearer tok-123");
    expect(headers[0]?.["Content-Type"]).toBe("application/json");
    expect(headers[1]?.["Content-Type"]).toBe("application/octet-stream");
  });

  it("begins, PUTs every chunk, finishes, and returns the handle", async () => {
    const calls: string[] = [];
    global.fetch = vi.fn(async (url: any, opts: any) => {
      const u = String(url);
      calls.push(`${opts?.method || "GET"} ${u.split("staged-uploads")[1] ?? u}`);
      if (u.endsWith("staged-uploads/")) {
        return { ok: true, json: async () => ({ upload_id: "abc", chunk_bytes: 8 }) } as any;
      }
      if (u.includes("/finish/")) {
        return { ok: true, json: async () => ({ state: "ready" }) } as any;
      }
      return { ok: true } as any; // chunk PUT -> 204
    }) as any;

    // 18 bytes -> ceil(18/8) = 3 chunks
    const handle = await stageFile(file("0123456789abcdefXY"), cap);

    expect(handle).toBe("abc");
    expect(calls.filter((c) => c.includes("/chunks/")).length).toBe(3);
    expect(calls.some((c) => c.includes("/finish/"))).toBe(true);
  });

  it("reports progress to 1 across all chunks", async () => {
    global.fetch = vi.fn(async (url: any) => {
      const u = String(url);
      if (u.endsWith("staged-uploads/")) return { ok: true, json: async () => ({ upload_id: "x", chunk_bytes: 8 }) } as any;
      if (u.includes("/finish/")) return { ok: true, json: async () => ({ state: "ready" }) } as any;
      return { ok: true } as any;
    }) as any;
    const fractions: number[] = [];
    await stageFile(file("0123456789abcdefXY"), cap, { onProgress: (f) => fractions.push(f) });
    expect(Math.max(...fractions)).toBeCloseTo(1, 5);
  });

  it("maps a 413 begin to a friendly error", async () => {
    global.fetch = vi.fn(async () => ({ ok: false, status: 413, json: async () => ({}) }) as any) as any;
    await expect(stageFile(file("big"), cap)).rejects.toThrow(/larger than/i);
  });
});
