/**
 * What the CCP4i2 proxy does with a response that has no body.
 *
 * Django answers a successful DELETE with 204 No Content, and DRF still
 * labels it application/json. The proxy's JSON branch rebuilds the response
 * from an arrayBuffer, and the Response constructor rejects any body -- even
 * a zero-length one -- alongside a 204. That threw, was caught as a proxy
 * error, and every successful DELETE reached the browser as a 500: the
 * campaign panel reported "Failed to withdraw verdict" over a verdict the
 * server had just withdrawn.
 */

import { describe, it, expect, vi, afterEach } from "vitest";
import { NextRequest } from "next/server";
import { DELETE } from "../app/api/proxy/ccp4i2/[...path]/route";

const path = [
  "projectgroups",
  "7",
  "sites",
  "3",
  "evaluation",
  "12",
];

const request = () =>
  new NextRequest(
    new URL(`http://localhost:3000/api/proxy/ccp4i2/${path.join("/")}`),
    { method: "DELETE" }
  );

afterEach(() => {
  vi.restoreAllMocks();
});

describe("a no-content response through the proxy", () => {
  it("passes a 204 through rather than failing on it", async () => {
    vi.spyOn(globalThis, "fetch").mockResolvedValue(
      // What DRF sends for Response(status=204): no body, but still typed.
      new Response(null, {
        status: 204,
        headers: { "Content-Type": "application/json" },
      })
    );

    const response = await DELETE(request(), {
      params: Promise.resolve({ path }),
    });

    expect(response.status).toBe(204);
    expect(await response.text()).toBe("");
  });

  it("still carries a body-bearing response through", async () => {
    vi.spyOn(globalThis, "fetch").mockResolvedValue(
      new Response(JSON.stringify({ verdict: "hit" }), {
        status: 200,
        headers: { "Content-Type": "application/json" },
      })
    );

    const response = await DELETE(request(), {
      params: Promise.resolve({ path }),
    });

    expect(response.status).toBe(200);
    expect(await response.json()).toEqual({ verdict: "hit" });
  });
});
