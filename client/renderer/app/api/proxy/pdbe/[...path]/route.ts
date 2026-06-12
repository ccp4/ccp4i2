import { NextRequest, NextResponse } from "next/server";

/**
 * PDBe Proxy Route
 *
 * /api/proxy/pdbe/* -> https://www.ebi.ac.uk/pdbe/*
 *
 * Why this exists: the Moorhen viewer needs SharedArrayBuffer, so the
 * pages serving it carry `Cross-Origin-Embedder-Policy: require-corp`.
 * Under COEP, the browser blocks any cross-origin fetch whose response
 * does not advertise `Cross-Origin-Resource-Policy: cross-origin`. The
 * PDBe API and its `entry-files` download host both send permissive
 * CORS but neither sends CORP — so direct browser fetches fail with
 * the famously opaque `TypeError: Failed to fetch`.
 *
 * This route re-emits the upstream response with `Cross-Origin-Resource-
 * Policy: cross-origin` added, which satisfies COEP.
 *
 * Covers everything under /pdbe/ — `/pdbe/api/pdb/entry/files/{pdbId}`
 * (the JSON metadata) and `/pdbe/entry-files/download/{pdbId}.cif`
 * (the actual coords) both flow through here.
 *
 * Constraint: this proxy never sends anything from the user's local
 * structures to PDBe. It only forwards the URL path and query string.
 */

const UPSTREAM = "https://www.ebi.ac.uk/pdbe";

async function proxy(
  req: NextRequest,
  params: Promise<{ path: string[] }>,
  method: "GET" | "HEAD",
) {
  const { path } = await params;
  const subPath = path ? path.join("/") : "";
  const search = req.nextUrl.searchParams.toString();
  const queryString = search ? `?${search}` : "";
  const targetUrl = `${UPSTREAM}/${subPath}${queryString}`;

  try {
    const res = await fetch(targetUrl, { method, redirect: "follow" });
    const body = method === "HEAD" ? null : await res.arrayBuffer();

    const headers = new Headers();
    const contentType = res.headers.get("Content-Type");
    if (contentType) headers.set("Content-Type", contentType);
    const contentLength = res.headers.get("Content-Length");
    if (contentLength) headers.set("Content-Length", contentLength);
    headers.set("Cross-Origin-Resource-Policy", "cross-origin");

    return new NextResponse(body, { status: res.status, headers });
  } catch (err: any) {
    console.error("[PDBE PROXY] Error:", err);
    return NextResponse.json(
      { error: `PDBe proxy error: ${err.message}` },
      { status: 502 },
    );
  }
}

export async function GET(
  req: NextRequest,
  { params }: { params: Promise<{ path: string[] }> },
) {
  return proxy(req, params, "GET");
}

export async function HEAD(
  req: NextRequest,
  { params }: { params: Promise<{ path: string[] }> },
) {
  return proxy(req, params, "HEAD");
}
