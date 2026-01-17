import { NextRequest, NextResponse } from "next/server";

/**
 * UniProt Proxy Route
 *
 * Proxies requests to UniProt API to avoid CORS issues.
 * /api/proxy/uniprot/* -> https://www.uniprot.org/*
 */

export async function GET(
  req: NextRequest,
  { params }: { params: Promise<{ path: string[] }> }
) {
  const { path } = await params;
  const uniprotPath = path ? path.join("/") : "";
  const searchParams = req.nextUrl.searchParams.toString();
  const queryString = searchParams ? `?${searchParams}` : "";
  const targetUrl = `https://www.uniprot.org/${uniprotPath}${queryString}`;

  console.log("[UNIPROT PROXY] Forwarding to:", targetUrl);

  try {
    const res = await fetch(targetUrl);
    if (!res.ok) {
      return NextResponse.json(
        { error: "Failed to fetch from UniProt" },
        { status: res.status }
      );
    }

    const data = await res.text();
    return new NextResponse(data, {
      status: 200,
      headers: {
        "Content-Type": res.headers.get("Content-Type") || "text/plain",
        // Allow loading from COEP-enabled pages (Moorhen)
        "Cross-Origin-Resource-Policy": "cross-origin",
      },
    });
  } catch (err: any) {
    console.error("[UNIPROT PROXY] Error:", err);
    return NextResponse.json(
      { error: `UniProt proxy error: ${err.message}` },
      { status: 500 }
    );
  }
}
