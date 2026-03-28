/*
 * Copyright (C) 2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
import { NextRequest, NextResponse } from "next/server";

/**
 * UniProt Proxy Route
 *
 * Proxies requests to UniProt REST API to avoid CORS issues.
 * /api/proxy/uniprot/* -> https://rest.uniprot.org/*
 *
 * Uses rest.uniprot.org (the REST API) rather than www.uniprot.org (the website).
 * Follows redirects (UniProt redirects mnemonic IDs like CDK2_HUMAN to accession IDs).
 */

export async function GET(
  req: NextRequest,
  { params }: { params: Promise<{ path: string[] }> }
) {
  const { path } = await params;
  const uniprotPath = path ? path.join("/") : "";
  const searchParams = req.nextUrl.searchParams.toString();
  const queryString = searchParams ? `?${searchParams}` : "";
  const targetUrl = `https://rest.uniprot.org/${uniprotPath}${queryString}`;

  try {
    // follow redirects (UniProt redirects mnemonic IDs to accession IDs)
    const res = await fetch(targetUrl, { redirect: "follow" });
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
