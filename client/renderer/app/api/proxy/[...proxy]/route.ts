import { NextRequest, NextResponse } from "next/server";

export async function GET(
  req: NextRequest,
  { params }: { params: Promise<{ proxy: string[] }> }
) {
  return await handleProxy(req, await params);
}

export async function POST(
  req: NextRequest,
  { params }: { params: Promise<{ proxy: string[] }> }
) {
  return await handleProxy(req, await params);
}

export async function PUT(
  req: NextRequest,
  { params }: { params: Promise<{ proxy: string[] }> }
) {
  return await handleProxy(req, await params);
}

export async function PATCH(
  req: NextRequest,
  { params }: { params: Promise<{ proxy: string[] }> }
) {
  return await handleProxy(req, await params);
}

export async function DELETE(
  req: NextRequest,
  { params }: { params: Promise<{ proxy: string[] }> }
) {
  return await handleProxy(req, await params);
}

interface RequestInitWithDuplex extends RequestInit {
  duplex?: string;
}
// Common handler for all HTTP methods
async function handleProxy(req: NextRequest, params: { proxy: string[] }) {
  // Priority: runtime API_BASE_URL (set by Electron) > build-time NEXT_PUBLIC > default
  // Note: NEXT_PUBLIC_* vars are embedded at build time, so we check API_BASE_URL first
  // for runtime configuration in packaged Electron apps
  let backendBaseUrl =
    process.env.API_BASE_URL ||
    process.env.NEXT_PUBLIC_API_BASE_URL ||
    "http://localhost:8000"; // Default backend URL

  if (req.headers.get("x-backend-url")) {
    backendBaseUrl = req.headers.get("x-backend-url") as string;
  }

  if (!backendBaseUrl) {
    return NextResponse.json(
      { error: "Backend URL is not configured" },
      { status: 500 }
    );
  }

  if (!backendBaseUrl.endsWith("/")) {
    backendBaseUrl += "/";
  }

  // Ensure params.proxy is properly handled
  const path = params.proxy ? params.proxy.join("/") : "";

  if (params?.proxy?.at(0) === "uniprot") {
    let targetUrl = `https://www.uniprot.org/${path}`;
    //console.log("uniprot targetUrl", targetUrl);
    try {
      const res = await fetch(targetUrl);
      if (!res.ok) {
        return NextResponse.json(
          { error: "Failed to fetch from UniProt" },
          { status: res.status }
        );
      }

      const data = await res.text();
      //console.log({data})

      return new NextResponse(data, {
        status: 200,
        headers: {
          "Content-Type": "text/plain",
        },
      });
    } catch (err) {
      return NextResponse.json({ error: "Server error" }, { status: 500 });
    }
  }

  let targetUrl = `${backendBaseUrl}${path}`;

  // Ensure the backend URL ends with a slash
  //console.log("req_url", req.url);
  if (!targetUrl.includes("/djangostatic") && !targetUrl.endsWith("/")) {
    targetUrl += "/";
  }

  // Append query parameters if any
  const searchParams = req.nextUrl.searchParams.toString();
  if (searchParams) {
    targetUrl += `?${searchParams}`;
  }
  if (["PATCH", "DELETE", "PUT", "POST"].includes(req.method)) {
    if (!targetUrl.endsWith("/")) {
      targetUrl += "/";
    }
  }
  //console.log("targetUrl", targetUrl, "req_url", req.url);
  try {
    // Clone the headers from the incoming request
    const headers = new Headers(req.headers);

    // Detect if the request is multipart/form-data (file upload)
    const isMultipart = headers
      .get("content-type")
      ?.startsWith("multipart/form-data");

    // Set up fetch options
    const fetchOptions: RequestInitWithDuplex = {
      method: req.method,
      headers,
      body: isMultipart || req.method !== "GET" ? req.body : undefined, // Forward body directly
      duplex: "half",
    };

    const response = await fetch(targetUrl, fetchOptions);

    return new NextResponse(response.body, {
      status: response.status,
      headers: response.headers,
    });
  } catch (error) {
    console.error("Error forwarding request:", error);
    return NextResponse.json(
      { error: `Proxy error: ${error.message}` },
      { status: 500 }
    );
  }
}
