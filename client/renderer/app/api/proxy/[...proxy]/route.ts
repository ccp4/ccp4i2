import { NextRequest, NextResponse } from "next/server";

/**
 * Check if authentication should be required.
 *
 * Authentication is required when:
 * - REQUIRE_PROXY_AUTH=true is set (explicit opt-in)
 * - We're running in Azure (WEBSITE_INSTANCE_ID or similar Azure indicator is present)
 *
 * This ensures local development works without auth while Azure is protected.
 */
function isAuthRequired(): boolean {
  // Explicit opt-in
  if (process.env.REQUIRE_PROXY_AUTH?.toLowerCase() === "true") {
    return true;
  }

  // Running in Azure Container Apps (has these environment variables)
  if (process.env.CONTAINER_APP_NAME || process.env.CONTAINER_APP_ENV_DNS_SUFFIX) {
    return true;
  }

  return false;
}

/**
 * Extract authentication token from request.
 * Supports:
 * - Authorization: Bearer <token>
 * - X-MS-TOKEN-AAD-ACCESS-TOKEN header (Azure Easy Auth)
 * - Cookie-based session from Easy Auth
 */
function extractAuthToken(req: NextRequest): string | null {
  // Check Authorization header first
  const authHeader = req.headers.get("Authorization");
  if (authHeader?.startsWith("Bearer ")) {
    return authHeader.substring(7);
  }

  // Check Azure Easy Auth token header
  const easyAuthToken = req.headers.get("X-MS-TOKEN-AAD-ACCESS-TOKEN");
  if (easyAuthToken) {
    return easyAuthToken;
  }

  return null;
}

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
  // Check authentication if required
  if (isAuthRequired()) {
    const token = extractAuthToken(req);
    if (!token) {
      // Check for Azure Easy Auth principal header (indicates authenticated user)
      const principalId = req.headers.get("X-MS-CLIENT-PRINCIPAL-ID");
      if (!principalId) {
        return NextResponse.json(
          {
            success: false,
            error: "Authentication required. Please sign in.",
          },
          { status: 401 }
        );
      }
    }
  }

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

  // Construct target URL with /api/ccp4i2/ prefix for multi-app integration
  let targetUrl = `${backendBaseUrl}api/ccp4i2/${path}`;

  // Ensure trailing slash for Django REST Framework endpoints
  if (!targetUrl.endsWith("/")) {
    targetUrl += "/";
  }

  // Append query parameters if any (after the trailing slash)
  const searchParams = req.nextUrl.searchParams.toString();
  if (searchParams) {
    targetUrl += `?${searchParams}`;
  }
  //console.log("targetUrl", targetUrl, "req_url", req.url);
  try {
    // Clone the headers from the incoming request
    const headers = new Headers(req.headers);

    // Forward authentication headers to backend
    // This allows Django's auth middleware to validate the token
    const token = extractAuthToken(req);
    if (token) {
      // Forward as Authorization header for Django middleware
      headers.set("Authorization", `Bearer ${token}`);
    }

    // Forward Azure Easy Auth headers if present
    const easyAuthHeaders = [
      "X-MS-TOKEN-AAD-ACCESS-TOKEN",
      "X-MS-TOKEN-AAD-ID-TOKEN",
      "X-MS-CLIENT-PRINCIPAL",
      "X-MS-CLIENT-PRINCIPAL-ID",
      "X-MS-CLIENT-PRINCIPAL-NAME",
    ];

    for (const headerName of easyAuthHeaders) {
      const headerValue = req.headers.get(headerName);
      if (headerValue) {
        headers.set(headerName, headerValue);
      }
    }

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
