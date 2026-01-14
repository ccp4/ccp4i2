import { NextRequest, NextResponse } from "next/server";

/**
 * CCP4i2 API Proxy Route
 *
 * Handles /api/proxy/ccp4i2/* requests and forwards them to the Django backend
 * at /api/ccp4i2/*. This explicit namespacing keeps ccp4i2 API calls separate
 * from other apps (like compounds).
 */

/**
 * Check if authentication should be required.
 */
function isAuthRequired(): boolean {
  if (process.env.REQUIRE_PROXY_AUTH?.toLowerCase() === "true") {
    return true;
  }
  if (process.env.CONTAINER_APP_NAME || process.env.CONTAINER_APP_ENV_DNS_SUFFIX) {
    return true;
  }
  return false;
}

/**
 * Extract authentication token from request.
 */
function extractAuthToken(req: NextRequest): string | null {
  const authHeader = req.headers.get("Authorization");
  if (authHeader?.startsWith("Bearer ")) {
    return authHeader.substring(7);
  }
  const easyAuthToken = req.headers.get("X-MS-TOKEN-AAD-ACCESS-TOKEN");
  if (easyAuthToken) {
    return easyAuthToken;
  }
  return null;
}

export async function GET(
  req: NextRequest,
  { params }: { params: Promise<{ path: string[] }> }
) {
  return await handleProxy(req, await params);
}

export async function POST(
  req: NextRequest,
  { params }: { params: Promise<{ path: string[] }> }
) {
  return await handleProxy(req, await params);
}

export async function PUT(
  req: NextRequest,
  { params }: { params: Promise<{ path: string[] }> }
) {
  return await handleProxy(req, await params);
}

export async function PATCH(
  req: NextRequest,
  { params }: { params: Promise<{ path: string[] }> }
) {
  return await handleProxy(req, await params);
}

export async function DELETE(
  req: NextRequest,
  { params }: { params: Promise<{ path: string[] }> }
) {
  return await handleProxy(req, await params);
}

interface RequestInitWithDuplex extends RequestInit {
  duplex?: string;
}

async function handleProxy(req: NextRequest, params: { path: string[] }) {
  const path = params.path ? params.path.join("/") : "";
  console.log("[CCP4I2 PROXY] Handling request:", req.method, path);

  // Check authentication if required
  if (isAuthRequired()) {
    const token = extractAuthToken(req);
    if (!token) {
      const principalId = req.headers.get("X-MS-CLIENT-PRINCIPAL-ID");
      if (!principalId) {
        return NextResponse.json(
          { success: false, error: "Authentication required. Please sign in." },
          { status: 401 }
        );
      }
    }
  }

  // Get backend URL
  let backendBaseUrl =
    process.env.API_BASE_URL ||
    process.env.NEXT_PUBLIC_API_BASE_URL ||
    "http://localhost:8000";

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

  // Construct target URL - Django serves ccp4i2 API at /api/ccp4i2/
  let targetUrl = `${backendBaseUrl}api/ccp4i2/${path}`;

  // Ensure trailing slash for Django REST Framework endpoints
  if (!targetUrl.endsWith("/")) {
    targetUrl += "/";
  }

  // Append query parameters
  const searchParams = req.nextUrl.searchParams.toString();
  if (searchParams) {
    targetUrl += `?${searchParams}`;
  }

  console.log("[CCP4I2 PROXY] Forwarding to:", targetUrl);

  try {
    const headers = new Headers(req.headers);

    // Forward authentication headers
    const token = extractAuthToken(req);
    if (token) {
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

    const isMultipart = headers
      .get("content-type")
      ?.startsWith("multipart/form-data");

    const fetchOptions: RequestInitWithDuplex = {
      method: req.method,
      headers,
      body: isMultipart || req.method !== "GET" ? req.body : undefined,
      duplex: "half",
    };

    const response = await fetch(targetUrl, fetchOptions);

    return new NextResponse(response.body, {
      status: response.status,
      headers: response.headers,
    });
  } catch (error: unknown) {
    const errorMessage = error instanceof Error ? error.message : String(error);
    console.error("[CCP4I2 PROXY] Error forwarding request:", error);
    return NextResponse.json(
      { error: `Proxy error: ${errorMessage}`, targetUrl },
      { status: 500 }
    );
  }
}
