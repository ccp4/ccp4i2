import { NextRequest, NextResponse } from "next/server";

/**
 * CCP4i2 API Proxy Route
 *
 * Handles /api/proxy/ccp4i2/* requests and forwards them to the Django backend
 * at /api/ccp4i2/*. This explicit namespacing keeps ccp4i2 API calls separate
 * from other apps (like compounds).
 */

/**
 * Public endpoints that don't require authentication.
 * - version: fetched from the app-selector page before login.
 * - health:  the launch-time readiness probe (LaunchGate / useServerReady) polls
 *   this before the app is up and before any token exists, so it must not require
 *   auth. The Django health_check view is itself unauthenticated (a plain view,
 *   no DRF permissions), so exposing it here matches its intent.
 */
const PUBLIC_ENDPOINTS = ["version", "health"];

/**
 * Check if the path is a public endpoint that doesn't require auth.
 */
function isPublicEndpoint(path: string): boolean {
  const normalizedPath = path.replace(/\/$/, ""); // Remove trailing slash
  return PUBLIC_ENDPOINTS.some(
    (endpoint) => normalizedPath === endpoint || normalizedPath.endsWith(`/${endpoint}`)
  );
}

/**
 * Check if authentication should be required.
 */
function isAuthRequired(path: string): boolean {
  // Public endpoints never require auth
  if (isPublicEndpoint(path)) {
    return false;
  }
  if (process.env.REQUIRE_PROXY_AUTH?.toLowerCase() === "true") {
    return true;
  }
  if (process.env.CONTAINER_APP_NAME || process.env.CONTAINER_APP_ENV_DNS_SUFFIX) {
    return true;
  }
  return false;
}

/**
 * File-grant plumbing.
 *
 * A grant is a signed, expiring, directory-scoped read capability minted by
 * the backend (see ccp4i2_api.file_grants). It exists for requests the
 * browser issues on its own behalf — the images, stylesheets and relative
 * fetches of a task's HTML report opened in a tab — which can never carry an
 * Authorization header.
 *
 * The initial navigation carries the grant as a query parameter; we hand it
 * back as a cookie scoped to that report's directory, so the page's own
 * subsequent requests carry it and nothing else on the origin does. The
 * backend re-checks the grant's scope on every request, so the cookie's path
 * is a convenience, not the security boundary.
 */
const GRANT_QUERY_PARAM = "file_grant";
const GRANT_COOKIE = "ccp4i2_file_grant";
const GRANT_HEADER = "X-CCP4I2-File-Grant";
// Matches the backend's default grant lifetime; an expired cookie is
// harmless (the backend rejects the grant and the tab must be reopened).
const GRANT_COOKIE_MAX_AGE = 60 * 60;

function extractGrant(req: NextRequest): string | null {
  return (
    req.nextUrl.searchParams.get(GRANT_QUERY_PARAM) ||
    req.cookies.get(GRANT_COOKIE)?.value ||
    null
  );
}

/**
 * Whether a grant may stand in for a token at this edge.
 *
 * The backend validates the grant properly — signature, expiry, and that the
 * request path lies inside the subtree it was minted for. This check mirrors
 * the two limits that are cheap to enforce here, so that presenting any
 * grant-shaped string cannot widen what reaches the backend beyond read
 * requests for the path-based file endpoint. It matters when the backend's
 * own auth is misconfigured: without it, a deployment that requires proxy
 * auth but has no auth middleware behind it would be fully exposed by an
 * arbitrary ?file_grant= value.
 */
function grantMayAuthenticate(req: NextRequest, path: string): boolean {
  if (req.method !== "GET" && req.method !== "HEAD") return false;
  if (!path.includes("files_by_path/")) return false;
  return Boolean(extractGrant(req));
}

/**
 * Extract authentication token from request.
 * Checks (in order): Authorization header, Azure Easy Auth header, query parameter.
 * Query parameter is needed for file downloads where anchor clicks don't send headers.
 */
function extractAuthToken(req: NextRequest): string | null {
  // Check Authorization header first
  const authHeader = req.headers.get("Authorization");
  if (authHeader?.startsWith("Bearer ")) {
    return authHeader.substring(7);
  }
  // Check Azure Easy Auth header
  const easyAuthToken = req.headers.get("X-MS-TOKEN-AAD-ACCESS-TOKEN");
  if (easyAuthToken) {
    return easyAuthToken;
  }
  // Check query parameter (used by file downloads via anchor clicks)
  const queryToken = req.nextUrl.searchParams.get("access_token");
  if (queryToken) {
    return queryToken;
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

  // Check authentication if required
  if (isAuthRequired(path) && !grantMayAuthenticate(req, path)) {
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

  // Ensure trailing slash for Django REST Framework endpoints,
  // but not for path-based file serving where the path IS the resource.
  const isFilePath = path.includes("files_by_path/");
  if (!isFilePath && !targetUrl.endsWith("/")) {
    targetUrl += "/";
  }

  // Append query parameters
  const searchParams = req.nextUrl.searchParams.toString();
  if (searchParams) {
    targetUrl += `?${searchParams}`;
  }

  try {
    const headers = new Headers(req.headers);

    // Forward authentication headers
    const token = extractAuthToken(req);
    if (token) {
      headers.set("Authorization", `Bearer ${token}`);
    }

    const grant = extractGrant(req);
    if (grant) {
      headers.set(GRANT_HEADER, grant);
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

    // Check if this is a JSON response - if so, use arrayBuffer to handle gzip decompression
    const contentType = response.headers.get('Content-Type') || '';
    if (contentType.includes('application/json')) {
      const data = await response.arrayBuffer();
      const headers = new Headers();
      headers.set('Content-Type', 'application/json');
      headers.set('Content-Length', String(data.byteLength));
      // CORP header allows loading from pages with COEP: require-corp (Moorhen)
      headers.set('Cross-Origin-Resource-Policy', 'cross-origin');
      return new NextResponse(data, {
        status: response.status,
        headers,
      });
    }

    // For non-JSON responses (file downloads, etc.), stream directly
    // Create new headers, excluding Content-Encoding since fetch() auto-decompresses
    // the response body but the original header would tell the browser to decompress again
    const responseHeaders = new Headers();
    response.headers.forEach((value, key) => {
      // Skip Content-Encoding as fetch() already decompressed the body
      // Also skip Content-Length as it may be wrong after decompression
      if (key.toLowerCase() !== 'content-encoding' && key.toLowerCase() !== 'content-length') {
        responseHeaders.set(key, value);
      }
    });
    // CORP header allows loading from pages with COEP: require-corp (Moorhen)
    responseHeaders.set('Cross-Origin-Resource-Policy', 'cross-origin');

    const fileResponse = new NextResponse(response.body, {
      status: response.status,
      headers: responseHeaders,
    });
    // The navigation that carried ?file_grant= is the report page itself:
    // hand the grant back as a cookie scoped to its directory so the page's
    // relative requests inherit it, and nothing above that directory does.
    const seededGrant = req.nextUrl.searchParams.get(GRANT_QUERY_PARAM);
    if (seededGrant && response.ok) {
      fileResponse.cookies.set(GRANT_COOKIE, seededGrant, {
        httpOnly: true,
        secure: process.env.NODE_ENV === "production",
        sameSite: "lax",
        path: req.nextUrl.pathname.replace(/[^/]*$/, ""),
        maxAge: GRANT_COOKIE_MAX_AGE,
      });
    }
    return fileResponse;
  } catch (error: any) {
    // Only log unexpected errors, not connection-refused during startup.
    // fetch (undici) reports that as a bare "fetch failed" with the refusal
    // in .cause, which is why the launcher's health polls were logged as
    // errors on every start-up.
    const code = error.code ?? error.cause?.code;
    if (code !== "ECONNREFUSED") {
      // undici's "fetch failed" carries the real reason in .cause; without it
      // a truncated upload and a refused connection log the same line.
      const cause = error.cause?.message ? ` (${error.cause.message})` : "";
      console.error("[CCP4I2 PROXY] Error:", error.message + cause, "→", targetUrl);
    }
    return NextResponse.json(
      { error: `Proxy error: ${error.message}`, targetUrl, code: error.code },
      { status: 500 }
    );
  }
}
