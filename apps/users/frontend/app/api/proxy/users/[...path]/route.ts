import { NextRequest, NextResponse } from 'next/server';

/**
 * Users API Proxy Route
 *
 * Handles /api/proxy/users/* requests and forwards them to the Django backend
 * at /api/users/*. Mirrors the auth forwarding pattern from compounds proxy.
 */

// Use same env var pattern as main ccp4i2 proxy for consistency
// Priority: runtime API_BASE_URL > build-time NEXT_PUBLIC > default
const DJANGO_URL = process.env.API_BASE_URL || process.env.NEXT_PUBLIC_API_BASE_URL || 'http://localhost:8000';

// API path on Django backend
const API_PATH = '/api/users';

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

/**
 * Build headers for forwarding to Django, including auth headers.
 */
function buildForwardHeaders(req: NextRequest): Headers {
  const headers = new Headers();
  headers.set('Accept', 'application/json');

  // Forward authentication token
  const token = extractAuthToken(req);
  if (token) {
    headers.set("Authorization", `Bearer ${token}`);
  }

  // Forward user email header (fallback for when access token lacks email claim)
  const userEmail = req.headers.get("X-User-Email");
  if (userEmail) {
    headers.set("X-User-Email", userEmail);
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

  return headers;
}

/**
 * Check auth and return 401 response if required but missing.
 */
function checkAuth(req: NextRequest): NextResponse | null {
  if (isAuthRequired()) {
    const token = extractAuthToken(req);
    if (!token) {
      const principalId = req.headers.get("X-MS-CLIENT-PRINCIPAL-ID");
      if (!principalId) {
        return NextResponse.json(
          { success: false, error: "Authentication required. Provide Authorization: Bearer <token>" },
          { status: 401 }
        );
      }
    }
  }
  return null;
}

export async function GET(
  request: NextRequest,
  { params }: { params: Promise<{ path: string[] }> }
) {
  const authError = checkAuth(request);
  if (authError) return authError;

  const { path } = await params;
  const pathString = path.join('/');
  const searchParams = request.nextUrl.searchParams.toString();
  const queryString = searchParams ? `?${searchParams}` : '';
  const url = `${DJANGO_URL}${API_PATH}/${pathString}/${queryString}`;

  console.log(`[Users Proxy] GET ${url}`);

  try {
    const response = await fetch(url, {
      headers: buildForwardHeaders(request),
    });

    const data = await response.json();
    return NextResponse.json(data, { status: response.status });
  } catch (error) {
    console.error('[Users Proxy] Error:', error);
    return NextResponse.json(
      { error: 'Failed to fetch from Django' },
      { status: 500 }
    );
  }
}

export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ path: string[] }> }
) {
  const authError = checkAuth(request);
  if (authError) return authError;

  const { path } = await params;
  const pathString = path.join('/');
  const url = `${DJANGO_URL}${API_PATH}/${pathString}/`;

  console.log(`[Users Proxy] POST ${url}`);

  try {
    const headers = buildForwardHeaders(request);
    let body: any = null;
    const contentLength = request.headers.get('content-length');
    if (contentLength && parseInt(contentLength) > 0) {
      try {
        body = await request.json();
      } catch (e) {
        // Ignore JSON parse errors for empty bodies
      }
    }

    if (body !== null) {
      headers.set('Content-Type', 'application/json');
    }

    const response = await fetch(url, {
      method: 'POST',
      headers,
      body: body !== null ? JSON.stringify(body) : undefined,
    });

    const data = await response.json();
    return NextResponse.json(data, { status: response.status });
  } catch (error) {
    console.error('[Users Proxy] Error:', error);
    return NextResponse.json(
      { error: 'Failed to fetch from Django' },
      { status: 500 }
    );
  }
}

export async function PATCH(
  request: NextRequest,
  { params }: { params: Promise<{ path: string[] }> }
) {
  const authError = checkAuth(request);
  if (authError) return authError;

  const { path } = await params;
  const pathString = path.join('/');
  const url = `${DJANGO_URL}${API_PATH}/${pathString}/`;
  const body = await request.json();

  console.log(`[Users Proxy] PATCH ${url}`);

  try {
    const headers = buildForwardHeaders(request);
    headers.set('Content-Type', 'application/json');

    const response = await fetch(url, {
      method: 'PATCH',
      headers,
      body: JSON.stringify(body),
    });

    const data = await response.json();
    return NextResponse.json(data, { status: response.status });
  } catch (error) {
    console.error('[Users Proxy] Error:', error);
    return NextResponse.json(
      { error: 'Failed to fetch from Django' },
      { status: 500 }
    );
  }
}

export async function DELETE(
  request: NextRequest,
  { params }: { params: Promise<{ path: string[] }> }
) {
  const authError = checkAuth(request);
  if (authError) return authError;

  const { path } = await params;
  const pathString = path.join('/');
  const url = `${DJANGO_URL}${API_PATH}/${pathString}/`;

  console.log(`[Users Proxy] DELETE ${url}`);

  try {
    const response = await fetch(url, {
      method: 'DELETE',
      headers: buildForwardHeaders(request),
    });

    if (response.status === 204) {
      return new NextResponse(null, { status: 204 });
    }

    const data = await response.json();
    return NextResponse.json(data, { status: response.status });
  } catch (error) {
    console.error('[Users Proxy] Error:', error);
    return NextResponse.json(
      { error: 'Failed to fetch from Django' },
      { status: 500 }
    );
  }
}
