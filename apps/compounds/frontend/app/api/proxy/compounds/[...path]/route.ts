import { NextRequest, NextResponse } from 'next/server';

/**
 * Compounds API Proxy Route
 *
 * Handles /api/proxy/compounds/* requests and forwards them to the Django backend
 * at /api/compounds/*. Mirrors the auth forwarding pattern from ccp4i2 proxy.
 */

// Use same env var pattern as main ccp4i2 proxy for consistency
// Priority: runtime API_BASE_URL > build-time NEXT_PUBLIC > default
const DJANGO_URL = process.env.API_BASE_URL || process.env.NEXT_PUBLIC_API_BASE_URL || 'http://localhost:8000';

// API path on Django backend - matches ccp4i2 multi-app URL structure
const API_PATH = '/api/compounds';

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

  console.log(`[Compounds Proxy] GET ${url}`);

  try {
    const response = await fetch(url, {
      headers: buildForwardHeaders(request),
      redirect: 'manual', // Handle redirects manually for file downloads
    });

    // Handle redirects (e.g., to Azure SAS URLs for file downloads)
    if (response.status >= 300 && response.status < 400) {
      const location = response.headers.get('Location');
      if (location) {
        return NextResponse.redirect(location, response.status);
      }
    }

    // Check content type to determine how to handle response
    const contentType = response.headers.get('Content-Type') || '';

    // For file downloads (non-JSON responses), stream the response
    if (!contentType.includes('application/json')) {
      const headers = new Headers();
      // Forward relevant headers
      const forwardHeaders = ['Content-Type', 'Content-Disposition', 'Content-Length'];
      for (const header of forwardHeaders) {
        const value = response.headers.get(header);
        if (value) headers.set(header, value);
      }

      return new NextResponse(response.body, {
        status: response.status,
        headers,
      });
    }

    // For JSON responses, stream directly to avoid memory issues with large payloads
    // Previous approach of parsing and re-serializing caused issues with ~12MB compound lists
    const headers = new Headers();
    headers.set('Content-Type', 'application/json');
    // Forward content-length if available
    const contentLength = response.headers.get('Content-Length');
    if (contentLength) {
      headers.set('Content-Length', contentLength);
    }

    return new NextResponse(response.body, {
      status: response.status,
      headers,
    });
  } catch (error) {
    console.error('[Compounds Proxy] Error:', error);
    return NextResponse.json(
      { error: 'Failed to fetch from Django', detail: String(error) },
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

  console.log(`[Compounds Proxy] POST ${url}`);

  try {
    const contentType = request.headers.get('content-type') || '';
    const headers = buildForwardHeaders(request);
    let response: Response;

    if (contentType.includes('multipart/form-data')) {
      // Handle file uploads - pass through the FormData
      // Don't set Content-Type - let fetch set it with boundary for multipart
      headers.delete('Content-Type');
      const formData = await request.formData();

      // Log form data fields for debugging
      const formDataFields: string[] = [];
      for (const [key, value] of formData.entries()) {
        if (value instanceof File || (value && typeof value === 'object' && 'name' in value)) {
          formDataFields.push(`${key}: File(${(value as File).name}, ${(value as File).size} bytes)`);
        } else {
          formDataFields.push(`${key}: ${String(value).substring(0, 50)}`);
        }
      }
      console.log(`[Compounds Proxy] FormData fields: ${formDataFields.join(', ')}`);

      response = await fetch(url, {
        method: 'POST',
        headers,
        body: formData,
      });
    } else {
      // Handle JSON requests (or empty body for action endpoints)
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

      response = await fetch(url, {
        method: 'POST',
        headers,
        body: body !== null ? JSON.stringify(body) : undefined,
      });
    }

    // Try to parse response as JSON, handling non-JSON error responses
    const responseText = await response.text();
    try {
      const data = JSON.parse(responseText);
      return NextResponse.json(data, { status: response.status });
    } catch (jsonError) {
      // Django returned non-JSON response (likely HTML error page)
      console.error('[Compounds Proxy] Non-JSON response:', response.status, responseText.substring(0, 500));
      return NextResponse.json(
        {
          error: 'Server error',
          detail: responseText.substring(0, 1000),
          status: response.status
        },
        { status: response.status || 500 }
      );
    }
  } catch (error) {
    console.error('[Compounds Proxy] Error:', error);
    return NextResponse.json(
      { error: 'Failed to fetch from Django', detail: String(error) },
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

  console.log(`[Compounds Proxy] PATCH ${url}`);

  try {
    const headers = buildForwardHeaders(request);
    headers.set('Content-Type', 'application/json');

    const response = await fetch(url, {
      method: 'PATCH',
      headers,
      body: JSON.stringify(body),
    });

    // Try to parse response as JSON, handling non-JSON error responses
    const responseText = await response.text();
    try {
      const data = JSON.parse(responseText);
      return NextResponse.json(data, { status: response.status });
    } catch (jsonError) {
      console.error('[Compounds Proxy] Non-JSON response:', response.status, responseText.substring(0, 500));
      return NextResponse.json(
        {
          error: 'Server error',
          detail: responseText.substring(0, 1000),
          status: response.status
        },
        { status: response.status || 500 }
      );
    }
  } catch (error) {
    console.error('[Compounds Proxy] Error:', error);
    return NextResponse.json(
      { error: 'Failed to fetch from Django', detail: String(error) },
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

  console.log(`[Compounds Proxy] DELETE ${url}`);

  try {
    const response = await fetch(url, {
      method: 'DELETE',
      headers: buildForwardHeaders(request),
    });

    if (response.status === 204) {
      return new NextResponse(null, { status: 204 });
    }

    // Try to parse response as JSON, handling non-JSON error responses
    const responseText = await response.text();
    try {
      const data = JSON.parse(responseText);
      return NextResponse.json(data, { status: response.status });
    } catch (jsonError) {
      console.error('[Compounds Proxy] Non-JSON response:', response.status, responseText.substring(0, 500));
      return NextResponse.json(
        {
          error: 'Server error',
          detail: responseText.substring(0, 1000),
          status: response.status
        },
        { status: response.status || 500 }
      );
    }
  } catch (error) {
    console.error('[Compounds Proxy] Error:', error);
    return NextResponse.json(
      { error: 'Failed to fetch from Django', detail: String(error) },
      { status: 500 }
    );
  }
}
