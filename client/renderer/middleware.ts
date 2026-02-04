import { NextResponse } from "next/server";
import type { NextRequest } from "next/server";

/**
 * Next.js Middleware for authentication gating and Cross-Origin headers.
 *
 * Authentication (when NEXT_PUBLIC_REQUIRE_AUTH=true):
 * - Checks for auth-session cookie on all requests
 * - Unauthenticated users are redirected to /auth/login
 * - /auth/login triggers MSAL Azure AD authentication
 * - After login, auth-provider sets the cookie
 *
 * This ensures NO application code loads for unauthenticated users,
 * reducing attack surface and preventing information disclosure.
 *
 * Cross-Origin headers (for Moorhen):
 * - Adds CORP headers for static files
 * - Adds COEP/COOP headers for Moorhen pages (SharedArrayBuffer support)
 */

// Paths that don't require authentication
const AUTH_EXEMPT_PATHS = [
  "/api/health", // Health check for Azure Container Apps probes
  "/api/auth/", // Auth session API (cookie management)
  "/auth/login", // Login page that triggers MSAL
  "/auth/callback", // MSAL redirect callback (completes auth flow)
  "/auth/teams-start", // Teams auth start (redirects to Azure AD)
  "/auth/teams-callback", // Teams auth callback (receives token)
  "/_next", // Next.js internals (static assets, webpack HMR)
  "/favicon.ico",
];

// Paths that accept Bearer token authentication (for CLI tools like i2remote)
// These paths handle their own auth validation in the route handler
const BEARER_AUTH_PATHS = [
  "/api/proxy/ccp4i2/", // CCP4i2 API proxy - validates Bearer token in route handler
  "/api/proxy/compounds/", // Compounds API proxy - validates Bearer token in route handler
];

// Check if path is exempt from authentication
function isAuthExempt(pathname: string): boolean {
  return AUTH_EXEMPT_PATHS.some((path) => pathname.startsWith(path));
}

// Check if this is a static file request
function isStaticFile(pathname: string): boolean {
  return (
    pathname.endsWith(".js") ||
    pathname.endsWith(".wasm") ||
    pathname.endsWith(".css") ||
    pathname.endsWith(".json") ||
    pathname.endsWith(".png") ||
    pathname.endsWith(".svg") ||
    pathname.endsWith(".woff") ||
    pathname.endsWith(".woff2") ||
    pathname.endsWith(".data") ||
    pathname.endsWith(".gz") ||
    pathname.endsWith(".html")
  );
}

export function middleware(request: NextRequest) {
  const { pathname } = request.nextUrl;

  // Check if authentication is required (set by environment variable)
  const requireAuth = process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true";

  // Authentication check (if enabled)
  if (requireAuth && !isAuthExempt(pathname) && !isStaticFile(pathname)) {
    const authSession = request.cookies.get("auth-session");

    // Allow Bearer token auth for API proxy paths (CLI tools like i2remote)
    const hasBearerToken = request.headers.get("Authorization")?.startsWith("Bearer ");
    const isBearerAuthPath = BEARER_AUTH_PATHS.some((path) => pathname.startsWith(path));

    if (!authSession && !(hasBearerToken && isBearerAuthPath)) {
      // Redirect to login page with return URL
      const loginUrl = new URL("/auth/login", request.url);
      loginUrl.searchParams.set("returnUrl", pathname + request.nextUrl.search);
      return NextResponse.redirect(loginUrl);
    }
  }

  // Cross-Origin headers for static files (needed for Moorhen)
  if (isStaticFile(pathname)) {
    const response = NextResponse.next();
    response.headers.set("Cross-Origin-Resource-Policy", "cross-origin");
    return response;
  }

  // Full cross-origin isolation for Moorhen pages (SharedArrayBuffer support)
  if (pathname.startsWith("/ccp4i2/moorhen-page")) {
    const response = NextResponse.next();
    response.headers.set("Cross-Origin-Opener-Policy", "same-origin");
    response.headers.set("Cross-Origin-Embedder-Policy", "require-corp");
    response.headers.set("Cross-Origin-Resource-Policy", "same-origin");
    return response;
  }

  return NextResponse.next();
}

// Configure which paths the middleware runs on
export const config = {
  matcher: [
    // Match all paths except static assets (which are handled separately)
    "/((?!_next/static|_next/image).*)",
  ],
};
