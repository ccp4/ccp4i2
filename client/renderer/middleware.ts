import { NextResponse } from "next/server";
import type { NextRequest } from "next/server";

/**
 * Middleware to add Cross-Origin headers for static files.
 *
 * Next.js's headers() config in next.config.ts doesn't apply to files in the
 * public/ directory - they're served directly without going through routing.
 * This middleware intercepts requests and adds the necessary headers.
 *
 * Key headers:
 * - Cross-Origin-Resource-Policy: cross-origin
 *   Allows resources to be loaded from pages with COEP: require-corp
 *   (needed for Moorhen's SharedArrayBuffer support)
 */
export function middleware(request: NextRequest) {
  const { pathname } = request.nextUrl;

  // Check if this is a static file that needs CORP headers
  const isStaticFile =
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
    pathname.endsWith(".html");

  // For static files, add CORP header
  if (isStaticFile) {
    const response = NextResponse.next();
    response.headers.set("Cross-Origin-Resource-Policy", "cross-origin");
    return response;
  }

  // For Moorhen pages, add full cross-origin isolation headers
  if (pathname.startsWith("/ccp4i2/moorhen-page")) {
    const response = NextResponse.next();
    response.headers.set("Cross-Origin-Opener-Policy", "same-origin");
    response.headers.set("Cross-Origin-Embedder-Policy", "require-corp");
    return response;
  }

  return NextResponse.next();
}

// Configure which paths the middleware runs on
export const config = {
  matcher: [
    // Moorhen pages need COEP/COOP headers
    "/ccp4i2/moorhen-page",
    "/ccp4i2/moorhen-page/:path*",
    // Static files that need CORP headers (explicit extensions)
    "/:path*.js",
    "/:path*.wasm",
    "/:path*.css",
    "/:path*.json",
    "/:path*.png",
    "/:path*.svg",
    "/:path*.woff",
    "/:path*.woff2",
    "/:path*.data",
    "/:path*.gz",
  ],
};
