import { NextResponse } from "next/server";

/**
 * Auth session cookie management.
 *
 * POST: Set the auth-session cookie after successful MSAL login
 * DELETE: Clear the auth-session cookie on logout
 *
 * This cookie allows the middleware to gate all requests server-side,
 * preventing any app code from loading for unauthenticated users.
 */

export async function POST() {
  const response = NextResponse.json({ success: true });

  // Set HTTP-only cookie that middleware can verify
  // maxAge matches typical Azure AD token lifetime (8 hours)
  //
  // IMPORTANT: Use sameSite: "none" in production for Teams iframe support.
  // With "lax", cookies aren't sent in third-party iframe contexts, causing login loops.
  // sameSite: "none" requires secure: true (enforced by browsers).
  const isProduction = process.env.NODE_ENV === "production";

  response.cookies.set("auth-session", "authenticated", {
    httpOnly: true,
    secure: isProduction,
    sameSite: isProduction ? "none" : "lax",
    path: "/",
    maxAge: 60 * 60 * 8, // 8 hours
  });

  return response;
}

export async function DELETE() {
  const response = NextResponse.json({ success: true });
  response.cookies.delete("auth-session");
  return response;
}
