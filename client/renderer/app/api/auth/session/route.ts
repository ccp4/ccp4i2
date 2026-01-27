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
  response.cookies.set("auth-session", "authenticated", {
    httpOnly: true,
    secure: process.env.NODE_ENV === "production",
    sameSite: "lax",
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
