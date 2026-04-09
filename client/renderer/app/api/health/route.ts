import { NextResponse } from "next/server";

/**
 * Health check endpoint for Azure Container Apps probes.
 * This endpoint must be accessible without authentication.
 */
export async function GET() {
  return NextResponse.json(
    {
      status: "healthy",
      timestamp: new Date().toISOString(),
    },
    { status: 200 }
  );
}
