import { NextResponse } from "next/server";

/**
 * Version endpoint returning build information.
 * Used by the app-selector to display deployment info.
 *
 * Environment variables (set at build time):
 * - NEXT_PUBLIC_BUILD_TIMESTAMP: When the image was built
 * - NEXT_PUBLIC_GIT_COMMIT: Git commit hash
 */
export async function GET() {
  return NextResponse.json(
    {
      web: {
        buildTimestamp: process.env.NEXT_PUBLIC_BUILD_TIMESTAMP || "dev",
        gitCommit: process.env.NEXT_PUBLIC_GIT_COMMIT || "unknown",
      },
    },
    { status: 200 }
  );
}
