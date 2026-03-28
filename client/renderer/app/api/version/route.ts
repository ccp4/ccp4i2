/*
 * Copyright (C) 2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
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
