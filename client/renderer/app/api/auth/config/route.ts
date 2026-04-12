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
 * Auth configuration endpoint.
 *
 * Serves Azure AD configuration at runtime so it doesn't need to be
 * baked into the JS bundle at build time. This allows a single Docker
 * image to be deployed to multiple instances with different app
 * registrations — the container's environment variables determine the
 * auth config at request time.
 *
 * Environment variables (set at container runtime, NOT build time):
 * - NEXT_PUBLIC_AAD_CLIENT_ID: Azure AD application client ID
 * - NEXT_PUBLIC_AAD_TENANT_ID: Azure AD tenant ID
 * - NEXT_PUBLIC_REQUIRE_AUTH: Whether auth is required ("true"/"false")
 */
export async function GET() {
  return NextResponse.json(
    {
      clientId: process.env.NEXT_PUBLIC_AAD_CLIENT_ID || "",
      tenantId: process.env.NEXT_PUBLIC_AAD_TENANT_ID || "",
      requireAuth: process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true",
    },
    {
      headers: {
        "Cache-Control": "public, max-age=3600, s-maxage=3600",
      },
    }
  );
}
