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

/**
 * Runtime auth configuration.
 *
 * Fetches Azure AD config from /api/auth/config at runtime instead of
 * reading NEXT_PUBLIC_* env vars that are baked in at build time. This
 * allows a single Docker image to serve multiple instances with
 * different app registrations.
 *
 * Falls back to process.env.NEXT_PUBLIC_* for local dev and Electron.
 */

export interface AuthConfig {
  clientId: string;
  tenantId: string;
  requireAuth: boolean;
}

let cachedConfig: AuthConfig | null = null;

/**
 * Fetch auth configuration from the server.
 * Caches the result for the lifetime of the page.
 * Falls back to NEXT_PUBLIC_ env vars for local development / Electron.
 */
export async function getAuthConfig(): Promise<AuthConfig> {
  if (cachedConfig) return cachedConfig;

  try {
    const res = await fetch("/api/auth/config");
    if (res.ok) {
      cachedConfig = await res.json();
      return cachedConfig!;
    }
  } catch {
    // Fetch failed (e.g., during SSR or offline) — fall through to env vars
  }

  // Fallback to build-time env vars (local dev, Electron)
  cachedConfig = {
    clientId: process.env.NEXT_PUBLIC_AAD_CLIENT_ID || "",
    tenantId: process.env.NEXT_PUBLIC_AAD_TENANT_ID || "",
    requireAuth: process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true",
  };
  return cachedConfig;
}

/**
 * Synchronous access to the cached auth config.
 * Only valid after getAuthConfig() has been called (i.e., after
 * AuthProvider has initialised). Falls back to env vars if the
 * cache hasn't been populated yet.
 */
export function getAuthConfigSync(): AuthConfig {
  if (cachedConfig) return cachedConfig;

  return {
    clientId: process.env.NEXT_PUBLIC_AAD_CLIENT_ID || "",
    tenantId: process.env.NEXT_PUBLIC_AAD_TENANT_ID || "",
    requireAuth: process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true",
  };
}
