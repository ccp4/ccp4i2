'use client';
import { PropsWithChildren } from 'react';

// RequireAuth is available from main ccp4i2 client when overlaid via Docker
// AuthProvider is in root layout so all routes share one MSAL instance
import RequireAuth from '../../components/require-auth';

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === 'true';

/**
 * Admin layout - wraps admin pages with authentication
 * when NEXT_PUBLIC_REQUIRE_AUTH is enabled.
 *
 * AuthProvider is at root level; this layout only adds RequireAuth gate.
 */
export default function AdminLayout({ children }: PropsWithChildren) {
  if (REQUIRE_AUTH) {
    return <RequireAuth>{children}</RequireAuth>;
  }

  return <>{children}</>;
}
