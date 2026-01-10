'use client';
import { PropsWithChildren } from 'react';

// These components exist in the main ccp4i2 client and will be available
// when this layout is overlaid via Docker
import AuthProvider from '../../components/auth-provider';
import RequireAuth from '../../components/require-auth';

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === 'true';

/**
 * Registry layout - wraps all registry pages with authentication
 * when NEXT_PUBLIC_REQUIRE_AUTH is enabled.
 *
 * This uses the same auth components as ccp4i2 to ensure
 * the compound registry is protected with the same Azure AD + Teams
 * membership checks.
 */
export default function RegistryLayout({ children }: PropsWithChildren) {
  if (REQUIRE_AUTH) {
    return (
      <AuthProvider>
        <RequireAuth>
          {children}
        </RequireAuth>
      </AuthProvider>
    );
  }

  return <>{children}</>;
}
