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
'use client';
import { PropsWithChildren } from 'react';

// RequireAuth is available from main ccp4i2 client when overlaid via Docker
// AuthProvider is in root layout so all routes share one MSAL instance
import RequireAuth from '../../components/require-auth';

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === 'true';

/**
 * Registry layout - wraps all registry pages with authentication
 * when NEXT_PUBLIC_REQUIRE_AUTH is enabled.
 *
 * AuthProvider is at root level; this layout only adds RequireAuth gate.
 */
export default function RegistryLayout({ children }: PropsWithChildren) {
  if (REQUIRE_AUTH) {
    return <RequireAuth>{children}</RequireAuth>;
  }

  return <>{children}</>;
}
