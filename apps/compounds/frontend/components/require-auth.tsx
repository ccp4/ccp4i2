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

import { ReactNode } from 'react';

interface RequireAuthProps {
  children: ReactNode;
}

/**
 * Stub RequireAuth component for standalone compounds app development.
 *
 * This is a pass-through component that allows the app to build and run
 * without authentication. During Docker overlay for web deployment,
 * this file is overwritten by the real RequireAuth component from the
 * main ccp4i2 client (client/renderer/components/require-auth.tsx).
 *
 * The real component integrates with Azure AD MSAL and Teams membership checks.
 */
export default function RequireAuth({ children }: RequireAuthProps) {
  // In standalone mode, just render children without auth checks
  return <>{children}</>;
}
