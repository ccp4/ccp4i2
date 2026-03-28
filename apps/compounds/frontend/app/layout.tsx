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

import './globals.css';
import { PropsWithChildren } from 'react';
import { AppRouterCacheProvider } from '@mui/material-nextjs/v15-appRouter';
import { RDKitProvider } from '@/lib/compounds/rdkit-context';
import { CompoundsThemeProvider } from '@/lib/compounds/theme-provider';
import { AuthProvider } from '@/lib/compounds/auth-context';

export default function RootLayout({ children }: PropsWithChildren) {
  return (
    <html lang="en">
      <body>
        <AppRouterCacheProvider>
          <CompoundsThemeProvider>
            <AuthProvider>
              <RDKitProvider>
                {children}
              </RDKitProvider>
            </AuthProvider>
          </CompoundsThemeProvider>
        </AppRouterCacheProvider>
      </body>
    </html>
  );
}
