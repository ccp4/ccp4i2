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
