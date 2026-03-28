/*
 * Copyright (C) 2025-2026 Newcastle University
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
import "./globals.css";
import { PropsWithChildren, Suspense } from "react";
import { CCP4i2ThemeProvider } from "../theme/theme-provider";
import MsalAuthProvider from "../components/auth-provider";
import { RDKitProvider } from "../providers/rdkit-provider";
import { AuthProvider as CompoundsAuthProvider } from "../lib/compounds/auth-context";
import { TeamsRoutePersistence } from "../components/teams-route-persistence";

export const metadata = {
  title: "ccp4i2",
  description: "Software for Macromolecular X-Ray Crystallography",
};

// Skip static prerendering of all pages. This works around a Next.js 15.5 bug
// where prerendering fails with "Expected workUnitAsyncStorage to have a store"
// (tracked in vercel/next.js#84026, #87719). This is safe for our use case:
// Electron runs a local Next.js server, and Azure serves pages dynamically.
export const dynamic = "force-dynamic";

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true";

/**
 * Root layout - wraps all apps with theme and auth providers.
 *
 * Two auth providers are used:
 * - MsalAuthProvider: Azure AD authentication (gets access tokens)
 * - CompoundsAuthProvider: User info from Django API (role, operating level)
 *
 * Both are at root level so all routes share the same auth state,
 * enabling seamless navigation between ccp4i2, registry, and assays.
 * RDKitProvider enables molecule rendering in compounds app.
 */
export default function RootLayout(props: PropsWithChildren) {
  return (
    <html lang="en">
      <body>
        <CCP4i2ThemeProvider>
          <RDKitProvider>
            {REQUIRE_AUTH ? (
              <MsalAuthProvider>
                <CompoundsAuthProvider>
                  <Suspense fallback={null}>
                    <TeamsRoutePersistence />
                  </Suspense>
                  {props.children}
                </CompoundsAuthProvider>
              </MsalAuthProvider>
            ) : (
              props.children
            )}
          </RDKitProvider>
        </CCP4i2ThemeProvider>
      </body>
    </html>
  );
}
