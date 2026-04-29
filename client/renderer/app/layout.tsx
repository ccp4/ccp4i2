import "./globals.css";
import { PropsWithChildren, Suspense } from "react";
import { CCP4i2ThemeProvider } from "../theme/theme-provider";
import MsalAuthProvider from "../components/auth-provider";
import { RDKitProvider } from "../providers/rdkit-provider";
import { AuthProvider as CompoundsAuthProvider } from "../lib/compounds/auth-context";
import { TeamsRoutePersistence } from "../components/teams-route-persistence";

// Per-instance title read at request time from INSTANCE_TITLE (server-only env
// var — NEXT_PUBLIC_* would bake in at build). Lets a single Docker image serve
// differently-branded deployments. Falls back to "CCP4i2".
export async function generateMetadata() {
  return {
    title: process.env.INSTANCE_TITLE || "CCP4i2",
    description: "Software for Macromolecular X-Ray Crystallography",
  };
}

// Skip static prerendering of all pages. This works around a Next.js 15.5 bug
// where prerendering fails with "Expected workUnitAsyncStorage to have a store"
// (tracked in vercel/next.js#84026, #87719). This is safe for our use case:
// Electron runs a local Next.js server, and Azure serves pages dynamically.
export const dynamic = "force-dynamic";

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true";
// Electron desktop sets CCP4I2_LOCAL_SESSION_TOKEN before spawning the
// Next.js server, so the env var is visible during SSR. When present,
// the renderer needs MsalAuthProvider mounted so it can wire up the
// LocalSession token-getter (Django's LocalSessionAuthMiddleware will
// reject every API call without a Bearer header).
const HAS_LOCAL_SESSION = !!process.env.CCP4I2_LOCAL_SESSION_TOKEN;
const NEEDS_AUTH_PROVIDER = REQUIRE_AUTH || HAS_LOCAL_SESSION;

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
            {NEEDS_AUTH_PROVIDER ? (
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
