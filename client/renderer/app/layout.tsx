import "./globals.css";
import { PropsWithChildren } from "react";
import { CCP4i2ThemeProvider } from "../theme/theme-provider";
import AuthProvider from "../components/auth-provider";
import { RDKitProvider } from "../providers/rdkit-provider";

export const metadata = {
  title: "ccp4i2",
  description: "Software for Macromolecular X-Ray Crystallography",
};

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true";

/**
 * Root layout - wraps all apps with theme and auth providers.
 * AuthProvider is at root level so all routes share one MSAL instance,
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
              <AuthProvider>{props.children}</AuthProvider>
            ) : (
              props.children
            )}
          </RDKitProvider>
        </CCP4i2ThemeProvider>
      </body>
    </html>
  );
}
