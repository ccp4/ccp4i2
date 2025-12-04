import "./globals.css";
import { PropsWithChildren } from "react";
import { DeleteDialogProvider } from "../providers/delete-dialog";
import { CCP4i2ThemeProvider } from "../theme/theme-provider";
import { CCP4i2App } from "../providers/ccp4i2-app";
import AuthProvider from "../components/auth-provider";
import RequireAuth from "../components/require-auth";
export const metadata = {
  title: "CCP4",
  description: "Software for Macromolecular X-Ray Crystallography",
};

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true";

export default function RootLayout(props: PropsWithChildren) {
  return (
    <html lang="en">
      <body>
        {REQUIRE_AUTH ? (
          <AuthProvider>
            <CCP4i2ThemeProvider>
              <DeleteDialogProvider>
                <RequireAuth>
                  <CCP4i2App>{props.children}</CCP4i2App>
                </RequireAuth>
              </DeleteDialogProvider>
            </CCP4i2ThemeProvider>
          </AuthProvider>
        ) : (
          <CCP4i2ThemeProvider>
            <DeleteDialogProvider>
              <CCP4i2App>{props.children}</CCP4i2App>
            </DeleteDialogProvider>
          </CCP4i2ThemeProvider>
        )}
      </body>
    </html>
  );
}
