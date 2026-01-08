import "./globals.css";
import { PropsWithChildren } from "react";
import { CCP4i2ThemeProvider } from "../theme/theme-provider";

export const metadata = {
  title: "CCP4 Cloud",
  description: "Software for Macromolecular X-Ray Crystallography",
};

/**
 * Root layout - minimal wrapper for all apps
 * App-specific providers (auth, etc.) are in each app's layout
 */
export default function RootLayout(props: PropsWithChildren) {
  return (
    <html lang="en">
      <body>
        <CCP4i2ThemeProvider>
          {props.children}
        </CCP4i2ThemeProvider>
      </body>
    </html>
  );
}
