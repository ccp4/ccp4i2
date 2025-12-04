"use client";
import { ReactNode, useEffect, useState } from "react";
import { MsalProvider } from "@azure/msal-react";
import { PublicClientApplication } from "@azure/msal-browser";

const msalConfig = {
  auth: {
    clientId: process.env.NEXT_PUBLIC_AAD_CLIENT_ID || "",
    authority: `https://login.microsoftonline.com/${process.env.NEXT_PUBLIC_AAD_TENANT_ID}`,
    redirectUri: "/",
  },
};

const pca = new PublicClientApplication(msalConfig);

interface AuthProviderProps {
  children: ReactNode;
}

export default function AuthProvider({ children }: AuthProviderProps) {
  const [initialized, setInitialized] = useState(false);

  useEffect(() => {
    pca
      .initialize()
      .then(() => {
        // Handle redirect promises when app loads (for redirect-based auth flows)
        return pca.handleRedirectPromise();
      })
      .then(() => {
        setInitialized(true);
      })
      .catch((error) => {
        console.error(
          "MSAL initialization or redirect handling failed:",
          error
        );
        setInitialized(true); // Initialize anyway to prevent blocking
      });
  }, []);

  if (!initialized) return null; // or a loading spinner

  return <MsalProvider instance={pca}>{children}</MsalProvider>;
}
