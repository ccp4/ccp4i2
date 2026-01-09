import type { NextConfig } from "next";

const isElectron = process.env.BUILD_TARGET === "electron";
const isWeb = process.env.BUILD_TARGET === "web";
const isDevelopment = process.env.NODE_ENV === "development";

// Set the Content Security Policy header
const csp = {
  defaultSrc: "'self'",
  imgSrc: "'self' data: blob:",
  connectSrc:
    "'self' https://www.ebi.ac.uk https://www.uniprot.org https://pubmed.ncbi.nlm.nih.gov https://raw.githubusercontent.com/MonomerLibrary/monomers/master/ " +
    "https://login.microsoftonline.com https://graph.microsoft.com https://*.microsoftonline.com https://*.microsoft.com " +
    "https://graph.windows.net https://management.azure.com " +
    process.env.NEXT_PUBLIC_API_BASE_URL,
  styleSrc:
    "'self' https://cdn.jsdelivr.net 'unsafe-inline' https://fonts.googleapis.com/css2",
  fontSrc:
    "'self' https://cdn.jsdelivr.net 'unsafe-inline' https://fonts.gstatic.com",
  scriptSrc: "'self' https://cdn.jsdelivr.net 'unsafe-inline' 'unsafe-eval'",
  workerSrc: "'self' blob:",
  frameSrc:
    "'self' https://login.microsoftonline.com https://*.microsoftonline.com",
  frameAncestors: "'self'",
};

const cspString = [
  `default-src ${csp.defaultSrc}`,
  `img-src ${csp.imgSrc}`,
  `connect-src ${csp.connectSrc}`,
  `style-src ${csp.styleSrc}`,
  `font-src ${csp.fontSrc}`,
  `script-src ${csp.scriptSrc}`,
  `worker-src ${csp.workerSrc}`,
  `frame-src ${csp.frameSrc}`,
  `frame-ancestors ${csp.frameAncestors}`,
]
  .join("; ")
  .trim(); // Clean up whitespace

const nextConfig: NextConfig = {
  trailingSlash: isElectron,
  images: {
    unoptimized: isElectron || isWeb,
  },

  // No basePath - routes are organized at app level (/ccp4i2/*, /compounds/*)
  // This allows multiple apps in the same Next.js deployment

  rewrites: isWeb
    ? async () => [
        {
          // Proxy API requests to Django backend
          // /api/ccp4i2/* -> Django /api/ccp4i2/*
          source: "/api/ccp4i2/:path*",
          destination:
            (process.env.NEXT_PUBLIC_API_BASE_URL || "http://server:8000") +
            "/api/ccp4i2/:path*",
        },
        // Future: /api/compounds/* -> compounds backend
      ]
    : undefined,

  webpack: (config, { isServer }) => {
    config.module.parser.javascript.importMeta = false;
    if (isServer) {
      config.externals = [...config.externals, "moorhen"];
    }
    if (isElectron && !isServer) {
      config.target = "electron-renderer";
    }
    return config;
  },

  async headers() {
    // Development headers (more permissive)
    if (isDevelopment) {
      return [
        {
          source: "/:path*",
          headers: [
            { key: "Content-Security-Policy", value: cspString },
            { key: "Access-Control-Allow-Origin", value: "*" },
            {
              key: "Access-Control-Allow-Headers",
              value:
                "Origin, X-Requested-With, Content, Accept, Content-Type, Authorization",
            },
            {
              key: "Access-Control-Allow-Methods",
              value: "GET, POST, PUT, DELETE, PATCH, OPTIONS",
            },
          ],
        },
      ];
    }

    // Production headers (WebAssembly-compatible)
    return [
      {
        source: "/:path*",
        headers: [
          { key: "Content-Security-Policy", value: cspString },
          // Restrictive CORS for production
          {
            key: "Access-Control-Allow-Origin",
            value: process.env.CORS_ALLOWED_ORIGINS || "https://yourdomain.com",
          },
          {
            key: "Access-Control-Allow-Headers",
            value:
              "Origin, X-Requested-With, Content, Accept, Content-Type, Authorization",
          },
          {
            key: "Access-Control-Allow-Methods",
            value: "GET, POST, PUT, DELETE, PATCH, OPTIONS",
          },
          // Cross-origin isolation for WebAssembly
          {
            key: "Cross-Origin-Opener-Policy",
            value: "same-origin",
          },
          {
            key: "Cross-Origin-Embedder-Policy",
            value: "require-corp",
          },
        ],
      },
    ];
  },

  distDir: ".next",
};

export default nextConfig;
