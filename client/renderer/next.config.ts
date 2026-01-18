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
  scriptSrc: "'self' https://cdn.jsdelivr.net 'unsafe-inline' 'unsafe-eval' blob:",
  workerSrc: "'self' blob:",
  frameSrc:
    "'self' https://login.microsoftonline.com https://*.microsoftonline.com",
  // Allow Teams to embed the app
  frameAncestors: "'self' https://teams.microsoft.com https://*.teams.microsoft.com https://*.office.com https://*.microsoft.com",
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
    // Cross-origin isolation headers required for SharedArrayBuffer (Moorhen WebAssembly)
    // Using "credentialless" instead of "require-corp" - it still enables SharedArrayBuffer
    // but doesn't require every sub-resource to have CORP headers
    const webAssemblyHeaders = [
      {
        key: "Cross-Origin-Opener-Policy",
        value: "same-origin",
      },
      {
        key: "Cross-Origin-Embedder-Policy",
        value: "credentialless",
      },
    ];

    // Development headers (more permissive)
    if (isDevelopment) {
      const devCommonHeaders = [
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
      ];

      // Cross-Origin-Resource-Policy: cross-origin allows loading from any context
      // This is safe and enables both COEP pages and cross-origin embedding
      const corpHeader = {
        key: "Cross-Origin-Resource-Policy",
        value: "cross-origin",
      };

      return [
        // All static assets get CORP: cross-origin to work with both COEP and Teams embedding
        {
          source: "/_next/:path*",
          headers: [corpHeader],
        },
        // Root-level static files (moorhen.js, CootWorker.js, etc.)
        {
          source: "/:file.js",
          headers: [corpHeader],
        },
        {
          source: "/:file.wasm",
          headers: [corpHeader],
        },
        // Nested static files
        {
          source: "/:path*/:file.js",
          headers: [corpHeader],
        },
        {
          source: "/:path*/:file.wasm",
          headers: [corpHeader],
        },
        {
          source: "/:path*.css",
          headers: [corpHeader],
        },
        {
          source: "/:path*.json",
          headers: [corpHeader],
        },
        {
          source: "/:path*.png",
          headers: [corpHeader],
        },
        {
          source: "/:path*.svg",
          headers: [corpHeader],
        },
        {
          source: "/:path*.woff",
          headers: [corpHeader],
        },
        {
          source: "/:path*.woff2",
          headers: [corpHeader],
        },
        // Data files for Moorhen (rotamer data, etc.)
        {
          source: "/:path*.data",
          headers: [corpHeader],
        },
        {
          source: "/:path*.gz",
          headers: [corpHeader],
        },
        // Moorhen pages need strict cross-origin isolation for WebAssembly
        {
          source: "/ccp4i2/moorhen-page",
          headers: [...devCommonHeaders, ...webAssemblyHeaders],
        },
        {
          source: "/ccp4i2/moorhen-page/:path*",
          headers: [...devCommonHeaders, ...webAssemblyHeaders],
        },
        // All other routes
        {
          source: "/:path*",
          headers: devCommonHeaders,
        },
      ];
    }

    // Common headers for all routes
    const commonHeaders = [
      { key: "Content-Security-Policy", value: cspString },
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
    ];

    // Cross-Origin-Resource-Policy: cross-origin allows loading from any context
    // This is safe and enables both COEP pages and Teams embedding
    const corpHeader = {
      key: "Cross-Origin-Resource-Policy",
      value: "cross-origin",
    };

    // Production headers with route-specific policies
    return [
      // Static assets need CORP header to be loadable from COEP-enabled pages
      // Apply to all common static asset patterns
      {
        source: "/_next/:path*",
        headers: [corpHeader],
      },
      // Root-level static files (moorhen.js, CootWorker.js, etc.)
      {
        source: "/:file.js",
        headers: [corpHeader],
      },
      {
        source: "/:file.wasm",
        headers: [corpHeader],
      },
      // Nested static files
      {
        source: "/:path*/:file.js",
        headers: [corpHeader],
      },
      {
        source: "/:path*/:file.wasm",
        headers: [corpHeader],
      },
      {
        source: "/:path*.css",
        headers: [corpHeader],
      },
      {
        source: "/:path*.json",
        headers: [corpHeader],
      },
      {
        source: "/:path*.png",
        headers: [corpHeader],
      },
      {
        source: "/:path*.svg",
        headers: [corpHeader],
      },
      {
        source: "/:path*.woff",
        headers: [corpHeader],
      },
      {
        source: "/:path*.woff2",
        headers: [corpHeader],
      },
      // Data files for Moorhen (rotamer data, etc.)
      {
        source: "/:path*.data",
        headers: [corpHeader],
      },
      {
        source: "/:path*.gz",
        headers: [corpHeader],
      },
      // Moorhen pages need strict cross-origin isolation for WebAssembly/SharedArrayBuffer
      // Need both patterns: one for /moorhen-page itself, one for subroutes
      {
        source: "/ccp4i2/moorhen-page",
        headers: [...commonHeaders, ...webAssemblyHeaders],
      },
      {
        source: "/ccp4i2/moorhen-page/:path*",
        headers: [...commonHeaders, ...webAssemblyHeaders],
      },
      // All other routes - no cross-origin isolation, allows Teams embedding
      {
        source: "/:path*",
        headers: commonHeaders,
      },
    ];
  },

  distDir: ".next",
};

export default nextConfig;
