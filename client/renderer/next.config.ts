import type { NextConfig } from "next";

const isElectron = process.env.BUILD_TARGET === "electron";
const isWeb = process.env.BUILD_TARGET === "web";
const isDevelopment = process.env.NODE_ENV === "development";

// Set the Content Security Policy header
const csp = {
  defaultSrc: "'self'",
  imgSrc: "'self' data: blob: https://*.blob.core.windows.net",
  connectSrc:
    "'self' https://www.ebi.ac.uk https://www.uniprot.org https://pubmed.ncbi.nlm.nih.gov https://raw.githubusercontent.com/MonomerLibrary/monomers/master/ " +
    "https://files.rcsb.org https://www.rcsb.org https://alphafold.ebi.ac.uk https://pdb-redo.eu https://cdn.jsdelivr.net " +
    "https://login.microsoftonline.com https://graph.microsoft.com https://*.microsoftonline.com https://*.microsoft.com " +
    "https://graph.windows.net https://management.azure.com " +
    "https://res.cdn.office.net " + // Teams SDK validDomains fetch
    process.env.NEXT_PUBLIC_API_BASE_URL,
  styleSrc:
    "'self' https://cdn.jsdelivr.net 'unsafe-inline' https://fonts.googleapis.com/css2",
  fontSrc:
    "'self' https://cdn.jsdelivr.net 'unsafe-inline' https://fonts.gstatic.com data:",
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

  // Increase body size limit for large fixture file imports (AssayCompounds can be >10MB)
  experimental: {
    serverActions: {
      bodySizeLimit: '100mb',
    },
    // For API routes (route handlers), use middlewareClientMaxBodySize
    middlewareClientMaxBodySize: '100mb',
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

  webpack: (config, { isServer, webpack }) => {
    config.module.parser.javascript.importMeta = false;
    if (isServer) {
      config.externals = [...config.externals, "moorhen"];
    } else {
      // pptxgenjs's ESM bundle imports node:fs and node:https at the top
      // for its Node.js execution mode. Those code paths are gated at
      // runtime and never run in the browser, but webpack walks the tree
      // at build time and fails on the unresolvable node: scheme. Stub
      // both the prefixed and bare forms so the client bundle compiles.
      config.resolve.fallback = {
        ...config.resolve.fallback,
        fs: false,
        https: false,
        http: false,
      };
      config.plugins.push(
        new webpack.NormalModuleReplacementPlugin(
          /^node:(fs|https|http)$/,
          (resource: { request: string }) => {
            resource.request = resource.request.replace(/^node:/, "");
          },
        ),
      );
    }
    // Note: we intentionally do NOT set config.target = "electron-renderer".
    // Our Electron uses contextIsolation:true + nodeIntegration:false, so the
    // renderer is a standard browser environment. The electron-renderer target
    // causes webpack to use the CommonJS branch of UMD modules (like moorhen),
    // which triggers deep require() parsing that corrupts export detection.
    return config;
  },

  async headers() {
    // Cross-origin isolation headers required for SharedArrayBuffer (Moorhen WebAssembly)
    // require-corp is needed for full Web Worker SharedArrayBuffer support
    // All resources must have Cross-Origin-Resource-Policy: cross-origin header
    const webAssemblyHeaders = [
      {
        key: "Cross-Origin-Opener-Policy",
        value: "same-origin",
      },
      {
        key: "Cross-Origin-Embedder-Policy",
        value: "require-corp",
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
        // Reinspect SPA is path-mounted at /reinspect/* and instantiates
        // Moorhen for event review; same COOP/COEP requirement as
        // /ccp4i2/moorhen-page above. Mirror it here so dev-mode parity
        // with production is maintained. See production block for full
        // rationale.
        {
          source: "/reinspect",
          headers: [...devCommonHeaders, ...webAssemblyHeaders],
        },
        {
          source: "/reinspect/:path*",
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
      // Reinspect SPA is path-mounted at /reinspect/* and instantiates Moorhen
      // (via CootWorker.js) when it shows event review pages. The SPA page
      // therefore needs the same cross-origin-isolation pair (COOP same-origin
      // + COEP require-corp) that /ccp4i2/moorhen-page gets, so the browser
      // accepts SharedArrayBuffer and Worker construction. The Reinspect
      // upstream sets these too (CLOUD_DEPLOYMENT.md "Cross-Origin headers")
      // and our /reinspect/* proxy preserves them, but emitting them at the
      // Materia edge as well guarantees consistency across the redirect chain
      // (the SPA index.html, the Reinspect API responses, asset 404 / error
      // pages, ...) -- without this, the WASM Worker load was blocked with
      // "Refused to load worker because of Cross-Origin-Embedder-Policy"
      // (BAZ2B smoke, 2026-06-07).
      {
        source: "/reinspect",
        headers: [...commonHeaders, ...webAssemblyHeaders],
      },
      {
        source: "/reinspect/:path*",
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
