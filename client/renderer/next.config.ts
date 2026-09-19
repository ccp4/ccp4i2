import type { NextConfig } from "next";
import { existsSync, readFileSync } from "fs";
import * as nodePath from "path";

const isElectron = process.env.BUILD_TARGET === "electron";
const isWeb = process.env.BUILD_TARGET === "web";
const isDevelopment = process.env.NODE_ENV === "development";

// The Moorhen package version, read from its package.json at build time and
// handed to the client as NEXT_PUBLIC_MOORHEN_VERSION.
//
// It becomes a path segment in the URL the web build uses for Moorhen's runtime
// assets (see lib/moorhen-asset-path.ts). Those responses are served
// `immutable` for a year -- necessary, because one window loads moorhen.js 33
// times for its pthread workers -- so the URL has to change when the bytes do,
// or a returning browser keeps last release's worker and WASM against the
// current React library. Resolved defensively: a build that cannot read the
// package falls back to the unversioned path rather than failing.
const moorhenVersion = (() => {
  // Read package.json off disk rather than `require("moorhen/package.json")`:
  // Moorhen's package has an `exports` map that does not expose
  // ./package.json, so both require and require.resolve throw
  // ERR_PACKAGE_PATH_NOT_EXPORTED. Walk up from this file to find the
  // node_modules that holds it (it is installed in client/, while this config
  // lives in client/renderer/), so the lookup does not depend on cwd.
  let dir = __dirname;
  for (let up = 0; up < 5; up++) {
    const candidate = nodePath.join(dir, "node_modules", "moorhen", "package.json");
    if (existsSync(candidate)) {
      const version = JSON.parse(readFileSync(candidate, "utf8")).version;
      if (typeof version === "string" && version) return version;
      break;
    }
    const parent = nodePath.dirname(dir);
    if (parent === dir) break;
    dir = parent;
  }
  // Unversioned assets still work; they just lose the cache-busting, so say so
  // loudly rather than shipping a silently broken cache policy.
  console.warn(
    "[next.config] Moorhen package.json not found: asset URLs will be unversioned " +
      "and the immutable cache will not bust on upgrade."
  );
  return "";
})();

// Set the Content Security Policy header
const csp = {
  defaultSrc: "'self'",
  // pdb-redo.eu: the PDB-REDO results page (served from the job directory)
  // renders status icons straight from the databank host.
  imgSrc: "'self' data: blob: https://*.blob.core.windows.net https://pdb-redo.eu",
  connectSrc:
    "'self' https://www.ebi.ac.uk https://www.uniprot.org https://pubmed.ncbi.nlm.nih.gov https://raw.githubusercontent.com/MonomerLibrary/monomers/master/ " +
    "https://files.rcsb.org https://www.rcsb.org https://alphafold.ebi.ac.uk https://pdb-redo.eu https://cdn.jsdelivr.net " +
    "https://login.microsoftonline.com https://graph.microsoft.com https://*.microsoftonline.com https://*.microsoft.com " +
    "https://graph.windows.net https://management.azure.com " +
    "https://res.cdn.office.net " + // Teams SDK validDomains fetch
    process.env.NEXT_PUBLIC_API_BASE_URL,
  styleSrc:
    "'self' https://cdn.jsdelivr.net 'unsafe-inline' https://fonts.googleapis.com",
  // Host, not a path: a path-scoped source only matches URLs starting
  // with it, so "/css2" silently blocked the older "/css?family="
  // endpoint that embedded task reports (MrParse) still use.
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
  env: {
    NEXT_PUBLIC_MOORHEN_VERSION: moorhenVersion,
  },

  trailingSlash: isElectron,
  images: {
    unoptimized: isElectron || isWeb,
  },

  // Body-size caps for uploads. Two of them: middlewareClientMaxBodySize (every
  // request through middleware has its body cloned with this cap; over-cap bodies
  // are silently truncated, Next's default 10 MB) and serverActions.bodySizeLimit.
  // Only the WEB (cloud) build is held at 100 MB -- its cloned body sits in memory
  // and a larger cap widens the denial-of-service surface. The desktop app (and
  // dev) is a single local user importing their own project zips and cryo-EM maps,
  // which run to hundreds of MB (a cryoSPARC volume is ~0.25 GB), so it gets a
  // generous 1900 MB ceiling. Verified: a 250 MB map uploads through both caps.
  //
  // The gate is `isWeb`, not `isElectron`: only the cloud is built BUILD_TARGET=web,
  // and `start:electron` doesn't set BUILD_TARGET at runtime, so keying off isWeb
  // makes packaged AND dev generous while leaving only the web build capped. (1900,
  // not 2 GB: '2gb' is exactly 2^31 bytes -- best avoided. Short-term fix; the
  // proper one is the desktop local-path bypass, #512.)
  //
  // Note for the desktop app: this file is not shipped in the package, so the
  // runtime gets these values through .next/required-server-files.json (see
  // client/main/ccp4i2-next-config.ts) - editing here is still the right place.
  experimental: {
    serverActions: {
      bodySizeLimit: isWeb ? '100mb' : '1900mb',
    },
    middlewareClientMaxBodySize: isWeb ? '100mb' : '1900mb',
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
    // Import .md files as raw source strings (e.g. the scene grammar embedded
    // into the LLM "Copy prompt" scaffold). No .md import is treated as a module
    // elsewhere, so this is safe across the whole tree.
    config.module.rules.push({
      test: /\.md$/,
      type: "asset/source",
    });
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
