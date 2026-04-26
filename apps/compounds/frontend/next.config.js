/** @type {import('next').NextConfig} */
const nextConfig = {
  // API proxy is handled by app/api/proxy/compounds/[...path]/route.ts
  // This gives us more control and better error handling than rewrites

  // Increase body size limit for large fixture file imports (AssayCompounds can be >10MB)
  experimental: {
    serverActions: {
      bodySizeLimit: '100mb',
    },
    // For API routes (route handlers), use middlewareClientMaxBodySize
    middlewareClientMaxBodySize: '100mb',
  },

  webpack: (config, { isServer, webpack }) => {
    if (!isServer) {
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
          (resource) => {
            resource.request = resource.request.replace(/^node:/, '');
          },
        ),
      );
    }
    return config;
  },

  // Allow images from Django backend
  images: {
    remotePatterns: [
      {
        protocol: 'http',
        hostname: 'localhost',
        port: '8000',
        pathname: '/media/**',
      },
    ],
  },

  // Proxy /drf/ to Django for the DRF browsable API
  async rewrites() {
    return [
      {
        source: '/drf/:path*',
        destination: 'http://localhost:8000/api/compounds/:path*',
      },
    ];
  },
};

module.exports = nextConfig;
