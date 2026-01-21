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
