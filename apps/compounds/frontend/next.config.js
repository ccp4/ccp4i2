/** @type {import('next').NextConfig} */
const nextConfig = {
  // API proxy is handled by app/api/proxy/compounds/[...path]/route.ts
  // This gives us more control and better error handling than rewrites

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
};

module.exports = nextConfig;
