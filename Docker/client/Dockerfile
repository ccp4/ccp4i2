# CCP4i2 Web Client Dockerfile
#
# This Dockerfile builds the Next.js web frontend for CCP4i2.
# It creates a standalone web deployment (not Electron).
#
# Build:
#   docker build -t ccp4i2-client -f Docker/client/Dockerfile .
#
# Run:
#   docker run -p 3000:3000 \
#     -e NEXT_PUBLIC_API_URL=http://localhost:8000 \
#     ccp4i2-client

# Stage 1: Dependencies
FROM node:20-alpine AS deps
WORKDIR /app

# Copy package files
COPY client/package.json client/package-lock.json ./

# Install dependencies
RUN npm ci --legacy-peer-deps

# Stage 2: Builder
FROM node:20-alpine AS builder
WORKDIR /app

# Copy dependencies from deps stage
COPY --from=deps /app/node_modules ./node_modules

# Copy client source
COPY client/ .

# Set build-time environment
ENV BUILD_TARGET=web \
    NEXT_TELEMETRY_DISABLED=1

# Copy assets (moorhen, rdkit, report files)
RUN npm run copy-assets 2>/dev/null || true

# Build the Next.js application
RUN npm run build:web

# Stage 3: Production runner
FROM node:20-alpine AS runner
WORKDIR /app

ENV NODE_ENV=production \
    NEXT_TELEMETRY_DISABLED=1 \
    PORT=3000

# Create non-root user
RUN addgroup --system --gid 1001 nodejs && \
    adduser --system --uid 1001 nextjs

# Copy built application
COPY --from=builder /app/renderer/.next/standalone ./
COPY --from=builder /app/renderer/.next/static ./renderer/.next/static
COPY --from=builder /app/renderer/public ./renderer/public

# Set ownership
RUN chown -R nextjs:nodejs /app

USER nextjs

EXPOSE 3000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=30s --retries=3 \
    CMD wget --no-verbose --tries=1 --spider http://localhost:3000/ || exit 1

CMD ["node", "renderer/server.js"]
