# CCP4i2 Web Client Dockerfile
#
# This Dockerfile builds the Next.js web frontend for CCP4i2.
# It creates a web deployment (not Electron).
#
# Build (from repo root):
#   docker build -t ccp4i2-client -f Docker/client/Dockerfile .
#
# Run:
#   docker run -p 3000:3000 \
#     -e NEXT_PUBLIC_API_URL=http://localhost:8000 \
#     ccp4i2-client

FROM node:20-alpine3.20 AS builder

ARG NEXT_PUBLIC_AAD_CLIENT_ID
ARG NEXT_PUBLIC_AAD_TENANT_ID
ARG NEXT_PUBLIC_REQUIRE_AUTH=true
ENV NEXT_PUBLIC_AAD_CLIENT_ID=$NEXT_PUBLIC_AAD_CLIENT_ID
ENV NEXT_PUBLIC_AAD_TENANT_ID=$NEXT_PUBLIC_AAD_TENANT_ID
ENV NEXT_PUBLIC_REQUIRE_AUTH=$NEXT_PUBLIC_REQUIRE_AUTH

WORKDIR /app

# Copy package files (build context is repo root)
COPY client/package*.json ./

# Remove Electron-specific packages from package.json for web-only build
# This significantly reduces install time and avoids downloading ~300MB of Electron binaries
RUN apk add --no-cache jq && \
    jq 'del(.devDependencies.electron) | del(.devDependencies["electron-builder"]) | del(.devDependencies["@electron/packager"]) | del(.dependencies["electron-store"])' package.json > package.web.json && \
    mv package.web.json package.json

# Install dependencies (excluding Electron) with retry logic
# Using npm install instead of npm ci to auto-update lock file when package.json changes
RUN npm config set fetch-retry-mintimeout 20000 && \
    npm config set fetch-retry-maxtimeout 120000 && \
    npm config set fetch-retries 5 && \
    npm install --legacy-peer-deps --verbose || (npm cache clean --force && npm install --legacy-peer-deps --verbose)

# Copy client source code
COPY client/ .

# Copy icons to Next.js public directory for direct serving
# This avoids the Django proxy hop for static assets
COPY server/ccp4i2/qticons/ ./renderer/public/qticons/
COPY server/ccp4i2/svgicons/ ./renderer/public/svgicons/

# =============================================================================
# Overlay compounds frontend (registry, assays & constructs sub-apps)
# =============================================================================
# Copy app routes
COPY apps/compounds/frontend/app/registry/ ./renderer/app/registry/
COPY apps/compounds/frontend/app/assays/ ./renderer/app/assays/
COPY apps/compounds/frontend/app/constructs/ ./renderer/app/constructs/
COPY apps/compounds/frontend/app/api/proxy/compounds/ ./renderer/app/api/proxy/compounds/

# Copy app selector as the root page (replaces the redirect-only page)
COPY apps/compounds/frontend/app/app-selector/page.tsx ./renderer/app/page.tsx

# Copy components, lib, and types (already in compounds/ subfolders)
COPY apps/compounds/frontend/components/compounds/ ./renderer/components/compounds/
COPY apps/compounds/frontend/lib/compounds/ ./renderer/lib/compounds/
COPY apps/compounds/frontend/types/compounds/ ./renderer/types/compounds/

# =============================================================================
# Overlay users frontend (platform admin)
# =============================================================================
COPY apps/users/frontend/app/admin/ ./renderer/app/admin/
COPY apps/users/frontend/app/api/proxy/users/ ./renderer/app/api/proxy/users/
COPY apps/users/frontend/lib/users/ ./renderer/lib/users/

# Add path aliases to tsconfig.json for @/ imports used by compounds
# The ccp4i2 client uses relative imports, but compounds uses @/ aliases
RUN apk add --no-cache jq && \
    jq '.compilerOptions += {"baseUrl": ".", "paths": {"@/*": ["./*"]}}' \
    ./renderer/tsconfig.json > ./renderer/tsconfig.new.json && \
    mv ./renderer/tsconfig.new.json ./renderer/tsconfig.json

# Build for web (from client root directory as originally intended)
RUN npx cross-env BUILD_TARGET=web npm run copy-assets

RUN cp package.json renderer/

RUN cd renderer && npx next build

FROM node:20-alpine3.20 AS runner

WORKDIR /app

# Add curl for health checks
RUN apk add --no-cache curl

# Copy the built application
COPY --from=builder /app/renderer ./renderer
COPY --from=builder /app/package*.json ./

# Install production dependencies only
RUN npm install --omit=dev --legacy-peer-deps && npm cache clean --force

# Create a non-root user
RUN addgroup -g 1001 -S nodejs
RUN adduser -S nextjs -u 1001

# Change ownership of the app directory
RUN chown -R nextjs:nodejs /app
USER nextjs

EXPOSE 3000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=30s --retries=3 \
    CMD curl -f http://localhost:3000/ || exit 1

# Use ENTRYPOINT with CMD for better Azure Container Apps compatibility
ENTRYPOINT ["sh", "-c"]
CMD ["cd renderer && exec npx next start"]
