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

ARG NEXT_PUBLIC_REQUIRE_AUTH=true
ARG NEXT_PUBLIC_BUILD_TIMESTAMP
ARG NEXT_PUBLIC_GIT_COMMIT
ENV NEXT_PUBLIC_REQUIRE_AUTH=$NEXT_PUBLIC_REQUIRE_AUTH
ENV NEXT_PUBLIC_BUILD_TIMESTAMP=$NEXT_PUBLIC_BUILD_TIMESTAMP
ENV NEXT_PUBLIC_GIT_COMMIT=$NEXT_PUBLIC_GIT_COMMIT

WORKDIR /app

# Copy package files and local dependency tarballs (build context is repo root)
COPY client/package*.json ./
COPY client/moorhen-0.23.1-alpha.0.tgz ./

# Remove Electron-specific packages and add markdown rendering packages for web-only build
# This significantly reduces install time and avoids downloading ~300MB of Electron binaries
# git is required by some npm packages during install (e.g., for resolving git dependencies)
RUN apk add --no-cache jq git && \
    jq 'del(.devDependencies.electron) | del(.devDependencies["electron-builder"]) | del(.devDependencies["@electron/packager"]) | del(.dependencies["electron-store"]) | .dependencies["react-markdown"] = "^9.0.1" | .dependencies["remark-gfm"] = "^4.0.0"' package.json > package.web.json && \
    mv package.web.json package.json

# Build the @ccp4/ccp4i2-auth shared package before consumer install.
# client/package.json declares it as `file:../packages/ccp4i2-auth`; npm
# symlinks the file: dep but does not run its build, so the package's
# lib/ must already exist before Next.js resolves @ccp4/ccp4i2-auth.
COPY packages/ccp4i2-auth /packages/ccp4i2-auth
RUN cd /packages/ccp4i2-auth && \
    npm install && \
    npm run build

# Install dependencies (excluding Electron) with retry logic
# Using npm install instead of npm ci to auto-update lock file when package.json changes
# --mount=type=cache persists npm cache between builds for faster installs
RUN --mount=type=cache,target=/root/.npm \
    npm config set fetch-retry-mintimeout 20000 && \
    npm config set fetch-retry-maxtimeout 120000 && \
    npm config set fetch-retries 5 && \
    npm install --legacy-peer-deps || npm install --legacy-peer-deps

# Copy client source code
COPY client/ .

# Copy icons to Next.js public directory for direct serving
COPY client/assets/qticons/ ./renderer/public/qticons/
COPY client/assets/svgicons/ ./renderer/public/svgicons/

# Build for web (from client root directory as originally intended)
RUN npx cross-env BUILD_TARGET=web npm run copy-assets

RUN cp package.json renderer/

RUN cd renderer && npx next build

FROM node:20-alpine3.20 AS runner

WORKDIR /app

# Add curl for health checks
RUN apk add --no-cache curl

# Copy the built application and local dependency tarballs
COPY --from=builder /app/renderer ./renderer
COPY --from=builder /app/package*.json ./
COPY --from=builder /app/moorhen-0.23.1-alpha.0.tgz ./

# Install production dependencies only
# --mount=type=cache persists npm cache between builds
RUN --mount=type=cache,target=/root/.npm \
    npm install --omit=dev --legacy-peer-deps

# Create a non-root user
RUN addgroup -g 1001 -S nodejs
RUN adduser -S nextjs -u 1001

# Change ownership of the app directory
RUN chown -R nextjs:nodejs /app
USER nextjs

EXPOSE 3000

# Health check - uses /api/health which bypasses auth middleware
HEALTHCHECK --interval=30s --timeout=10s --start-period=30s --retries=3 \
    CMD curl -f http://localhost:3000/api/health || exit 1

# Use ENTRYPOINT with CMD for better Azure Container Apps compatibility
ENTRYPOINT ["sh", "-c"]
CMD ["cd renderer && exec npx next start"]
