# CCP4 Distribution Curation Guide

This document describes the strategy and procedures for managing CCP4 software distributions in the CCP4i2 Docker deployment infrastructure.

## Overview

CCP4i2 uses a **layered Docker image strategy** to efficiently manage the large CCP4 software suite (~10-15GB). Rather than rebuilding the entire stack for each code change, CCP4 and ARP/wARP are packaged as base images that change infrequently, while application code builds as thin layers on top.

### Image Hierarchy

```
┌─────────────────────────────────────────────────────────┐
│  Layer 4: Application Code (ccp4i2/server:TIMESTAMP)    │  ← Frequent rebuilds
│  Django app, Azure extensions, Compounds overlay        │     (~2-3GB total)
├─────────────────────────────────────────────────────────┤
│  Layer 3: ARP/wARP (ccp4i2/base-arpwarp:ccp4-YYYYMMDD) │  ← Optional layer
│  ARP/wARP molecular replacement suite                   │     (~1-2GB)
├─────────────────────────────────────────────────────────┤
│  Layer 2: CCP4 Suite (ccp4i2/base:ccp4-YYYYMMDD)       │  ← Rarely changes
│  Full CCP4 distribution, dereferenced                   │     (~10-15GB)
├─────────────────────────────────────────────────────────┤
│  Layer 1: OS Base (python:3.11-slim + system libs)      │  ← Stable
└─────────────────────────────────────────────────────────┘
```

## Directory Structure

```
Docker/
├── base/
│   ├── Dockerfile.ccp4-base      # CCP4 base image definition
│   └── Dockerfile.arpwarp        # ARP/wARP layer definition
│
├── azure/
│   └── scripts/
│       ├── build-base-image.sh   # Build CCP4 base image
│       └── upload-ccp4data.sh    # Upload CCP4 to Azure Files
│
├── azure-uksouth/
│   └── scripts/
│       ├── build-base-image.sh   # Regional variant
│       └── build-arpwarp-image.sh # ARP/wARP layer build
│
└── server/
    ├── Dockerfile                # Server (CCP4 mounted at runtime)
    └── Dockerfile.with-ccp4      # Server (CCP4 bundled in image)
```

## Preparing Dereferenced Tarballs

### Why Dereference?

CCP4 installations contain numerous symbolic links. These links:
- Won't work across different platforms (e.g., macOS → Linux)
- May break when extracted in Docker build contexts
- Reference relative paths that don't exist in containers

**Dereferencing** converts all symlinks to actual file copies, creating a self-contained, portable distribution.

### Creating a Dereferenced CCP4 Tarball

```bash
# Navigate to the directory containing the CCP4 installation
cd ~/Developer/ccp4data

# Create dereferenced tarball (preserves directory structure)
tar --dereference -czf ccp4-20251105-dereferenced.tgz ccp4-20251105/

# Verify the tarball
tar -tzf ccp4-20251105-dereferenced.tgz | head -20
```

**Expected size:** 10-15GB depending on CCP4 version

### Creating a Dereferenced ARP/wARP Tarball

```bash
cd ~/Developer/ccp4data

# Create dereferenced tarball
tar --dereference -czf arp-warp-8.0-dereferenced.tgz arp_warp_8.0/

# Verify
tar -tzf arp-warp-8.0-dereferenced.tgz | head -20
```

**Expected size:** 1-2GB

## Building Base Images

### Step 1: Build CCP4 Base Image

The base image contains only the CCP4 suite. It's built once per CCP4 version and reused by all server images.

```bash
# Source the deployment environment
. ./Docker/azure/.env.deployment

# Build the base image
./Docker/azure/scripts/build-base-image.sh \
    ~/Developer/ccp4data/ccp4-20251105-dereferenced.tgz \
    ccp4-20251105
```

**What this does:**
1. Creates a temporary build context directory
2. Extracts the dereferenced tarball into the context
3. Copies `Docker/base/Dockerfile.ccp4-base` to the context
4. Builds the Docker image locally
5. Tags as `${ACR_LOGIN_SERVER}/ccp4i2/base:ccp4-20251105`
6. Pushes to Azure Container Registry

**Build time:** 10-30 minutes (first push), subsequent builds use cached layers

### Step 2: Build ARP/wARP Layer (Optional)

If ARP/wARP functionality is required, build it as a separate layer on top of the CCP4 base:

```bash
./Docker/azure-uksouth/scripts/build-arpwarp-image.sh \
    ~/Developer/ccp4data/arp-warp-8.0-dereferenced.tgz \
    arp_warp_8.0
```

**Output image:** `${ACR_LOGIN_SERVER}/ccp4i2/base-arpwarp:ccp4-20251105`

### Step 3: Update Server Dockerfile

Edit `Docker/server/Dockerfile.with-ccp4` to use the appropriate base:

```dockerfile
# Without ARP/wARP
FROM ${ACR_LOGIN_SERVER}/ccp4i2/base:ccp4-20251105

# With ARP/wARP
FROM ${ACR_LOGIN_SERVER}/ccp4i2/base-arpwarp:ccp4-20251105
```

## Version Naming Conventions

### CCP4 Versions

| Format | Example | Usage |
|--------|---------|-------|
| Date-based | `ccp4-20251105` | Current standard |
| Major version | `ccp4-9`, `ccp4-10` | Legacy |

The version string must match:
- The directory name inside the tarball
- The `CCP4_VERSION` environment variable
- The Docker image tag suffix

### Image Tags

| Component | Format | Example |
|-----------|--------|---------|
| CCP4 base | `base:ccp4-YYYYMMDD` | `base:ccp4-20251105` |
| ARP/wARP layer | `base-arpwarp:ccp4-YYYYMMDD` | `base-arpwarp:ccp4-20251105` |
| Server | `server:YYYYMMDD-HHMMSS` | `server:20260116-095050` |
| Web | `web:YYYYMMDD-HHMMSS` | `web:20260116-095050` |

## Azure Storage Upload

For runtime-mounted deployments (as opposed to bundled images), CCP4 must be uploaded to Azure Files.

### Upload CCP4 to Azure Files

```bash
# Source environment
. ./Docker/azure/.env.deployment

# Upload (uses azcopy for efficient transfer)
./Docker/azure/scripts/upload-ccp4data.sh \
    ~/Developer/ccp4data/ccp4-20251105 \
    ccp4-20251105
```

**Process:**
1. Creates a compressed tarball locally
2. Uploads to Azure Blob storage via `azcopy`
3. Extracts on the Azure Files share
4. Available at `/mnt/ccp4data/ccp4-20251105` in containers

### Storage Architecture

```
Azure Storage Account (storpub...)
└── Azure Files Share: ccp4data
    ├── ccp4-20251105/          # CCP4 distribution
    │   ├── bin/
    │   │   └── ccp4.setup-sh
    │   ├── lib/
    │   └── py-packages/
    ├── arp_warp_8.0/           # Optional ARP/wARP
    └── ccp4i2-projects/        # User project data
```

## Regional Deployments

CCP4i2 supports multiple Azure regions with independent deployments:

| Region | Directory | ACR |
|--------|-----------|-----|
| North Europe | `Docker/azure/` | `ccp4acrnekmay.azurecr.io` |
| UK South | `Docker/azure-uksouth/` | `ccp4acrukbwmx.azurecr.io` |

Each region maintains:
- Its own Container Registry
- Independent base images
- Region-specific deployment scripts
- Shared Bicep templates (infrastructure as code)

### Cross-Region Image Sharing

Base images can be replicated between regions:

```bash
# Pull from source region
docker pull ccp4acrnekmay.azurecr.io/ccp4i2/base:ccp4-20251105

# Tag for destination region
docker tag ccp4acrnekmay.azurecr.io/ccp4i2/base:ccp4-20251105 \
           ccp4acrukbwmx.azurecr.io/ccp4i2/base:ccp4-20251105

# Push to destination
docker push ccp4acrukbwmx.azurecr.io/ccp4i2/base:ccp4-20251105
```

## Runtime Environment Configuration

### Environment Variables

| Variable | Purpose | Example |
|----------|---------|---------|
| `CCP4_VERSION` | CCP4 directory name | `ccp4-20251105` |
| `CCP4_DATA_PATH` | Mount point for CCP4 data | `/mnt/ccp4data` |
| `CCP4` | Full CCP4 installation path | `/mnt/ccp4data/ccp4-20251105` |

### Startup Detection (startup-server.sh)

The server container detects CCP4 location in priority order:

1. **`CCP4` environment variable** - Explicit path (set by bundled Dockerfile)
2. **`/opt/ccp4`** - Bundled in image at build time
3. **`${CCP4_DATA_PATH}/${CCP4_VERSION}`** - Mounted at runtime
4. **Auto-scan** - Search `CCP4_DATA_PATH` for `ccp4-*` directories

### PYTHONPATH Configuration

The startup script configures Python paths to ensure proper module resolution:

```
1. /usr/src/app                              # Application code (highest priority)
2. /usr/local/lib/python3.11/site-packages   # Container-installed packages
3. /mnt/ccp4data/ccp4-20251105/py-packages   # CCP4 Python packages (lowest)
```

This ordering ensures container-specific packages (e.g., `psycopg2-binary` for PostgreSQL) take precedence over CCP4's bundled versions.

## Complete Deployment Checklist

### New CCP4 Version

```
□ 1. Download/install new CCP4 version locally
□ 2. Create dereferenced tarball
      tar --dereference -czf ccp4-YYYYMMDD-dereferenced.tgz ccp4-YYYYMMDD/
□ 3. Build base image
      ./Docker/azure/scripts/build-base-image.sh <tarball> <version>
□ 4. (Optional) Build ARP/wARP layer
      ./Docker/azure-uksouth/scripts/build-arpwarp-image.sh <tarball> <version>
□ 5. Update Dockerfile.with-ccp4 base image reference
□ 6. Update .env files with new CCP4_VERSION
□ 7. Build application images
      ./Docker/azure/scripts/build-and-push.sh server
      ./Docker/azure/scripts/build-and-push.sh web
□ 8. Deploy to Container Apps
      ./Docker/azure/scripts/deploy-applications.sh
□ 9. Verify deployment
      az containerapp logs show --name ccp4i2-bicep-server ...
```

### New ARP/wARP Version

```
□ 1. Download/install new ARP/wARP locally
□ 2. Create dereferenced tarball
      tar --dereference -czf arp-warp-X.Y-dereferenced.tgz arp_warp_X.Y/
□ 3. Build ARP/wARP layer
      ./Docker/azure-uksouth/scripts/build-arpwarp-image.sh <tarball> arp_warp_X.Y
□ 4. Update Dockerfile.with-ccp4 base image reference
□ 5. Rebuild and deploy server images
```

## Troubleshooting

### Image Build Fails

**Symptom:** `tar: Cannot open: No such file or directory`

**Solution:** Ensure the tarball path is absolute and the file exists:
```bash
ls -la ~/Developer/ccp4data/ccp4-20251105-dereferenced.tgz
```

### CCP4 Not Found at Runtime

**Symptom:** `CCP4 setup script not found`

**Check:**
1. Verify `CCP4_VERSION` environment variable
2. Check mount point: `ls /mnt/ccp4data/`
3. Verify Azure Files share is mounted
4. Check container logs for startup script output

### Import Errors for CCP4 Modules

**Symptom:** `ModuleNotFoundError: No module named 'gemmi'`

**Solution:** Ensure CCP4 setup script is sourced:
```bash
source /mnt/ccp4data/ccp4-20251105/bin/ccp4.setup-sh
```

### Large Push Times

**Symptom:** `docker push` takes hours

**Cause:** First push of 10-15GB base image

**Mitigation:**
- Use wired network connection
- Build during off-peak hours
- Subsequent pushes use layer caching

## Dockerfile Reference

### Dockerfile.ccp4-base

```dockerfile
FROM python:3.11-slim

# System dependencies for CCP4
RUN apt-get update && apt-get install -y \
    libgl1 libglu1-mesa libfreetype6 \
    libxcb-* libxkbcommon0 \
    && rm -rf /var/lib/apt/lists/*

# Copy dereferenced CCP4 distribution
ARG CCP4_VERSION=ccp4-20251105
COPY ${CCP4_VERSION} /opt/ccp4

# Set environment
ENV CCP4=/opt/ccp4
ENV PATH="${CCP4}/bin:${PATH}"
```

### Dockerfile.arpwarp

```dockerfile
ARG BASE_IMAGE
FROM ${BASE_IMAGE}

# Copy ARP/wARP distribution
ARG ARPWARP_VERSION=arp_warp_8.0
COPY ${ARPWARP_VERSION} /mnt/ccp4data/${ARPWARP_VERSION}

# Update PATH
ENV PATH="/mnt/ccp4data/${ARPWARP_VERSION}/bin:${PATH}"
```

## Related Documentation

- [Docker/README.md](README.md) - General Docker setup
- [Docker/azure/README.md](azure/README.md) - Azure deployment guide
- [CLAUDE.md](../CLAUDE.md) - Project-wide development guide
