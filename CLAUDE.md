# CCP4i2 - Django Branch

CCP4i2 provides an environment for crystallographic computing. This branch (ccp4i2-django) implements a Qt-free architecture using Django for the backend and React/Electron for the frontend.

## Architecture

### Backend
- **Django REST API** (`server/`) - Database interactions via ORM, REST API with Django REST Framework
- **Core Modules** (`core/`) - Business logic, data containers, task management

### Frontend
- **Electron/React App** (`client/`) - Next.js 15, React 19, Moorhen structure viewer

### Compatibility Layer
- **BaseLayer** (`baselayer/`) - Qt-free implementations of PySide2 APIs

### Compounds App (Optional)
- **Compounds Registry & Assays** (`apps/compounds/`) - Compound registration, batch tracking, and assay management
- **Compounds Frontend** (`apps/compounds/frontend/`) - Next.js pages/components for compounds features

## Key Directories

| Directory | Purpose |
|-----------|---------|
| `server/` | Django backend (models, REST API, job execution) |
| `client/` | Electron/React frontend for desktop |
| `apps/` | Optional apps (compounds registry, assays) |
| `apps/compounds/frontend/` | Compounds frontend components (overlaid during Docker build) |
| `baselayer/` | Qt-free compatibility layer (Signal, Slot, QObject stubs) |
| `core/` | Core Python modules and business logic |
| `cli/` | Command-line tools (i2run, utilities) |
| `wrappers/` | Task wrappers for crystallographic programs |
| `wrappers2/` | Additional task wrappers |
| `pipelines/` | Multi-step crystallographic pipelines |
| `pimple/` | Matplotlib-based graph generation |
| `smartie/` | Log parsing utilities |
| `report/` | Report generation and parsing |
| `tests/` | Test suite |
| `Docker/` | Docker configuration for web/cloud deployment |

## Environment Detection

The system automatically detects which backend to use:

1. **Explicit**: Set `CCP4I2_BACKEND=django` or `CCP4I2_BACKEND=qt`
2. **Django context**: Presence of `DJANGO_SETTINGS_MODULE`
3. **Auto-detection**: Try importing PySide2; if unavailable, use Django mode

## CCP4 Environment Setup

**IMPORTANT**: CCP4i2 requires the CCP4 suite with `ccp4-python` to run. Before running any Python commands:

```bash
# Source the CCP4 setup script (adjust path as needed)
source ../ccp4-20251105/bin/ccp4.setup-sh

# This sets up:
# - ccp4-python interpreter with all dependencies (gemmi, clipper, etc.)
# - CCP4 environment variables
# - CCP4 binaries in PATH
```

All Python commands below should use `ccp4-python` instead of `python`.

## Running

### Django Mode
```bash
export CCP4I2_BACKEND=django
export DJANGO_SETTINGS_MODULE=ccp4x.config.settings
ccp4-python -m cli.i2run <task> [options]
```

### Development Server
```bash
cd server
ccp4-python manage.py runserver
```

### Frontend
```bash
cd client
npm install
npm run dev
```

### Tests

Run tests with proper CCP4 environment configuration:
```bash
# Run all tests
./run_test.sh

# Run a specific test file
./run_test.sh tests/i2run/test_substitute_ligand.py

# Run tests matching a pattern
./run_test.sh -k "test_aimless"
```

The `run_test.sh` script sets up the required environment variables and paths.

## Import Pattern

Code in wrappers/pipelines uses the baselayer compatibility module:
```python
from ccp4i2.baselayer import QtCore, Signal, Slot
```

## Plugin Registry

The plugin registry (`core/task_manager/plugin_registry.py`) provides lazy loading of task plugins:

```python
from core import CCP4Modules
TASKMANAGER = CCP4Modules.TASKMANAGER()

# List all available plugins
plugins = TASKMANAGER.list_plugins()  # Returns 176 plugins

# Get a specific plugin class
acorn_class = TASKMANAGER.get_plugin_class('acorn')
```

### Regenerating the Registry

If plugins are added or modified, regenerate the registry:

```bash
source ../ccp4-20251105/bin/ccp4.setup-sh
export CCP4I2_ROOT=$(pwd)
ccp4-python core/task_manager/plugin_lookup.py
```

This scans `wrappers/`, `wrappers2/`, and `pipelines/` directories and generates:
- `plugin_registry.py` - Explicit imports for lazy loading
- `plugin_lookup.json` - Plugin metadata cache

## Files Preserved from Legacy

- `wrappers/`, `wrappers2/`, `pipelines/` - Crystallographic logic
- `pimple/`, `smartie/` - Utility modules
- `qticons/`, `svgicons/` - Task icons
- `tipsOfTheDay/` - User tips
- `docs/` - Documentation

## Docker Build Process

### Web Image Build (Compounds Overlay)

The web Docker image (`ccp4i2/web`) is built using `Docker/client/Dockerfile`. This builds a **unified frontend** by overlaying the compounds app onto the ccp4i2 base:

**Base**: `client/renderer/` (CCP4i2 Electron/Next.js app)

**Overlay process** (performed in Dockerfile):
1. Copy compounds routes into the renderer app directory:
   - `apps/compounds/frontend/app/registry/` → `renderer/app/registry/`
   - `apps/compounds/frontend/app/assays/` → `renderer/app/assays/`
   - `apps/compounds/frontend/app/api/proxy/compounds/` → `renderer/app/api/proxy/compounds/`

2. **Replace the root page** with the app selector:
   - `apps/compounds/frontend/app/app-selector/page.tsx` → `renderer/app/page.tsx`
   - This provides the landing page with links to both CCP4i2 and Compounds features

3. Copy shared components/lib/types:
   - `apps/compounds/frontend/components/compounds/` → `renderer/components/compounds/`
   - `apps/compounds/frontend/lib/compounds/` → `renderer/lib/compounds/`
   - `apps/compounds/frontend/types/compounds/` → `renderer/types/compounds/`

4. Add `@/*` path aliases to tsconfig.json (compounds uses `@/` imports)

**Key file**: `apps/compounds/frontend/app/app-selector/page.tsx` - The front page of the deployed web app

### Local Development vs Docker Build

| Mode | Frontend Source | Root Page |
|------|-----------------|-----------|
| Desktop (Electron) | `client/renderer/` only | CCP4i2 projects list |
| Docker/Cloud | `client/renderer/` + compounds overlay | App selector (`app-selector/page.tsx`) |

**Important**: Changes to the compounds frontend won't appear in the deployed app until you rebuild the web image.

## Azure Deployment

### Configuration

- **Region**: UK South
- **Environment file**: `Docker/azure-uksouth/.env.deployment` - Contains ACR name, resource group, and deployment settings
- **ACR Registry**: `ccp4acrukbwmx.azurecr.io`
- **Resource Group**: `ccp4i2-bicep-rg-uksouth`

### Layered Container Build Architecture

The deployment uses a **3-layer approach** to optimize image sizes and rebuild times:

#### Layer 1: CCP4 Base Image
- **Dockerfile**: `Docker/base/Dockerfile.ccp4-base`
- **Image**: `ccp4i2/base:ccp4-20251105`
- **Purpose**: Contains the large, static CCP4 installation (~10GB)
- **Rebuild**: Only when CCP4 version updates

#### Layer 2: ARP/wARP Layer
- **Dockerfile**: `Docker/base/Dockerfile.arpwarp`
- **Image**: `ccp4i2/base-arpwarp:ccp4-20251105`
- **Purpose**: Adds ARP/wARP tools on top of CCP4 base
- **Rebuild**: Only when ARP/wARP updates

#### Layer 3: Application Images
- **Server Dockerfile**: `Docker/server/Dockerfile.with-ccp4`
- **Web Dockerfile**: `Docker/client/Dockerfile`
- **Images**: `ccp4i2/server:timestamp`, `ccp4i2/web:timestamp`
- **Purpose**: Application code only (~500MB)
- **Rebuild**: On every code change (fast)

### Container Images

| Repository | Description |
|------------|-------------|
| `ccp4i2/base` | CCP4 installation base layer |
| `ccp4i2/base-arpwarp` | CCP4 + ARP/wARP layer |
| `ccp4i2/server` | Django backend (built on base-arpwarp) |
| `ccp4i2/web` | Next.js frontend (CCP4i2 + Compounds overlay) |

### Container Apps

| App Name | Purpose | Scaling |
|----------|---------|---------|
| `ccp4i2-bicep-web` | Next.js frontend | 1-5 replicas (HTTP) |
| `ccp4i2-bicep-server` | Django REST API | 1-10 replicas (CPU/HTTP) |
| `ccp4i2-bicep-worker` | Background job processor | 0-20 replicas (queue depth) |

### Deployment Commands

```bash
# Source environment variables
. ./Docker/azure-uksouth/.env.deployment

# Build base layers (one-time or on CCP4/ARP updates)
./Docker/azure-uksouth/scripts/build-base-image.sh ~/ccp4data/ccp4-20251105-dereferenced.tgz
./Docker/azure-uksouth/scripts/build-arpwarp-image.sh ~/ccp4data/arp-warp-8.0-dereferenced.tgz

# Build and push application images (frequent - on code changes)
./Docker/azure-uksouth/scripts/build-and-push.sh web
./Docker/azure-uksouth/scripts/build-and-push.sh server

# Deploy to container apps
./Docker/azure-uksouth/scripts/deploy-applications.sh web
./Docker/azure-uksouth/scripts/deploy-applications.sh server
./Docker/azure-uksouth/scripts/deploy-applications.sh worker

# Check current deployed image
az containerapp show --name ccp4i2-bicep-web --resource-group "$RESOURCE_GROUP" --query "properties.template.containers[0].image" -o tsv

# Check available image tags
az acr repository show-tags --name "$ACR_NAME" --repository ccp4i2/web --orderby time_desc --top 5
```

### Key Differences from Previous Architecture

| Aspect | Previous (North Europe) | Current (UK South) |
|--------|------------------------|-------------------|
| Region | North Europe | UK South |
| ACR | `ccp4acrnekmay` | `ccp4acrukbwmx` |
| Config directory | `Docker/azure/` | `Docker/azure-uksouth/` |
| CCP4 deployment | Mounted from Azure Files at runtime | Bundled in container images |
| Image build | Single monolithic build | 3-layer approach |
| Server image size | ~15GB (full CCP4) | ~500MB (app only) |
| Rebuild speed | Slow (rebuilds entire CCP4) | Fast (only app layer) |
