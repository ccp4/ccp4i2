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

- **Environment file**: `Docker/azure/.env.deployment` - Contains ACR name, resource group, and other deployment settings
- **ACR Registry**: `ccp4acrnekmay.azurecr.io`

### Container Images

| Repository | Description |
|------------|-------------|
| `ccp4i2/web` | Unified Next.js frontend (CCP4i2 + Compounds overlay) |
| `ccp4i2/server` | Django backend |

### Container Apps

| App Name | Purpose |
|----------|---------|
| `ccp4i2-bicep-web` | Unified frontend (CCP4i2 + Compounds) |
| `ccp4i2-bicep-server` | Unified backend API (CCP4i2 + Compounds) |
| `ccp4i2-bicep-worker` | Background job worker |

### Deployment Commands

```bash
# Source environment variables
. ./Docker/azure/.env.deployment

# Build and push an image (web, server, or worker)
./Docker/azure/scripts/build-and-push.sh web

# Deploy to container app
./Docker/azure/scripts/deploy-applications.sh web

# Check current deployed image
az containerapp show --name ccp4i2-bicep-web --resource-group "$RESOURCE_GROUP" --query "properties.template.containers[0].image" -o tsv

# Check available image tags
az acr repository show-tags --name "$ACR_NAME" --repository ccp4i2/web --orderby time_desc --top 5
```
