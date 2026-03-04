# CCP4i2 Frontend

Modern Electron + Next.js frontend for CCP4i2.

## Prerequisites

- **Node.js 18+** and **npm 9+**
- **CCP4 installation** with `ccp4-python` available
- **ccp4i2 pip-installed** into ccp4-python (`cd server && ccp4-python -m pip install -e .`)

The Electron app does not bundle any Python code. It relies on `ccp4-python` having ccp4i2 installed as a package, and launches the Django backend via `ccp4-python -m uvicorn ccp4i2.config.asgi:application`.

## Quick Start

```bash
# Install dependencies
npm install

# Development (Electron mode)
npm run start:electron

# Development (web mode — requires Django server running separately)
npm run start:web

# Package for distribution
npm run package-mac        # macOS .dmg
npm run package-win        # Windows .exe
npm run package-linux-x64  # Linux .AppImage
```

Packaged apps are written to `client/release/`.

## Architecture Overview

```
client/
├── main/               # Electron main process
│   ├── ccp4i2-django-server.ts  # Django/Uvicorn lifecycle
│   ├── ccp4i2-master.ts         # CCP4 detection, app init
│   └── ccp4i2-ipc.ts            # IPC handlers
├── renderer/           # Next.js application (React)
│   ├── app/           # Next.js App Router pages
│   ├── components/    # React components
│   ├── providers/     # React Context providers
│   ├── types/         # TypeScript type definitions
│   └── api.ts         # API layer with SWR hooks
├── preload/           # Electron preload scripts
└── assets/            # Static assets
```

## Technology Stack

| Technology | Purpose |
|------------|---------|
| **Electron** | Desktop application wrapper |
| **Next.js 15** | React framework with App Router |
| **React 19** | UI library |
| **MUI v6** | Component library |
| **SWR** | Data fetching and caching |
| **Redux Toolkit** | Global state management |
| **Moorhen** | 3D molecular visualization |
| **RDKit** | Chemistry/2D structure rendering |
| **Chart.js** | Data visualization |

## How the Backend is Started

In both development and packaged modes, the Electron main process:
1. Detects the CCP4 installation (scans sibling directories, then standard paths)
2. Finds `ccp4-python` within CCP4
3. Runs `ccp4-python -m django migrate` to set up the database
4. Starts `ccp4-python -m uvicorn ccp4i2.config.asgi:application --workers 2`

No `PYTHONPATH` manipulation is needed — `ccp4-python` has ccp4i2 in its site-packages.

## Development Guide

See [FRONTEND_DEVELOPMENT.md](FRONTEND_DEVELOPMENT.md) for detailed documentation.

## Key Directories

- `renderer/components/task/` - Task interface components
- `renderer/components/report/` - Job report visualization
- `renderer/providers/` - React Context providers
- `main/` - Electron IPC and server management
