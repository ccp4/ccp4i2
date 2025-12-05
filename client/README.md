# CCP4i2 Frontend

Modern Electron + Next.js frontend for CCP4i2.

## Quick Start

```bash
# Install dependencies
npm install

# Development (web mode)
npm run start:web

# Development (Electron mode)
npm run start:electron

# Build for production
npm run build:electron

# Package for distribution
npm run package-mac      # macOS
npm run package-win      # Windows
npm run package-linux-x64   # Linux
```

## Architecture Overview

```
client/
├── main/               # Electron main process
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

## Development Guide

See [FRONTEND_DEVELOPMENT.md](FRONTEND_DEVELOPMENT.md) for detailed documentation.

## Key Directories

- `renderer/components/task/` - Task interface components
- `renderer/components/report/` - Job report visualization
- `renderer/providers/` - React Context providers
- `main/` - Electron IPC and server management
