# Multi-App Integration Strategy

## Overview

This document describes the strategy for integrating multiple applications (CCP4i2, Compounds, etc.) into a single deployment while maintaining separation of concerns and independent development.

## Architecture Decision

### The Problem

CCP4i2 needs to coexist with other applications (e.g., a Compounds app) under a single deployment on Azure Container Apps. The requirements are:

1. **Single domain deployment** - All apps accessible from one URL (e.g., `https://myapp.azurecontainerapps.io`)
2. **Shared infrastructure** - One Django backend, one Next.js frontend container
3. **Independent routing** - Each app has its own URL namespace
4. **Future extensibility** - Easy to add new apps without restructuring

### Chosen Approach: Route-Based Multi-App (Path B)

We chose a **route-based architecture** over Next.js's `basePath` configuration:

| Approach | Pros | Cons |
|----------|------|------|
| **basePath** | Simple for single app | Cannot serve multiple apps; hard-coded at build time |
| **Route Groups (chosen)** | Multiple apps under one deployment; flexible | Requires reorganizing directory structure |

## URL Structure

### Frontend Routes

```
/                           → Landing page (redirects to /ccp4i2)
/ccp4i2                     → CCP4i2 home/projects list
/ccp4i2/project/[id]        → Project view
/ccp4i2/project/[id]/job/[jobid] → Job view
/ccp4i2/config              → Settings page
/ccp4i2/import-project      → Import project wizard
/ccp4i2/moorhen-page/...    → Moorhen molecular viewer
/ccp4i2/graph-viewer/...    → Graph viewer

# Future apps:
/compounds                  → Compounds app home
/compounds/...              → Compounds app pages
```

### API Routes

```
/api/ccp4i2/projects/       → CCP4i2 projects API
/api/ccp4i2/jobs/           → CCP4i2 jobs API
/api/ccp4i2/files/          → CCP4i2 files API
/api/ccp4i2/...             → All CCP4i2 API endpoints

# Future apps:
/api/compounds/...          → Compounds API endpoints
```

### Static Assets

```
/djangostatic/              → Django static files (admin, etc.)
/files/                     → User-uploaded media files
/_next/                     → Next.js static assets (no prefix needed)
```

## Directory Structure

### Frontend (`client/renderer/`)

```
app/
├── layout.tsx              # Root layout (minimal - just html/body/ThemeProvider)
├── page.tsx                # Root page (redirects to /ccp4i2)
├── globals.css             # Global styles
├── api/
│   └── proxy/[...proxy]/   # API proxy to Django
│       └── route.ts
├── ccp4i2/                 # CCP4i2 app namespace
│   ├── layout.tsx          # CCP4i2 layout (AuthProvider, DeleteDialog, etc.)
│   ├── page.tsx            # CCP4i2 home (projects list)
│   ├── config/
│   ├── import-project/
│   ├── new-project/
│   ├── project/
│   │   └── [id]/
│   │       ├── layout.tsx
│   │       ├── page.tsx
│   │       ├── job/[jobid]/
│   │       ├── files/
│   │       ├── jobs/
│   │       └── network/
│   ├── moorhen-page/
│   └── graph-viewer/
│
├── compounds/              # Future: Compounds app namespace
│   ├── layout.tsx
│   └── ...
│
components/                 # Shared components (outside app/)
providers/                  # Shared providers
api/                        # API client utilities
utils/                      # Shared utilities
types/                      # TypeScript types
```

### Backend (`server/`)

```
ccp4i2/
├── api/
│   ├── urls.py             # Routes under /api/ccp4i2/
│   └── views/
├── config/
│   └── settings.py         # STATIC_URL, MEDIA_URL config
└── ...

# Future:
compounds/
├── api/
│   ├── urls.py             # Routes under /api/compounds/
│   └── views/
└── ...
```

## Implementation Status

### ✅ Completed

1. **Next.js Configuration** (`client/renderer/next.config.ts`)
   - Removed `basePath` and `assetPrefix`
   - Updated rewrites to proxy `/api/ccp4i2/*` to Django

2. **Root Layout** (`client/renderer/app/layout.tsx`)
   - Minimal layout with just `<html>`, `<body>`, `<ThemeProvider>`
   - No app-specific providers here

3. **Root Landing Page** (`client/renderer/app/page.tsx`)
   - Redirects to `/ccp4i2`

4. **CCP4i2 Layout** (`client/renderer/app/ccp4i2/layout.tsx`)
   - Contains `AuthProvider`, `DeleteDialogProvider`, `CCP4i2App`
   - App-specific context and providers

5. **Django URL Routing** (`server/ccp4i2/api/urls.py`)
   - All API routes wrapped under `path("api/ccp4i2/", ...)`

6. **Django Settings** (`server/ccp4i2/config/settings.py`)
   - `STATIC_URL = "/djangostatic/"`
   - `MEDIA_URL = "/files/"`

7. **API Proxy** (`client/renderer/app/api/proxy/[...proxy]/route.ts`)
   - Updated to target `api/ccp4i2/${path}`

8. **Navigation Updates**
   - All `router.push()` calls updated to include `/ccp4i2` prefix
   - All `<Link href>` attributes updated

9. **Docker/Azure Configuration**
   - Health checks updated to `/ccp4i2`
   - Nginx proxy configuration updated

10. **Import Path Fixes**
    - All relative imports corrected for new directory depth
    - Build verified successful

### 🔮 Future Work

1. **Add Compounds App**
   - Create `client/renderer/app/compounds/` directory
   - Create `server/compounds/` Django app
   - Add `/api/compounds/` URL routing
   - Update Next.js rewrites for new API namespace

2. **Shared Components**
   - Identify components useful across apps
   - Move to shared location outside app namespaces

3. **Authentication**
   - Consider shared auth at root layout level
   - Or per-app auth providers as currently implemented

## Key Configuration Files

### `client/renderer/next.config.ts`

```typescript
const nextConfig: NextConfig = {
  output: "standalone",
  // NO basePath - we use route-based multi-app
  async rewrites() {
    return [
      {
        source: "/api/ccp4i2/:path*",
        destination: `${serverUrl}/api/ccp4i2/:path*`,
      },
      // Future: /api/compounds/:path*
    ];
  },
};
```

### `server/ccp4i2/api/urls.py`

```python
urlpatterns = [
    path("api/ccp4i2/", include([
        path("projects/", include(projects_router.urls)),
        path("jobs/", include(jobs_router.urls)),
        # ... other routes
    ])),
    path("djangostatic/", ...),
    path("files/", ...),
]
```

## Adding a New App

To add a new app (e.g., "compounds"):

### 1. Frontend

```bash
# Create app directory
mkdir -p client/renderer/app/compounds

# Create layout with app-specific providers
# Create pages under compounds/
```

### 2. Backend

```bash
# Create Django app
cd server
python manage.py startapp compounds

# Add to INSTALLED_APPS
# Create api/urls.py with /api/compounds/ routes
```

### 3. Configuration

```typescript
// next.config.ts - add rewrite
{
  source: "/api/compounds/:path*",
  destination: `${serverUrl}/api/compounds/:path*`,
}
```

### 4. Update Root Landing Page

Consider updating `/app/page.tsx` to show app selection instead of redirecting to ccp4i2.

## Electron Considerations

The Electron desktop app loads specific routes:

- Main window: `/ccp4i2`
- Config window: `/ccp4i2/config`

These are configured in `client/main/index.ts`. When adding new apps, consider whether they should be accessible in Electron mode.

## Testing

After changes, verify:

1. `npm run build` - No module resolution errors
2. `npm run start` - App loads at `/ccp4i2`
3. API calls work - Check browser network tab
4. Navigation works - All links include `/ccp4i2` prefix
5. Electron mode - Desktop app launches correctly

## References

- [Next.js App Router](https://nextjs.org/docs/app)
- [Next.js Route Groups](https://nextjs.org/docs/app/building-your-application/routing/route-groups)
- [Azure Container Apps](https://docs.microsoft.com/en-us/azure/container-apps/)
