# Compounds Frontend

Standalone Next.js frontend for the Compounds Registry and Assays management.

## Prerequisites

- Node.js 18+
- The Django backend running on port 8000

## Quick Start

### 1. Start the Django Backend

In one terminal, start the compounds Django server:

```bash
# From the ccp4i2 root directory
cd /Users/martinnoble/Developer/ccp4i2

# Set up environment
export PYTHONPATH="$PWD:$PWD/apps:$PWD/server"
export DJANGO_SETTINGS_MODULE=compounds.settings

# Run the development server
ccp4-python server/manage.py runserver
```

### 2. Start the Frontend

In another terminal:

```bash
cd /Users/martinnoble/Developer/ccp4i2/apps/compounds/frontend

# Install dependencies (first time only)
npm install

# Start the development server
npm run dev
```

The frontend will be available at **http://localhost:3001**

## Architecture

```
frontend/
├── app/                      # Next.js App Router pages
│   ├── page.tsx              # Home page
│   └── registry/
│       ├── targets/
│       │   ├── page.tsx      # Targets list
│       │   └── [id]/
│       │       └── page.tsx  # Target detail + compounds
│       ├── compounds/
│       │   └── [id]/
│       │       └── page.tsx  # Compound detail + batches
│       └── batches/
│           └── [id]/
│               └── page.tsx  # Batch detail + QC files
├── components/               # Reusable components
│   ├── Breadcrumbs.tsx       # Hierarchical navigation
│   └── DataTable.tsx         # Sortable, searchable table
├── lib/
│   └── api.ts                # SWR-based API hooks
└── types/
    └── models.ts             # TypeScript types
```

## Data Hierarchy

```
Targets
  └── Compounds (filtered by target)
        └── Batches (filtered by compound)
              └── QC Files (filtered by batch)
```

## API Endpoints

The frontend proxies requests to Django via Next.js rewrites:

- `GET /api/compounds/targets/` - List all targets
- `GET /api/compounds/targets/{id}/` - Target detail
- `GET /api/compounds/compounds/?target={id}` - Compounds for a target
- `GET /api/compounds/compounds/{id}/` - Compound detail
- `GET /api/compounds/batches/?compound={id}` - Batches for a compound
- `GET /api/compounds/batches/{id}/` - Batch detail
- `GET /api/compounds/batch-qc-files/?batch={id}` - QC files for a batch

## Development

This is a standalone app that doesn't affect the main CCP4i2 client. You can:

- Run it independently on port 3001
- Make changes without impacting core CCP4i2
- Test API integrations in isolation

## Integration with CCP4i2 Client

Once the frontend is mature, it can be integrated into the main CCP4i2 client by:
1. Copying components to `client/renderer/components/compounds/`
2. Adding routes under `client/renderer/app/compounds/`
3. Adapting to use the shared API hooks and theme
