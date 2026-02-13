# Frontend Development Guide

Comprehensive guide for developing the CCP4i2 frontend.

## Table of Contents

- [Architecture](#architecture)
- [Directory Structure](#directory-structure)
- [Components](#components)
- [Data Flow](#data-flow)
- [API Layer](#api-layer)
- [State Management](#state-management)
- [Task Interface System](#task-interface-system)
- [Report System](#report-system)
- [Styling](#styling)
- [Adding New Features](#adding-new-features)

---

## Architecture

### Dual-Mode Operation

The frontend runs in two modes:

1. **Electron Mode** (`npm run start:electron`)
   - Desktop application
   - Main process manages Django server lifecycle
   - IPC communication between main and renderer

2. **Web Mode** (`npm run start:web`)
   - Browser-based development
   - Connects to external Django server
   - Faster development iteration

### Process Model (Electron)

```
┌─────────────────────────────────────────────────────────────┐
│                    Electron Main Process                     │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────────┐   │
│  │ Django Server │  │ Next.js      │  │ IPC Handlers     │   │
│  │ Manager       │  │ Server       │  │ (ccp4i2-ipc.ts)  │   │
│  └──────────────┘  └──────────────┘  └──────────────────┘   │
└─────────────────────────────────────────────────────────────┘
                              │
                              │ IPC
                              ▼
┌─────────────────────────────────────────────────────────────┐
│                  Electron Renderer Process                   │
│                    (Next.js Application)                     │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────────┐   │
│  │ React        │  │ SWR Data     │  │ Redux Store      │   │
│  │ Components   │  │ Fetching     │  │ (global state)   │   │
│  └──────────────┘  └──────────────┘  └──────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

---

## Directory Structure

### Main Process (`main/`)

```
main/
├── ccp4i2-master.ts        # Application entry point
├── ccp4i2-django-server.ts # Django server lifecycle management
├── ccp4i2-next-server.ts   # Next.js server management
├── ccp4i2-ipc.ts           # IPC handlers (file dialogs, system info)
├── ccp4i2-menu.ts          # Application menu
├── ccp4i2-create-window.ts # Window creation
├── ccp4i2-session.ts       # Session management
├── ccp4i2-setup-sh.ts      # Unix environment setup
└── ccp4i2-setup-windows.ts # Windows environment setup
```

### Renderer Process (`renderer/`)

```
renderer/
├── app/                    # Next.js App Router
│   ├── layout.tsx         # Root layout with providers
│   ├── page.tsx           # Home page (project list)
│   ├── project/[id]/      # Project pages
│   │   ├── layout.tsx     # Project layout (sidebar + content)
│   │   ├── page.tsx       # Project overview
│   │   ├── job/[jobid]/   # Job view
│   │   ├── jobs/          # Jobs list
│   │   └── files/         # Files list
│   ├── moorhen-page/      # Molecular visualization
│   └── api/proxy/         # API proxy route
│
├── components/             # React components
│   ├── task/              # Task interface system
│   │   ├── task-elements/ # Data type UI components
│   │   └── task-interfaces/ # Task-specific UIs
│   ├── report/            # Job report components
│   └── ...                # Other components
│
├── providers/             # React Context providers
│   ├── task-provider.tsx  # Task interface state
│   ├── coot-provider.tsx  # Moorhen integration
│   └── ...
│
├── types/                 # TypeScript definitions
│   └── models.ts          # Data model interfaces
│
├── api.ts                 # API hooks (SWR)
├── api-fetch.ts           # Low-level fetch utilities
├── utils.ts               # Utility functions
└── store.ts               # Redux store configuration
```

---

## Components

### Component Hierarchy

```
App (layout.tsx)
└── CCP4i2App (providers/ccp4i2-app.tsx)
    └── ProjectLayout (project/[id]/layout.tsx)
        ├── Sidebar (jobs list, files)
        └── Content Area
            └── JobView (components/job-view.tsx)
                ├── JobHeader
                ├── Tabs
                │   ├── TaskContainer (task interface)
                │   ├── Report (CCP4i2ReportXMLView)
                │   ├── Diagnostics
                │   └── ...
                └── JobMenu (context menu)
```

### Task Element Components

Located in `components/task/task-elements/`, these render UI for CData types:

| Component | CData Type | Description |
|-----------|------------|-------------|
| `cstring.tsx` | CString | Text input |
| `cint.tsx` | CInt | Integer input |
| `cfloat.tsx` | CFloat | Float input |
| `cboolean.tsx` | CBoolean | Checkbox/toggle |
| `cdatafile.tsx` | CDataFile | File selector |
| `cpdbdatafile.tsx` | CPdbDataFile | PDB file with selection |
| `clist.tsx` | CList | List of items |
| `ccontainer.tsx` | CContainer | Container rendering |
| `csimple.tsx` | Various | Enum/choice selector |

### Report Components

Located in `components/report/`, these render job output:

| Component | Purpose |
|-----------|---------|
| `CCP4i2ReportXMLView.tsx` | Main report container |
| `CCP4i2ReportElement.tsx` | Generic element renderer |
| `CCP4i2ReportTable.tsx` | Table rendering |
| `CCP4i2ReportGrid.tsx` | Grid layout |
| `CCP4i2ReportVerdict.tsx` | Pass/fail verdict display |
| `ChartLib.tsx` | Chart.js integration |

---

## Data Flow

### API Request Flow

```
Component
    │
    ▼
useApi() / useSWR() hook (api.ts)
    │
    ▼
apiJson() / apiPost() (api-fetch.ts)
    │
    ▼
fetch() to /api/proxy/...
    │
    ▼
Next.js API Route (app/api/proxy/[...proxy]/route.ts)
    │
    ▼
Django Backend (localhost:8001)
```

### State Management

1. **Server State** (SWR)
   - Projects, jobs, files
   - Auto-revalidation
   - Cache management

2. **UI State** (React Context)
   - Task interface state (TaskProvider)
   - Modal dialogs
   - Navigation

3. **Global State** (Redux)
   - User preferences
   - Application configuration

---

## API Layer

### SWR Hooks (`api.ts`)

```typescript
// Fetch a project
const { data: project } = useApi<Project>({
  type: "projects",
  id: projectId,
  endpoint: ""
});

// Fetch jobs for a project
const { data: jobs } = useApi<Job[]>({
  type: "projects",
  id: projectId,
  endpoint: "jobs"
});

// Custom hook for job data
const { job, params_xml, validation } = useJob(jobId);
```

### API Functions (`api-fetch.ts`)

```typescript
// GET request
const data = await apiJson<T>("/endpoint/");

// POST request
const result = await apiPost("/endpoint/", { key: "value" });

// PATCH request
const updated = await apiPatch("/endpoint/", { field: "newValue" });

// DELETE request
await apiDelete("/endpoint/");
```

### Response Format

All API responses follow this structure:

```typescript
interface ApiResponse<T> {
  success: boolean;
  data?: T;
  error?: string;
  status: number;
}
```

---

## State Management

### TaskProvider

Manages task interface state:

```typescript
const {
  container,      // Current container data
  setContainer,   // Update container
  validation,     // Validation errors
  setParameter,   // Set parameter value
} = useTask();
```

### Key Providers

| Provider | Purpose | File |
|----------|---------|------|
| `TaskProvider` | Task interface state | `task-provider.tsx` |
| `CootProvider` | Moorhen 3D viewer | `coot-provider.tsx` |
| `RunCheckProvider` | Job validation/run | `run-check-provider.tsx` |
| `FilePreviewContext` | File preview modal | `file-preview-context.tsx` |
| `JobTabProvider` | Job view tab state | `job-tab-provider.tsx` |

---

## Task Interface System

For the **definitive guide** to building task interfaces — including layout, conditional visibility, reactive behaviour, pitfalls, and worked examples — see:

**[Task Interface Implementation Guide](renderer/components/task/task-elements/TASK_INTERFACE_IMPLEMENTATION_GUIDE.md)**

### How Task Interfaces Work (Summary)

1. **Container Loading**: Job container JSON is fetched from the Django backend via SWR
2. **Element Rendering**: `CCP4i2TaskElement` dispatches to type-specific widgets based on `_class`
3. **Parameter Changes**: User edits trigger `setParameter` API calls with optimistic SWR cache patching
4. **Validation**: Changes trigger re-validation; borders turn red/green accordingly
5. **Custom interfaces** live in `components/task/task-interfaces/<taskname>.tsx` and are registered in `task-container.tsx`

---

## Report System

### Report XML Structure

Reports are XML documents with elements like:

```xml
<ccp4i2_body>
  <fold title="Results">
    <table>...</table>
    <graph>...</graph>
  </fold>
  <verdict>...</verdict>
</ccp4i2_body>
```

### Report Components

The `CCP4i2ReportXMLView` recursively renders report XML:

```typescript
// Maps XML element names to React components
const elementMap = {
  "fold": CCP4i2ReportFold,
  "table": CCP4i2ReportTable,
  "graph": CCP4i2ReportFlotGraphGroup,
  "verdict": CCP4i2ReportVerdict,
  // ...
};
```

### Adding a New Report Element

1. Create component in `components/report/`:

```typescript
// CCP4i2ReportMyElement.tsx
export const CCP4i2ReportMyElement: React.FC<{element: Element}> = ({
  element
}) => {
  const data = element.getAttribute("data");
  return <Box>{/* render data */}</Box>;
};
```

2. Register in `CCP4i2ReportElement.tsx`:

```typescript
case "myelement":
  return <CCP4i2ReportMyElement element={element} />;
```

---

## Styling

### MUI Theme

Theme configuration in `theme/`:

```typescript
// theme/theme-provider.tsx
const theme = createTheme({
  palette: {
    mode: "light", // or "dark"
    primary: { main: "#1976d2" },
  },
});
```

### Component Styling

Use MUI's `sx` prop or styled components:

```typescript
<Box sx={{
  display: "flex",
  gap: 2,
  p: 2,
  backgroundColor: "background.paper",
}}>
  {children}
</Box>
```

### Theme-Aware Scrollbars

```typescript
<Box sx={(theme) => ({
  overflowY: "auto",
  scrollbarColor: `${theme.palette.action.disabled} transparent`,
  scrollbarWidth: "thin",
  "&::-webkit-scrollbar": { width: 8 },
  "&::-webkit-scrollbar-thumb": {
    backgroundColor: theme.palette.action.disabled,
    borderRadius: 4,
  },
})}>
```

---

## Adding New Features

### Adding a New Page

1. Create page in `app/`:

```typescript
// app/my-feature/page.tsx
export default function MyFeaturePage() {
  return <MyFeatureComponent />;
}
```

2. Add to navigation if needed

### Adding a New Provider

1. Create provider:

```typescript
// providers/my-provider.tsx
const MyContext = createContext<MyContextType | null>(null);

export const MyProvider: React.FC<{children: ReactNode}> = ({children}) => {
  const [state, setState] = useState(initialState);
  return (
    <MyContext.Provider value={{ state, setState }}>
      {children}
    </MyContext.Provider>
  );
};

export const useMyContext = () => {
  const ctx = useContext(MyContext);
  if (!ctx) throw new Error("useMyContext must be used within MyProvider");
  return ctx;
};
```

2. Add to provider hierarchy in `app/layout.tsx`

### Adding an IPC Handler (Electron)

1. Add handler in `main/ccp4i2-ipc.ts`:

```typescript
ipcMain.handle("my-action", async (event, args) => {
  // Handle action
  return result;
});
```

2. Call from renderer via preload:

```typescript
const result = await window.electron.invoke("my-action", args);
```

---

## Testing

### Running in Development

```bash
# Web mode (faster iteration)
npm run start:web

# Electron mode (full desktop experience)
npm run start:electron
```

### Environment Variables

- `.env.web` - Web mode configuration
- `.env.electron` - Electron mode configuration
- `.env.local` - Local overrides (not committed)

---

## Troubleshooting

### Common Issues

| Issue | Solution |
|-------|----------|
| API calls fail | Check Django server is running on port 8001 |
| Moorhen doesn't load | Ensure `copy-assets` has run |
| Type errors | Run `npm install` to update types |
| Build fails | Clear `.next` cache: `rm -rf renderer/.next` |

### Debug Logging

```typescript
// Enable in browser console
localStorage.setItem("debug", "ccp4i2:*");
```

---

## See Also

- [API Overview](../mddocs/api/API_OVERVIEW.md) - REST API documentation
- [QUICK_REFERENCE.md](../mddocs/QUICK_REFERENCE.md) - Backend plugin patterns
