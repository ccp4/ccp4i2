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

## Key Directories

| Directory | Purpose |
|-----------|---------|
| `baselayer/` | Qt-free compatibility layer (Signal, Slot, QObject stubs) |
| `server/` | Django backend (models, REST API, job execution) |
| `client/` | Electron/React frontend |
| `core/` | Core Python modules and business logic |
| `cli/` | Command-line tools (i2run, utilities) |
| `wrappers/` | Task wrappers for crystallographic programs |
| `wrappers2/` | Additional task wrappers |
| `pipelines/` | Multi-step crystallographic pipelines |
| `pimple/` | Matplotlib-based graph generation |
| `smartie/` | Log parsing utilities |
| `report/` | Report generation and parsing |
| `tests/` | Test suite |

## Running

### Django Mode
```bash
export DJANGO_SETTINGS_MODULE=ccp4i2.config.settings
python -m cli.i2run <task> [options]
```

### Development Server
```bash
cd server
python manage.py runserver
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

## Files Preserved from Legacy

- `wrappers/`, `wrappers2/`, `pipelines/` - Crystallographic logic
- `pimple/`, `smartie/` - Utility modules
- `qticons/`, `svgicons/` - Task icons
- `tipsOfTheDay/` - User tips
- `docs/` - Documentation
