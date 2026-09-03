# CCP4i2

**CCP4i2** is the graphical user interface and scripting environment for the [CCP4 software suite](https://www.ccp4.ac.uk/), providing tools for macromolecular crystallography.

## Overview

CCP4i2 provides:

- **Task wrappers** for crystallographic programs (Refmac, Phaser, Servalcat, etc.)
- **Pipelines** that chain multiple programs together
- **Data management** for crystallographic projects
- **Report generation** with interactive visualizations
- **Desktop application** (Electron + React); the Django server also runs standalone for development

## Architecture

CCP4i2 uses a modern web-based architecture:

| Component | Technology |
|-----------|-----------|
| **Database** | Django ORM + SQLite |
| **Backend** | Django REST Framework + Uvicorn (ASGI) |
| **Frontend** | Next.js 15 / React 19 + Electron |
| **Job execution** | asyncio subprocess with crash-safe wrappers |
| **3D Viewer** | Moorhen (web-based Coot) |

## Quick Start

### For Testers (Packaged Electron App)

**See [Give it a try](docs/give-it-a-try.md)** for the full walkthrough. In short:

1. **Install CCP4 10** from [CCP4 Downloads](https://ccp4serv6.rc-harwell.ac.uk/10/downloads/) (provides `ccp4-python`).
2. **Download the desktop app** for your OS from [Releases](https://github.com/ccp4/ccp4i2/releases) — during the alpha, pick the top entry marked **Pre-release** and grab the installer for your OS (macOS `.dmg`, Windows `.exe`, Linux `.deb` **recommended** or `.AppImage`). On macOS, clear quarantine: `xattr -cr "/Applications/ccp4i2-django.app"`. On Ubuntu/Debian, `sudo apt install ./ccp4i2-django_*.deb`.
3. **Launch it**, point it at your CCP4 installation, and click **Install** when prompted — the app installs the matching `ccp4i2` backend into `ccp4-python` for you.

> Do **not** `pip install ccp4i2` against system Python — the backend must live in CCP4's `ccp4-python`, which the app handles automatically. Each app build is pinned to one exact backend version (shown on the launch screen), so you never match versions by hand.

### For Developers

See [Development Setup](mddocs/setup/DEVELOPMENT_SETUP.md) for the full guide. In brief:

```bash
# 1. Install CCP4 (development build)
#    Download from https://ccp4serv6.rc-harwell.ac.uk/10/downloads/packages/

# 2. Clone and enter the repository
git clone https://github.com/ccp4/ccp4i2.git
cd ccp4i2

# 3. Source CCP4 environment (use whichever CCP4 10 build you installed)
source /path/to/ccp4-<build>/bin/ccp4.setup-sh

# 4. Install ccp4i2 into ccp4-python (editable mode)
cd server
ccp4-python -m pip install -e .
cd ..

# 5. Verify
ccp4-python -c "import ccp4i2; print('ccp4i2 OK')"
ccp4-python -c "import django; print(f'Django {django.__version__}')"

# 6. Run tests (from server/, paths are relative to the ccp4i2 package)
cd server
ccp4-python -m pytest ccp4i2/tests/unit/ -v               # fast, no CCP4 binaries
ccp4-python -m pytest ccp4i2/tests/i2run/test_parrot.py -v  # end-to-end task test
cd ..

# 7. Run the Electron client (optional)
cd client
npm install
npm run start:electron
```

### Running Tests

Run from `server/`; test paths are relative to the `ccp4i2` package:

```bash
cd server

# Fast unit tests (no CCP4 binaries)
ccp4-python -m pytest ccp4i2/tests/unit/ -v

# A single end-to-end task test (needs CCP4 environment sourced)
ccp4-python -m pytest ccp4i2/tests/i2run/test_parrot.py -v

# Multiple tests in parallel
ccp4-python -m pytest ccp4i2/tests/i2run/ -n 4
```

(`./run_test.sh` remains as a mac/linux convenience wrapper that sources the
environment for you; the cross-platform way is `ccp4-python -m pytest`.)

Test results are stored in `~/.cache/ccp4i2-tests/`. See cleanup instructions printed after each test run.

## Documentation

**[Documentation Map](docs/README.md)** — the single index to all developer docs.

### Setup & Testing

- [Development Setup](mddocs/setup/DEVELOPMENT_SETUP.md) - **Start here** - Complete environment setup guide
- [Testing Guide](mddocs/setup/TESTING.md) - Running and writing tests

### Authoring a Task

- [Authoring a Task](docs/authoring-a-task.md) - **Start here** - end-to-end guide to adding a wrapper
- [def.xml Reference](docs/def-xml-reference.md) - declaring a task's data model
- [PHIL Task Guide](server/ccp4i2/wrappers/PHIL_TASK_GUIDE.md) - the PHIL alternative for Phenix/PhaserTNG/DIALS-style tools

### Frontend Development

- [Frontend README](client/README.md) - Quick start for the Electron/React frontend
- [Frontend Development Guide](client/FRONTEND_DEVELOPMENT.md) - Comprehensive developer documentation
- [Task Interface Implementation Guide](client/renderer/components/task/task-elements/TASK_INTERFACE_IMPLEMENTATION_GUIDE.md) - building a bespoke task UI

### Architecture & API

- [Quick Reference](mddocs/QUICK_REFERENCE.md) - Common operations and examples
- [Task Registry](server/ccp4i2/core/tasks.py) - tasks are registered in the `TASKS` dict (one entry per task)
- [API Overview](mddocs/api/API_OVERVIEW.md) - REST API endpoints and data models
- [Architecture Overview](mddocs/architecture/) - System design documentation

## Directory Structure

```
ccp4i2/
├── server/                 # Django backend (pip-installable package)
│   ├── ccp4i2/
│   │   ├── api/            # REST API ViewSets
│   │   ├── config/         # Django settings, ASGI entry point
│   │   ├── core/           # Core Python modules (CCP4Data, CCP4File, tasks.py registry)
│   │   ├── db/             # Django ORM models, job execution
│   │   ├── cli/            # Command-line tools (i2, i2run)
│   │   ├── lib/            # Core library (parameters, containers, jobs)
│   │   ├── wrappers/       # Single-program task wrappers
│   │   ├── pipelines/      # Multi-program workflows
│   │   ├── pimple/         # Matplotlib-based graph generation
│   │   ├── smartie/        # Log parsing utilities
│   │   ├── report/         # Report generation
│   │   ├── scripts/        # Shell scripts (run_job_safe.sh)
│   │   ├── demo_data/      # Sample data for testing
│   │   └── tests/          # Test suite (unit/, i2run/, api/, db/, async/)
│   ├── run_test.sh         # End-to-end test runner (run from server/)
│   └── pyproject.toml      # Package definition and dependencies
├── client/                 # Electron + Next.js/React frontend
│   ├── main/               # Electron main process
│   ├── renderer/           # Next.js application
│   └── preload/            # Electron preload scripts
├── packages/               # Shared API contract (npm @ccp4/ccp4i2-api, PyPI ccp4i2-api)
├── docs/                   # Developer & user documentation (see docs/README.md)
├── mddocs/                 # Setup, testing, architecture, API reference
└── deploy/                 # Deployment tooling (test VMs, infra scripts)
```

## Deployment Modes

| Mode | Description |
|------|-------------|
| **Electron (Desktop)** | The first-class route: the packaged app finds your CCP4 installation and launches Django via `ccp4-python` |
| **Development server** | `cd server && ccp4-python manage.py runserver` plus `cd client && npm run dev` |

In both modes, `ccp4i2` is installed as a Python package within `ccp4-python`. The Electron app does not bundle any Python code.

> The Docker/Azure deployment tooling that used to live here moved to the
> `newcastleuniversity/materia` repository with the compounds-app extraction; a
> self-contained lab-wide Docker deployment for CCP4i2 itself is a planned
> future addition (see CLAUDE.md).

## Contributing

See [Development Setup](mddocs/setup/DEVELOPMENT_SETUP.md) for environment configuration.

Key guidelines:
- All changes should include tests
- Run the test suite before submitting PRs

## License

CCP4i2 is part of the CCP4 Software Suite. See [CCP4 License](https://www.ccp4.ac.uk/ccp4license.php) for details.

## Links

- [CCP4 Homepage](https://www.ccp4.ac.uk/)
- [CCP4 Downloads](https://ccp4serv6.rc-harwell.ac.uk/10/)
- [CCP4 Development Builds](https://ccp4serv6.rc-harwell.ac.uk/10/downloads/packages/)
- [CCP4 Documentation](https://www.ccp4.ac.uk/html/INDEX.html)
