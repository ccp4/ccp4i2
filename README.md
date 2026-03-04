# CCP4i2

**CCP4i2** is the graphical user interface and scripting environment for the [CCP4 software suite](https://www.ccp4.ac.uk/), providing tools for macromolecular crystallography.

## Overview

CCP4i2 provides:

- **Task wrappers** for crystallographic programs (Refmac, Phaser, Servalcat, etc.)
- **Pipelines** that chain multiple programs together
- **Data management** for crystallographic projects
- **Report generation** with interactive visualizations
- **Desktop application** (Electron + React) and **web deployment** (Docker/Azure)

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

1. **Install CCP4** with a development build from [CCP4 Downloads](https://ccp4serv6.rc-harwell.ac.uk/10/downloads/packages/)
2. **Install ccp4i2 into ccp4-python**:
   ```bash
   source /path/to/ccp4-20251105/bin/ccp4.setup-sh
   cd server
   ccp4-python -m pip install -e .
   ```
3. **Download and run** the packaged Electron app (macOS `.dmg`, Windows `.exe`, or Linux `.AppImage` from [GitHub Actions](../../actions))

The Electron app automatically detects your CCP4 installation and uses `ccp4-python` (with ccp4i2 pip-installed) to run the Django backend. No Python code is bundled in the Electron app itself.

### For Developers

See [Development Setup](mddocs/setup/DEVELOPMENT_SETUP.md) for the full guide. In brief:

```bash
# 1. Install CCP4 (development build)
#    Download from https://ccp4serv6.rc-harwell.ac.uk/10/downloads/packages/

# 2. Clone and enter the repository
git clone https://github.com/ccp4/ccp4i2.git
cd ccp4i2

# 3. Source CCP4 environment
source /path/to/ccp4-20251105/bin/ccp4.setup-sh

# 4. Install ccp4i2 into ccp4-python (editable mode)
cd server
ccp4-python -m pip install -e .
cd ..

# 5. Verify
ccp4-python -c "import ccp4i2; print('ccp4i2 OK')"
ccp4-python -c "import django; print(f'Django {django.__version__}')"

# 6. Run tests
./run_test.sh tests/i2run/test_parrot.py -v

# 7. Run the Electron client (optional)
cd client
npm install
npm run start:electron
```

### Running Tests

```bash
# Run a single test
./run_test.sh tests/i2run/test_parrot.py -v

# Run multiple tests in parallel
./run_test.sh tests/i2run/ -n 4
```

Test results are stored in `~/.cache/ccp4i2-tests/`. See cleanup instructions printed after each test run.

## Documentation

### Setup & Testing

- [Development Setup](mddocs/setup/DEVELOPMENT_SETUP.md) - **Start here** - Complete environment setup guide
- [Testing Guide](mddocs/setup/TESTING.md) - Running and writing tests

### Frontend Development

- [Frontend README](client/README.md) - Quick start for the Electron/React frontend
- [Frontend Development Guide](client/FRONTEND_DEVELOPMENT.md) - Comprehensive developer documentation

### Architecture & API

- [Quick Reference](mddocs/QUICK_REFERENCE.md) - Common operations and examples
- [Plugin Registry](mddocs/PLUGIN_REGISTRY_README.md) - Plugin discovery and registration
- [API Reference](mddocs/api/) - REST API endpoints and data models
- [Architecture Overview](mddocs/architecture/) - System design documentation

## Directory Structure

```
ccp4i2/
├── server/             # Django backend (pip-installable package)
│   ├── ccp4i2/
│   │   ├── api/        # REST API ViewSets
│   │   ├── config/     # Django settings, ASGI entry point
│   │   ├── db/         # Django ORM models
│   │   ├── i2run/      # Job runner
│   │   ├── lib/        # Core library (parameters, containers, jobs)
│   │   └── scripts/    # Shell scripts (run_job_safe.sh)
│   └── pyproject.toml  # Package definition and dependencies
├── client/             # Electron + Next.js/React frontend
│   ├── main/           # Electron main process
│   ├── renderer/       # Next.js application
│   └── preload/        # Electron preload scripts
├── core/               # Core Python modules (CCP4Data, CCP4File, etc.)
├── wrappers/           # Single-program task wrappers
├── wrappers2/          # Additional wrapper implementations
├── pipelines/          # Multi-program workflows
├── tests/              # Test suite
├── Docker/             # Docker configuration for web/cloud deployment
├── mddocs/             # Documentation
└── demo_data/          # Sample data for testing
```

## Deployment Modes

| Mode | Description |
|------|-------------|
| **Electron (Desktop)** | Packaged app finds CCP4 installation, launches Django via `ccp4-python` |
| **Docker Compose** | Server + web containers for local/dev deployment |
| **Azure** | Container Apps with auto-scaling workers for production |

In all modes, `ccp4i2` is installed as a Python package within `ccp4-python`. The Electron app does not bundle any Python code.

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
