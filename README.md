# CCP4i2

**CCP4i2** is the graphical user interface and scripting environment for the [CCP4 software suite](https://www.ccp4.ac.uk/), providing tools for macromolecular crystallography.

## Overview

CCP4i2 provides:

- **Task wrappers** for crystallographic programs (Refmac, Phaser, Servalcat, etc.)
- **Pipelines** that chain multiple programs together
- **Data management** for crystallographic projects
- **Report generation** with interactive visualizations

## Architecture Migration (ccp4i2-django branch)

This branch represents a major architectural migration from Qt/PySide2 to a modern web-based architecture:

| Component | Legacy (Qt) | Modern (Django) |
|-----------|-------------|-----------------|
| **Database** | QtSql + SQLite | Django ORM + SQLite |
| **Backend** | Qt event loop | Django + asyncio |
| **API** | Qt signals/slots | REST API (DRF) |
| **Frontend** | Qt Widgets | React/Electron (planned) |
| **Job execution** | Qt subprocess | asyncio subprocess |

### Key Features of the Migration

- **Dual-mode compatibility**: Code works in both Qt and Django environments via the `baselayer` module
- **Preserved plugins**: All wrappers and pipelines remain functional with minimal changes
- **Environment detection**: Automatic detection of Qt vs Django mode
- **Modern testing**: pytest-based test suite with Django integration

## Documentation

### For Users

- [Quick Reference](mddocs/QUICK_REFERENCE.md) - Common operations and examples

### For Developers

- [Development Setup](mddocs/setup/DEVELOPMENT_SETUP.md) - **Start here** - Complete environment setup guide
- [Testing Guide](mddocs/setup/TESTING.md) - Running and writing tests
- [Migration Strategy](mddocs/MIGRATION_STRATEGY.md) - Technical details of the Qt→Django migration
- [Architecture Overview](mddocs/architecture/) - System design documentation

### API Documentation

- [Plugin Registry](mddocs/PLUGIN_REGISTRY_README.md) - Plugin discovery and registration
- [API Reference](mddocs/api/) - REST API endpoints and data models

## Quick Start

### Prerequisites

1. **CCP4 Suite 2024/2025** - Download from [CCP4 Downloads](https://ccp4serv6.rc-harwell.ac.uk/10/)
2. **Python 3.11** (included with CCP4)
3. **Git**

### Setup

```bash
# Clone the repository
git clone https://github.com/ccp4/ccp4i2.git
cd ccp4i2
git checkout ccp4i2-django

# Source CCP4 environment
source /path/to/ccp4-20251105/bin/ccp4.setup-sh

# Install additional Django dependencies into ccp4-python
ccp4-python -m pip install django-cors-headers django-filter djangorestframework whitenoise

# Verify setup
ccp4-python -c "import django; print(f'Django {django.__version__}')"
```

### Running Tests

```bash
# Run a single test
./run_test.sh tests/i2run/test_parrot.py -v

# Run multiple tests
./run_test.sh tests/i2run/ -n 4  # parallel execution
```

Test results are stored in `~/.cache/ccp4i2-tests/`. See cleanup instructions printed after each test run.

## Directory Structure

```
ccp4i2/
├── baselayer/          # Qt-free compatibility layer (PySide2 stubs)
├── core/               # Core Python modules (CCP4Data, CCP4File, etc.)
├── pipelines/          # Multi-program workflows
├── wrappers/           # Single-program task wrappers
├── wrappers2/          # Modern wrapper implementations
├── server/             # Django backend
│   └── ccp4x/
│       ├── api/        # REST API ViewSets
│       ├── db/         # Django ORM models
│       ├── i2run/      # Job runner
│       └── config/     # Django settings
├── client/             # React/Electron frontend (planned)
├── tests/              # Test suite
│   └── i2run/          # Integration tests
├── mddocs/             # Documentation
├── demo_data/          # Sample data for testing
└── report/             # Report templates and generation
```

## Contributing

See [Development Setup](mddocs/setup/DEVELOPMENT_SETUP.md) for environment configuration.

Key guidelines:
- All changes should include tests
- Use `baselayer` imports instead of direct PySide2 imports
- Run the test suite before submitting PRs

## License

CCP4i2 is part of the CCP4 Software Suite. See [CCP4 License](https://www.ccp4.ac.uk/ccp4license.php) for details.

## Links

- [CCP4 Homepage](https://www.ccp4.ac.uk/)
- [CCP4 Downloads](https://ccp4serv6.rc-harwell.ac.uk/10/)
- [CCP4 Documentation](https://www.ccp4.ac.uk/html/INDEX.html)
