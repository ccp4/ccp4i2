# CCP4i2 Documentation

Reference documentation for the ccp4i2-django project. This folder is for live reference material only; in-flight planning notes, progress reports, and one-time refactor logs are not kept here (they live in commits and PR descriptions).

## Getting Started

- **[Development Setup](setup/DEVELOPMENT_SETUP.md)** — Environment setup for development and testing
- **[Testing Guide](setup/TESTING.md)** — Running and writing tests

## Architecture

- **[Quick Reference](QUICK_REFERENCE.md)** — Async execution infrastructure reference
- **[Plugin Registry](PLUGIN_REGISTRY_README.md)** — Plugin discovery and registration system
- **[Stub Modules for Plugin Discovery](STUBS_README.md)** — PySide2/Qt stub packages used by `plugin_lookup.py`
- **[Stub/Implementation Inheritance Pattern](STUB_IMPLEMENTATION_INHERITANCE_PATTERN.md)** — How full-fat classes inherit from both their generated stub and their full-fat parent (multiple-inheritance MRO concern)
- **[Qt-Based Task GUI Guide](qt_task_gui_guide.md)** — Comprehensive guide to the legacy Qt task GUI system
- [architecture/](architecture/) — System design documentation
- [migration/](migration/) — Migration notes (subdirectory; legacy)

## REST API

- **[API Overview](api/API_OVERVIEW.md)** — REST API endpoints and usage
- [api/](api/) — Detailed API documentation

## Pipeline Development

- **[Error Handling Patterns](pipeline/ERROR_HANDLING_PATTERNS.md)** — CErrorReport usage, try/except, ERROR_CODES
- **[Validity Patterns](pipeline/VALIDITY_PATTERNS.md)** — Content-aware validation via `validity()` overrides

## Frontend Development

- **[Frontend README](../client/README.md)** — Quick start and overview
- **[Frontend Development Guide](../client/FRONTEND_DEVELOPMENT.md)** — Comprehensive developer documentation

## CLI Documentation

- **[CLI Reference](cli/CLI.md)** — Full command-line interface documentation
- **[CLI Quick Reference](cli/CLI_QUICK_REFERENCE.md)** — One-page CLI cheat sheet
- **[i2run Guide](cli/I2RUN_GUIDE.md)** — Task runner with complex data object construction
