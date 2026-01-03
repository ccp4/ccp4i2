# CCP4i2 Documentation

This directory contains documentation for the ccp4i2-django project.

## Getting Started

- **[Development Setup](setup/DEVELOPMENT_SETUP.md)** - Environment setup for development and testing
- **[Testing Guide](setup/TESTING.md)** - Running and writing tests

## Architecture

- **[Migration Strategy](MIGRATION_STRATEGY.md)** - Technical details of the Qt to Django migration
- **[Quick Reference](QUICK_REFERENCE.md)** - Common operations and code patterns
- **[Plugin Registry](PLUGIN_REGISTRY_README.md)** - Plugin discovery and registration system

### Detailed Architecture Docs

- [architecture/](architecture/) - System design documentation

### REST API

- **[API Overview](api/API_OVERVIEW.md)** - REST API endpoints and usage
- [api/](api/) - Detailed API documentation and analysis

## Migration Documentation

Historical documentation from the Qt-to-Django migration:

- [migration/](migration/) - Detailed migration notes and progress tracking
- [API_HARMONIZATION_*.md](.) - API harmonization across CData classes
- [MILESTONE_ASYNC_EXECUTION.md](MILESTONE_ASYNC_EXECUTION.md) - Async execution system milestone

## Pipeline Development

Best practices and patterns for developing CCP4i2 pipelines:

- **[Error Handling Patterns](pipeline/ERROR_HANDLING_PATTERNS.md)** - CErrorReport usage, try/except patterns, ERROR_CODES
- **[Validity Patterns](pipeline/VALIDITY_PATTERNS.md)** - Content-aware validation using validity() overrides

## Frontend Development

- **[Frontend README](../client/README.md)** - Quick start and overview
- **[Frontend Development Guide](../client/FRONTEND_DEVELOPMENT.md)** - Comprehensive developer documentation
  - Architecture, components, data flow
  - Task interface system
  - Report system
  - Adding new features

## CLI Documentation

- **[CLI Reference](cli/CLI.md)** - Full command-line interface documentation
- **[CLI Quick Reference](cli/CLI_QUICK_REFERENCE.md)** - One-page CLI cheat sheet
- **[i2run Guide](cli/I2RUN_GUIDE.md)** - Task runner with complex data object construction

## Historical Notes

These files document specific implementation decisions and completed work:

| File | Description |
|------|-------------|
| `BASE_CLASS_DECISION.md` | Base class hierarchy decisions |
| `CDATAFILE_STUB_ANALYSIS.md` | CDataFile stub generation analysis |
| `CONTENT_FLAG*.md` | Content flag introspection implementation |
| `MAKEHKLIN_*.md` | Makehklin refactoring documentation |
| `MTZ_*.md` | MTZ conversion system documentation |
| `STUB_*.md` | Stub generation patterns |
