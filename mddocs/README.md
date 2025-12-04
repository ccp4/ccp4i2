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
- [api/](api/) - REST API documentation

## Migration Documentation

Historical documentation from the Qt-to-Django migration:

- [migration/](migration/) - Detailed migration notes and progress tracking
- [API_HARMONIZATION_*.md](.) - API harmonization across CData classes
- [MILESTONE_ASYNC_EXECUTION.md](MILESTONE_ASYNC_EXECUTION.md) - Async execution system milestone

## CLI Documentation

- [cli/](cli/) - Command-line interface documentation

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
