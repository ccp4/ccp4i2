# Legacy Code

This directory contains legacy Qt-dependent code that has been superseded by the Django-based ccp4i2-django system.

## Contents

### utils/
Legacy utility scripts that depend on the old Qt-based `dbapi` and `qtgui` modules:

- `startup.py` - Qt application startup
- `exportI2Project.py` - Project export (Qt-dependent)
- `ImportAllProjects.py` - Batch project import
- `ExportAllProjects.py` - Batch project export
- `importDir.py` - Directory import

## Migration Notes

These utilities have been replaced by:

- **Project export/import**: `ccp4i2 projects export` and `ccp4i2 projects import` CLI commands
- **Database access**: Django models in `server/ccp4x/db/models.py`
- **Startup**: Django management commands

## Stub Module

The `stubs/dbapi.py` module provides stub classes that will raise `NotImplementedError` if any legacy code attempts to use the old database API, with guidance to use the new Django-based system instead.
