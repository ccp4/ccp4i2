# Stub Modules for Plugin Discovery

This directory contains minimal stub implementations of external dependencies
to allow `plugin_lookup.py` to import and introspect CCP4i2 plugin files without
requiring the full dependency stack.

## Purpose

The `core/task_manager/plugin_lookup.py` script dynamically imports Python files
from CCP4i2 to discover plugin classes and extract their metadata (TASKNAME,
TASKVERSION, ERROR_CODES, etc.). However, these files import many external
libraries.

## Approach

We use a hybrid approach:

1. **Real pip packages** for libraries that can be easily installed:
   - lxml, gemmi, networkx, requests, cctbx-base
   - See `requirements-plugin-discovery.txt`

2. **Stub packages** for what we're replacing or can't easily install:
   - PySide2 - We're replacing Qt with our own signal system
   - Other complex/proprietary packages as needed

## Stub Packages

### PySide2/ - Qt GUI Framework Stubs
Provides minimal Qt classes and decorators:
- **QtCore.py**: Slot decorator, Signal, QObject, Qt enums, Property
  - Note: Our `Slot` decorator is compatible with `core/base_object/signal_system.py`

**Why stub instead of install?** We're replacing PySide2/Qt with our own
modern Python signal/slot system (see `core/base_object/signal_system.py`).
The stub allows imports to work during plugin discovery without requiring Qt.

## Results

With real packages + PySide2 stubs, plugin_lookup.py successfully discovers:
- **145 plugins** from CCP4i2
- Full metadata including TASKNAME, TASKVERSION, ERROR_CODES, etc.
- Module and class information for dynamic loading

## Installation

To set up the environment for plugin discovery:

```bash
# Install required packages
pip install -r requirements-plugin-discovery.txt

# Or install individually
pip install lxml gemmi networkx requests cctbx-base
```

Note: The PySide2 stub is already in this repository and doesn't need installation.

## Adding More Stubs

If you encounter import errors for additional libraries, follow this pattern:

1. Create a directory with the package name (e.g., `requests/`)
2. Add `__init__.py` with basic exports
3. Add submodules as needed (e.g., `requests/api.py`)
4. Provide minimal stub classes/functions that return dummy values

The goal is **not** to reimplement the library, just to prevent ImportError and
allow the introspection code to run.

## Example: Adding a new stub

```python
# mystub/__init__.py
"""Stub for mystub library."""

class MyClass:
    """Stub class."""
    def __init__(self, *args, **kwargs):
        pass

    def method(self):
        return None

def my_function(*args, **kwargs):
    """Stub function."""
    pass
```

## Notes

- These stubs are ONLY used during plugin discovery
- They are NOT used at runtime when actually running plugins
- They have no real functionality - just enough to make imports work
