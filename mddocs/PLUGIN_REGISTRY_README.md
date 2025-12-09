# Plugin Registry - Lazy Loading System

## Overview

The plugin registry provides **lazy loading** of CCP4i2 plugins, allowing you to:
- Browse all available plugins without importing them
- Get metadata instantly (no imports needed)
- Load plugin classes only when actually needed
- Automatic caching of loaded plugins

## Key Benefits

✅ **No startup penalty** - Plugins only imported when used
✅ **Fast metadata access** - Get TASKNAME, VERSION, ERROR_CODES without imports
✅ **Automatic caching** - Plugins loaded once, reused forever
✅ **148 plugins discovered** - Full CCP4i2 plugin catalog

## Architecture

### Discovery Phase (One-time)

```bash
export CCP4I2_ROOT=/path/to/ccp4i2
python3 core/task_manager/plugin_lookup.py
```

Scans CCP4i2 directory and generates:
- `plugin_lookup.json` - JSON format (for compatibility/debugging)
- **`plugin_registry.py`** - Python module with lazy loading

### Runtime Usage

```python
from ccp4i2.core.CCP4TaskManager import TASKMANAGER

tm = TASKMANAGER()

# List all plugins (no imports!)
plugins = tm.list_plugins()  # Returns list of 148 plugin names

# Get metadata without importing
meta = tm.get_plugin_metadata('pointless')
# Returns: {TASKNAME, TASKVERSION, TASKTITLE, ERROR_CODES, ...}

# Lazy load a plugin (imports only when called)
PointlessClass = tm.get_plugin_class('pointless')
# Returns the class, ready to instantiate

# Second call uses cache (no import)
PointlessClass2 = tm.get_plugin_class('pointless')
# PointlessClass is PointlessClass2 == True
```

## How It Works

### 1. Metadata Storage
All plugin metadata is stored as a Python dict in `plugin_registry.py`:
```python
PLUGIN_METADATA = {
    'pointless': {
        'TASKNAME': 'pointless',
        'TASKVERSION': 0.0,
        'TASKTITLE': 'Analyse unmerged dataset (POINTLESS)',
        'ERROR_CODES': {...},
        '_import_module': 'wrappers.pointless.script.pointless',
        '_import_class': 'pointless'
    },
    # ... 147 more plugins
}
```

### 2. Lazy Import on Demand
```python
class PluginRegistry:
    def get_plugin_class(self, task_name: str):
        # Check cache first
        if task_name in self._cache:
            return self._cache[task_name]

        # Get import info
        metadata = PLUGIN_METADATA[task_name]
        module_name = metadata['_import_module']
        class_name = metadata['_import_class']

        # Import only now
        module = __import__(module_name, fromlist=[class_name])
        plugin_class = getattr(module, class_name)

        # Cache it
        self._cache[task_name] = plugin_class
        return plugin_class
```

### 3. CTaskManager Integration

`CTaskManager` provides the high-level API:
```python
class CTaskManager:
    def get_plugin_class(self, task_name, version=None):
        """Get a plugin class with lazy loading."""
        return self.plugin_registry.get_plugin_class(task_name, version)

    def get_plugin_metadata(self, task_name):
        """Get metadata without importing."""
        return self.plugin_registry.get_plugin_metadata(task_name)

    def list_plugins(self):
        """List all available plugins."""
        return self.plugin_registry.list_plugins()
```

## Example: Building a Plugin Launcher

```python
from ccp4i2.core.CCP4TaskManager import TASKMANAGER

tm = TASKMANAGER()

# Show user a list of all plugins
for plugin_name in tm.list_plugins():
    meta = tm.get_plugin_metadata(plugin_name)
    print(f"{plugin_name}: {meta['TASTTITLE']}")

# User selects "pointless"
# Only NOW do we import it
PointlessClass = tm.get_plugin_class('pointless')

# Instantiate and use
plugin = PointlessClass()
plugin.TASKNAME  # 'pointless'
```

## Regenerating the Registry

When plugins are added or modified in CCP4i2:

```bash
export CCP4I2_ROOT=/path/to/ccp4i2
python3 core/task_manager/plugin_lookup.py
```

This regenerates both:
- `core/task_manager/plugin_lookup.json`
- `core/task_manager/plugin_registry.py` ← **Auto-generated, do not edit**

## Testing

```bash
export CCP4I2_ROOT=/path/to/ccp4i2
pytest tests/test_plugin_registry.py -v
```

Tests cover:
- Listing plugins without imports
- Getting metadata without imports
- Lazy loading individual plugins
- Plugin caching
- Loading multiple plugins
- Error handling for non-existent plugins

## Performance

| Operation | Time | Imports |
|-----------|------|---------|
| List 148 plugins | < 1ms | None |
| Get metadata | < 1ms | None |
| First `get_plugin_class()` | ~50ms | 1 plugin |
| Subsequent `get_plugin_class()` | < 1ms | Cached |

**Startup time: 0ms** - No plugins loaded until requested!

## Files Generated

- **`core/task_manager/plugin_registry.py`** - Generated Python module (DO NOT EDIT)
  - ~500KB of pure Python
  - No runtime JSON parsing
  - Type hints included
  - Singleton registry pattern

- `core/task_manager/plugin_lookup.json` - Generated JSON (for compatibility)
  - ~484KB
  - Backward compatible
  - Useful for debugging

## Dependencies

See `requirements-plugin-discovery.txt`:
- lxml, gemmi, networkx, requests (real packages)
- PySide2 (stub in repo - uses our signal system)
