# Plugin Registry - Lazy Loading System

## Overview

The plugin registry provides **lazy loading** of CCP4i2 plugins using a pure Python approach with explicit imports. This enables:
- Fast O(1) plugin name lookup
- Plugin classes loaded only when actually needed
- IDE-friendly explicit import statements
- Automatic caching of loaded plugins

## Key Benefits

- **No startup penalty** - Plugins only imported when used
- **O(1) membership testing** - `PLUGIN_NAMES` set for instant lookup
- **Explicit imports** - IDE navigation and traceability
- **Automatic caching** - Plugins loaded once, reused forever
- **~140 plugins discovered** - Full CCP4i2 plugin catalog

## Architecture

### Generated Registry (`plugin_registry.py`)

The registry uses **explicit import statements** rather than dynamic `__import__`:

```python
def _get_plugin_class(plugin_name: str) -> Optional[Type]:
    """Get a plugin class by name using explicit imports."""
    if plugin_name == 'pointless':
        from wrappers.pointless.script.pointless import pointless
        return pointless
    if plugin_name == 'refmac':
        from wrappers.refmac_i2.script.refmac_i2 import refmac_i2
        return refmac_i2
    # ... ~140 more plugins
    return None
```

This approach provides:
- **IDE support** - Click-through navigation works
- **Static analysis** - Tools can trace imports
- **Clear traceability** - Explicit module paths

### Fast Lookup Set

```python
PLUGIN_NAMES: set[str] = {
    'pointless',
    'refmac',
    'aimless',
    # ... all plugin names
}
```

Enables O(1) membership testing before attempting import.

### PluginRegistry Class

```python
class PluginRegistry:
    def __init__(self):
        self._cache: Dict[str, Type] = {}

    def get_plugin_class(self, task_name: str, version: Optional[str] = None):
        # Check cache first
        if task_name in self._cache:
            return self._cache[task_name]

        # Fast rejection for unknown plugins
        if task_name not in PLUGIN_NAMES:
            return None

        # Lazy import via explicit function
        plugin_class = _get_plugin_class(task_name)
        if plugin_class:
            self._cache[task_name] = plugin_class
        return plugin_class
```

## Runtime Usage

```python
from ccp4i2.core.CCP4TaskManager import TASKMANAGER

tm = TASKMANAGER()

# List all plugins (no imports!)
plugins = tm.list_plugins()  # Returns sorted list of ~140 plugin names

# Get metadata without importing (from JSON file)
meta = tm.get_plugin_metadata('pointless')
# Returns: {TASKNAME, TASKVERSION, TASKTITLE, ERROR_CODES, ...}

# Lazy load a plugin (imports only when called)
PointlessClass = tm.get_plugin_class('pointless')

# Second call uses cache (no import)
PointlessClass2 = tm.get_plugin_class('pointless')
assert PointlessClass is PointlessClass2  # Same object
```

## Discovery & Generation

### Running the Generator

```bash
export CCP4I2_ROOT=/path/to/ccp4i2
python3 core/task_manager/plugin_lookup.py
```

### What It Does

1. **Scans directories**: `wrappers/`, `wrappers2/`, `pipelines/`
2. **Imports each module** and finds `CPluginScript` subclasses
3. **Extracts metadata**: `TASKNAME`, `TASKVERSION`, `TASKTITLE`, `ERROR_CODES`, etc.
4. **Generates two files**:
   - `plugin_registry.py` - Python module with explicit imports
   - `plugin_lookup.json` - JSON metadata (for debugging/optional use)

### Generated Files

| File | Purpose | Usage |
|------|---------|-------|
| `plugin_registry.py` | Explicit imports + registry class | **Primary** - used at runtime |
| `plugin_lookup.json` | Plugin metadata in JSON | Optional - for debugging/metadata |

## CTaskManager Integration

`CTaskManager` provides the high-level API via `TASKMANAGER()` singleton:

```python
class CTaskManager:
    @property
    def plugin_registry(self):
        """Lazy load registry on first access."""
        if self._plugin_registry is None:
            from .task_manager.plugin_registry import get_registry
            self._plugin_registry = get_registry()
        return self._plugin_registry

    def get_plugin_class(self, task_name, version=None):
        return self.plugin_registry.get_plugin_class(task_name, version)

    def get_plugin_metadata(self, task_name):
        return self.plugin_registry.get_plugin_metadata(task_name)

    def list_plugins(self):
        return self.plugin_registry.list_plugins()
```

## Example: Building a Plugin Launcher

```python
from ccp4i2.core.CCP4TaskManager import TASKMANAGER

tm = TASKMANAGER()

# Show user all plugins (no imports yet!)
for plugin_name in tm.list_plugins():
    meta = tm.get_plugin_metadata(plugin_name)
    if meta:
        print(f"{plugin_name}: {meta.get('TASKTITLE', 'No title')}")

# User selects "pointless"
# Only NOW do we import it
PointlessClass = tm.get_plugin_class('pointless')

# Instantiate and use
plugin = PointlessClass()
print(plugin.TASKNAME)  # 'pointless'
```

## Regenerating the Registry

When plugins are added or modified:

```bash
export CCP4I2_ROOT=/path/to/ccp4i2
python3 core/task_manager/plugin_lookup.py
```

Output:
```
Building plugin lookup from: /path/to/ccp4i2
Scanning wrappers...
  Found 95 plugins in wrappers
Scanning wrappers2...
  Found 2 plugins in wrappers2
Scanning pipelines...
  Found 43 plugins in pipelines
Finished scanning, found 140 plugins
Writing JSON to: .../plugin_lookup.json
Writing Python registry to: .../plugin_registry.py
```

## Performance

| Operation | Time | Imports |
|-----------|------|---------|
| `list_plugins()` | < 1ms | None |
| `task_name in PLUGIN_NAMES` | O(1) | None |
| `get_plugin_metadata()` | < 1ms | None (reads JSON) |
| First `get_plugin_class()` | ~50ms | 1 plugin |
| Subsequent `get_plugin_class()` | < 1ms | Cached |

**Startup time: 0ms** - No plugins loaded until requested!

## Key Files

| File | Type | Purpose |
|------|------|---------|
| `core/task_manager/plugin_lookup.py` | Generator | Scans codebase, generates registry |
| `core/task_manager/plugin_registry.py` | Generated | Explicit imports + PluginRegistry class |
| `core/task_manager/plugin_lookup.json` | Generated | JSON metadata (optional) |
| `core/CCP4TaskManager.py` | Core | High-level API via TASKMANAGER() |

## Why Explicit Imports?

The previous approach used dynamic imports:
```python
# Old approach - no IDE support
module = __import__(module_name, fromlist=[class_name])
plugin_class = getattr(module, class_name)
```

The new approach uses explicit imports:
```python
# New approach - IDE-friendly
if plugin_name == 'pointless':
    from wrappers.pointless.script.pointless import pointless
    return pointless
```

Benefits:
- **IDE navigation** - Ctrl+click works
- **Static analysis** - Linters can check imports
- **Debugging** - Stack traces show real import paths
- **Refactoring** - Tools can track usage
