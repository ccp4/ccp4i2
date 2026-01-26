import sys
import os
import importlib.util
import inspect
import logging
import json
from typing import Dict, Any, Type

# Get CCP4I2_ROOT from environment variable
CCP4I2_ROOT = os.environ.get("CCP4I2_ROOT")
if not CCP4I2_ROOT:
    raise ValueError("CCP4I2_ROOT environment variable is not set")

CCP4I2_ROOT = os.path.abspath(CCP4I2_ROOT)

# Configure Django before importing any plugins that might use Django models
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "ccp4i2.config.settings")
try:
    import django
    django.setup()
except Exception as e:
    logging.warning(f"Failed to configure Django: {e}")

TASKATTRIBUTES = [
    "TASKMODULE",
    "TASKTITLE",
    "TASKNAME",
    "TASKVERSION",
    "WHATNEXT",
    "ASYNCHRONOUS",
    "PERFORMANCECLASS",
    "PURGESEARCHLIST",
    "ERROR_CODES",
]


def setup_logger() -> logging.Logger:
    logger = logging.getLogger("build_lookup")
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter("%(levelname)s %(asctime)s %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    logger.setLevel(logging.WARNING)
    return logger


logger = setup_logger()


def is_plugin_script_subclass(obj: Type) -> bool:
    """Return True if obj is a subclass of CPluginScript (by name, robust to import path)."""
    return any(base.__name__.endswith("CPluginScript") for base in obj.__mro__)


def discover_python_files(root_dir: str):
    """Yield (dirpath, filename, fullpath, module_name) for each plugin candidate .py file."""
    for dirpath, dirnames, filenames in os.walk(root_dir):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        for fname in filenames:
            if (
                not fname.endswith(".py")
                or fname.startswith("__")
                or fname == "setup.py"
            ):
                continue
            fpath = os.path.join(dirpath, fname)
            # Compute path relative to CCP4I2_ROOT (not root_dir) for proper module naming
            # e.g., "ccp4i2.wrappers.phaser_analysis.script.phaser_analysis"
            rel_path = os.path.relpath(fpath, CCP4I2_ROOT)
            module_name = "ccp4i2." + rel_path[:-3].replace(os.sep, ".")
            yield dirpath, fname, fpath, module_name


def import_module_from_file(module_name: str, fpath: str):
    """Import a module by name, returning the module object or None.

    Since ccp4i2 is installed as an editable package, we can use standard
    importlib.import_module which properly handles relative imports.
    """
    try:
        # Use standard import_module - this properly handles relative imports
        # because ccp4i2 is installed as an editable package
        mod = importlib.import_module(module_name)
        return mod
    except SystemExit as e:
        # Some ccp4i2 modules call sys.exit() when dependencies are missing
        logger.warning(f"Module {fpath} called sys.exit({e.code})")
    except Exception as e:
        logger.warning(f"Failed to import {fpath}: {e}")
    return None


def extract_plugin_classes(mod, module_name: str) -> Dict[str, Any]:
    """Return a dict of {task_name: plugin_metadata} for all plugin classes in a module."""
    plugins = {}
    # Get the directory of this script for relative path calculation
    script_dir = os.path.dirname(os.path.abspath(__file__))
    for name, obj in inspect.getmembers(mod, inspect.isclass):
        if obj.__name__ == "CPluginScript":
            continue
        if is_plugin_script_subclass(obj):
            entry = {"module": module_name, "class": name}
            # If the class has a __module__ attribute, try to get its file path
            mod_file = getattr(obj, "__module__", None)
            if mod_file:
                try:
                    # Try to get the file path of the module
                    mod_obj = sys.modules.get(mod_file)
                    if mod_obj and hasattr(mod_obj, "__file__"):
                        abs_path = os.path.abspath(mod_obj.__file__)
                        rel_path = os.path.relpath(abs_path, script_dir)
                        entry["module_file"] = rel_path
                except Exception:
                    pass
            for attr in TASKATTRIBUTES:
                if hasattr(obj, attr):
                    entry[attr] = getattr(obj, attr)
            task_name = getattr(obj, "TASKNAME", None)
            if task_name:
                plugins[task_name] = entry
                logger.debug(f"Registered plugin: {task_name} ({module_name}.{name})")
            else:
                logger.warning(
                    f"Class {name} in {module_name} is a CPluginScript subclass but has no TASKNAME"
                )
    return plugins


def build_lookup_from_dir(root_dir: str) -> Dict[str, Any]:
    """
    Crawl a directory tree, import .py files, and build a lookup mapping
    taskName to plugin class metadata.
    Returns a dict: {taskName: {module, class, attributes...}}
    """
    lookup = {}
    for _, fname, fpath, module_name in discover_python_files(root_dir):
        logger.debug(f"Scanning file: {fpath}")
        mod = import_module_from_file(module_name, fpath)
        if mod:
            plugins = extract_plugin_classes(mod, module_name)
            lookup.update(plugins)
    return lookup


def sanitize_for_repr(obj):
    """Convert objects to JSON-compatible representations."""
    if isinstance(obj, dict):
        return {k: sanitize_for_repr(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [sanitize_for_repr(v) for v in obj]
    elif hasattr(obj, 'value'):  # Enum objects
        return obj.value
    elif hasattr(obj, '__dict__') and not isinstance(obj, (str, int, float, bool, type(None))):
        # Custom objects - convert to string
        return str(obj)
    else:
        return obj


def generate_plugin_registry(lookup: Dict[str, Any], output_path: str):
    """
    Generate a Python module with lazy-loading plugin registry.

    Creates a module that imports plugins only when requested, using
    explicit import statements for clarity and IDE support.
    """
    lines = [
        '"""',
        'Auto-generated plugin registry for CCP4i2 plugins.',
        '',
        'This file is generated by plugin_lookup.py and provides lazy loading',
        'of plugin classes using explicit import statements.',
        '',
        'DO NOT EDIT THIS FILE MANUALLY - it will be regenerated.',
        '"""',
        '',
        'from typing import Optional, Type, Dict, Any',
        '',
        '',
        'def _get_plugin_class(plugin_name: str) -> Optional[Type]:',
        '    """',
        '    Get a plugin class by name using explicit imports.',
        '',
        '    This function uses explicit import statements for each plugin,',
        '    providing clear traceability and IDE support.',
        '    """',
    ]

    # Generate explicit import statements for each plugin
    for task_name, plugin_data in sorted(lookup.items()):
        module_name = plugin_data['module']
        class_name = plugin_data['class']
        # Ensure module path is rooted at ccp4i2 package
        if not module_name.startswith('ccp4i2.'):
            module_name = f'ccp4i2.{module_name}'
        lines.extend([
            f'    if plugin_name == {repr(task_name)}:',
            f'        from {module_name} import {class_name}',
            f'        return {class_name}',
        ])

    lines.extend([
        '    return None',
        '',
        '',
        '# Plugin names for fast lookup without loading metadata',
        'PLUGIN_NAMES: set[str] = {',
    ])

    # Add plugin names as a set for fast membership testing
    for task_name in sorted(lookup.keys()):
        lines.append(f'    {repr(task_name)},')

    lines.extend([
        '}',
        '',
        '',
        '# Plugin metadata loaded lazily from JSON',
        '_PLUGIN_METADATA: Optional[Dict[str, Dict[str, Any]]] = None',
        '',
        '',
        'def _load_metadata() -> Dict[str, Dict[str, Any]]:',
        '    """Load plugin metadata from JSON file."""',
        '    global _PLUGIN_METADATA',
        '    if _PLUGIN_METADATA is None:',
        '        import json',
        '        import os',
        '        script_dir = os.path.dirname(os.path.abspath(__file__))',
        '        json_path = os.path.join(script_dir, "plugin_lookup.json")',
        '        try:',
        '            with open(json_path, "r") as f:',
        '                _PLUGIN_METADATA = json.load(f)',
        '        except Exception:',
        '            _PLUGIN_METADATA = {}',
        '    return _PLUGIN_METADATA',
        '',
        '',
        'class PluginRegistry:',
        '    """Registry for lazy-loading plugin classes."""',
        '',
        '    def __init__(self):',
        '        self._cache: Dict[str, Type] = {}',
        '',
        '    def get_plugin_class(self, task_name: str, version: Optional[str] = None) -> Optional[Type]:',
        '        """',
        '        Get a plugin class by name and optional version.',
        '',
        '        The plugin is imported lazily on first access and cached.',
        '',
        '        Args:',
        '            task_name: Name of the task/plugin',
        '            version: Optional version (currently ignored - uses latest)',
        '',
        '        Returns:',
        '            Plugin class, or None if not found',
        '        """',
        '        cache_key = f"{task_name}:{version}" if version else task_name',
        '        if cache_key in self._cache:',
        '            return self._cache[cache_key]',
        '',
        '        if task_name not in PLUGIN_NAMES:',
        '            return None',
        '',
        '        try:',
        '            plugin_class = _get_plugin_class(task_name)',
        '            if plugin_class is not None:',
        '                self._cache[cache_key] = plugin_class',
        '            return plugin_class',
        '        except Exception as e:',
        '            import warnings',
        '            warnings.warn(f"Failed to import plugin {task_name}: {e}")',
        '            return None',
        '',
        '    def get_plugin_metadata(self, task_name: str) -> Optional[Dict[str, Any]]:',
        '        """Get plugin metadata without importing the plugin."""',
        '        metadata = _load_metadata()',
        '        return metadata.get(task_name)',
        '',
        '    def list_plugins(self) -> list[str]:',
        '        """Get list of all available plugin names."""',
        '        return sorted(PLUGIN_NAMES)',
        '',
        '    def clear_cache(self):',
        '        """Clear the plugin cache (for testing/reloading)."""',
        '        self._cache.clear()',
        '',
        '',
        '# Singleton instance',
        '_registry = None',
        '',
        '',
        'def get_registry() -> PluginRegistry:',
        '    """Get the singleton plugin registry instance."""',
        '    global _registry',
        '    if _registry is None:',
        '        _registry = PluginRegistry()',
        '    return _registry',
    ])

    # Write the generated module
    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))
        f.write('\n')


if __name__ == "__main__":
    try:
        root_directory = CCP4I2_ROOT
        print(f"Building plugin lookup from: {root_directory}")
        logger.info(f"Building plugin lookup from: {root_directory}")

        # Only scan plugin directories, not the entire project
        plugin_dirs = ['wrappers', 'wrappers2', 'pipelines']
        result = {}

        for plugin_dir in plugin_dirs:
            dir_path = os.path.join(root_directory, plugin_dir)
            if os.path.exists(dir_path):
                print(f"Scanning {plugin_dir}...")
                plugins = build_lookup_from_dir(dir_path)
                # Module paths are already computed relative to CCP4I2_ROOT
                # so they include the plugin_dir (e.g., "ccp4i2.wrappers.pyphaser_mr...")
                result.update(plugins)
                print(f"  Found {len(plugins)} plugins in {plugin_dir}")

        print(f"Finished scanning, found {len(result)} plugins")

        # Write to script's own directory
        script_dir = os.path.dirname(os.path.abspath(__file__))

        # Write JSON for backward compatibility / debugging
        json_output_path = os.path.join(script_dir, "plugin_lookup.json")
        print(f"Writing JSON to: {json_output_path}")
        with open(json_output_path, "w") as f:
            json.dump(result, f, indent=2)

        # Write Python module with lazy loading
        py_output_path = os.path.join(script_dir, "plugin_registry.py")
        print(f"Writing Python registry to: {py_output_path}")
        generate_plugin_registry(result, py_output_path)

        print(f"Plugin lookup written to: {json_output_path}")
        print(f"Plugin registry written to: {py_output_path}")
        print(f"Found {len(result)} plugins")

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
