import sys
import pathlib
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

# Add project root to sys.path for stub modules (PySide2, lxml, etc.)
# The project root is two levels up from this script
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(script_dir))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# Add server directory to sys.path for Django imports
server_path = os.path.join(project_root, "server")
if server_path not in sys.path:
    sys.path.insert(0, server_path)

# Add CCP4I2_ROOT to sys.path so plugin imports can work
if CCP4I2_ROOT not in sys.path:
    sys.path.insert(0, CCP4I2_ROOT)

# Configure Django before importing any plugins that might use Django models
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "ccp4x.config.settings")
try:
    import django
    django.setup()
except Exception as e:
    logging.warning(f"Failed to configure Django: {e}")

TASKATTRIBUTES = [
    "COMTEMPLATE",
    "COMTEMPLATEFILE",
    "TASKMODULE",
    "TASKTITLE",
    "TASKNAME",
    "TASKVERSION",
    "WHATNEXT",
    "ASYNCHRONOUS",
    "TIMEOUT_PERIOD",
    "MAXNJOBS",
    "PERFORMANCECLASS",
    "SUBTASKS",
    "RUNEXTERNALPROCESS",
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
            rel_path = os.path.relpath(fpath, root_dir)
            # Module name is relative to ccp4i2 root (no "ccp4i2." prefix)
            # since plugins use imports like "from core.CCP4PluginScript import"
            module_name = rel_path[:-3].replace(os.sep, ".")
            yield dirpath, fname, fpath, module_name


def import_module_from_file(module_name: str, fpath: str):
    """Import a module from a file path, returning the module object or None."""
    try:
        spec = importlib.util.spec_from_file_location(module_name, fpath)
        if spec and spec.loader:
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            return mod
        else:
            logger.warning(f"Could not create import spec for {fpath}")
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
    if root_dir not in sys.path:
        sys.path.insert(0, root_dir)
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

    Creates a module that imports plugins only when requested, avoiding
    the overhead of importing all plugins at startup.
    """
    lines = [
        '"""',
        'Auto-generated plugin registry for CCP4i2 plugins.',
        '',
        'This file is generated by plugin_lookup.py and provides lazy loading',
        'of plugin classes. Plugins are only imported when actually requested.',
        '',
        'DO NOT EDIT THIS FILE MANUALLY - it will be regenerated.',
        '"""',
        '',
        'from typing import Optional, Type, Dict, Any',
        'import sys',
        '',
        '',
        '# Plugin metadata (task attributes without requiring imports)',
        'PLUGIN_METADATA: Dict[str, Dict[str, Any]] = {',
    ]

    # Add metadata for each plugin
    for task_name, plugin_data in sorted(lookup.items()):
        # Extract metadata (everything except module/class which are for import)
        metadata = {k: v for k, v in plugin_data.items()
                   if k not in ('module', 'class', 'module_file')}
        metadata['_import_module'] = plugin_data['module']
        metadata['_import_class'] = plugin_data['class']

        # Sanitize metadata to remove non-JSON-compatible objects
        metadata = sanitize_for_repr(metadata)

        lines.append(f'    {repr(task_name)}: {repr(metadata)},')

    lines.extend([
        '}',
        '',
        '',
        'class PluginRegistry:',
        '    """Registry for lazy-loading plugin classes."""',
        '    ',
        '    def __init__(self):',
        '        self._cache: Dict[str, Type] = {}',
        '    ',
        '    def get_plugin_class(self, task_name: str, version: Optional[str] = None) -> Optional[Type]:',
        '        """',
        '        Get a plugin class by name and optional version.',
        '        ',
        '        The plugin is imported lazily on first access and cached.',
        '        ',
        '        Args:',
        '            task_name: Name of the task/plugin',
        '            version: Optional version (currently ignored - uses latest)',
        '        ',
        '        Returns:',
        '            Plugin class, or None if not found',
        '        """',
        '        # Check cache first',
        '        cache_key = f"{task_name}:{version}" if version else task_name',
        '        if cache_key in self._cache:',
        '            return self._cache[cache_key]',
        '        ',
        '        # Get metadata',
        '        if task_name not in PLUGIN_METADATA:',
        '            return None',
        '        ',
        '        metadata = PLUGIN_METADATA[task_name]',
        '        ',
        '        # TODO: Version filtering if needed',
        '        # For now, just return the registered version',
        '        ',
        '        # Lazy import',
        '        try:',
        '            module_name = metadata["_import_module"]',
        '            class_name = metadata["_import_class"]',
        '            ',
        '            # Import the module',
        '            module = __import__(module_name, fromlist=[class_name])',
        '            ',
        '            # Get the class',
        '            plugin_class = getattr(module, class_name)',
        '            ',
        '            # Cache it',
        '            self._cache[cache_key] = plugin_class',
        '            ',
        '            return plugin_class',
        '        except Exception as e:',
        '            import warnings',
        '            warnings.warn(f"Failed to import plugin {task_name}: {e}")',
        '            return None',
        '    ',
        '    def get_plugin_metadata(self, task_name: str) -> Optional[Dict[str, Any]]:',
        '        """Get plugin metadata without importing the plugin."""',
        '        return PLUGIN_METADATA.get(task_name)',
        '    ',
        '    def list_plugins(self) -> list[str]:',
        '        """Get list of all available plugin names."""',
        '        return list(PLUGIN_METADATA.keys())',
        '    ',
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

                # Fix module paths to be relative to CCP4I2_ROOT, not the plugin dir
                for task_name, plugin_data in plugins.items():
                    if 'module' in plugin_data:
                        # Prepend the plugin directory name
                        plugin_data['module'] = f"{plugin_dir}.{plugin_data['module']}"

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
