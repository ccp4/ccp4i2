"""
Report class discovery script for CCP4i2.

This script discovers report classes (subclasses of Report from CCP4ReportParser)
in the wrappers, wrappers2, and pipelines directories and generates:
- report_lookup.json: JSON metadata for all discovered reports
- report_registry.py: Python module with lazy loading report registry

Usage:
    python report_lookup.py

Environment:
    CCP4I2_ROOT must be set to the project root directory.
"""

import sys
import os
import importlib.util
import inspect
import logging
import json
from typing import Dict, Any, Type, Optional

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

# Report class attributes to extract
REPORT_ATTRIBUTES = [
    "TASKNAME",
    "TASKTITLE",
    "RUNNING",
    "WATCHED_FILE",
    "FAILED",
    "USEPROGRAMXML",
    "SEPARATEDATA",
    "ERROR_CODES",
]


def setup_logger() -> logging.Logger:
    logger = logging.getLogger("build_report_lookup")
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter("%(levelname)s %(asctime)s %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    logger.setLevel(logging.WARNING)
    return logger


logger = setup_logger()


def is_report_subclass(obj: Type) -> bool:
    """Return True if obj is a subclass of Report (by name, robust to import path)."""
    if not inspect.isclass(obj):
        return False
    # Check by name in MRO to handle different import paths
    return any(base.__name__ == "Report" for base in obj.__mro__ if base is not obj)


def discover_report_files(root_dir: str):
    """Yield (dirpath, filename, fullpath, module_name) for each *_report.py file."""
    for dirpath, dirnames, filenames in os.walk(root_dir):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        for fname in filenames:
            # Only look at *_report.py files
            if not fname.endswith("_report.py"):
                continue
            if fname.startswith("__"):
                continue
            fpath = os.path.join(dirpath, fname)
            rel_path = os.path.relpath(fpath, root_dir)
            # Module name is relative to the plugin dir root
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
        logger.warning(f"Module {fpath} called sys.exit({e.code})")
    except Exception as e:
        logger.warning(f"Failed to import {fpath}: {e}")
    return None


def extract_report_classes(mod, module_name: str) -> Dict[str, Any]:
    """Return a dict of {task_name: report_metadata} for all Report subclasses in a module."""
    reports = {}
    script_dir = os.path.dirname(os.path.abspath(__file__))

    for name, obj in inspect.getmembers(mod, inspect.isclass):
        # Skip the base Report class itself
        if name == "Report":
            continue
        if not is_report_subclass(obj):
            continue

        # Skip classes that were imported from other modules (not defined in this file)
        # This prevents picking up GenericReport etc. from wildcard imports
        obj_module = getattr(obj, "__module__", None)
        if obj_module and obj_module != mod.__name__:
            # Check if it's a class from CCP4ReportParser (common wildcard import)
            if "CCP4ReportParser" in obj_module or "report.CCP4ReportParser" in obj_module:
                continue

        entry = {"module": module_name, "class": name}

        # Get module file path for reference
        mod_file = getattr(obj, "__module__", None)
        if mod_file:
            try:
                mod_obj = sys.modules.get(mod_file)
                if mod_obj and hasattr(mod_obj, "__file__"):
                    abs_path = os.path.abspath(mod_obj.__file__)
                    rel_path = os.path.relpath(abs_path, script_dir)
                    entry["module_file"] = rel_path
            except Exception:
                pass

        # Extract report attributes
        for attr in REPORT_ATTRIBUTES:
            if hasattr(obj, attr):
                value = getattr(obj, attr)
                # Handle special cases
                if value is None or isinstance(value, (str, int, float, bool, list, dict)):
                    entry[attr] = value
                else:
                    # Convert non-serializable to string
                    entry[attr] = str(value)

        # Use TASKNAME as key if available
        task_name = getattr(obj, "TASKNAME", None)
        if task_name:
            reports[task_name] = entry
            logger.debug(f"Registered report: {task_name} ({module_name}.{name})")
        else:
            # Fall back to class name (strip _report suffix if present)
            fallback_name = name
            if fallback_name.endswith("_report"):
                fallback_name = fallback_name[:-7]
            reports[fallback_name] = entry
            logger.warning(
                f"Report class {name} in {module_name} has no TASKNAME, using '{fallback_name}'"
            )

    return reports


def build_report_lookup_from_dir(root_dir: str) -> Dict[str, Any]:
    """
    Crawl a directory tree, import *_report.py files, and build a lookup mapping
    taskName to report class metadata.
    """
    lookup = {}
    for _, fname, fpath, module_name in discover_report_files(root_dir):
        logger.debug(f"Scanning report file: {fpath}")
        mod = import_module_from_file(module_name, fpath)
        if mod:
            reports = extract_report_classes(mod, module_name)
            lookup.update(reports)
    return lookup


def sanitize_for_json(obj):
    """Convert objects to JSON-compatible representations."""
    if isinstance(obj, dict):
        return {k: sanitize_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [sanitize_for_json(v) for v in obj]
    elif hasattr(obj, 'value'):  # Enum objects
        return obj.value
    elif hasattr(obj, '__dict__') and not isinstance(obj, (str, int, float, bool, type(None))):
        return str(obj)
    else:
        return obj


def generate_report_registry(lookup: Dict[str, Any], output_path: str):
    """
    Generate a Python module with lazy-loading report registry.

    Creates a module that imports reports only when requested, using
    explicit import statements for clarity and IDE support.
    """
    lines = [
        '"""',
        'Auto-generated report registry for CCP4i2 reports.',
        '',
        'This file is generated by report_lookup.py and provides lazy loading',
        'of report classes using explicit import statements.',
        '',
        'DO NOT EDIT THIS FILE MANUALLY - it will be regenerated.',
        '"""',
        '',
        'from typing import Optional, Type, Dict, Any',
        '',
        '',
        'def _get_report_class(task_name: str) -> Optional[Type]:',
        '    """',
        '    Get a report class by task name using explicit imports.',
        '',
        '    This function uses explicit import statements for each report,',
        '    providing clear traceability and IDE support.',
        '    """',
    ]

    # Generate explicit import statements for each report
    for task_name, report_data in sorted(lookup.items()):
        module_name = report_data['module']
        class_name = report_data['class']
        # Ensure module path is rooted at ccp4i2 package
        if not module_name.startswith('ccp4i2.'):
            module_name = f'ccp4i2.{module_name}'
        lines.extend([
            f'    if task_name == {repr(task_name)}:',
            f'        from {module_name} import {class_name}',
            f'        return {class_name}',
        ])

    lines.extend([
        '    return None',
        '',
        '',
        '# Report names for fast lookup without loading metadata',
        'REPORT_NAMES: set[str] = {',
    ])

    # Add report names as a set for fast membership testing
    for task_name in sorted(lookup.keys()):
        lines.append(f'    {repr(task_name)},')

    lines.extend([
        '}',
        '',
        '',
        '# Report metadata loaded lazily from JSON',
        '_REPORT_METADATA: Optional[Dict[str, Dict[str, Any]]] = None',
        '',
        '',
        'def _load_metadata() -> Dict[str, Dict[str, Any]]:',
        '    """Load report metadata from JSON file."""',
        '    global _REPORT_METADATA',
        '    if _REPORT_METADATA is None:',
        '        import json',
        '        import os',
        '        script_dir = os.path.dirname(os.path.abspath(__file__))',
        '        json_path = os.path.join(script_dir, "report_lookup.json")',
        '        try:',
        '            with open(json_path, "r") as f:',
        '                _REPORT_METADATA = json.load(f)',
        '        except Exception:',
        '            _REPORT_METADATA = {}',
        '    return _REPORT_METADATA',
        '',
        '',
        'class ReportRegistry:',
        '    """Registry for lazy-loading report classes."""',
        '',
        '    def __init__(self):',
        '        self._cache: Dict[str, Type] = {}',
        '',
        '    def get_report_class(self, task_name: str, version: Optional[str] = None) -> Optional[Type]:',
        '        """',
        '        Get a report class by task name.',
        '',
        '        The report is imported lazily on first access and cached.',
        '',
        '        Args:',
        '            task_name: Name of the task (e.g., "refmac", "pointless")',
        '            version: Optional version (currently ignored)',
        '',
        '        Returns:',
        '            Report class, or None if not found',
        '        """',
        '        cache_key = f"{task_name}:{version}" if version else task_name',
        '        if cache_key in self._cache:',
        '            return self._cache[cache_key]',
        '',
        '        if task_name not in REPORT_NAMES:',
        '            return None',
        '',
        '        try:',
        '            report_class = _get_report_class(task_name)',
        '            if report_class is not None:',
        '                self._cache[cache_key] = report_class',
        '            return report_class',
        '        except Exception as e:',
        '            import warnings',
        '            warnings.warn(f"Failed to import report {task_name}: {e}")',
        '            return None',
        '',
        '    def get_report_metadata(self, task_name: str) -> Optional[Dict[str, Any]]:',
        '        """Get report metadata without importing the report class."""',
        '        metadata = _load_metadata()',
        '        return metadata.get(task_name)',
        '',
        '    def list_reports(self) -> list[str]:',
        '        """Get list of all available report task names."""',
        '        return sorted(REPORT_NAMES)',
        '',
        '    def has_report(self, task_name: str) -> bool:',
        '        """Check if a report exists for a task name."""',
        '        return task_name in REPORT_NAMES',
        '',
        '    def clear_cache(self):',
        '        """Clear the report cache (for testing/reloading)."""',
        '        self._cache.clear()',
        '',
        '',
        '# Singleton instance',
        '_registry = None',
        '',
        '',
        'def get_report_registry() -> ReportRegistry:',
        '    """Get the singleton report registry instance."""',
        '    global _registry',
        '    if _registry is None:',
        '        _registry = ReportRegistry()',
        '    return _registry',
    ])

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))
        f.write('\n')


if __name__ == "__main__":
    try:
        root_directory = CCP4I2_ROOT
        print(f"Building report lookup from: {root_directory}")

        # Only scan plugin directories
        plugin_dirs = ['wrappers', 'wrappers2', 'pipelines']
        result = {}

        for plugin_dir in plugin_dirs:
            dir_path = os.path.join(root_directory, plugin_dir)
            if os.path.exists(dir_path):
                print(f"Scanning {plugin_dir}...")
                reports = build_report_lookup_from_dir(dir_path)

                # Fix module paths to be relative to CCP4I2_ROOT
                for task_name, report_data in reports.items():
                    if 'module' in report_data:
                        report_data['module'] = f"{plugin_dir}.{report_data['module']}"

                result.update(reports)
                print(f"  Found {len(reports)} reports in {plugin_dir}")

        print(f"Finished scanning, found {len(result)} reports")

        # Write to script's own directory
        script_dir = os.path.dirname(os.path.abspath(__file__))

        # Write JSON for debugging/reference
        json_output_path = os.path.join(script_dir, "report_lookup.json")
        print(f"Writing JSON to: {json_output_path}")
        with open(json_output_path, "w") as f:
            json.dump(result, f, indent=2)

        # Write Python module with lazy loading
        py_output_path = os.path.join(script_dir, "report_registry.py")
        print(f"Writing Python registry to: {py_output_path}")
        generate_report_registry(result, py_output_path)

        print(f"Report lookup written to: {json_output_path}")
        print(f"Report registry written to: {py_output_path}")
        print(f"Found {len(result)} reports")

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
