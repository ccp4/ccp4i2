"""The CData classes a def.xml may name, by bare class name.

One registry, two kinds of source:

* the core implementation modules CCP4i2 ships (``CORE_MODULES``), and
* the modules a task declares in its ``Task.dataTypes`` entry in
  ``core/tasks.py``: a task directory that needs a composed type of its own
  registers it there and never touches ``core/``.

Both ``def_xml_handler`` (building a task's container from its def.xml) and
``lib/utils/files/digest`` (instantiating a file class by name) resolve
through :func:`cdata_classes`. The registry is built once per process, on
first use; a module that fails to import is reported and skipped, so one
broken plugin cannot take every task down with it.
"""
import importlib
import logging
from functools import lru_cache
from typing import Dict, Iterable, Type

from ccp4i2.core.base_object.base_classes import CData, CContainer
from ccp4i2.core.base_object.fundamental_types import (
    CInt, CFloat, CBoolean, CString, CList,
)

logger = logging.getLogger(f"ccp4i2:{__name__}")

FUNDAMENTAL: Dict[str, Type[CData]] = {
    "CInt": CInt, "CFloat": CFloat, "CBoolean": CBoolean, "CString": CString,
    "CContainer": CContainer, "CList": CList,
}

# The implementation modules under ccp4i2.core. A type used by one task only
# does not belong here: declare it in that task's ``dataTypes`` instead.
CORE_MODULES = (
    "ccp4i2.core.CCP4Annotation",
    "ccp4i2.core.CCP4ComFilePatchManager",
    "ccp4i2.core.CCP4CootData",
    "ccp4i2.core.CCP4CustomTaskManager",
    "ccp4i2.core.CCP4Data",
    "ccp4i2.core.CCP4File",
    "ccp4i2.core.CCP4ImportedJobManager",
    "ccp4i2.core.CCP4MathsData",
    "ccp4i2.core.CCP4ModelData",
    "ccp4i2.core.CCP4PerformanceData",
    "ccp4i2.core.CCP4Preferences",
    "ccp4i2.core.CCP4RefmacData",
    "ccp4i2.core.CCP4XtalData",
    "ccp4i2.core.CDmDomain",
    "ccp4i2.core.CMoorhenSceneDataFile",
)


def classes_in(module_name: str) -> Dict[str, Type[CData]]:
    """Every CData subclass a module defines or re-exports, by class name."""
    found: Dict[str, Type[CData]] = {}
    try:
        module = importlib.import_module(module_name)
    except Exception as err:  # noqa: BLE001 -- one bad module must not take the rest down
        logger.warning("CData module %s could not be imported: %s: %s",
                       module_name, type(err).__name__, err)
        return found
    for attr_name in dir(module):
        if attr_name.startswith("_"):
            continue
        attr = getattr(module, attr_name)
        if (isinstance(attr, type) and issubclass(attr, CData)
                and attr is not CData and attr is not CContainer):
            found[attr.__name__] = attr
    return found


def task_data_type_modules() -> Iterable[str]:
    """The ``dataTypes`` modules every registered task declares, in order."""
    from ccp4i2.core.tasks import TASKS
    seen = []
    for task in TASKS.values():
        for module_name in getattr(task, "dataTypes", ()) or ():
            if module_name not in seen:
                seen.append(module_name)
    return seen


@lru_cache(maxsize=1)
def cdata_classes() -> Dict[str, Type[CData]]:
    """Bare class name -> class, for everything a def.xml may name."""
    registry: Dict[str, Type[CData]] = dict(FUNDAMENTAL)
    for module_name in CORE_MODULES:
        registry.update(classes_in(module_name))
    for module_name in task_data_type_modules():
        for name, cls in classes_in(module_name).items():
            if name in registry and registry[name] is not cls:
                logger.warning("CData class %s from %s shadows the one already registered "
                               "from %s", name, module_name, registry[name].__module__)
            registry[name] = cls
    return registry


def reset_cache() -> None:
    """For tests that register a task after the registry was first built."""
    cdata_classes.cache_clear()
