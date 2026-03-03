"""
PhilPluginScript - Base class for CCP4i2 wrappers around PHIL-based tools.

This extends CPluginScript to support tools that use libtbx.phil for parameter
definitions (e.g., Phenix, PhaserTNG, DIALS). The key differences from standard
CPluginScript:

1. controlParameters are populated at RUNTIME from the tool's master_phil,
   not from a static .def.xml definition.
2. inputData/outputData still use CCP4i2's rich file types (from .def.xml).
3. At execution time, all parameters (including rich file type shims) are
   assembled into a working_phil via master_phil.fetch().
4. Parameter persistence uses the standard input_params.xml mechanism.

Subclass contract:
    get_master_phil()          -> libtbx.phil scope
    get_phil_exclude_scopes()  -> list[str]
    get_shim_definitions()     -> list[PhilShim]
    get_command_target()       -> str or list

Usage:
    class my_phenix_task(PhilPluginScript):
        TASKNAME = "my_phenix_task"
        TASKCOMMAND = "ccp4-python"

        def get_master_phil(self):
            from my_tool import master_phil
            return master_phil

        def get_shim_definitions(self):
            return [MtzFileShim("HKLIN", "input.hklin")]

        def get_command_target(self):
            return "my_tool.run"
"""

import os
import logging
from typing import Optional
from pathlib import Path

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.base_object.base_classes import CContainer

logger = logging.getLogger(__name__)


class PhilPluginScript(CPluginScript):
    """Base class for CCP4i2 wrappers around PHIL-based tools."""

    # Subclasses may set this as a class-level list of excluded scope paths
    PHIL_EXCLUDE_SCOPES = []

    def __init__(self,
                 parent=None,
                 name=None,
                 xmlFile=None,
                 workDirectory=None,
                 dummy=False,
                 **kwargs):
        """Initialize PhilPluginScript.

        The init order is critical:
        1. super().__init__() loads .def.xml (inputData/outputData) but NOT params
        2. _merge_phil_parameters() populates controlParameters from master_phil
        3. loadDataFromXml() overlays saved user values onto the complete container
        """
        # Pop xmlFile so super().__init__() doesn't load params before
        # PHIL parameters are merged into the container
        super().__init__(
            parent=parent,
            name=name,
            xmlFile=None,  # Defer params loading
            workDirectory=workDirectory,
            dummy=dummy,
            **kwargs
        )

        # Merge PHIL-derived parameters into controlParameters
        if not dummy:
            self._merge_phil_parameters()

        # Now load saved params (if any) onto the complete container
        if xmlFile:
            self.loadDataFromXml(xmlFile)

    def _merge_phil_parameters(self):
        """Import master_phil and merge into controlParameters."""
        try:
            master_phil = self.get_master_phil()
            if master_phil is None:
                logger.warning("No master_phil available for %s", self.TASKNAME)
                return

            from ccp4i2.utils.phil_to_cdata import Phil2CData
            converter = Phil2CData(
                master_phil,
                exclude_scopes=self.get_phil_exclude_scopes(),
            )
            phil_container = converter.convert(root_name="controlParameters")

            # Merge children into existing controlParameters.
            # We must detach children from the temporary Phil2CData root first,
            # otherwise _setup_hierarchy_for_value sees they already have a parent
            # and won't re-parent them to the real controlParameters.
            existing_cp = self.container.controlParameters
            for child in list(phil_container.children()):
                child_name = child.objectName()
                # Don't overwrite existing parameters from .def.xml
                try:
                    getattr(existing_cp, child_name)
                    logger.debug("Skipping PHIL param %s (already in def.xml)", child_name)
                except AttributeError:
                    child.set_parent(None)
                    setattr(existing_cp, child_name, child)

        except ImportError as e:
            logger.warning("Cannot import PHIL module for %s: %s", self.TASKNAME, e)
        except Exception as e:
            logger.error("Error merging PHIL parameters for %s: %s",
                         self.TASKNAME, e, exc_info=True)

    # --- Subclass contract (override these) ---

    def get_master_phil(self):
        """Return the tool's master PHIL scope.

        Override in subclass to import and return the libtbx.phil scope.
        """
        raise NotImplementedError(
            f"{self.__class__.__name__} must implement get_master_phil()"
        )

    def get_phil_exclude_scopes(self):
        """Return list of PHIL scope paths to exclude from the GUI.

        These are scopes whose values are set by shims (e.g., file I/O paths
        that come from CCP4i2's rich file types in inputData).
        """
        return list(self.PHIL_EXCLUDE_SCOPES)

    def get_shim_definitions(self):
        """Return list of PhilShim instances for converting rich types.

        Override in subclass.
        """
        return []

    def get_command_target(self):
        """Return the command-line target for this tool.

        Override in subclass. Can return a string (single argument) or
        a list of strings.
        """
        raise NotImplementedError(
            f"{self.__class__.__name__} must implement get_command_target()"
        )

    # --- PHIL parameter extraction and working_phil assembly ---

    def extract_phil_parameters(self):
        """Walk controlParameters extracting user-set values with their PHIL paths.

        Returns:
            list of (phil_dotted_path, value_string) tuples
        """
        result = []
        self._extract_from_container(self.container.controlParameters, result)
        return result

    def _extract_from_container(self, container, result):
        """Recursively extract set parameters from a container."""
        for name in container.dataOrder():
            obj = getattr(container, name)
            if isinstance(obj, CContainer):
                self._extract_from_container(obj, result)
            elif obj.isSet(allowDefault=False):
                # Only extract parameters the user explicitly changed (not defaults)
                # Use stored philPath qualifier if available, else reverse the __ mapping
                phil_path = obj.get_qualifier("philPath")
                if phil_path is None:
                    phil_path = name.replace("__", ".")

                # Convert value to string, handling comma-separated lists
                # (PHIL uses whitespace-separated values)
                val = str(obj.get()).split()
                val = " ".join([v[:-1] if v.endswith(",") else v for v in val])
                result.append((phil_path, val))

    def build_working_phil(self):
        """Assemble a complete working_phil file using master_phil.fetch().

        This is the key improvement over the old name=value line approach.
        Using fetch() ensures proper PHIL syntax, type validation, and
        correct default propagation.

        Returns:
            str: Path to the written working.phil file
        """
        from libtbx.phil import parse

        master_phil = self.get_master_phil()

        # Collect user-set PHIL parameters from controlParameters
        user_params = self.extract_phil_parameters()

        # Run shims to convert rich CCP4i2 types to PHIL values
        shim_params = []
        work_dir = str(self.getWorkDirectory())
        for shim in self.get_shim_definitions():
            shim_params.extend(shim.convert(self.container, work_dir))

        # Build user PHIL string from all sources
        all_params = user_params + shim_params
        if all_params:
            user_lines = [f"{name}={val}" for name, val in all_params]
            user_phil = parse("\n".join(user_lines))
            working_phil = master_phil.fetch(sources=[user_phil])
        else:
            working_phil = master_phil

        # Write to working.phil in job directory
        phil_path = os.path.join(work_dir, "working.phil")
        from io import StringIO
        buf = StringIO()
        working_phil.show(out=buf)
        with open(phil_path, "w") as f:
            f.write(buf.getvalue())

        return phil_path

    def makeCommandAndScript(self):
        """Build working_phil and construct the command line.

        Subclasses can override for more complex command construction,
        but the default handles the common case.
        """
        phil_path = self.build_working_phil()

        command_target = self.get_command_target()
        if isinstance(command_target, list):
            self.appendCommandLine(command_target)
        else:
            self.appendCommandLine([command_target])

        self.appendCommandLine([phil_path])

        return CPluginScript.SUCCEEDED
