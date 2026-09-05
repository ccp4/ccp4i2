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
from ccp4i2.core.base_object.fundamental_types import CList

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

            existing_cp = self.container.controlParameters

            # Add expert level meta-parameter for controlling visibility
            # and serialization filtering. Matches PHIL convention:
            # 0 = basic user-facing, higher = increasingly expert.
            from ccp4i2.core.base_object.fundamental_types import CInt, CList
            from ccp4i2.core.base_object.base_classes import ValueState
            expert_level = CInt()
            expert_level._skip_validation = True
            expert_level._name = "PHIL_EXPERT_LEVEL"
            expert_level.set_qualifier("guiLabel", "Expert level")
            expert_level.set_qualifier("toolTip",
                "Controls which PHIL parameters are visible and written "
                "to working.phil. 0 = basic, higher = more expert.")
            expert_level.set_qualifier("min", 0)
            expert_level.set_qualifier("max", 10)
            expert_level.value = 0
            if hasattr(expert_level, "_value_states"):
                expert_level._value_states["value"] = ValueState.DEFAULT
            expert_level._skip_validation = False
            setattr(existing_cp, "PHIL_EXPERT_LEVEL", expert_level)

            # Merge children into existing controlParameters.
            # We must detach children from the temporary Phil2CData root first,
            # otherwise _setup_hierarchy_for_value sees they already have a parent
            # and won't re-parent them to the real controlParameters.
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

    # --- Declaring where the tool's PHIL comes from ---
    #
    # Every wrapper wrote the same few lines of import-and-return, in one of
    # three shapes. Declare which shape applies and the base class does it,
    # with one consistent failure when the tool is not installed. Override
    # get_master_phil() directly for anything these do not cover.

    #: "module.path:attribute" — a scope object a module already exposes,
    #: e.g. "xia2.Handlers.Phil:master_phil".
    PHIL_SCOPE = None

    #: "package:relative/path.params" — a PHIL file shipped inside a package,
    #: e.g. "phaser:phenix_interface/__init__.params".
    PHIL_PARAMS_FILE = None

    #: "module.path:ClassName" — a CCTBX Program template whose master_phil is
    #: assembled by CCTBXParser, e.g. "phasertng.programs.picard:Program".
    PHIL_PROGRAM = None

    def get_master_phil(self):
        """Return the tool's master PHIL scope.

        Resolves whichever of PHIL_SCOPE / PHIL_PARAMS_FILE / PHIL_PROGRAM the
        subclass declared. Returns None when the tool is not installed, so the
        task still opens in the interface and can say what it needs, rather
        than failing to load at all.
        """
        declared = {
            "PHIL_SCOPE": self.PHIL_SCOPE,
            "PHIL_PARAMS_FILE": self.PHIL_PARAMS_FILE,
            "PHIL_PROGRAM": self.PHIL_PROGRAM,
        }
        named = {k: v for k, v in declared.items() if v}

        if not named:
            raise NotImplementedError(
                f"{self.__class__.__name__} must declare one of "
                f"{', '.join(declared)}, or override get_master_phil()"
            )
        if len(named) > 1:
            raise ValueError(
                f"{self.__class__.__name__} declares more than one PHIL "
                f"source ({', '.join(sorted(named))}); exactly one is allowed"
            )

        kind, declaration = next(iter(named.items()))
        resolver = {
            "PHIL_SCOPE": self._phil_from_scope,
            "PHIL_PARAMS_FILE": self._phil_from_params_file,
            "PHIL_PROGRAM": self._phil_from_program,
        }[kind]
        try:
            return resolver(declaration)
        except (ImportError, AttributeError, OSError) as err:
            logger.warning(
                "%s: cannot resolve %s=%r (%s); the tool is probably not "
                "installed", self.TASKNAME, kind, declaration, err
            )
            return None

    @staticmethod
    def _phil_from_scope(declaration):
        """"module.path:attribute" -> that attribute."""
        import importlib

        module_path, _, attribute = declaration.partition(":")
        return getattr(importlib.import_module(module_path), attribute)

    @staticmethod
    def _phil_from_params_file(declaration):
        """"package:relative/file.params" -> the parsed contents of that file."""
        import importlib

        from iotbx import phil

        package, _, relative = declaration.partition(":")
        root = Path(importlib.import_module(package).__file__).parent
        return phil.parse((root / relative).read_text())

    @staticmethod
    def _phil_from_program(declaration):
        """"module.path:ClassName" -> the CCTBX program template's master_phil."""
        import importlib

        from iotbx.cli_parser import CCTBXParser

        module_path, _, class_name = declaration.partition(":")
        program = getattr(importlib.import_module(module_path), class_name)
        parser = CCTBXParser(
            program_class=program, logger=None, parse_phil=False
        )
        return parser.master_phil

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
        """The user-set scalar parameters as (phil_dotted_path, value_string).

        Respects PHIL_EXPERT_LEVEL: only parameters whose expertLevel
        qualifier is at or below the selected level are serialized, so
        expert defaults do not leak into working.phil unasked.

        Items of a repeated *scope* cannot be expressed as flat pairs and are
        left out here; extract_phil_lines() renders everything, blocks
        included, and is what build_working_phil() uses.

        Returns:
            list of (phil_dotted_path, value_string) tuples
        """
        return [(path, value) for kind, path, value
                in self._collect_phil_entries(self.container.controlParameters,
                                              self._phil_max_level(), "")
                if kind == "leaf"]

    def extract_phil_lines(self):
        """The user-set parameters rendered as PHIL text lines: `path = value`
        for scalars, one line per item for a repeated definition, and one
        `path { ... }` block per item of a repeated scope."""
        entries = self._collect_phil_entries(self.container.controlParameters,
                                             self._phil_max_level(), "")
        return self._render_phil_entries(entries)

    def _phil_max_level(self):
        try:
            max_level = self.container.controlParameters.PHIL_EXPERT_LEVEL.get()
            if max_level is None:
                max_level = 0
        except (AttributeError, Exception):
            max_level = 0
        return max_level

    def _collect_phil_entries(self, container, max_level, prefix):
        """Walk `container` collecting ("leaf", path, value) and
        ("block", path, [entries]) in dataOrder.

        Paths are relative to `prefix` (a scope path plus a dot, or "") so
        that entries inside a block read as PHIL requires. Only user-set
        values are written, inside blocks as at the top level: libtbx treats
        a repeated-scope instance identical to the master's template as the
        template, not an instance, so writing an item's defaults out could
        not make it count anyway.
        """
        entries = []
        for name in container.dataOrder():
            # Skip the meta-parameter itself -- not a PHIL parameter
            if name == "PHIL_EXPERT_LEVEL":
                continue
            obj = getattr(container, name)
            level = (obj.get_qualifier("expertLevel")
                     if hasattr(obj, "get_qualifier") else None)
            if level is not None and level > max_level:
                continue
            full_path, path = self._phil_path_of(obj, name, prefix)

            if isinstance(obj, CList):
                for item in obj:
                    if isinstance(item, CContainer):
                        # Children carry absolute philPaths, so the prefix
                        # to strip is the absolute scope path however deep
                        # this block is nested
                        inner = self._collect_phil_entries(
                            item, max_level, full_path + ".")
                        entries.append(("block", path, inner))
                    elif not hasattr(item, "isSet"):
                        # CList.append() keeps a plain str/int as it is
                        entries.append(("leaf", path, self._phil_value(item)))
                    elif item.isSet(allowUndefined=False, allowDefault=True):
                        # Membership is the explicit act; a default-valued
                        # item is still an item
                        entries.append(("leaf", path, self._phil_value(item)))
            elif isinstance(obj, CContainer):
                entries.extend(self._collect_phil_entries(
                    obj, max_level, prefix))
            elif obj.isSet(allowUndefined=False, allowDefault=False):
                entries.append(("leaf", path, self._phil_value(obj)))
        return entries

    @staticmethod
    def _phil_path_of(obj, name, prefix):
        """The PHIL path of `obj`: (absolute, relative to `prefix`)."""
        phil_path = (obj.get_qualifier("philPath")
                     if hasattr(obj, "get_qualifier") else None)
        if phil_path is None:
            phil_path = name.replace("__", ".")
        relative = phil_path
        if prefix and phil_path.startswith(prefix):
            relative = phil_path[len(prefix):]
        return phil_path, relative

    @staticmethod
    def _phil_value(obj):
        """The value as PHIL text: whitespace-separated, commas stripped."""
        raw = obj.get() if hasattr(obj, "get") else obj
        val = str(raw).split()
        return " ".join([v[:-1] if v.endswith(",") else v for v in val])

    @classmethod
    def _render_phil_entries(cls, entries, indent=0):
        pad = "  " * indent
        lines = []
        for kind, path, payload in entries:
            if kind == "leaf":
                lines.append(f"{pad}{path} = {payload}")
            else:
                lines.append(f"{pad}{path} {{")
                lines.extend(cls._render_phil_entries(payload, indent + 1))
                lines.append(f"{pad}}}")
        return lines

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

        # Collect user-set PHIL parameters from controlParameters, repeated
        # scopes rendered as blocks
        user_lines = self.extract_phil_lines()

        # Run shims to convert rich CCP4i2 types to PHIL values
        work_dir = str(self.getWorkDirectory())
        for shim in self.get_shim_definitions():
            user_lines.extend(f"{name}={val}" for name, val
                              in shim.convert(self.container, work_dir))

        # Build user PHIL string from all sources
        if user_lines:
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
