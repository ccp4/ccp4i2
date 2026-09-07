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


def _legacy_field(container, name):
    """`name` in any of a classic container's sections (inputData,
    controlParameters, keywords...), or None."""
    for section in container.children():
        if isinstance(section, CContainer) and not isinstance(section, CList):
            try:
                found = getattr(section, name)
            except AttributeError:
                continue
            if found is not None and not callable(found):
                return found
    return None


def _legacy_is_set(obj):
    if isinstance(obj, CList):
        return len(obj) > 0
    return bool(obj.isSet()) if hasattr(obj, "isSet") else obj is not None


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
            # What a shim writes is not offered in the tree as well
            excluded = list(self.get_phil_exclude_scopes())
            for shim in self.get_shim_definitions():
                excluded.extend(t for t in shim.phil_targets() if t not in excluded)
            if self.PHIL_MODE_PATH and self.PHIL_MODE_PATH not in excluded:
                excluded.append(self.PHIL_MODE_PATH)
            converter = Phil2CData(master_phil, exclude_scopes=excluded,
                                   mode=self.PHIL_MODE)
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

    #: A tool that runs in one of several modes, with parameters tagged by
    #: mode in .style (Phaser's `phaser:mode:`, phasertng's `tng:input:`),
    #: can fix the mode per task: only the parameters that apply are
    #: offered, the mode parameter itself leaves the tree, and the working
    #: phil carries `PHIL_MODE_PATH = PHIL_MODE`.
    PHIL_MODE = None
    PHIL_MODE_PATH = None

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

    # --- Handing parameters between PHIL-hosted tasks ---------------------
    #
    # A pipeline that hosts a tool's PHIL hands the tree to the sub-job that
    # runs the tool. Two containers converted from the same master by the
    # same mode have the same shape, so the copy walks them in parallel.

    @staticmethod
    def compose_master_phil(base_master_phil, extra_phil_text):
        """A master phil that is `base_master_phil` with the scopes of
        `extra_phil_text` adopted -- a pipeline's own parameters beside the
        tool's. libtbx's adopt_scope() mutates in place, so the base is
        copied first and the caller's object is left as it was."""
        import copy
        from libtbx.phil import parse
        # customized_copy() is shallow: adopt_scope() on it reaches the base
        composed = copy.deepcopy(base_master_phil)
        composed.adopt_scope(parse(extra_phil_text))
        return composed

    #: A classic task's values that are PHIL parameters here: old field name
    #: -> PHIL path. Adopted from a legacy job on clone.
    LEGACY_PHIL_VALUES = {}
    #: Typed inputs that changed name: this task's name -> the classic name.
    LEGACY_INPUT_RENAMES = {}

    def adopt_legacy_container(self, old):
        """Take the front page of a classic task's container: typed inputs
        of the same name (or a declared rename), and the few values that
        were parameters there and are PHIL here. Everything else is left to
        the PHIL defaults. Returns the names adopted."""
        from ccp4i2.core.base_object.fundamental_types import CInt, CFloat, CString, CBoolean
        adopted = []
        inp = self.container.inputData
        for name in inp.dataOrder():
            src = _legacy_field(old, self.LEGACY_INPUT_RENAMES.get(name, name))
            if src is None or not _legacy_is_set(src):
                continue
            dst = getattr(inp, name)
            enumerators = dst.get_qualifier("enumerators") if hasattr(dst, "get_qualifier") else None
            if isinstance(enumerators, str):
                enumerators = [e.strip() for e in enumerators.split(",")]
            if enumerators and str(src) not in [str(e) for e in enumerators]:
                continue
            try:
                dst.set(src.get() if isinstance(src, (CInt, CFloat, CString, CBoolean)) else src)
            except Exception as err:
                logger.warning("%s: could not adopt %s from the legacy job: %s", self.TASKNAME, name, err)
                continue
            adopted.append(name)
        for old_name, phil_path in self.LEGACY_PHIL_VALUES.items():
            src = _legacy_field(old, old_name)
            if src is None or not _legacy_is_set(src):
                continue
            try:
                self.set_phil(phil_path, src.get() if hasattr(src, "get") else src)
            except Exception as err:
                logger.warning("%s: could not adopt %s as %s: %s", self.TASKNAME, old_name, phil_path, err)
                continue
            adopted.append(f"{old_name} -> {phil_path}")
        return adopted

    def find_phil(self, phil_path, container=None):
        """The CData object in controlParameters whose philPath is
        `phil_path`, or None."""
        container = self.container.controlParameters if container is None else container
        for name in container.dataOrder():
            obj = getattr(container, name)
            if hasattr(obj, "get_qualifier") and obj.get_qualifier("philPath") == phil_path:
                return obj
            if isinstance(obj, CContainer) and not isinstance(obj, CList):
                found = self.find_phil(phil_path, obj)
                if found is not None:
                    return found
        return None

    def set_phil(self, phil_path, value):
        """Set one PHIL parameter by its dotted path: a scalar, or for a
        repeated definition a list of scalars."""
        obj = self.find_phil(phil_path)
        if obj is None:
            raise AttributeError(f"{self.TASKNAME}: no PHIL parameter {phil_path!r}")
        if isinstance(obj, CList):
            obj.set(list(value) if isinstance(value, (list, tuple)) else [value])
        else:
            obj.set(value)
        return obj

    @classmethod
    def copy_phil_tree(cls, source, target):
        """Copy every user-set value from `source` into `target`, two
        containers of the same shape. Defaults are left as the target has
        them; lists are rebuilt item by item."""
        for name in source.dataOrder():
            if name == "PHIL_EXPERT_LEVEL":
                continue
            src = getattr(source, name)
            try:
                dst = getattr(target, name)
            except AttributeError:
                continue
            if isinstance(src, CList):
                if not isinstance(dst, CList):
                    continue
                dst.clear()
                for item in src:
                    if isinstance(item, CContainer):
                        new = dst.makeItem()
                        cls.copy_phil_tree(item, new)
                        dst.append(new)
                    else:
                        dst.append(item.get() if hasattr(item, "get") else item)
            elif isinstance(src, CContainer):
                if isinstance(dst, CContainer):
                    cls.copy_phil_tree(src, dst)
            elif hasattr(src, "isSet") and src.isSet(allowUndefined=False, allowDefault=False):
                dst.set(src.get())

    def hand_phil_to(self, other):
        """Give `other`, a PhilPluginScript over the same master, this
        task's PHIL parameters and expert level."""
        self.copy_phil_tree(self.container.controlParameters, other.container.controlParameters)
        try:
            other.container.controlParameters.PHIL_EXPERT_LEVEL.set(
                self.container.controlParameters.PHIL_EXPERT_LEVEL.get())
        except AttributeError:
            pass

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
                in self._collect_phil_entries(self.container.controlParameters, "")
                if kind == "leaf"]

    def extract_phil_lines(self):
        """The user-set parameters rendered as PHIL text lines: `path = value`
        for scalars, one line per item for a repeated definition, and one
        `path { ... }` block per item of a repeated scope."""
        entries = self._collect_phil_entries(self.container.controlParameters, "")
        return self._render_phil_entries(entries)

    def _collect_phil_entries(self, container, prefix):
        """Walk `container` collecting ("leaf", path, value) and
        ("block", path, [entries]) in dataOrder.

        Paths are relative to `prefix` (a scope path plus a dot, or "") so
        that entries inside a block read as PHIL requires. Only user-set
        values are written, inside blocks as at the top level: libtbx treats
        a repeated-scope instance identical to the master's template as the
        template, not an instance, so writing an item's defaults out could
        not make it count anyway.

        PHIL_EXPERT_LEVEL plays no part here. It is a display choice -- the
        client hides parameters above it -- and a value the user set while
        looking at a deeper level is still theirs after they come back up.
        Filtering on it could only ever drop explicitly set values, since
        defaults are not written anyway.
        """
        entries = []
        for name in container.dataOrder():
            # Skip the meta-parameter itself -- not a PHIL parameter
            if name == "PHIL_EXPERT_LEVEL":
                continue
            obj = getattr(container, name)
            full_path, path = self._phil_path_of(obj, name, prefix)

            if isinstance(obj, CList):
                for item in obj:
                    if isinstance(item, CContainer):
                        # Children carry absolute philPaths, so the prefix
                        # to strip is the absolute scope path however deep
                        # this block is nested
                        inner = self._collect_phil_entries(
                            item, full_path + ".")
                        entries.append(("block", path, inner))
                    elif not hasattr(item, "isSet"):
                        # CList.append() keeps a plain str/int as it is
                        entries.append(("leaf", path, self._phil_value(item)))
                    elif item.isSet(allowUndefined=False, allowDefault=True):
                        # Membership is the explicit act; a default-valued
                        # item is still an item
                        entries.append(("leaf", path, self._phil_value(item)))
            elif isinstance(obj, CContainer):
                entries.extend(self._collect_phil_entries(obj, prefix))
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
    def _entries_from_pairs(cls, pairs):
        """Shim output -- (path, value) pairs, or (path, [pairs]) for one
        instance of a repeated scope -- as collected entries."""
        entries = []
        for path, payload in pairs:
            if isinstance(payload, (list, tuple)) and payload and all(
                    isinstance(e, tuple) for e in payload):
                entries.append(("block", path, cls._entries_from_pairs(payload)))
            elif isinstance(payload, (list, tuple)):
                entries.append(("leaf", path, " ".join(str(v) for v in payload)))
            else:
                entries.append(("leaf", path, str(payload)))
        return entries

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
        # scopes rendered as blocks; a fixed mode comes first
        user_lines = self.extract_phil_lines()
        if self.PHIL_MODE and self.PHIL_MODE_PATH:
            user_lines.insert(0, f"{self.PHIL_MODE_PATH} = {self.PHIL_MODE}")

        # Run shims to convert rich CCP4i2 types to PHIL values; a shim may
        # hand back blocks for a repeated scope as well as pairs
        work_dir = str(self.getWorkDirectory())
        for shim in self.get_shim_definitions():
            user_lines.extend(self._render_phil_entries(
                self._entries_from_pairs(shim.convert(self.container, work_dir))))

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
