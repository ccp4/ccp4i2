"""
Phil2CData - Convert libtbx.phil scopes directly to CData object hierarchies.

This module provides runtime conversion of PHIL parameter definitions into
CCP4i2's CData type system, enabling native PHIL support without intermediate
XML generation.

Usage:
    from libtbx.phil import parse
    from ccp4i2.utils.phil_to_cdata import Phil2CData

    master_phil = parse("...")
    converter = Phil2CData(master_phil, exclude_scopes=["output"])
    control_params = converter.convert(root_name="controlParameters")
"""

import re
import logging
from typing import Optional

from ccp4i2.core.base_object.base_classes import CData, CContainer, ValueState
from ccp4i2.core.base_object.fundamental_types import CInt, CFloat, CBoolean, CString, CList

logger = logging.getLogger(__name__)


def parse_phil_style(style):
    """Read the libtbx GUI conventions out of a `.style` string.

    `.style` is free text the Phenix GUI interprets; the tokens that say
    something a CCP4i2 GUI can act on are translated to qualifiers and the
    rest (bold, box, auto_align, noauto, OnChange:, renderer:) are left in
    the raw `style` qualifier. Mode tags come in two spellings: Phaser's
    `phaser:mode:MR_AUTO,EP*` and phasertng's `tng:input:+frf+ftf`.
    """
    result = {"min": None, "max": None, "hidden": False, "multiLine": False,
              "inputFile": False, "fileType": None, "directory": False,
              "modes": [], "ignored": False}
    for tok in (style or "").split():
        if tok == "hidden":
            result["hidden"] = True
        elif tok == "input_file":
            result["inputFile"] = True
        elif tok == "directory":
            result["directory"] = True
        elif tok == "phaser:ignore":
            result["ignored"] = True
        elif tok.startswith("file_type:"):
            result["fileType"] = tok.partition(":")[2]
        elif tok.startswith("height:"):
            result["multiLine"] = True
        elif tok.startswith(("min=", "max=")):
            key, _, number = tok.partition("=")
            try:
                result[key] = float(number) if "." in number else int(number)
            except ValueError:
                pass
        elif tok.startswith("phaser:mode:"):
            # A layout word can ride along after a comma ("MR_OCC,box"):
            # a mode is upper case, underscores and a trailing wildcard
            result["modes"].extend(
                m for m in tok[len("phaser:mode:"):].split(",")
                if re.fullmatch(r"[A-Z][A-Z0-9_]*\*?|\*", m))
        elif tok.startswith("tng:input:"):
            result["modes"].extend(
                m for m in tok[len("tng:input:"):].split("+") if m)
    return result


def match_modes(current_mode, keyword_modes):
    """Phaser's rule for whether a parameter tagged `keyword_modes` applies
    in `current_mode`: an exact name, `*` for every mode, or a trailing
    wildcard (`MR*`). An untagged parameter (None or empty) is decided by
    its nearest tagged ancestor, which the converter passes down; at the
    top it applies everywhere.
    """
    if not keyword_modes:
        return True
    for mode in keyword_modes:
        if mode == current_mode or mode == "*":
            return True
        if mode.endswith("*") and current_mode.startswith(mode[:-1]):
            return True
    return False


class Phil2CData:
    """Convert a libtbx.phil scope directly to a CData object hierarchy.

    Type mapping (PHIL -> CData):
        str        -> CString
        int        -> CInt
        float      -> CFloat
        bool       -> CBoolean
        choice     -> CString (with enumerators/onlyEnumerators qualifiers)
        path       -> CString (rich file types live in inputData, not here)
        ternary    -> CString (with enumerators for True/False/original)
        floats     -> CString (space-separated float values)
        ints       -> CString (space-separated int values)
        strings    -> CString (space-separated string values)
        choices    -> CString (multi-choice, space-separated)
        filesystem -> CString (file path with extension constraints)
        mtzcol     -> CString (MTZ column name)
        scatterer  -> CString (scatterer element symbol)
        scatterers -> CString (space-separated scatterer symbols)
        unit_cell  -> CString (space-separated unit cell parameters)
        uuid       -> CString (UUID string)

    PHIL attributes are mapped to CData qualifiers:
        short_caption  -> guiLabel
        help           -> toolTip
        expert_level   -> expertLevel
        style          -> style
        caption        -> caption (used for menuText on choices)
        multiple       -> multiple
        value_min      -> min
        value_max      -> max

    Within `.style`, the libtbx GUI conventions are read too (see
    parse_phil_style): `spinner min= max=` -> min/max where the type gave
    none; `hidden` -> hidden; `height:` -> guiMode multiLine; `input_file`
    and `file_type:` -> philInputFile/philFileType (tagged for shims, not
    turned into file objects); `directory` -> isDirectory; mode tags ->
    philModes; `phaser:ignore` -> philIgnored.

    A `.multiple = True` scope becomes a CList whose items are containers
    shaped like the scope (one generated CContainer subclass per scope, so
    makeItem(), params.xml import and the JSON _subItem all work unchanged);
    a `.multiple = True` definition becomes a CList of the leaf type. Both
    start empty: libtbx's fetch() yields no instance of a multiple scope
    unless the working phil supplies one, so an empty list is what the tool
    itself sees, and the master's single instance is a template of defaults
    rather than data. Each new item carries those defaults.
    """

    # Standard PHIL types
    PHIL_TYPE_MAP = {
        "str": CString,
        "int": CInt,
        "float": CFloat,
        "bool": CBoolean,
        "choice": CString,
        "path": CString,
        "ternary": CString,
    }

    # Custom/extended PHIL types (from phasertng, cctbx, etc.)
    # All map to CString since they represent values as text
    PHIL_CUSTOM_TYPE_MAP = {
        "floats": CString,       # e.g., "0.6 1.2"
        "ints": CString,         # e.g., "1 3"
        "strings": CString,      # e.g., "P1 P21"
        "choices": CString,      # multi-choice variant
        "filesystem": CString,   # file path with extension info
        "mtzcol": CString,       # MTZ column label
        "scatterer": CString,    # element symbol e.g. "S"
        "scatterers": CString,   # space-separated elements
        "unit_cell": CString,    # "a b c alpha beta gamma"
        "uuid": CString,         # UUID string
    }

    def __init__(self, phil_scope, exclude_scopes=None, max_expert_level=None,
                 mode=None):
        """
        Args:
            phil_scope: A libtbx.phil scope object (the master_phil).
            exclude_scopes: List of dotted PHIL paths to skip entirely
                (e.g., I/O scopes handled by shims).
            max_expert_level: If set, skip definitions and scopes with
                expert_level > this value. PHIL convention:
                0 = user-facing, 1+ = increasingly expert/internal.
                None means include all levels (no filtering).
            mode: If set, keep only the definitions and scopes whose mode
                tags (`phaser:mode:` / `tng:input:` in .style, own or
                inherited from the nearest tagged ancestor) match it --
                see match_modes(). None means no mode filtering.
        """
        self.phil_scope = phil_scope
        self.exclude_scopes = set(exclude_scopes or [])
        self.max_expert_level = max_expert_level
        self.mode = mode
        # One generated item class per multiple scope, keyed by PHIL path
        self._item_classes = {}

    def convert(self, root_name="controlParameters"):
        """Convert the full PHIL scope to a CContainer hierarchy.

        Returns:
            CContainer populated with CData children mirroring the PHIL tree.
        """
        root = CContainer()
        root._name = root_name
        self._convert_scope(self.phil_scope, root)
        # Re-enable validation on all created objects
        self._reenable_validation(root)
        return root

    def _is_excluded(self, full_path):
        """Check if a PHIL path or any of its ancestors is excluded."""
        for excluded in self.exclude_scopes:
            if full_path == excluded or full_path.startswith(excluded + "."):
                return True
        return False

    def _convert_scope(self, scope, container, inherited_modes=None):
        """Recursively convert scope children to CData children.

        `inherited_modes` are the mode tags of the nearest tagged ancestor,
        which an untagged child takes as its own.
        """
        for obj in scope.objects:
            full_path = obj.full_path()

            if self._is_excluded(full_path):
                continue

            modes = parse_phil_style(obj.style)["modes"] or inherited_modes
            if self.mode is not None and not match_modes(self.mode, modes):
                continue

            if obj.is_definition:
                if obj.multiple:
                    self._convert_multiple_definition(obj, container)
                else:
                    self._convert_definition(obj, container)
            elif obj.is_scope:
                if obj.multiple:
                    self._convert_multiple_scope(obj, container, modes)
                    continue
                sub = CContainer()
                sub._name = full_path.replace(".", "__")
                self._apply_scope_qualifiers(obj, sub)
                setattr(container, sub._name, sub)
                if hasattr(container, 'declare_content'):
                    container.declare_content(sub._name)
                self._convert_scope(obj, sub, modes)

    def _convert_multiple_scope(self, scope, container, modes=None):
        """A `.multiple` scope: a CList of containers shaped like the scope."""
        lst = CList()
        lst._name = scope.full_path().replace(".", "__")
        self._apply_scope_qualifiers(scope, lst)
        lst.set_qualifier("philPath", scope.full_path())
        lst.set_qualifier("multiple", True)
        lst.set_qualifier("listMinLength", 0)
        lst.set_qualifier("subItem", {
            "class": self._item_class_for_scope(scope, modes),
            "qualifiers": {},
        })
        setattr(container, lst._name, lst)
        if hasattr(container, 'declare_content'):
            container.declare_content(lst._name)

    def _item_class_for_scope(self, scope, modes=None):
        """A CContainer subclass whose every instance is one repetition of
        `scope`, children converted exactly as a single scope's would be.

        Generated rather than declared with content(), because the defaults
        and qualifiers of the children come from the PHIL at runtime and
        _convert_definition already knows how to apply them; the class only
        has to reproduce that for each makeItem(). `modes` are the scope's
        effective mode tags, inherited by its untagged children.
        """
        path = scope.full_path()
        cls = self._item_classes.get(path)
        if cls is None:
            converter = self

            def __init__(self, parent=None, name=None, **kwargs):
                CContainer.__init__(self, parent=parent, name=name, **kwargs)
                converter._convert_scope(scope, self, modes)
                converter._reenable_validation(self)

            cls = type(path.replace(".", "__"), (CContainer,), {
                "__init__": __init__,
                "__module__": __name__,
                "__doc__": f"One repetition of the PHIL scope {path}.",
                "PHIL_SCOPE_PATH": path,
            })
            self._item_classes[path] = cls
        return cls

    def _convert_multiple_definition(self, keyword, container):
        """A `.multiple` definition: a CList of the leaf type, each item
        carrying the definition's qualifiers."""
        phil_type = keyword.type.phil_type
        value = keyword.extract()
        if phil_type == "bool" and str(value) not in ("True", "False"):
            phil_type = "ternary"
        cls = self.PHIL_TYPE_MAP.get(phil_type)
        if cls is None:
            cls = self.PHIL_CUSTOM_TYPE_MAP.get(phil_type, CString)

        # Let the ordinary path work out the item qualifiers on a template
        template = cls()
        template._skip_validation = True
        self._apply_definition_qualifiers(keyword, template, phil_type, value)
        item_qualifiers = {k: v for k, v in template._qualifiers.items()
                           if k != "multiple"}

        lst = CList()
        lst._name = keyword.full_path().replace(".", "__")
        for key in ("guiLabel", "toolTip", "expertLevel"):
            if key in item_qualifiers:
                lst.set_qualifier(key, item_qualifiers[key])
        lst.set_qualifier("philPath", keyword.full_path())
        lst.set_qualifier("multiple", True)
        lst.set_qualifier("listMinLength", 0)
        lst.set_qualifier("subItem", {"class": cls, "qualifiers": item_qualifiers})
        setattr(container, lst._name, lst)
        if hasattr(container, 'declare_content'):
            container.declare_content(lst._name)

    def _convert_definition(self, keyword, container):
        """Convert a single PHIL definition to a CData leaf node."""
        phil_type = keyword.type.phil_type
        value = keyword.extract()

        # Detect ternary booleans (bool type but value is not True/False)
        if phil_type == "bool" and str(value) not in ("True", "False"):
            phil_type = "ternary"

        # Create the CData object — check standard types first, then custom
        cls = self.PHIL_TYPE_MAP.get(phil_type)
        if cls is None:
            cls = self.PHIL_CUSTOM_TYPE_MAP.get(phil_type, CString)
        obj = cls()
        obj._skip_validation = True
        obj._name = keyword.full_path().replace(".", "__")

        # Store original PHIL path for reverse mapping at execution time
        obj.set_qualifier("philPath", keyword.full_path())

        # Apply PHIL attributes as qualifiers
        self._apply_definition_qualifiers(keyword, obj, phil_type, value)

        # Set default value with correct state tracking
        self._apply_default_value(obj, phil_type, value, keyword)

        # Add to parent container
        setattr(container, obj._name, obj)
        if hasattr(container, 'declare_content'):
            container.declare_content(obj._name)

    def _apply_definition_qualifiers(self, keyword, obj, phil_type, value):
        """Map PHIL definition attributes to CData qualifiers."""
        # guiLabel
        if keyword.short_caption is not None:
            obj.set_qualifier("guiLabel", _sanitize(keyword.short_caption))
        else:
            obj.set_qualifier("guiLabel", keyword.name)

        # toolTip
        if keyword.help is not None:
            obj.set_qualifier("toolTip", _sanitize(keyword.help))

        # expertLevel
        if keyword.expert_level is not None:
            obj.set_qualifier("expertLevel", keyword.expert_level)

        # style
        if keyword.style is not None:
            obj.set_qualifier("style", _sanitize(str(keyword.style)))

        # multiple (for scopes that allow repeated definitions)
        if keyword.multiple is not None:
            obj.set_qualifier("multiple", keyword.multiple)

        # Type-specific qualifiers
        if phil_type in ("choice", "choices"):
            is_multi = getattr(keyword.type, "multi", False)
            if is_multi:
                # Multi-choice doesn't map well to enumerators; use plain string
                pass
            elif hasattr(keyword, "words") and keyword.words:
                enumerators = [re.sub(r"\*", "", w.value) for w in keyword.words]
                obj.set_qualifier("enumerators", enumerators)
                obj.set_qualifier("onlyEnumerators", True)
                if keyword.caption is not None:
                    menu_items = [re.sub("_", " ", item)
                                  for item in keyword.caption.split()]
                    obj.set_qualifier("menuText", menu_items)

        elif phil_type == "ternary":
            obj.set_qualifier("enumerators", ["True", "False", _sanitize(str(value))])
            obj.set_qualifier("onlyEnumerators", True)

        elif phil_type in ("int", "float"):
            if keyword.type.value_min is not None:
                obj.set_qualifier("min", keyword.type.value_min)
            if keyword.type.value_max is not None:
                obj.set_qualifier("max", keyword.type.value_max)

        # After the type: a bound the type declares wins over a spinner's
        self._apply_style_qualifiers(keyword.style, obj, phil_type)

    def _apply_style_qualifiers(self, style, obj, phil_type):
        """Qualifiers from the GUI conventions in `.style`."""
        if not style:
            return
        parsed = parse_phil_style(style)
        if phil_type in ("int", "float"):
            for key in ("min", "max"):
                if parsed[key] is not None and obj.get_qualifier(key) is None:
                    obj.set_qualifier(key, parsed[key])
        if parsed["hidden"]:
            obj.set_qualifier("hidden", True)
        if parsed["multiLine"] and phil_type == "str":
            obj.set_qualifier("guiMode", "multiLine")
        if parsed["inputFile"]:
            obj.set_qualifier("philInputFile", True)
        if parsed["fileType"]:
            obj.set_qualifier("philFileType", parsed["fileType"])
        if parsed["directory"]:
            obj.set_qualifier("isDirectory", True)
        if parsed["modes"]:
            obj.set_qualifier("philModes", parsed["modes"])
        if parsed["ignored"]:
            obj.set_qualifier("philIgnored", True)

    def _apply_default_value(self, obj, phil_type, value, keyword):
        """Set the default value on a CData object with DEFAULT state."""
        if value is None:
            return

        try:
            if phil_type == "choice":
                if keyword.type.multi:
                    # Multi-choice: join as string
                    default_str = keyword.as_str().split("=")[1].strip()
                    obj.value = default_str
                elif isinstance(value, list):
                    obj.value = str(value[0]) if value else ""
                else:
                    obj.value = _sanitize(str(value))
            elif phil_type == "ternary":
                obj.value = _sanitize(str(value))
            elif phil_type == "bool":
                obj.value = value
            elif phil_type in ("int", "float"):
                # In PHIL, a default of 0 with value_min > 0 means "auto/unset".
                # Skip setting the default to avoid CData min/max validation errors
                # during clone/import — the program treats unset as auto anyway.
                min_val = obj.get_qualifier("min") if hasattr(obj, "get_qualifier") else None
                max_val = obj.get_qualifier("max") if hasattr(obj, "get_qualifier") else None
                if min_val is not None and value < min_val:
                    return
                if max_val is not None and value > max_val:
                    return
                obj.value = value
            elif phil_type in ("floats", "ints", "strings", "scatterers"):
                # List types: convert to space-separated string
                if isinstance(value, (list, tuple)):
                    obj.value = " ".join(str(v) for v in value)
                else:
                    obj.value = str(value)
            elif phil_type in ("choices",):
                # Multi-choice: extract from keyword string representation
                if isinstance(value, (list, tuple)):
                    obj.value = " ".join(str(v) for v in value)
                else:
                    default_str = keyword.as_str().split("=")[1].strip()
                    obj.value = default_str
            else:
                # str, path, filesystem, mtzcol, scatterer, unit_cell, uuid, etc.
                obj.value = str(value)

            # Mark as DEFAULT, not EXPLICITLY_SET
            if hasattr(obj, "_value_states"):
                obj._value_states["value"] = ValueState.DEFAULT

        except Exception as e:
            logger.debug("Could not set default for %s: %s", keyword.full_path(), e)

    def _apply_scope_qualifiers(self, scope, container):
        """Apply scope-level qualifiers (label, tooltip)."""
        if scope.short_caption is not None:
            container.set_qualifier("guiLabel", _sanitize(scope.short_caption))
        elif scope.name:
            container.set_qualifier("guiLabel", scope.name)
        if scope.help is not None:
            container.set_qualifier("toolTip", _sanitize(scope.help))
        if scope.expert_level is not None:
            container.set_qualifier("expertLevel", scope.expert_level)
        if scope.style is not None:
            container.set_qualifier("style", _sanitize(str(scope.style)))
            self._apply_style_qualifiers(scope.style, container, None)

    def _reenable_validation(self, obj):
        """Recursively re-enable validation after construction."""
        if hasattr(obj, "_skip_validation"):
            obj._skip_validation = False
        if hasattr(obj, "children"):
            for child in obj.children():
                if isinstance(child, CData):
                    self._reenable_validation(child)


def _sanitize(s):
    """Sanitize text for XML-safe display."""
    if s is None:
        return s
    return s.replace("<", "&lt;").replace(">", "&gt;")
