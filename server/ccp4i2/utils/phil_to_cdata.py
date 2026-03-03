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
from ccp4i2.core.base_object.fundamental_types import CInt, CFloat, CBoolean, CString

logger = logging.getLogger(__name__)


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

    def __init__(self, phil_scope, exclude_scopes=None):
        """
        Args:
            phil_scope: A libtbx.phil scope object (the master_phil).
            exclude_scopes: List of dotted PHIL paths to skip entirely
                (e.g., I/O scopes handled by shims).
        """
        self.phil_scope = phil_scope
        self.exclude_scopes = set(exclude_scopes or [])

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

    def _convert_scope(self, scope, container):
        """Recursively convert scope children to CData children."""
        for obj in scope.objects:
            full_path = obj.full_path()

            if self._is_excluded(full_path):
                continue

            if obj.is_definition:
                self._convert_definition(obj, container)
            elif obj.is_scope:
                sub = CContainer()
                sub._name = full_path.replace(".", "__")
                self._apply_scope_qualifiers(obj, sub)
                setattr(container, sub._name, sub)
                self._convert_scope(obj, sub)

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
