# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
"""
Tests for Phil2CData — runtime conversion of libtbx.phil scopes to CData hierarchies.
"""

import pytest
from libtbx.phil import parse

from ccp4i2.utils.phil_to_cdata import Phil2CData
from ccp4i2.core.base_object.base_classes import CContainer, ValueState
from ccp4i2.core.base_object.fundamental_types import CInt, CFloat, CBoolean, CString


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

SIMPLE_PHIL = parse("""
    refinement {
      resolution = 2.0
        .short_caption = "Resolution limit"
        .help = "High resolution limit for refinement"
        .type = float(value_min=0.5, value_max=20.0)
      cycles = 10
        .type = int
        .expert_level = 1
      use_solvent = True
        .type = bool
      output_prefix = "refined"
        .type = str
      model_file = None
        .type = path
    }
""")

CHOICE_PHIL = parse("""
    settings {
      mode = *auto manual careful
        .type = choice
        .caption = "Automatic" "Manual" "Careful"
      format = pdb *mmcif
        .type = choice
    }
""")

NESTED_PHIL = parse("""
    top {
      middle {
        bottom {
          value = 42
            .type = int
        }
      }
    }
""")

EXCLUSION_PHIL = parse("""
    data {
      hklin = None
        .type = path
    }
    refinement {
      resolution = 2.0
        .type = float
    }
    output {
      prefix = "out"
        .type = str
      directory = None
        .type = path
    }
""")


# ---------------------------------------------------------------------------
# Type mapping tests
# ---------------------------------------------------------------------------

class TestTypeMapping:
    """Test PHIL type → CData type conversion."""

    def test_float_type(self):
        container = Phil2CData(SIMPLE_PHIL).convert()
        obj = container.refinement.refinement__resolution
        assert isinstance(obj, CFloat)
        assert obj.get() == 2.0

    def test_int_type(self):
        container = Phil2CData(SIMPLE_PHIL).convert()
        obj = container.refinement.refinement__cycles
        assert isinstance(obj, CInt)
        assert obj.get() == 10

    def test_bool_type(self):
        container = Phil2CData(SIMPLE_PHIL).convert()
        obj = container.refinement.refinement__use_solvent
        assert isinstance(obj, CBoolean)
        assert obj.get() is True

    def test_str_type(self):
        container = Phil2CData(SIMPLE_PHIL).convert()
        obj = container.refinement.refinement__output_prefix
        assert isinstance(obj, CString)
        assert obj.get() == "refined"

    def test_path_type(self):
        container = Phil2CData(SIMPLE_PHIL).convert()
        obj = container.refinement.refinement__model_file
        assert isinstance(obj, CString)

    def test_choice_type(self):
        container = Phil2CData(CHOICE_PHIL).convert()
        obj = container.settings.settings__mode
        assert isinstance(obj, CString)
        assert obj.get() == "auto"


# ---------------------------------------------------------------------------
# Qualifier tests
# ---------------------------------------------------------------------------

class TestQualifiers:
    """Test PHIL attributes → CData qualifier mapping."""

    def test_gui_label_from_short_caption(self):
        container = Phil2CData(SIMPLE_PHIL).convert()
        obj = container.refinement.refinement__resolution
        assert obj.get_qualifier("guiLabel") == "Resolution limit"

    def test_gui_label_falls_back_to_name(self):
        container = Phil2CData(SIMPLE_PHIL).convert()
        obj = container.refinement.refinement__cycles
        assert obj.get_qualifier("guiLabel") == "cycles"

    def test_tooltip(self):
        container = Phil2CData(SIMPLE_PHIL).convert()
        obj = container.refinement.refinement__resolution
        assert "High resolution limit" in obj.get_qualifier("toolTip")

    def test_expert_level(self):
        container = Phil2CData(SIMPLE_PHIL).convert()
        obj = container.refinement.refinement__cycles
        assert obj.get_qualifier("expertLevel") == 1

    def test_min_max(self):
        container = Phil2CData(SIMPLE_PHIL).convert()
        obj = container.refinement.refinement__resolution
        assert obj.get_qualifier("min") == 0.5
        assert obj.get_qualifier("max") == 20.0

    def test_choice_enumerators(self):
        container = Phil2CData(CHOICE_PHIL).convert()
        obj = container.settings.settings__mode
        enums = obj.get_qualifier("enumerators")
        assert "auto" in enums
        assert "manual" in enums
        assert "careful" in enums
        assert obj.get_qualifier("onlyEnumerators") is True

    def test_choice_menu_text(self):
        container = Phil2CData(CHOICE_PHIL).convert()
        obj = container.settings.settings__mode
        menu = obj.get_qualifier("menuText")
        assert "Automatic" in menu
        assert "Manual" in menu
        assert "Careful" in menu

    def test_choice_without_caption(self):
        container = Phil2CData(CHOICE_PHIL).convert()
        obj = container.settings.settings__format
        enums = obj.get_qualifier("enumerators")
        assert "pdb" in enums
        assert "mmcif" in enums
        # No caption → no menuText (may be None or empty list)
        menu = obj.get_qualifier("menuText")
        assert menu is None or menu == []


# ---------------------------------------------------------------------------
# philPath qualifier
# ---------------------------------------------------------------------------

class TestPhilPath:
    """Test that original PHIL dotted paths are stored for reverse mapping."""

    def test_leaf_has_phil_path(self):
        container = Phil2CData(SIMPLE_PHIL).convert()
        obj = container.refinement.refinement__resolution
        assert obj.get_qualifier("philPath") == "refinement.resolution"

    def test_nested_leaf_has_phil_path(self):
        container = Phil2CData(NESTED_PHIL).convert()
        obj = container.top.top__middle.top__middle__bottom.top__middle__bottom__value
        assert obj.get_qualifier("philPath") == "top.middle.bottom.value"

    def test_choice_has_phil_path(self):
        container = Phil2CData(CHOICE_PHIL).convert()
        obj = container.settings.settings__mode
        assert obj.get_qualifier("philPath") == "settings.mode"


# ---------------------------------------------------------------------------
# Scope structure
# ---------------------------------------------------------------------------

class TestScopeStructure:
    """Test PHIL scope → CContainer hierarchy."""

    def test_scope_becomes_container(self):
        container = Phil2CData(SIMPLE_PHIL).convert()
        assert isinstance(container.refinement, CContainer)

    def test_nested_scopes(self):
        container = Phil2CData(NESTED_PHIL).convert()
        assert isinstance(container.top, CContainer)
        assert isinstance(container.top.top__middle, CContainer)
        assert isinstance(container.top.top__middle.top__middle__bottom, CContainer)

    def test_scope_name_uses_double_underscore(self):
        container = Phil2CData(NESTED_PHIL).convert()
        middle = container.top.top__middle
        assert middle._name == "top__middle"

    def test_root_name(self):
        container = Phil2CData(SIMPLE_PHIL).convert(root_name="myParams")
        assert container._name == "myParams"


# ---------------------------------------------------------------------------
# Default state tracking
# ---------------------------------------------------------------------------

class TestDefaultState:
    """Test that defaults use ValueState.DEFAULT, not EXPLICITLY_SET."""

    def test_default_value_not_user_set(self):
        container = Phil2CData(SIMPLE_PHIL).convert()
        obj = container.refinement.refinement__resolution
        # Has a default value
        assert obj.get() == 2.0
        # But isSet(allowDefault=False) should return False
        assert obj.isSet(allowDefault=False) is False

    def test_user_change_marks_as_set(self):
        container = Phil2CData(SIMPLE_PHIL).convert()
        obj = container.refinement.refinement__resolution
        assert obj.isSet(allowDefault=False) is False

        # Simulate user changing the value
        obj.value = 3.0
        assert obj.get() == 3.0
        assert obj.isSet(allowDefault=False) is True

    def test_all_defaults_not_user_set(self):
        container = Phil2CData(SIMPLE_PHIL).convert()
        for child in container.refinement.children():
            if not isinstance(child, CContainer):
                assert child.isSet(allowDefault=False) is False, \
                    f"{child._name} should not be marked as user-set"


# ---------------------------------------------------------------------------
# Scope exclusion
# ---------------------------------------------------------------------------

class TestScopeExclusion:
    """Test exclude_scopes filtering."""

    def test_excluded_scope_not_present(self):
        converter = Phil2CData(EXCLUSION_PHIL, exclude_scopes=["output"])
        container = converter.convert()
        assert hasattr(container, "refinement")
        assert hasattr(container, "data")
        assert not hasattr(container, "output")

    def test_excluded_children_not_present(self):
        converter = Phil2CData(EXCLUSION_PHIL, exclude_scopes=["output"])
        container = converter.convert()
        # output.prefix and output.directory should not exist
        children_names = [c._name for c in container.children()]
        assert "output" not in children_names

    def test_non_excluded_scopes_present(self):
        converter = Phil2CData(EXCLUSION_PHIL, exclude_scopes=["output"])
        container = converter.convert()
        assert container.refinement.refinement__resolution.get() == 2.0

    def test_no_exclusions(self):
        converter = Phil2CData(EXCLUSION_PHIL)
        container = converter.convert()
        assert hasattr(container, "output")
        assert hasattr(container, "refinement")
        assert hasattr(container, "data")

    def test_child_of_excluded_scope(self):
        """Excluding a parent scope should also exclude all its children."""
        converter = Phil2CData(EXCLUSION_PHIL, exclude_scopes=["output"])
        container = converter.convert()
        # Should have no "output__prefix" anywhere
        all_names = [c._name for c in _all_leaves(container)]
        assert all("output" not in n for n in all_names)


# ---------------------------------------------------------------------------
# Ternary boolean
# ---------------------------------------------------------------------------

class TestTernary:
    """Test ternary bool detection (bool with non-True/False default)."""

    def test_ternary_becomes_string_with_enumerators(self):
        ternary_phil = parse("""
            setting {
              auto_flag = Auto
                .type = bool
            }
        """)
        container = Phil2CData(ternary_phil).convert()
        obj = container.setting.setting__auto_flag
        assert isinstance(obj, CString)
        enums = obj.get_qualifier("enumerators")
        assert "True" in enums
        assert "False" in enums


# ---------------------------------------------------------------------------
# Validation re-enabled
# ---------------------------------------------------------------------------

class TestValidation:
    """Test that validation is re-enabled after construction."""

    def test_skip_validation_is_false_after_convert(self):
        container = Phil2CData(SIMPLE_PHIL).convert()
        obj = container.refinement.refinement__resolution
        assert obj._skip_validation is False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _all_leaves(container):
    """Recursively yield all non-CContainer children."""
    for child in container.children():
        if isinstance(child, CContainer):
            yield from _all_leaves(child)
        else:
            yield child
