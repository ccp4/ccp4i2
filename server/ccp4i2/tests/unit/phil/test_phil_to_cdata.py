"""
Tests for Phil2CData — runtime conversion of libtbx.phil scopes to CData hierarchies.
"""

import pytest

# See test_phil_plugin_script.py: libtbx is CCP4/cctbx-only, no pip wheel.
parse = pytest.importorskip("libtbx.phil", reason="needs libtbx (CCP4/cctbx)").parse

from ccp4i2.utils.phil_to_cdata import Phil2CData, parse_phil_style, match_modes
from ccp4i2.core.base_object.base_classes import CContainer, ValueState
from ccp4i2.core.base_object.fundamental_types import CInt, CFloat, CBoolean, CString, CList


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


# ---------------------------------------------------------------------------
# .multiple scopes and definitions
# ---------------------------------------------------------------------------

MULTIPLE_PHIL = parse("""
    composition {
      solvent = None
        .type = float
      chain
        .short_caption = "Macromolecular chain"
        .optional = True
        .multiple = True
      {
        chain_type = *protein na
          .type = choice
        nres = None
          .short_caption = "Number of residues"
          .type = int
        num = 1
          .type = int(value_min=1)
        dataset
          .multiple = True
        {
          label = None
            .type = str
        }
      }
    }
    model = None
      .short_caption = "Model file"
      .type = path
      .multiple = True
""")


class TestMultiple:
    """A .multiple scope is a list of scope-shaped containers, and a .multiple
    definition a list of leaves. Both start empty, as libtbx's fetch() leaves
    them unless the working phil supplies instances."""

    def setup_method(self):
        self.root = Phil2CData(MULTIPLE_PHIL).convert()

    def test_multiple_scope_is_an_empty_clist(self):
        chain = self.root.composition.composition__chain
        assert isinstance(chain, CList)
        assert len(chain) == 0
        assert chain.get_qualifier("multiple") is True
        assert chain.get_qualifier("philPath") == "composition.chain"
        assert chain.get_qualifier("guiLabel") == "Macromolecular chain"
        assert chain.get_qualifier("listMinLength") == 0

    def test_item_is_a_container_shaped_like_the_scope(self):
        chain = self.root.composition.composition__chain
        item = chain.makeItem()
        assert isinstance(item, CContainer)
        assert item.dataOrder() == [
            "composition__chain__chain_type", "composition__chain__nres",
            "composition__chain__num", "composition__chain__dataset"]
        assert item.composition__chain__num.get_qualifier("philPath") == "composition.chain.num"

    def test_item_carries_the_scope_defaults(self):
        item = self.root.composition.composition__chain.makeItem()
        assert item.composition__chain__num.value == 1
        assert item.composition__chain__num.getValueState() == ValueState.DEFAULT
        assert item.composition__chain__chain_type.value == "protein"
        assert item.composition__chain__nres.getValueState() == ValueState.NOT_SET
        assert item.composition__chain__num.get_qualifier("min") == 1

    def test_items_are_independent(self):
        chain = self.root.composition.composition__chain
        chain.append(chain.makeItem())
        chain.append(chain.makeItem())
        chain[0].composition__chain__nres.value = 120
        assert chain[1].composition__chain__nres.getValueState() == ValueState.NOT_SET
        assert chain[0].objectPath().endswith("composition__chain[0]")

    def test_nested_multiple_scope_inside_an_item(self):
        item = self.root.composition.composition__chain.makeItem()
        dataset = item.composition__chain__dataset
        assert isinstance(dataset, CList)
        assert len(dataset) == 0
        inner = dataset.makeItem()
        assert inner.dataOrder() == ["composition__chain__dataset__label"]

    def test_set_from_dicts_builds_items(self):
        # The client adds an item by sending the list with a new dict in it
        chain = self.root.composition.composition__chain
        chain.set([{"composition__chain__nres": 120, "composition__chain__num": 2},
                   {"composition__chain__nres": 50}])
        assert len(chain) == 2
        assert chain[0].composition__chain__num.value == 2
        assert chain[1].composition__chain__nres.value == 50
        assert chain[1].composition__chain__num.value == 1   # default kept

    def test_item_class_is_reused(self):
        chain = self.root.composition.composition__chain
        assert type(chain.makeItem()) is type(chain.makeItem())
        assert type(chain.makeItem()).PHIL_SCOPE_PATH == "composition.chain"

    def test_multiple_definition_is_a_clist_of_the_leaf_type(self):
        model = self.root.model
        assert isinstance(model, CList)
        assert len(model) == 0
        assert model.get_qualifier("philPath") == "model"
        assert model.get_qualifier("multiple") is True
        item = model.makeItem()
        assert isinstance(item, CString)
        assert item.get_qualifier("guiLabel") == "Model file"

    def test_single_scope_is_still_a_container(self):
        assert isinstance(self.root.composition, CContainer)
        assert isinstance(self.root.composition.composition__solvent, CFloat)


# ---------------------------------------------------------------------------
# .style: the libtbx GUI conventions become qualifiers
# ---------------------------------------------------------------------------

STYLE_PHIL = parse("""
    top
      .style = "phaser:mode:EP_AUTO box"
    {
      copies = 1
        .type = int
        .style = "spinner max=1000 min=1 bold"
      bounded = 5
        .type = int(value_min=2)
        .style = "spinner min=1"
      mute = None
        .type = bool
        .style = "hidden tribool"
      sequence = None
        .type = str
        .style = "height:48"
      seq_file = None
        .type = path
        .style = "input_file file_type:seq phaser:mode:MR*"
      out_dir = None
        .type = path
        .style = "directory"
      mode = *a b
        .type = choice
        .style = "bold phaser:ignore OnChange:update"
      resolution = None
        .type = float
        .style = "tng:input:+brf+frf+ftf"
      occupancy = None
        .type = float
        .style = "phaser:mode:MR_AUTO,MR_OCC,box auto_align"
    }
""")


class TestParsePhilStyle:

    def test_reads_every_convention(self):
        parsed = parse_phil_style("spinner max=1000 min=1 hidden height:48 "
                                  "input_file file_type:seq directory phaser:ignore")
        assert parsed["min"] == 1 and parsed["max"] == 1000
        assert parsed["hidden"] and parsed["multiLine"] and parsed["directory"]
        assert parsed["inputFile"] and parsed["fileType"] == "seq"
        assert parsed["ignored"]

    def test_both_mode_spellings(self):
        assert parse_phil_style("phaser:mode:MR*,EP_AUTO")["modes"] == ["MR*", "EP_AUTO"]
        assert parse_phil_style("tng:input:+brf+frf")["modes"] == ["brf", "frf"]

    def test_layout_word_after_a_comma_is_not_a_mode(self):
        assert parse_phil_style("phaser:mode:MR_AUTO,MR_OCC,box auto_align")["modes"] == [
            "MR_AUTO", "MR_OCC"]

    def test_ignores_what_it_does_not_know(self):
        parsed = parse_phil_style("bold box noauto OnChange:x renderer:y")
        assert parsed == parse_phil_style("")


class TestStyleQualifiers:

    def setup_method(self):
        # Hold the root: a container garbage-collected destroys its children
        self.root = Phil2CData(STYLE_PHIL).convert()
        self.top = self.root.top

    def test_spinner_bounds_where_the_type_gave_none(self):
        assert self.top.top__copies.get_qualifier("min") == 1
        assert self.top.top__copies.get_qualifier("max") == 1000

    def test_the_type_bound_wins_over_the_spinner(self):
        assert self.top.top__bounded.get_qualifier("min") == 2

    def test_hidden(self):
        assert self.top.top__mute.get_qualifier("hidden") is True
        assert self.top.top__copies.get_qualifier("hidden") is None

    def test_height_is_multiline_for_strings(self):
        assert self.top.top__sequence.get_qualifier("guiMode") == "multiLine"

    def test_input_files_are_tagged_not_converted(self):
        f = self.top.top__seq_file
        assert isinstance(f, CString)
        assert f.get_qualifier("philInputFile") is True
        assert f.get_qualifier("philFileType") == "seq"
        assert self.top.top__out_dir.get_qualifier("isDirectory") is True

    def test_modes_on_definitions_and_scopes(self):
        assert self.top.top__seq_file.get_qualifier("philModes") == ["MR*"]
        assert self.top.top__resolution.get_qualifier("philModes") == ["brf", "frf", "ftf"]
        assert self.top.top__occupancy.get_qualifier("philModes") == ["MR_AUTO", "MR_OCC"]
        assert self.top.get_qualifier("philModes") == ["EP_AUTO"]

    def test_ignored_is_recorded_not_hidden(self):
        assert self.top.top__mode.get_qualifier("philIgnored") is True
        assert self.top.top__mode.get_qualifier("hidden") is None

    def test_raw_style_is_kept(self):
        assert "OnChange:update" in self.top.top__mode.get_qualifier("style")


# ---------------------------------------------------------------------------
# Mode filtering
# ---------------------------------------------------------------------------

MODE_PHIL = parse("""
    tool {
      mode = *MR_AUTO EP_AUTO
        .type = choice
      title = None
        .type = str
      hklin = None
        .type = path
        .style = "phaser:mode:ANO,MR*"
      ensemble
        .multiple = True
        .style = "phaser:mode:MR*"
      {
        pdb = None
          .type = path
        rms = None
          .type = float
          .style = "phaser:mode:MR_AUTO"
      }
      crystal
        .style = "phaser:mode:EP_AUTO"
      {
        wavelength = None
          .type = float
      }
      keywords {
        resolution = None
          .type = float
        macmr
          .style = "phaser:mode:MR_AUTO,MR_RNP"
        {
          cycles = 50
            .type = int
        }
        xyzout = None
          .type = bool
          .style = "phaser:mode:*"
      }
    }
""")


class TestMatchModes:

    def test_untagged_applies_everywhere(self):
        assert match_modes("EP_AUTO", None) and match_modes("EP_AUTO", [])

    def test_exact_star_and_prefix(self):
        assert match_modes("MR_AUTO", ["ANO", "MR_AUTO"])
        assert match_modes("MR_FRF", ["*"])
        assert match_modes("MR_FRF", ["MR*"])
        assert not match_modes("EP_AUTO", ["MR*"])
        assert not match_modes("MR_AUTO", ["EP_AUTO"])


class TestModeFiltering:

    def names(self, mode):
        self.root = Phil2CData(MODE_PHIL, mode=mode).convert()
        found = []
        def walk(c):
            for n in c.dataOrder():
                o = getattr(c, n)
                found.append(n)
                if isinstance(o, CList):
                    walk(o.makeItem())
                elif hasattr(o, "dataOrder"):
                    walk(o)
        walk(self.root)
        return found

    def test_no_mode_keeps_everything(self):
        assert "tool__crystal" in self.names(None) and "tool__ensemble" in self.names(None)

    def test_mr_auto_sees_mr_things_and_shared_things(self):
        names = self.names("MR_AUTO")
        assert "tool__hklin" in names and "tool__ensemble" in names
        assert "tool__ensemble__rms" in names
        assert "tool__keywords__macmr" in names
        assert "tool__crystal" not in names
        assert "tool__title" in names and "tool__keywords__resolution" in names
        assert "tool__keywords__xyzout" in names

    def test_ep_auto_sees_ep_things_and_shared_things(self):
        names = self.names("EP_AUTO")
        assert "tool__crystal" in names and "tool__crystal__wavelength" in names
        assert "tool__hklin" not in names and "tool__ensemble" not in names
        assert "tool__keywords__macmr" not in names
        assert "tool__title" in names and "tool__keywords__xyzout" in names

    def test_an_untagged_child_inherits_the_scope_tag(self):
        # ensemble.pdb has no tag: MR* from its scope; ensemble.rms narrows to MR_AUTO
        names = self.names("MR_FRF")
        assert "tool__ensemble__pdb" in names
        assert "tool__ensemble__rms" not in names


def test_real_phaser_modes_partition_the_tree():
    pi = pytest.importorskip("phaser.phenix_interface", reason="needs phaser")
    def paths(mode):
        root = Phil2CData(pi.master_phil(), mode=mode).convert()
        out = set()
        def walk(c):
            for n in c.dataOrder():
                o = getattr(c, n)
                out.add(o.get_qualifier("philPath") or n.replace("__", "."))
                if hasattr(o, "dataOrder") and not isinstance(o, CList):
                    walk(o)
        walk(root)
        return out
    mr, ep = paths("MR_AUTO"), paths("EP_AUTO")
    assert {"phaser.ensemble", "phaser.search", "phaser.keywords.macmr"} <= mr
    assert {"phaser.crystal", "phaser.keywords.macsad", "phaser.keywords.llgcompletion"} <= ep
    assert not {"phaser.crystal", "phaser.keywords.macsad"} & mr
    assert not {"phaser.ensemble", "phaser.search", "phaser.keywords.macmr"} & ep
    assert {"phaser.composition", "phaser.keywords.resolution", "phaser.keywords.general.root"} <= mr & ep
