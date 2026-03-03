"""
Integration tests for phasertng_picard — the PhaserTNG Picard wrapper.

These tests require phasertng to be installed. They are skipped if phasertng
is not available.
"""

import json
import os
import pytest

from .utils import demoData, i2run

# Check if phasertng is available
try:
    from phasertng.programs import picard
    _has_phasertng = True
except ImportError:
    _has_phasertng = False

pytestmark = pytest.mark.skipif(
    not _has_phasertng,
    reason="phasertng not installed"
)


# ---------------------------------------------------------------------------
# Plugin instantiation tests
# ---------------------------------------------------------------------------

@pytest.mark.order("first")
def test_instantiation():
    """Test that phasertng_picard plugin instantiates with PHIL parameters."""
    from ccp4i2.core.tasks import get_plugin_class

    cls = get_plugin_class("phasertng_picard")
    assert cls is not None, "phasertng_picard should be registered in tasks.py"

    plugin = cls()

    # Should have rich file types in inputData
    input_children = list(plugin.container.inputData.children())
    input_names = [c._name for c in input_children]
    assert "HKLIN" in input_names
    assert "XYZIN" in input_names
    assert "ASUIN" in input_names
    assert "DICT" in input_names

    # Should have PHIL scope containers in controlParameters
    # The PHIL tree has top-level scopes (picard, phasertng) with nested params
    cp_children = list(plugin.container.controlParameters.children())
    assert len(cp_children) >= 2, \
        f"Expected at least 2 PHIL scope containers, got {len(cp_children)}"

    # Count total leaf parameters recursively
    def count_leaves(container):
        total = 0
        for child in container.children():
            if hasattr(child, 'children') and list(child.children()):
                total += count_leaves(child)
            else:
                total += 1
        return total

    leaf_count = count_leaves(plugin.container.controlParameters)
    assert leaf_count > 100, \
        f"Expected hundreds of PHIL leaf parameters, got {leaf_count}"


@pytest.mark.order("first")
def test_json_serialization():
    """Test that the container serializes to JSON correctly for the frontend."""
    from ccp4i2.core.tasks import get_plugin_class
    from ccp4i2.lib.utils.containers.json_encoder import CCP4i2JsonEncoder

    cls = get_plugin_class("phasertng_picard")
    plugin = cls()

    # Serialize to JSON (CCP4i2JsonEncoder handles CData objects directly)
    json_str = json.dumps(plugin.container, cls=CCP4i2JsonEncoder)
    data = json.loads(json_str)

    assert "inputData" in str(data)
    assert "controlParameters" in str(data)

    # Count leaf parameters in JSON — a leaf is a _class/_value node
    # where _value is not a dict (i.e., not a container)
    def count_leaves(obj, depth=0):
        if depth > 100:
            return 0
        if isinstance(obj, dict):
            if "_class" in obj and "_value" in obj:
                val = obj["_value"]
                if isinstance(val, dict):
                    # Container — recurse into children
                    return sum(count_leaves(v, depth + 1) for v in val.values())
                elif isinstance(val, list):
                    return sum(count_leaves(v, depth + 1) for v in val)
                else:
                    # Leaf parameter (scalar value)
                    return 1
            return sum(count_leaves(v, depth + 1) for v in obj.values())
        if isinstance(obj, list):
            return sum(count_leaves(v, depth + 1) for v in obj)
        return 0

    leaf_count = count_leaves(data)
    assert leaf_count > 100, \
        f"Expected >100 leaf parameters in JSON, got {leaf_count}"


@pytest.mark.order("first")
def test_extract_no_user_changes():
    """Initially no parameters should be marked as user-modified."""
    from ccp4i2.core.tasks import get_plugin_class

    cls = get_plugin_class("phasertng_picard")
    plugin = cls()

    params = plugin.extract_phil_parameters()
    assert params == [], \
        f"Expected no user-modified params, got {len(params)}: {params[:5]}"


@pytest.mark.order("first")
def test_phil_exclude_scopes():
    """Verify excluded scopes are not in controlParameters."""
    from ccp4i2.core.tasks import get_plugin_class

    cls = get_plugin_class("phasertng_picard")
    plugin = cls()

    # Collect all philPath qualifiers
    phil_paths = []

    def _collect(container):
        for child in container.children():
            if isinstance(child, type(container)):
                _collect(child)
            else:
                path = child.get_qualifier("philPath")
                if path:
                    phil_paths.append(path)

    _collect(plugin.container.controlParameters)

    # None of the excluded scopes should appear
    excluded_prefixes = plugin.PHIL_EXCLUDE_SCOPES
    for path in phil_paths:
        for excl in excluded_prefixes:
            assert not path.startswith(excl), \
                f"Excluded path {excl} found in controlParameters: {path}"


# ---------------------------------------------------------------------------
# i2run execution test
# ---------------------------------------------------------------------------

@pytest.mark.order("first")
def test_i2run_gamma():
    """Test phasertng_picard molecular replacement with gamma demo data.

    Gamma pyrone is a simple MR case — one copy in the ASU, high-quality
    data with xenon SAD. Should solve in ~15 seconds with R-factor < 30%.
    """
    args = ["phasertng_picard"]
    args += ["--HKLIN", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--ASUIN", demoData("gamma", "gamma.asu.xml")]

    with i2run(args) as job:
        # Verify working.phil was generated
        assert (job / "working.phil").exists(), \
            "working.phil should be generated in the job directory"

        # Verify phasertng database directory was created
        db_dir = job / "phasertng_picard"
        assert db_dir.is_dir(), \
            "phasertng_picard database directory should exist"

        # Verify best solution PDB was produced
        import gemmi
        best_pdb = db_dir / "best.1.coordinates.pdb"
        assert best_pdb.exists(), \
            "best.1.coordinates.pdb should exist in database directory"
        structure = gemmi.read_pdb(str(best_pdb))
        assert len(structure) > 0, "Output PDB should contain at least one model"
