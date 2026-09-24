"""dm_multidomain must complain on the field, not 200 lines into a run.

The assembly (which chains are copies of each other) and the rigid bodies
(which residues move as one unit) are cross-referenced by a role name the user
invents, and a role that matches nothing used to surface only as a traceback
from processInputFiles. These tests pin the checks that catch it while the
parameters are still being edited, and the preview that gives the interface
the one number that says whether a body is really a body: the superposition
RMSD of its residues against each copy.

CCP4-free: gemmi reads the demo model, no binary is run.
"""
import os

import pytest

import ccp4i2
from ccp4i2.core.tasks import get_plugin_class

pytest.importorskip("gemmi", reason="needs gemmi")

_DEMO = os.path.join(os.path.dirname(ccp4i2.__file__), "demo_data")
_CDK = os.path.join(_DEMO, "CDK1CyclinBCKS2", "1jst.pdb")
_MTZ = os.path.join(_DEMO, "CDK1CyclinBCKS2", "CDK1CyclinBCKS2Compound23.mtz")

# CDK2/cyclin A, two copies of the hetero-dimer: A+B and C+D.
ASSEMBLY = ["A=A B=B", "A=C B=D"]


def _plugin(tmp_path, assembly=None, bodies=(), with_data=False):
    if not os.path.exists(_CDK):
        pytest.skip("CDK/cyclin demo model not present")
    plugin = get_plugin_class("dm_multidomain")(
        workDirectory=str(tmp_path), name="dm_multidomain")
    plugin.container.inputData.XYZIN.setFullPath(_CDK)
    if with_data:
        plugin.container.inputData.F_SIGF.setFullPath(_MTZ)
    ctrl = plugin.container.controlParameters
    ctrl.PHASE_SOURCE.set("model")      # so ABCD is not required
    if assembly is not None:
        ctrl.ASSEMBLY.set(assembly)
    # A fresh container arrives with one blank DOMAINS row (listMinLength=1),
    # which is exactly what the interface has to replace on a new job.
    while len(ctrl.DOMAINS):
        ctrl.DOMAINS.remove(ctrl.DOMAINS[0])
    for spec, mode in bodies:
        ctrl.DOMAINS.append(ctrl.DOMAINS.makeItem())
        ctrl.DOMAINS[-1].segments.set(spec)
        ctrl.DOMAINS[-1].mode.set(mode)
    return plugin


# Codes this task raises itself. 299 ("program not found") comes from the
# framework and is expected in a CCP4-free environment, so it is not one.
_OWN_CODES = range(210, 230)


def _codes(report):
    return [e.get("code") for e in report.getErrors()]


def _own_codes(report):
    return [c for c in _codes(report) if c in _OWN_CODES]


def _detail(report, code):
    return next(str(e["details"]) for e in report.getErrors()
                if e.get("code") == code)


def _names(report):
    return {e.get("name") for e in report.getErrors()}


# -- validity(): strings only, polled while editing --------------------------


def test_body_naming_an_unknown_role_is_an_error(tmp_path):
    """The failure this task is most prone to: the role in the body does not
    match any role in the assembly, so nothing lines up and nothing says so."""
    report = _plugin(tmp_path, ASSEMBLY, [("CDK:1-298", "average")]).validity()
    assert 216 in _codes(report)
    assert "CDK" in _detail(report, 216)


def test_error_names_point_at_the_field(tmp_path):
    """Field-level display needs {task}.container.{section}.{field}."""
    report = _plugin(tmp_path, ASSEMBLY, [("CDK:1-298", "average")]).validity()
    assert ("dm_multidomain.container.controlParameters.DOMAINS"
            in _names(report))


def test_one_copy_is_not_an_assembly(tmp_path):
    report = _plugin(tmp_path, ["A=A B=B"], [("A:1-298", "average")]).validity()
    assert 210 in _codes(report)


def test_a_chain_belongs_to_exactly_one_copy(tmp_path):
    report = _plugin(tmp_path, ["A=A B=B", "A=A B=D"],
                     [("A:1-298", "average")]).validity()
    assert 211 in _codes(report)


def test_reversed_range_is_an_error(tmp_path):
    report = _plugin(tmp_path, ASSEMBLY, [("A:298-1", "average")]).validity()
    assert 215 in _codes(report)


def test_excluding_everything_leaves_nothing_to_average(tmp_path):
    report = _plugin(tmp_path, ASSEMBLY, [("A:1-298", "exclude")]).validity()
    assert 217 in _codes(report)


def test_overlapping_bodies_warn_but_do_not_block(tmp_path):
    """The masks get split by nearest atom, so this runs -- but two bodies
    claiming the same residues is almost always a typo."""
    report = _plugin(tmp_path, ASSEMBLY,
                     [("A:1-298", "average"), ("A:1-100", "average")]
                     ).validity()
    warning = next(e for e in report.getErrors() if e.get("code") == 218)
    assert warning["severity"] == 2          # SEVERITY_WARNING, not blocking


def test_an_empty_assembly_means_detect_and_is_not_an_error(tmp_path):
    """Leaving ASSEMBLY empty is legal -- the model is asked instead -- so
    validity must not invent a complaint about it."""
    report = _plugin(tmp_path, None, [("1-298", "average")]).validity()
    assert not {210, 211, 216} & set(_codes(report))


def test_a_clean_setup_adds_nothing(tmp_path):
    report = _plugin(tmp_path, ASSEMBLY, [("A:1-298", "average")]).validity()
    assert not _own_codes(report)


# -- ncs_preview(): what the interface shows the user ------------------------


def test_preview_suggests_the_assembly_the_model_implies(tmp_path):
    preview = _plugin(tmp_path).ncs_preview()
    assert preview["ok"]
    assert preview["suggestion"]["assembly"] == ASSEMBLY
    assert preview["suggestion"]["segments"] == "A:1-298,B:175-432"
    assert preview["assembly"]["source"] == "detected"


def test_preview_reports_matched_atoms_and_rmsd_per_copy(tmp_path):
    """The number that answers 'do these residues really move as one unit'."""
    preview = _plugin(tmp_path, ASSEMBLY,
                      [("A:1-298", "average")]).ncs_preview()
    body = preview["bodies"][0]
    assert body["segments"] == [{"role": "A", "lo": 1, "hi": 298}]
    assert [c["label"] for c in body["copies"]] == ["C+D"]
    assert body["copies"][0]["nCA"] == 298
    assert body["copies"][0]["rmsd"] < 2.0


def test_preview_handles_a_cross_chain_body(tmp_path):
    """A body spanning two roles is the reason roles exist; its CA count is
    the sum of both segments."""
    preview = _plugin(tmp_path, ASSEMBLY,
                      [("A:1-100,B:175-200", "average")]).ncs_preview()
    body = preview["bodies"][0]
    assert body["nReferenceCA"] == body["copies"][0]["nCA"] == 126


def test_preview_never_raises_without_a_model(tmp_path):
    plugin = get_plugin_class("dm_multidomain")(
        workDirectory=str(tmp_path), name="dm_multidomain")
    preview = plugin.ncs_preview()
    assert preview["ok"] is False and preview["error"]


def test_preview_flags_a_chain_the_model_does_not_have(tmp_path):
    preview = _plugin(tmp_path, ["A=A B=B", "A=Z B=D"],
                      [("A:1-298", "average")]).ncs_preview()
    assert 221 in [m["code"] for m in preview["messages"]]


# -- runTimeValidity(): the expensive checks, at submission ------------------


def test_runtime_validity_rejects_a_range_with_no_atoms(tmp_path):
    report = _plugin(tmp_path, ASSEMBLY, [("A:900-950", "average")],
                     with_data=True).runTimeValidity()
    assert 222 in _codes(report)


def test_runtime_validity_rejects_a_body_with_nothing_to_average_against(
        tmp_path):
    report = _plugin(tmp_path, ["A=A B=B", "A=Z B=D"],
                     [("A:1-298", "average")],
                     with_data=True).runTimeValidity()
    assert 223 in _codes(report)


def test_runtime_validity_passes_a_real_setup(tmp_path):
    report = _plugin(tmp_path, ASSEMBLY,
                     [("A:1-298", "average"), ("B:175-432", "average")],
                     with_data=True).runTimeValidity()
    assert not _own_codes(report)


def test_a_fresh_job_starts_with_one_blank_body(tmp_path):
    """listMinLength=1 means a new job opens with an empty rigid body, which
    validity reports -- the interface fills it rather than leaving the user
    looking at an empty box."""
    plugin = get_plugin_class("dm_multidomain")(
        workDirectory=str(tmp_path), name="dm_multidomain")
    assert len(plugin.container.controlParameters.DOMAINS) == 1
    assert plugin.container.controlParameters.DOMAINS[0].segments_spec() == ""
    assert 213 in _codes(plugin.validity())
