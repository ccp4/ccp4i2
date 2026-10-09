"""A structure given to Phaser as already placed must be in THIS crystal.

Haiku, solving CDK4/cyclin D1, gave the cyclin chain of a homologous entry
(6p8e, cell 62.4 67.5 187.3) as the fixed structure for data of cell
57.8 64.7 186.1; Phaser held it in that frame and placed CDK4 against it
(LLG -892). The structure's cell and point group are now checked against
the data's before the job runs.
"""
from pathlib import Path

import gemmi
import pytest

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.tasks import get_plugin_class

DEMO = Path(__file__).resolve().parents[3] / "demo_data" / "beta_blip"


def _model(path, cell, space_group):
    atoms = [l for l in (DEMO / "beta.pdb").read_text().splitlines() if l.startswith("ATOM")]
    path.write_text(f"CRYST1{cell[0]:9.3f}{cell[1]:9.3f}{cell[2]:9.3f}"
                    f"{cell[3]:7.2f}{cell[4]:7.2f}{cell[5]:7.2f} {space_group:<11s}\n"
                    + "\n".join(atoms) + "\nEND\n")
    return path


def _plugin(tmp_path, fixed, input_fixed=True):
    plugin = get_plugin_class("phaser_simple_phil")(workDirectory=str(tmp_path), name="mr")
    inp = plugin.container.inputData
    inp.F_SIGF.setFullPath(str(DEMO / "beta_blip_P3221.mtz"))
    inp.XYZIN.setFullPath(str(DEMO / "blip.pdb"))
    inp.INPUT_FIXED.set(input_fixed)
    inp.XYZIN_FIXED.setFullPath(str(fixed))
    return plugin


def _blocking_on_fixed(error):
    return [r for r in error._reports
            if str(r.get("name", "")).endswith("XYZIN_FIXED")
            and r["severity"] >= CCP4ErrorHandling.SEVERITY_ERROR]


@pytest.fixture
def data_cell():
    mtz = gemmi.read_mtz_file(str(DEMO / "beta_blip_P3221.mtz"))
    c = mtz.cell
    return (c.a, c.b, c.c, c.alpha, c.beta, c.gamma)


def test_a_structure_from_another_crystal_is_refused(tmp_path, data_cell):
    other = (data_cell[0] * 1.08, data_cell[1] * 1.08, data_cell[2]) + data_cell[3:]
    fixed = _model(tmp_path / "homologue.pdb", other, "P 32 2 1")
    assert _blocking_on_fixed(_plugin(tmp_path, fixed).runTimeValidity())


def test_a_predicted_model_is_refused(tmp_path):
    # AlphaFold models carry a 1 A cell in P 1: placed in no crystal
    fixed = _model(tmp_path / "AF.pdb", (1, 1, 1, 90, 90, 90), "P 1")
    assert _blocking_on_fixed(_plugin(tmp_path, fixed).runTimeValidity())


def test_a_structure_placed_in_this_crystal_is_accepted(tmp_path, data_cell):
    fixed = _model(tmp_path / "placed.pdb", data_cell, "P 32 2 1")
    assert not _blocking_on_fixed(_plugin(tmp_path, fixed).runTimeValidity())


def test_a_file_left_set_without_input_fixed_is_ignored(tmp_path):
    fixed = _model(tmp_path / "AF.pdb", (1, 1, 1, 90, 90, 90), "P 1")
    assert not _blocking_on_fixed(_plugin(tmp_path, fixed, input_fixed=False).runTimeValidity())


def _seq(name):
    return "".join(l.strip() for l in (DEMO / f"{name}.seq").read_text().splitlines()
                   if not l.startswith(">"))


def _asu(path, names):
    """An AU contents file of one copy of each named demo protein (BETA, BLIP)."""
    items = "".join(
        f"<CAsuContentSeq><sequence>{_seq(n.lower())}</sequence><nCopies>1</nCopies>"
        f"<polymerType>PROTEIN</polymerType><name>{n}</name></CAsuContentSeq>" for n in names)
    path.write_text(
        '<?xml version="1.0"?><ns0:ccp4i2 xmlns:ns0="http://www.ccp4.ac.uk/ccp4ns">'
        "<ccp4i2_header><function>ASUCONTENT</function></ccp4i2_header>"
        f"<ccp4i2_body><seqList>{items}</seqList></ccp4i2_body></ns0:ccp4i2>")
    return path


def _complex(path):
    """beta-lactamase and BLIP as chains A and B of one file: a complex
    template, searched as one rigid body (CDK4/cyclin D1 was solved so)."""
    out = gemmi.Structure()
    model = gemmi.Model("1")
    for name, chain_id in (("beta", "A"), ("blip", "B")):
        chain = gemmi.read_structure(str(DEMO / f"{name}.pdb"))[0][0]
        chain.name = chain_id
        model.add_chain(chain)
    out.add_model(model)
    out.setup_entities()
    out.write_pdb(str(path))
    return path


def _advice(error):
    return [r for r in error._reports if r["code"] == 222]


def test_a_kind_no_model_accounts_for_gets_advice(tmp_path, data_cell):
    # Haiku's case: CDK4 and cyclin D1 in the contents, one searched for.
    # Here BLIP is searched for and beta-lactamase is in the contents too.
    plugin = _plugin(tmp_path, _model(tmp_path / "placed.pdb", data_cell, "P 32 2 1"),
                     input_fixed=False)
    plugin.container.inputData.ASUFILE.setFullPath(str(_asu(tmp_path / "two.asu.xml", ["BETA", "BLIP"])))
    advice = _advice(plugin.runTimeValidity())
    assert advice and advice[0]["severity"] == CCP4ErrorHandling.SEVERITY_WARNING
    details = advice[0]["details"]
    assert "search model 1 covers BLIP (100%)" in details
    assert "nothing accounts for BETA" in details
    assert "one rigid body" in details and "in turn" in details

    # one searched for and one already placed (here, in this crystal) covers both
    plugin.container.inputData.INPUT_FIXED.set(True)
    assert not _advice(plugin.runTimeValidity())


def test_a_complex_template_accounts_for_every_kind_it_holds(tmp_path, data_cell):
    # CDK4/cyclin D1 jobs 3 and 9: one search model holding both chains. The
    # earlier count of models told that correct job it searched for one of two.
    plugin = _plugin(tmp_path, _model(tmp_path / "unused.pdb", data_cell, "P 32 2 1"),
                     input_fixed=False)
    plugin.container.inputData.XYZIN.setFullPath(str(_complex(tmp_path / "beta_blip.pdb")))
    plugin.container.inputData.ASUFILE.setFullPath(str(_asu(tmp_path / "two.asu.xml", ["BETA", "BLIP"])))
    assert not _advice(plugin.runTimeValidity())


def test_contents_of_one_kind_get_no_advice(tmp_path, data_cell):
    plugin = _plugin(tmp_path, _model(tmp_path / "placed.pdb", data_cell, "P 32 2 1"),
                     input_fixed=False)
    plugin.container.inputData.ASUFILE.setFullPath(str(_asu(tmp_path / "one.asu.xml", ["BLIP"])))
    assert not _advice(plugin.runTimeValidity())


def test_the_fixed_structure_is_never_filled_from_the_context(tmp_path, data_cell):
    # Following on from mrparse put hit 1 in XYZIN and hit 2 in XYZIN_FIXED by
    # itself, and the agent then only had to tick INPUT_FIXED. That something
    # is already placed is a decision, so the slot stays empty until made.
    # keep the plugin: a container outlives its plugin only until collection
    plugin = _plugin(tmp_path, _model(tmp_path / "placed.pdb", data_cell, "P 32 2 1"))
    inp = plugin.container.inputData
    assert inp.XYZIN.qualifiers("fromPreviousJob") is True
    assert inp.XYZIN_FIXED.qualifiers("fromPreviousJob") is False


def test_the_rnp_pipeline_validates_without_ensembles(tmp_path):
    # phaser_rnp_pipeline_phil inherits the coverage check but builds its
    # ensembles at run time, so it has no ENSEMBLES to read: every job failed
    # validation with "'CContainer' object has no attribute 'FIXENSEMBLES'"
    # (Opus, 2026-10-09, project 30 job 6).
    plugin = get_plugin_class("phaser_rnp_pipeline_phil")(workDirectory=str(tmp_path), name="rnp")
    inp = plugin.container.inputData
    inp.F_SIGF.setFullPath(str(DEMO / "beta_blip_P3221.mtz"))
    inp.XYZIN_PARENT.setFullPath(str(DEMO / "beta.pdb"))
    inp.ASUFILE.setFullPath(str(_asu(tmp_path / "two.asu.xml", ["BETA", "BLIP"])))
    error = plugin.runTimeValidity()
    assert not _advice(error)
    assert not [r for r in error._reports if "FIXENSEMBLES" in str(r.get("details", ""))]
