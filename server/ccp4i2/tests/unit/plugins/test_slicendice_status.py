"""SliceNDice's job status and annotation come from the pieces' search TFZs.

BAD1330 (a two-lobe periplasmic binding protein, 6fth Chainsaw model at 51%
identity, 1.26 A, 2026-10-09): both lobes placed at search TFZ 15.8 and
17.6, and modelcraft built the structure from that placement to R-free
0.253 -- yet ten REFMAC cycles left R-free at 0.509, so SliceNDice's own
test (R and R-free below 0.45) said "no solution" and the job ended
Unsatisfactory with its model labelled so. The solution block below is
that run's Phaser output, verbatim.
"""
import json
import types

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.tasks import get_plugin_class

SOLUTION = """\
** SINGLE solution

   Solution annotation (history):
   SOLU SET  RFZ=10.0 TFZ={tfz1} PAK=0 LLG=180 TFZ==17.5 RFZ=4.3 TFZ={tfz2} PAK=1 LLG=343 TFZ==19.6 LLG=463 TFZ==20.1 PAK=1
    LLG=463 TFZ==20.2
   SOLU SPAC P 21 21 21
   SOLU 6DIM ENSE pdb_6fth_chainsaw_noh_cluster_0 EULER  221.2  110.2   28.4 FRAC -0.37  0.91  0.21 BFAC -1.91#TFZ==17.5
   SOLU 6DIM ENSE pdb_6fth_chainsaw_noh_cluster_1 EULER  230.6  114.3   38.9 FRAC -0.68  0.94  0.13 BFAC  0.89#TFZ==20.2

   $$
"""


def _job(tmp_path, tfz1, tfz2, rfree):
    work = tmp_path / "slicendice_0"
    (work / "split_2").mkdir(parents=True)
    log = work / "split_2" / "split_2_phaser.log"
    log.write_text(SOLUTION.format(tfz1=tfz1, tfz2=tfz2))
    xyz = work / "split_2" / "split_2_refmac.pdb"
    xyz.write_text("CRYST1   50.036   96.476  100.074  90.00  90.00  90.00 P 21 21 21\nEND\n")
    (work / "slicendice_results.json").write_text(json.dumps({
        "slice": {"split_2": {"residues_ranges": {
            "0": [["M1:A:127", "M1:A:308"], ["M1:A:414", "M1:A:501"]],
            "1": [["M1:A:2", "M1:A:126"], ["M1:A:309", "M1:A:413"]]}}},
        "dice": {"split_2": {"phaser_llg": 463.0, "phaser_tfz": 20.1, "phaser_logfile": str(log),
                             "final_r_fact": rfree - 0.006, "final_r_free": rfree,
                             "xyzout": str(xyz), "hklout": str(work / "none.mtz")}}}))
    plugin = get_plugin_class("slicendice")(workDirectory=str(tmp_path), name="snd")
    plugin.splitHklout = lambda *a, **k: types.SimpleNamespace(maxSeverity=lambda: 0)
    return plugin


def _best(plugin):
    from xml.etree import ElementTree as ET
    best = ET.parse(plugin.makeFileName("PROGRAMXML")).find(".//RunInfo/Best")
    return {c.tag: c.text for c in best}


def test_every_piece_placed_finishes_whatever_the_rfree(tmp_path):
    plugin = _job(tmp_path, 15.8, 17.6, 0.5088)
    assert plugin.processOutputFiles() == CPluginScript.SUCCEEDED
    assert plugin.container.outputData.XYZOUT.annotation.startswith(
        "SliceNDice placement, not yet a solution: 2 splits, piece TFZ 15.8, 17.6, R-free 0.509")
    best = _best(plugin)
    assert (best["Solved"], best["Placed"], best["Partial"]) == ("False", "True", "False")


def test_a_piece_unplaced_is_a_partial_placement_and_unsatisfactory(tmp_path):
    # Lck's shape: one lobe at 26.4, the other at 6.0 with the refinement failing the 0.45 test
    plugin = _job(tmp_path, 26.4, 6.0, 0.5088)
    assert plugin.processOutputFiles() == CPluginScript.UNSATISFACTORY
    assert plugin.container.outputData.XYZOUT.annotation.startswith("SliceNDice partial placement: 2 splits, piece TFZ 26.4, 6.0")
    best = _best(plugin)
    assert (best["Solved"], best["Placed"], best["Partial"]) == ("False", "False", "True")


def test_slicendices_own_solution_is_still_a_solution(tmp_path):
    plugin = _job(tmp_path, 15.8, 17.6, 0.40)
    assert plugin.processOutputFiles() == CPluginScript.SUCCEEDED
    assert plugin.container.outputData.XYZOUT.annotation.startswith("SliceNDice solution: 2 splits")
    assert _best(plugin)["Solved"] == "True"


def test_a_space_group_change_takes_phasers_model_and_relabels_the_data(tmp_path):
    # BAD1330 merged as P 21 2 21, solved by Phaser in P 21 21 21 (SGALTERNATIVE
    # all): SliceNDice 0.1.3 reindexed the data with an axis permutation for
    # REFMAC while Phaser's model kept the data's setting, so the refinement
    # (R-free 0.562) and the refined XYZOUT's CRYST1 were of mismatched frames.
    import gemmi
    from pathlib import Path
    demo = Path(__file__).resolve().parents[3] / "demo_data" / "beta_blip"
    plugin = _job(tmp_path, 15.8, 17.6, 0.5619)
    inp = plugin.container.inputData
    inp.F_SIGF.setFullPath(str(demo / "beta_blip_P3221.mtz"))      # P 32 2 1
    inp.FREERFLAG.setFullPath(str(demo / "beta_blip_P3221.mtz"))
    plugin.hklin = str(demo / "beta_blip_P3221.mtz")
    phaser_dir = tmp_path / "slicendice_0" / "split_2"
    cell = gemmi.read_mtz_file(plugin.hklin).cell
    (phaser_dir / "split_2_phaser.pdb").write_text(
        "CRYST1%9.3f%9.3f%9.3f  90.00  90.00 120.00 P 31 2 1     6\nATOM      1  N   GLY A 127"
        "      10.000  10.000  10.000  1.00 20.00           N\nEND\n" % (cell.a, cell.b, cell.c))
    relabel_target = phaser_dir / "phaser_mr_output.1.mtz"
    mtz = gemmi.read_mtz_file(plugin.hklin)
    mtz.spacegroup = gemmi.find_spacegroup_by_name("P 31 2 1")
    mtz.write_to_file(str(relabel_target))
    assert plugin.processOutputFiles() == CPluginScript.SUCCEEDED
    out = plugin.container.outputData
    assert out.XYZOUT.annotation.startswith("SliceNDice placement in P 31 2 1, not the data's group: 2 splits, piece TFZ 15.8, 17.6")
    assert Path(str(out.XYZOUT.fullPath)).name == "split_2_phaser.pdb"
    assert gemmi.read_mtz_file(str(out.F_SIGF_OUT.fullPath)).spacegroup.hm == "P 31 2 1"
    assert gemmi.read_mtz_file(str(out.FREERFLAG_OUT.fullPath)).spacegroup.hm == "P 31 2 1"
    best = _best(plugin)
    assert (best["SpaceGroupInput"], best["SpaceGroup"], best["SpaceGroupChanged"], best["Solved"]) == \
        ("P 32 2 1", "P 31 2 1", "True", "False")
    assert not out.PERFORMANCEINDICATOR.RFree.isSet()
