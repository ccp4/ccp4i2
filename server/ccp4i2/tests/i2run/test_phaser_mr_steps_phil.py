"""The molecular-replacement steps as separate Phaser jobs over PHIL:
rotation function, then translation function on its list, then the
packing test on those solutions. Each nests inside the last: the
harness's test database belongs to the outermost context."""
import pickle
import xml.etree.ElementTree as ET

import pytest

from .utils import demoData, i2run


def _args(task, copies):
    # beta of beta-blip: a search model that is not already in place (the
    # gamma model is, and Phaser then skips the rotation function for an
    # identity trial)
    args = [task]
    args += ["--F_SIGF", f"fullPath={demoData('beta_blip', 'beta_blip_P3221.mtz')}",
             "columnLabels=/*/*/[Fobs,Sigma]"]
    args += ["--ENSEMBLES", "label=beta", "use=True", f"number={copies}",
             "pdbItemList/identity_to_target=0.9",
             f"pdbItemList/structure={demoData('beta_blip', 'beta.pdb')}"]
    args += ["--COMP_BY", "DEFAULT"]
    return args


@pytest.mark.order("first")
def test_rotation_then_translation_then_packing():
    with i2run(_args("phaser_mr_frf_phil", copies=1)) as frf:
        rlist = frf / "RFILEOUT.phaser_rlist.pkl"
        assert rlist.is_file()
        with open(rlist, "rb") as handle:
            assert pickle.load(handle).is_rlist()
        xml = ET.parse(frf / "program.xml")
        rotations = xml.findall("Rotations/Rotation")
        assert rotations and all(r.findtext("Name") == "beta" for r in rotations)
        assert max(float(r.findtext("RFZ")) for r in rotations) > 3
        assert "*MR_FRF" in (frf / "working.phil").read_text()
        assert not (frf / "PHASER.1.pdb").exists()

        with i2run(_args("phaser_mr_ftf_phil", copies=0) + ["--RFILEIN", str(rlist)]) as ftf:
            solutions = ftf / "SOLOUT.phaser_sol.pkl"
            assert solutions.is_file()
            with open(solutions, "rb") as handle:
                assert not pickle.load(handle).is_rlist()
            xml = ET.parse(ftf / "program.xml")
            # A translation-function solution is scored by its TFZ; the LLG
            # field of the solution set is filled by refinement, later
            tfzs = [float(e.text) for e in xml.findall("Solutions/Solution/TFZ")]
            assert tfzs and max(tfzs) > 5
            assert "*MR_FTF" in (ftf / "working.phil").read_text()

            with i2run(_args("phaser_mr_pak_phil", copies=0) + ["--SOLIN", str(solutions)]) as pak:
                assert (pak / "SOLOUT.phaser_sol.pkl").is_file()
                xml = ET.parse(pak / "program.xml")
                packed = xml.findall("Solutions/Solution")
                assert packed and len(packed) <= len(tfzs)
                assert "*MR_PAK" in (pak / "working.phil").read_text()
