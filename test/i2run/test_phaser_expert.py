from pathlib import Path
import xml.etree.ElementTree as ET
import gemmi
from .utils import demoData, i2run


def _beta_blip_args():
    args = ["phaser_pipeline"]
    args += [
        "--F_SIGF",
        f"fullPath={demoData('beta_blip', 'beta_blip_P3221.mtz')}",
        "columnLabels=/*/*/[Fobs,Sigma]",
    ]
    args += [
        "--ENSEMBLES",
        "use=True",
        "pdbItemList/identity_to_target=0.9",
        f"pdbItemList/structure={demoData('beta_blip', 'beta.pdb')}",
    ]
    args += [
        "--ENSEMBLES",
        "use=True",
        "pdbItemList/identity_to_target=0.9",
        f"pdbItemList/structure={demoData('beta_blip', 'blip.pdb')}",
    ]
    args += ["--F_OR_I", "F"]
    return args


def _check_output(job: Path, include_refinement: bool=True):
    gemmi.read_pdb(str(job / f"PHASER.1.pdb"))
    xml = ET.parse(job / "program.xml")
    if include_refinement:
        gemmi.read_pdb(str(job / "XYZOUT_REFMAC.pdb"))
        gemmi.read_pdb(str(job / "XYZOUT_SHEETBEND.pdb"))
        rworks = [float(e.text) for e in xml.iter("r_factor")]
        assert min(rworks) < 0.3
    for mtz in ["DIFMAPOUT_1", "MAPOUT_1", "PHASEOUT_1"]:
        gemmi.read_mtz_file(str(job / f"{mtz}.mtz"))
    llgs = [float(e.text) for e in xml.findall(".//Solution/LLG")]
    assert max(llgs) > 1000


def test_beta_blip_default():
    args = _beta_blip_args()
    args += ["--COMP_BY", "DEFAULT"]
    with i2run(args) as job:
        _check_output(job)


def test_beta_blip_asu():
    args = _beta_blip_args()
    args += [
        "--ASUFILE",
        f"seqFile={demoData('beta_blip', 'beta.seq')}",
        f"seqFile={demoData('beta_blip', 'blip.seq')}",
    ]
    args += ["--RUNSHEETBEND", "False"]
    args += ["--RUNREFMAC", "False"]
    with i2run(args) as job:
        _check_output(job, include_refinement=False)
