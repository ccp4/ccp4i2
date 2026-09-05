import xml.etree.ElementTree as ET

import gemmi
import pytest

from .utils import demoData, i2run


def _args(parrot="False"):
    args = ["phaser_ep_phil"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN_HA", demoData("gamma", "heavy_atoms.pdb")]
    args += ["--ASUFILE", demoData("gamma", "gamma.asu.xml")]
    args += ["--COMP_BY", "ASU"]
    args += ["--WAVELENGTH", "1.542"]
    args += ["--LLGC_CYCLES", "20"]
    args += ["--ELEMENTS", "Xe"]
    args += ["--RUNPARROT", parrot]
    return args


@pytest.mark.order("first")
def test_gamma_both_hands_no_parrot():
    with i2run(_args()) as job:
        for name in ("ABCDOUT_1", "ABCDOUT_2", "MAPOUT_1", "MAPOUT_2"):
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        for name in ("XYZOUT_1", "XYZOUT_2"):          # original and inverted hands
            gemmi.read_pdb(str(job / f"{name}.pdb"))
        model = gemmi.read_structure(str(job / "XYZOUT_1.pdb"))[0]
        occs = [cra.atom.occ for cra in model.all()]
        assert sum(occ > 0.11 for occ in occs) == 3
        record = ET.parse(job / "program.xml").find("PhaserEpResults")
        assert len(record.findall("Hands/Hand")) == 2
        assert float(record.findtext("Overall/fom")) > 0.3


@pytest.mark.order("first")
def test_gamma_with_parrot():
    with i2run(_args(parrot="True")) as job:
        xml = ET.parse(job / "program.xml")
        assert xml.find("original/ParrotResult") is not None
        assert xml.find("inverted/ParrotResult") is not None
        assert sorted(f.name for f in job.glob("FPHIOUT_*.mtz")) == ["FPHIOUT_1.mtz", "FPHIOUT_2.mtz"]
        assert len(list(job.glob("ABCDOUT_*.mtz"))) == 4      # two from Phaser, two from parrot
