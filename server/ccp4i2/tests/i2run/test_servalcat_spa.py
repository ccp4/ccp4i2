"""servalcat SPA (single-particle) refinement, end-to-end.

Every other servalcat test drives DATA_METHOD=xtal, so the SPA output harvest
went untested -- and quietly broke when servalcat 0.4 renamed its SPA outputs
(refined_diffmap.mtz -> refined_maps.mtz, and the normalised maps lost their
_diffmap infix). A *successful* SPA refinement then failed to find the file it
wanted to split, recorded a warning, and returned FAILED anyway (CryoMapMR
job 15). This refines the deposited human CAK model (PDB 7B5O) against its own
EMD-12042 half maps and asserts the job completes with its outputs harvested --
i2run() already fails the test if diagnostic.xml carries any error report.

Slow tier: needs the servalcat binary and downloads the EMD-12042 half maps
(cached between runs).
"""

import shutil
from xml.etree import ElementTree as ET

import pytest

from .urls import emdb_half_map, rcsb_pdb
from .utils import download, download_map, i2run

pytestmark = pytest.mark.skipif(
    shutil.which("servalcat") is None, reason="servalcat binary not on PATH")


def test_spa_emd12042_harvests_outputs():
    with download_map(emdb_half_map("12042", 1)) as hm1, \
         download_map(emdb_half_map("12042", 2)) as hm2, \
         download(rcsb_pdb("7b5o")) as model:
        args = [
            "servalcat",
            "--DATA_METHOD", "spa",
            "--XYZIN", model,
            "--MAPIN1", hm1,
            "--MAPIN2", hm2,
            "--RES_MIN", "2.5",
            "--NCYCLES", "1",
        ]
        # i2run() itself asserts diagnostic.xml carries no error report, so a
        # green run is the core regression proof: before the fix the harvest
        # tripped 207/991 on the renamed outputs and the job came back FAILED.
        with i2run(args) as job:
            # The reflection file was resolved under its servalcat-0.4 name...
            assert (job / "refined_maps.mtz").exists()
            # ...and split into the map-coefficient files (the step that failed).
            assert (job / "FPHIOUT.mtz").exists()
            assert (job / "DIFFPHIOUT.mtz").exists()
            xml = ET.parse(job / "program.xml")
            assert xml.getroot() is not None
