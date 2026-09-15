"""End-to-end molrep_map on real cryo-EM data (EMDB map + PDB model).

Places the deposited human CDK-activating kinase model (PDB 7B5O) into its own
EMDB map (EMD-12042, 2.5 A) -- a positive control: the model fits its map, so the
task should place both hands, score them by real-space map-model CC, and emit the
trimmed refinement package. This is also the first test to pull an EMDB fixture
through the harness (download_map + emdb_map).

Half maps are deliberately not supplied here: EMD-12042's primary map is re-boxed
(128^3) while its half maps are the full 256^3 reconstruction, so they do not share
a frame -- the half-map path needs box-independent trimming before that can be
tested (see the task's _prepare_half_maps).

Slow tier: needs the molrep binary and downloads ~10 MB (cached between runs).
"""

import shutil
from xml.etree import ElementTree as ET

import pytest

from .urls import emdb_map, rcsb_pdb
from .utils import download, download_map, i2run

pytestmark = pytest.mark.skipif(
    shutil.which("molrep") is None, reason="molrep binary not on PATH")


def test_cak_emd12042():
    with download_map(emdb_map("12042")) as mapin, \
         download(rcsb_pdb("7b5o")) as model:
        args = [
            "molrep_map",
            "--MAPIN", mapin,
            "--XYZIN", model,
            "--SEARCH_RESOLUTION", "4.0",   # triage speed
            "--DOWNSAMPLE", "2",            # coarse search grid
        ]
        with i2run(args) as job:
            # Both hands' refinement packages are emitted.
            for name in ("ORIGINALMODEL.pdb", "FLIPPEDMODEL.pdb",
                         "ORIGINALTRIMMEDMAP.map", "FLIPPEDTRIMMEDMAP.map",
                         "ORIGINALMASK.map", "FLIPPEDMASK.map"):
                assert (job / name).exists(), f"missing {name}: {list(job.iterdir())}"

            # The report records a hand recommendation, a confidence verdict, and
            # a real-space map-model CC for at least one hand.
            root = ET.parse(str(job / "program.xml")).getroot()
            rec = root.find("recommendation")
            assert rec is not None, "no <recommendation> in program.xml"
            assert rec.get("hand") in {"Original", "Flipped"}
            assert rec.get("confidence") in {"confident", "ambiguous", "weak", "single"}
            assert rec.get("cc_original") or rec.get("cc_flipped"), \
                "no map-model CC recorded for either hand"
