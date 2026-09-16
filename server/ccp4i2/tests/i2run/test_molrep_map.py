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

import gemmi
import numpy as np
import pytest

from .urls import emdb_half_map, emdb_map, rcsb_pdb
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
            for name in ("ORIGINALMODEL.pdb", "INVERTEDMODEL.pdb",
                         "ORIGINALTRIMMEDMAP.map", "INVERTEDTRIMMEDMAP.map",
                         "ORIGINALMASK.map", "INVERTEDMASK.map"):
                assert (job / name).exists(), f"missing {name}: {list(job.iterdir())}"

            # The report records a hand recommendation, a confidence verdict, and
            # a real-space map-model CC for at least one hand.
            root = ET.parse(str(job / "program.xml")).getroot()
            rec = root.find("recommendation")
            assert rec is not None, "no <recommendation> in program.xml"
            assert rec.get("hand") in {"Original", "Inverted"}
            assert rec.get("confidence") in {"confident", "ambiguous", "weak", "single"}
            assert rec.get("cc_original") or rec.get("cc_inverted"), \
                "no map-model CC recorded for either hand"


def test_cak_emd12042_halfmaps():
    """The half-map path with a real box mismatch.

    EMD-12042's primary map is a re-boxed 128^3 but its half maps are the full
    256^3 reconstruction (both origin-reset), so the emitted half maps only land
    on the model if the frames are reconciled by the central-crop offset. Assert
    the recommended model, shifted into the half-map frame, falls inside the
    cropped half map -- which it cannot if the reconciliation is wrong.
    """
    with download_map(emdb_map("12042")) as mapin, \
         download_map(emdb_half_map("12042", 1)) as half1, \
         download_map(emdb_half_map("12042", 2)) as half2, \
         download(rcsb_pdb("7b5o")) as model:
        args = [
            "molrep_map",
            "--MAPIN", mapin,
            "--XYZIN", model,
            "--HALFMAP1", half1,
            "--HALFMAP2", half2,
            "--SEARCH_RESOLUTION", "4.0",
            "--DOWNSAMPLE", "2",
        ]
        with i2run(args) as job:
            for name in ("HALFMAPOUT1.map", "HALFMAPOUT2.map"):
                assert (job / name).exists(), f"missing {name}: {list(job.iterdir())}"

            root = ET.parse(str(job / "program.xml")).getroot()
            hand = root.find("recommendation").get("hand")
            model_out = job / ("ORIGINALMODEL.pdb" if hand == "Original"
                               else "INVERTEDMODEL.pdb")

            primary = gemmi.read_ccp4_map(mapin)
            hmout = gemmi.read_ccp4_map(str(job / "HALFMAPOUT1.map"))
            hmout.setup(float("nan"), gemmi.MapSetup.Full)
            # central-crop offset primary(131.7) -> half(263.4) = 65.856 A
            off = tuple((hmout.grid.unit_cell.parameters[i]
                         - primary.grid.unit_cell.parameters[i]) / 2.0
                        for i in range(3))
            st = gemmi.read_structure(str(model_out))
            vals = [hmout.grid.interpolate_value(
                        gemmi.Position(a.pos.x + off[0], a.pos.y + off[1],
                                       a.pos.z + off[2]))
                    for m in st for c in m for r in c for a in r]
            finite = np.isfinite(np.array(vals)).mean()
            # The reconciled crop must cover the molecule: nearly every atom
            # (shifted into the half frame) falls inside the stored sub-block.
            assert finite > 0.9, f"only {finite:.0%} of atoms inside the cropped half map"
