# Copyright (C) 2026 University of York
import xml.etree.ElementTree as ET
from gemmi import read_mtz_file
from .utils import demoData, i2run


def test_gamma_native_vs_xe():
    """Scale Xe-derivative data against native using the gamma demo dataset.

    Both input files are merged anomalous intensities (I+/I-), so the wrapper
    exercises the IPAIR -> FMEAN conversion path inside myMakeHklInput before
    running scaleit.
    """
    args = ["scaleit"]
    args += ["--MERGEDFILES", demoData("gamma", "merged_intensities_native.mtz")]
    args += ["--MERGEDFILES", demoData("gamma", "merged_intensities_Xe.mtz")]

    with i2run(args) as job:
        # Output MTZ must exist and be readable
        read_mtz_file(str(job / "hklout.mtz"))

        # Check program.xml for expected scaleit output structure
        tree = ET.parse(job / "program.xml")
        root = tree.getroot()

        # Should have processed exactly 1 native + 1 derivative
        nderiv = root.find(".//SCALEITLOG/Nderivatives")
        assert nderiv is not None, "Nderivatives element missing from program.xml"
        assert int(nderiv.text) == 1, f"Expected 1 derivative, got {nderiv.text}"
