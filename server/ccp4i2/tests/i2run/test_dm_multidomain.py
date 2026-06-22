import xml.etree.ElementTree as ET

import gemmi

from .utils import demoData, i2run


def test_dm_multidomain():
    """Multi-domain NCS averaging on the AHIR driver case (6 NCS copies, T sym).

    The DARPIN / NADP-knob / helical-shell domains follow different NCS
    transformations, so each gets its own operator set and mask. Observed
    intensities are French-Wilsoned to F by makeHklin0; starting phases are
    supplied (ABCD).
    """
    args = ["dm_multidomain"]
    args += ["--F_SIGF", demoData("ahir", "ahir_observed.mtz")]
    args += ["--ABCD", demoData("ahir", "ahir_phases.mtz")]
    args += ["--FREERFLAG", demoData("ahir", "ahir_freer.mtz")]
    args += ["--XYZIN", demoData("ahir", "ahir_model.cif")]
    # homomer: ASSEMBLY auto-detected (chains A-F); segments use the implicit
    # role, so they are bare residue ranges on the reference copy.
    args += ["--DOMAINS", "segments=340-485", "mode=average"]
    args += ["--DOMAINS", "segments=140-339", "mode=refine"]
    args += ["--DOMAINS", "segments=13-139", "mode=exclude"]
    with i2run(args) as job:
        for name in ["ABCDOUT", "FPHIOUT"]:
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        ET.parse(job / "program.xml")
        # the two averaged bodies each yield a captured mask; the excluded
        # 13-139 yields none. The masks must be DISJOINT (nearest-atom
        # partition) -- dm needs non-overlapping NCS masks.
        import numpy as np
        masks = sorted(job.glob("mask_*.map"))
        assert [m.name for m in masks] == ["mask_140_339.map", "mask_340_485.map"]
        g0 = np.array(gemmi.read_ccp4_mask(str(masks[0])).grid, copy=False) > 0
        g1 = np.array(gemmi.read_ccp4_mask(str(masks[1])).grid, copy=False) > 0
        assert int(np.count_nonzero(g0 & g1)) == 0


def test_dm_multidomain_phases_from_model():
    """Same case but starting phases are CALCULATED from the model with
    servalcat sigmaa (bulk solvent + sigmaA weighting) -- no ABCD supplied.
    """
    args = ["dm_multidomain"]
    args += ["--F_SIGF", demoData("ahir", "ahir_observed.mtz")]
    args += ["--XYZIN", demoData("ahir", "ahir_model.cif")]
    args += ["--PHASE_SOURCE", "model"]
    args += ["--DOMAINS", "segments=340-485", "mode=average"]
    args += ["--DOMAINS", "segments=140-339", "mode=average"]
    args += ["--DOMAINS", "segments=13-139", "mode=average"]
    with i2run(args) as job:
        for name in ["ABCDOUT", "FPHIOUT"]:
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        ET.parse(job / "program.xml")


def test_dm_multidomain_explicit_assembly():
    """Same homomer case but with an EXPLICIT ASSEMBLY (named role 'm' over
    chains A-F) and role-qualified segments -- exercises the ASSEMBLY parsing
    and the (role, seqid) correspondence path end-to-end, rather than the
    auto-detect / implicit-role shortcut.
    """
    args = ["dm_multidomain"]
    args += ["--F_SIGF", demoData("ahir", "ahir_observed.mtz")]
    args += ["--ABCD", demoData("ahir", "ahir_phases.mtz")]
    args += ["--XYZIN", demoData("ahir", "ahir_model.cif")]
    for ch in ["A", "B", "C", "D", "E", "F"]:
        args += ["--ASSEMBLY", f"m={ch}"]
    args += ["--DOMAINS", "segments=m:340-485", "mode=average"]
    args += ["--DOMAINS", "segments=m:140-339", "mode=average"]
    with i2run(args) as job:
        for name in ["ABCDOUT", "FPHIOUT"]:
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        ET.parse(job / "program.xml")
