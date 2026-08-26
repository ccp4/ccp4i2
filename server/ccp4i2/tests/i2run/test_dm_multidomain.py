import xml.etree.ElementTree as ET

import gemmi
import numpy as np

from .utils import demoData, i2run


def _voxels_in_crystal_frame(mask_path):
    """A mask's set voxels placed on the parent grid it was cut from.

    Each mask is written as a tight sub-volume --- its own body's bounding box
    on the sampling all of them share --- so two masks are different shapes and
    cannot be compared element by element. NSTART (header words 5-7) says where
    the block sits on the parent grid and MX/MY/MZ (8-10) give that grid's
    sampling, exactly as tests/unit/plugins/test_dm_ncs_lib.py describes. Place
    both blocks on the parent grid and the question "do these overlap" can be
    asked at all.

    Comparing the raw arrays instead --- which this test used to do --- worked
    only while every body happened to produce the same box, and became a
    ValueError about mismatched shapes the moment they did not.
    """
    ccp4 = gemmi.read_ccp4_mask(str(mask_path))
    assert [ccp4.header_i32(w) for w in (17, 18, 19)] == [1, 2, 3], \
        "written in x,y,z axis order"
    block = np.array(ccp4.grid, copy=False) > 0
    start = [ccp4.header_i32(w) for w in (5, 6, 7)]
    sampling = [ccp4.header_i32(w) for w in (8, 9, 10)]
    parent = np.zeros(sampling, dtype=bool)
    axes = [(np.arange(s, s + n) % m)
            for s, n, m in zip(start, block.shape, sampling)]
    parent[np.ix_(*axes)] = block
    return parent


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
        masks = sorted(job.glob("mask_*.map"))
        assert [m.name for m in masks] == ["mask_140_339.map", "mask_340_485.map"]
        occupied = [_voxels_in_crystal_frame(m) for m in masks]
        assert int(np.count_nonzero(occupied[0] & occupied[1])) == 0


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
