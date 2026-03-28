# Copyright (C) 2025-2026 Newcastle University
# Copyright (C) 2025-2026 University of York
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
import xml.etree.ElementTree as ET
import gemmi
import pytest
from .utils import demoData, i2run


# TODO: Test long ligand names (e.g. 8xfm)


# NOTE: All phaser tests run FIRST (order=1) to avoid RDKit pickle contamination
# RDKit (imported by acedrg tests) modifies pickle module's dispatch table,
# causing phaser's pickle.dump() to fail. Running phaser tests before acedrg
# ensures pickle module is clean when phaser needs it.

@pytest.mark.order("first")
def test_substitute_ligand_no_ligand():
    """Test SubstituteLigand with LIGANDAS=NONE (no ligand fitting, just refinement)."""
    args = ["SubstituteLigand"]
    args += ["--XYZIN", demoData("mdm2", "4hg7.cif")]
    args += ["--UNMERGEDFILES", "file=" + demoData("mdm2", "mdm2_unmerged.mtz")]
    args += ["--LIGANDAS", "NONE"]
    args += ["--PIPELINE", "DIMPLE"]
    with i2run(args) as job:
        # Check expected MTZ outputs
        for name in (
            "DIFFPHIOUT",
            "F_SIGF_OUT",
            "FREERFLAG_OUT",
            "ANOMFPHIOUT",
        ):
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        gemmi.read_structure(str(job / "selected_atoms.pdb"), format=gemmi.CoorFormat.Mmcif)
        gemmi.read_structure(str(job / "XYZOUT.pdb"))
        xml = ET.parse(job / "program.xml")
        rworks = [float(e.text) for e in xml.iter("r_factor")]
        rfrees = [float(e.text) for e in xml.iter("r_free")]
        assert rworks[-1] < 0.23
        assert rfrees[-1] < 0.25
        # Check aimless_pipe sub-job wrote PERFORMANCE indicators to params.xml
        _check_aimless_pipe_performance(job)


@pytest.mark.order("first")
def test_substitute_ligand_with_smiles():
    """Test SubstituteLigand with SMILES input (full pipeline with ligand fitting)."""
    args = ["SubstituteLigand"]
    args += ["--XYZIN", demoData("mdm2", "4hg7.cif")]
    args += ["--UNMERGEDFILES", "file=" + demoData("mdm2", "mdm2_unmerged.mtz")]
    args += ["--SMILESIN", "CC(C)OC1=C(C=CC(=C1)OC)C2=NC(C(N2C(=O)N3CCNC(=O)C3)C4=CC=C(C=C4)Cl)C5=CC=C(C=C5)C"]
    args += ["--PIPELINE", "DIMPLE"]
    with i2run(args) as job:
        # Check ligand dictionary was generated
        doc = gemmi.cif.read(str(job / "DICTOUT.cif"))
        gemmi.make_chemcomp_from_block(doc[-1])
        # Check expected MTZ outputs
        for name in (
            "DIFFPHIOUT",
            "F_SIGF_OUT",
            "FREERFLAG_OUT",
        ):
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        gemmi.read_structure(str(job / "selected_atoms.pdb"), format=gemmi.CoorFormat.Mmcif)
        gemmi.read_structure(str(job / "XYZOUT.pdb"))
        xml = ET.parse(job / "program.xml")
        rworks = [float(e.text) for e in xml.iter("r_factor")]
        rfrees = [float(e.text) for e in xml.iter("r_free")]
        assert rworks[-1] < 0.23
        assert rfrees[-1] < 0.25
        # Check aimless_pipe sub-job wrote PERFORMANCE indicators to params.xml
        _check_aimless_pipe_performance(job)


def _find_aimless_pipe_dir(job_dir):
    """Find the aimless_pipe sub-job directory by checking params.xml pluginName."""
    for sub in sorted(job_dir.iterdir()):
        params = sub / "params.xml"
        if params.exists():
            tree = ET.parse(params)
            plugin = tree.find('.//pluginName')
            if plugin is not None and plugin.text == 'aimless_pipe':
                return sub
    return None


def _check_aimless_pipe_performance(job_dir):
    """Verify aimless_pipe wrote non-zero PERFORMANCE KPIs to its params.xml."""
    aimless_pipe_dir = _find_aimless_pipe_dir(job_dir)
    assert aimless_pipe_dir is not None, f"No aimless_pipe sub-job found in {job_dir}"
    params = aimless_pipe_dir / "params.xml"
    tree = ET.parse(params)
    perf = tree.find('.//outputData/PERFORMANCE')
    assert perf is not None, "No PERFORMANCE element in aimless_pipe params.xml"
    high_res = perf.find('highResLimit')
    assert high_res is not None and high_res.text, "highResLimit missing from PERFORMANCE"
    assert float(high_res.text) > 0, f"highResLimit is zero: {high_res.text}"
    r_meas = perf.find('rMeas')
    assert r_meas is not None and r_meas.text, "rMeas missing from PERFORMANCE"
    assert float(r_meas.text) > 0, f"rMeas is zero: {r_meas.text}"