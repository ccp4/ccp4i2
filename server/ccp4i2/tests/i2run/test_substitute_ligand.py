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
        _check_aimless_pipe_performance(job / "job_1")


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
        _check_aimless_pipe_performance(job / "job_1")


def _check_aimless_pipe_performance(aimless_pipe_dir):
    """Verify aimless_pipe wrote non-zero PERFORMANCE KPIs to its params.xml."""
    params = aimless_pipe_dir / "params.xml"
    assert params.exists(), f"params.xml not found at {params}"
    tree = ET.parse(params)
    perf = tree.find('.//outputData/PERFORMANCE')
    assert perf is not None, "No PERFORMANCE element in aimless_pipe params.xml"
    high_res = perf.find('highResLimit')
    assert high_res is not None and high_res.text, "highResLimit missing from PERFORMANCE"
    assert float(high_res.text) > 0, f"highResLimit is zero: {high_res.text}"
    r_meas = perf.find('rMeas')
    assert r_meas is not None and r_meas.text, "rMeas missing from PERFORMANCE"
    assert float(r_meas.text) > 0, f"rMeas is zero: {r_meas.text}"