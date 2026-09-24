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
        # The input is mmCIF, so the extracted selection stays mmCIF and is
        # named accordingly. This used to be selected_atoms.pdb holding mmCIF,
        # and the test asserted exactly that by forcing the reader past the
        # extension --- which is how the format lie survived. Read it by
        # content, and require the name to match.
        gemmi.read_structure(str(job / "selected_atoms.cif"),
                             format=gemmi.CoorFormat.Detect)
        assert not (job / "selected_atoms.pdb").exists()
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
        # The input is mmCIF, so the extracted selection stays mmCIF and is
        # named accordingly. This used to be selected_atoms.pdb holding mmCIF,
        # and the test asserted exactly that by forcing the reader past the
        # extension --- which is how the format lie survived. Read it by
        # content, and require the name to match.
        gemmi.read_structure(str(job / "selected_atoms.cif"),
                             format=gemmi.CoorFormat.Detect)
        assert not (job / "selected_atoms.pdb").exists()
        gemmi.read_structure(str(job / "XYZOUT.pdb"))
        xml = ET.parse(job / "program.xml")
        rworks = [float(e.text) for e in xml.iter("r_factor")]
        rfrees = [float(e.text) for e in xml.iter("r_free")]
        assert rworks[-1] < 0.23
        assert rfrees[-1] < 0.25
        # Check aimless_pipe sub-job wrote PERFORMANCE indicators to params.xml
        _check_aimless_pipe_performance(job)


@pytest.mark.order("first")
def test_substitute_ligand_phaser_rnp():
    """The Phaser route: rigid-body refinement of the whole model by
    phaser_rnp_pipeline_phil, then servalcat. Until this test the Phaser
    route had no coverage at all."""
    args = ["SubstituteLigand"]
    args += ["--XYZIN", demoData("mdm2", "4hg7.cif")]
    args += ["--UNMERGEDFILES", "file=" + demoData("mdm2", "mdm2_unmerged.mtz")]
    args += ["--LIGANDAS", "NONE"]
    args += ["--PIPELINE", "PHASER_RNP"]
    with i2run(args) as job:
        for name in ("DIFFPHIOUT", "F_SIGF_OUT", "FREERFLAG_OUT"):
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        gemmi.read_structure(str(job / "XYZOUT.pdb"))
        xml = ET.parse(job / "program.xml")
        # Rigid-body refinement (MR_RNP) gives refined solutions, no verdict
        record = xml.find(".//PhaserMrResults")
        assert record is not None
        llgs = [float(e.text) for e in record.findall("Solutions/Solution/LLG")]
        assert llgs and max(llgs) > 0
        assert xml.find(".//POINTLESS") is not None
        # The record's r_factor elements are the RNP pipeline's refmac (ten
        # jelly-body cycles from the rigid-body model); the final figures
        # are servalcat's, on its performance indicator. This route refines
        # less than DIMPLE's (which reaches 0.23): the refmac step is the
        # classic pipeline's, unchanged.
        r_work, r_free = _servalcat_r_factors(job)
        assert r_work < 0.30
        assert r_free < 0.32
        _check_aimless_pipe_performance(job)


def _find_subjob_dir(job_dir, plugin_name):
    """Find a sub-job directory by its params.xml pluginName."""
    for sub in sorted(job_dir.iterdir()):
        params = sub / "params.xml"
        if params.exists():
            tree = ET.parse(params)
            plugin = tree.find('.//pluginName')
            if plugin is not None and plugin.text == plugin_name:
                return sub
    return None


def _find_aimless_pipe_dir(job_dir):
    return _find_subjob_dir(job_dir, "aimless_pipe")


def _servalcat_r_factors(job_dir):
    sub = _find_subjob_dir(job_dir, "servalcat_pipe")
    assert sub is not None, f"No servalcat_pipe sub-job found in {job_dir}"
    perf = ET.parse(sub / "params.xml").find(".//outputData/PERFORMANCEINDICATOR")
    assert perf is not None, "No PERFORMANCEINDICATOR in servalcat_pipe params.xml"
    return float(perf.findtext("R1Factor")), float(perf.findtext("R1Free"))


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

# ---------------------------------------------------------------------------
# Merged data (OBSAS=MERGED)
#
# Until these, every SubstituteLigand test supplied UNMERGEDFILES, so the
# merged route -- the one a fragment campaign uses, where aimless never runs
# and the free-R set is whatever the user hands in -- had no end-to-end
# coverage at all. That is where a run of real failures turned up: the ligand
# written to the wrong field, the reference's own ligand left in the model,
# and a free-R set from another crystal refused by a cell check part-way
# through the job.
# ---------------------------------------------------------------------------


def _recelled_free_r(source_mtz, target_cell, destination):
    """A copy of a free-R set stamped with a different unit cell.

    Stands in for a campaign's shared free set: the same reflections and the
    same flags, but the cell of the crystal it was measured on rather than the
    one being refined. Built here rather than shipped as demo data so the
    difference is visible in the test, and adjustable.
    """
    import gemmi

    mtz = gemmi.read_mtz_file(str(source_mtz))
    mtz.cell = gemmi.UnitCell(*target_cell)
    for dataset in mtz.datasets:
        dataset.cell = gemmi.UnitCell(*target_cell)
    mtz.write_to_file(str(destination))
    return destination


def test_substitute_ligand_merged_data():
    """The merged route runs at all: no unmerged files, no aimless.

    OBSAS=MERGED sends F_SIGF_IN straight to refinement. This is the shape a
    fragment campaign uses, and what make_demo_campaign configures.
    """
    args = ["SubstituteLigand"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--OBSAS", "MERGED"]
    args += ["--F_SIGF_IN", demoData("gamma", "merged_intensities_native.mtz")]
    args += ["--FREERFLAG_IN", demoData("gamma", "freeR.mtz")]
    args += ["--LIGANDAS", "NONE"]
    args += ["--PIPELINE", "DIMPLE"]
    with i2run(args) as job:
        # No F_SIGF_OUT on this route: the observations were already merged
        # when they came in and pass through unchanged, so there is nothing to
        # re-export. FREERFLAG_OUT *is* written, because the supplied free-R
        # set is reconciled with the data before refinement.
        for name in ("FREERFLAG_OUT", "DIFFPHIOUT", "FPHIOUT"):
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        gemmi.read_structure(str(job / "XYZOUT.pdb"))


def test_merged_free_r_other_crystal(tmp_path):
    """A free-R set whose cell disagrees with the data still refines.

    This is the campaign case: one free set shared across a series of soaks,
    so its cell is the reference crystal's and differs from every member's by
    a percent or more. Every merge in the pipeline compares the two and would
    otherwise refuse -- the job used to stop inside i2Dimple with
    "Incompatible unit cells", part-way through, having already built the
    ligand.

    The pipeline reconciles the set up front instead: freerflag in COMPLETE
    mode joins by reflection index, so the existing flags keep the reflections
    they were assigned to, stamps the data's cell, and extends the set to the
    data's resolution. Re-stamping the cell alone would not do that last part.
    """
    from ccp4i2.core.CCP4XtalData import cells_are_compatible

    observations = demoData("gamma", "merged_intensities_native.mtz")
    shifted = _recelled_free_r(
        demoData("gamma", "freeR.mtz"),
        # ~2 A out on a, the scale of a real soak-to-soak drift and well past
        # Clipper's 1 A test.
        (36.15, 54.81, 68.00, 90.0, 90.0, 90.0),
        tmp_path / "free_from_another_crystal.mtz",
    )

    data_cell = gemmi.read_mtz_file(str(observations)).cell.parameters
    free_cell = gemmi.read_mtz_file(str(shifted)).cell.parameters
    assert not cells_are_compatible(data_cell, free_cell, 1.0)["validity"], (
        "this test is pointless unless the cells really do fail the strict check"
    )

    args = ["SubstituteLigand"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--OBSAS", "MERGED"]
    args += ["--F_SIGF_IN", observations]
    args += ["--FREERFLAG_IN", str(shifted)]
    args += ["--LIGANDAS", "NONE"]
    args += ["--PIPELINE", "DIMPLE"]
    with i2run(args) as job:
        gemmi.read_structure(str(job / "XYZOUT.pdb"))

        # The published free-R set is the reconciled one: it carries the
        # DATA's cell, not the one it came in with.
        reconciled = gemmi.read_mtz_file(str(job / "FREERFLAG_OUT.mtz"))
        assert cells_are_compatible(
            data_cell, reconciled.cell.parameters, 1.0
        )["validity"], (
            "FREERFLAG_OUT should carry the data's cell, so later jobs on this "
            "dataset can use it without repeating the reconciliation"
        )
