from .utils import demoData, i2run


def test_mdm2_two_models():
    """Test phaser.ensembler with two different MDM2 structures."""
    args = ["phaser_ensembler"]
    args += ["--XYZIN_LIST", demoData("mdm2", "4qo4.pdb")]
    args += ["--XYZIN_LIST", demoData("mdm2", "4hg7.pdb")]
    with i2run(args) as job:
        pdb_files = list(job.glob("*.pdb"))
        assert len(pdb_files) > 0, f"No PDB output: {list(job.iterdir())}"


def test_mdm2_models_with_selections():
    """A model with an atom selection is cut to it before ensembling. The
    selection added a residue-name list the selection language no longer
    parses, so the selected file was never written and Phaser failed."""
    args = ["phaser_ensembler"]
    args += ["--XYZIN_LIST", f"fullPath={demoData('mdm2', '4qo4.pdb')}", "selection/text=A/"]
    args += ["--XYZIN_LIST", f"fullPath={demoData('mdm2', '4hg7.pdb')}", "selection/text=A/"]
    with i2run(args) as job:
        assert (job / "selected_0.pdb").exists() and (job / "selected_1.pdb").exists(), \
            list(job.iterdir())
        assert list(job.glob("*.pdb")), list(job.iterdir())
