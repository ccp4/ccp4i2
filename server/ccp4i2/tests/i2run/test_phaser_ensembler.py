from .utils import demoData, i2run


def test_mdm2_two_models():
    """Test phaser.ensembler with two different MDM2 structures."""
    args = ["phaser_ensembler"]
    args += ["--XYZIN_LIST", demoData("mdm2", "4qo4.pdb")]
    args += ["--XYZIN_LIST", demoData("mdm2", "4hg7.pdb")]
    with i2run(args) as job:
        pdb_files = list(job.glob("*.pdb"))
        assert len(pdb_files) > 0, f"No PDB output: {list(job.iterdir())}"
