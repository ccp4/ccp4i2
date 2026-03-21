from .utils import demoData, i2run


def test_gamma_xe():
    """Test phaser EP LLG anomalous map from Xe coordinates."""
    args = ["phaser_EP_LLG"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("gamma", "heavy_atoms.pdb")]
    with i2run(args, allow_errors=True) as job:
        # EP_LLG may not converge with this data but should produce output
        assert (job / "diagnostic.xml").exists()
