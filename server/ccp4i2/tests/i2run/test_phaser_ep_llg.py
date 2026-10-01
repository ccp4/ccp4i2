import gemmi

from .utils import demoData, i2run


def test_gamma_xe():
    """Anomalous LLG map phased by the protein model, which shows the Xe.

    The map is written by Phaser as FLLG_AX,PHLLG_AX. Left under those names
    it failed the output check for map coefficients (F,PHI), and every run of
    the task ended in error. This test allowed errors and so never saw it; it
    also gave only the substructure, from which Phaser has no atoms to phase.
    """
    args = ["phaser_EP_LLG"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--PARTIALMODELORMAP", "MODEL"]
    args += ["--XYZIN_PARTIAL", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        mtz = gemmi.read_mtz_file(str(job / "LLGMAPOUT_1.mtz"))
        assert [c.label for c in mtz.columns][3:] == ["F", "PHI"]
        params = (job / "params.xml").read_text()
        # One hand, and no empty "sites" file offered as a structure.
        assert "original hand" not in params
        assert "PHASER.1.pdb" not in params
