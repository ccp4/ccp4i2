import gemmi

from .utils import i2run


def test_comit(mtz8xfm):
    args = ["comit"]
    args += ["--F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--F_PHI_IN", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FWT,PHWT]"]
    with i2run(args) as job:
        gemmi.read_mtz_file(str(job / "F_PHI_OUT.mtz"))


def test_comit_from_intensities(mtz8xfm, tmp_path):
    """comit reads amplitudes: given intensities (as from aimless) it died in
    clipper with "Missing column ... F_SIGF_F". The input is now converted."""
    import numpy as np
    mtz = gemmi.read_mtz_file(str(mtz8xfm))
    data = np.array(mtz, copy=False)
    labels = [c.label for c in mtz.columns]
    f, sigf = data[:, labels.index("FP")], data[:, labels.index("SIGFP")]
    out = gemmi.Mtz(with_base=True)
    out.spacegroup, out.cell = mtz.spacegroup, mtz.cell
    out.add_dataset("x")
    out.add_column("I", "J")
    out.add_column("SIGI", "Q")
    hkl = data[:, :3]
    out.set_data(np.column_stack([hkl, f * f, 2 * f * sigf]).astype(np.float32))
    intensities = tmp_path / "intensities.mtz"
    out.write_to_file(str(intensities))

    args = ["comit"]
    args += ["--F_SIGF", f"fullPath={intensities}", "columnLabels=/*/*/[I,SIGI]"]
    args += ["--F_PHI_IN", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FWT,PHWT]"]
    with i2run(args) as job:
        gemmi.read_mtz_file(str(job / "F_PHI_OUT.mtz"))
