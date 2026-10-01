from .utils import hasLongLigandName, i2run


def test_8xfm(cif8xfm, mtz8xfm):
    args = ["coot_rsr_morph"]
    args += ["--XYZIN", cif8xfm]
    args += ["--FPHIIN", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FWT,PHWT]"]
    with i2run(args) as job:
        assert hasLongLigandName(job / "XYZOUT.cif")
        # How far it moved the model is measured and reported (it said only
        # "finished"); a deposited model in its own map moves little.
        import xml.etree.ElementTree as ET
        moved = ET.parse(job / "program.xml").find(".//Shifts")
        assert int(moved.get("atoms")) > 1000 and float(moved.get("rms")) < 1.0, moved.attrib
