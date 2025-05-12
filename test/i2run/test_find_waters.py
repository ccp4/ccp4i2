from pathlib import Path
from tempfile import NamedTemporaryFile
import xml.etree.ElementTree as ET
import gemmi
from .utils import hasLongLigandName, i2run


def test_8xfm(cif8xfm, mtz8xfm):
    structure = gemmi.read_structure(cif8xfm)
    for chain in structure[0]:
        for i, residue in reversed(list(enumerate(chain))):
            if residue.name == "HOH":
                del chain[i]
    with NamedTemporaryFile(suffix=".cif", delete=False) as temp:
        xyzin = temp.name
    doc = structure.make_mmcif_document()
    doc.write_file(xyzin)
    assert not has_water(xyzin)
    args = ["coot_find_waters"]
    args += ["--XYZIN", xyzin]
    args += ["--FPHIIN", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FWT,PHWT]"]
    try:
        with i2run(args) as job:
            assert hasLongLigandName(job / "XYZOUT.cif")
            assert has_water(job / "XYZOUT.cif")
    finally:
        Path(xyzin).unlink(missing_ok=True)


def has_water(path):
    structure = gemmi.read_structure(str(path))
    for chain in structure[0]:
        for residue in chain:
            if residue.name == "HOH":
                return True
    return False
