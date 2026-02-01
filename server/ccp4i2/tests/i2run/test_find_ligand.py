from .utils import hasLongLigandName, i2run
from os import environ
from pathlib import Path
from tempfile import NamedTemporaryFile
import gemmi
import xml.etree.ElementTree as ET


@mark.skip(reason="Skipping temporarily until task is in registry")
def test_8xfm(cif8xfm, mtz8xfm):
    structure = gemmi.read_structure(cif8xfm)
    for chain in structure[0]:
        for i, residue in reversed(list(enumerate(chain))):
            if residue.name == "A1LU6":
                del chain[i]
    with NamedTemporaryFile(suffix=".cif", delete=False) as temp:
        xyzin = temp.name
    doc = structure.make_mmcif_document()
    doc.write_file(xyzin)
    args = ["coot_find_ligand"]
    args += ["--XYZIN", xyzin]
    args += ["--FPHI", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FWT,PHWT]"]
    args += ["--COMPID", "A1LU6"]
    args += ["--DICT", str(Path(environ["CLIBD_MON"], "a", "A1LU6.cif"))]
    try:
        with i2run(args) as job:
            tree = ET.parse(job / "program.xml")
            ligands = tree.findall(".//AddedLigand")
            assert len(ligands) > 0
            assert hasLongLigandName(job / "XYZOUT.cif")
    finally:
        Path(xyzin).unlink(missing_ok=True)
