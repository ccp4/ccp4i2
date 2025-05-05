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
    with NamedTemporaryFile(suffix=".cif") as temp:
        doc = structure.make_mmcif_document()
        doc.write_file(temp.name)
        args = ["coot_find_waters"]
        args += ["--XYZIN", temp.name]
        args += ["--FPHIIN", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FWT,PHWT]"]
        with i2run(args) as job:
            assert hasLongLigandName(job / "XYZOUT.cif")
            xml = ET.parse(job / "program.xml")
            waters_found = int(xml.find("WatersFound").text)
            assert waters_found > 0
            waters_in_structure = 0
            structure = gemmi.read_structure(str(job / "XYZOUT.cif"))
            for chain in structure[0]:
                for residue in chain:
                    if residue.name == "HOH":
                        waters_in_structure += 1
            assert waters_in_structure == waters_found
