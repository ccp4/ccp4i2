import xml.etree.ElementTree as ET
from os import environ
from pathlib import Path
from tempfile import NamedTemporaryFile
import gemmi
from .utils import demoData, hasLongLigandName, i2run


def test_substitute_ligand():
    args = ["SubstituteLigand"]
    args += ["--XYZIN", demoData("mdm2", "4hg7.cif")]
    args += ["--UNMERGEDFILES", "file=" + demoData("mdm2", "mdm2_unmerged.mtz")]
    args += ["--SMILESIN", '"CC(C)OC1=C(C=CC(=C1)OC)C2=NC(C(N2C(=O)N3CCNC(=O)C3)C4=CC=C(C=C4)Cl)C5=CC=C(C=C5)C"']
    args += ["--OBSAS", "UNMERGED"]
    args += ["--LIGANDAS", "SMILES"]
    args += ["--PIPELINE", "DIMPLE"]
    with i2run(args) as job:
        doc = gemmi.cif.read(str(job / "DICTOUT.cif"))
        gemmi.make_chemcomp_from_block(doc[-1])
        for name in (
            "DIFFPHIOUT",
            "F_SIGF_OUT_asFMEAN",
            "F_SIGF_OUT",
            "FREERFLAG_OUT",
        ):
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        gemmi.read_structure(str(job / "selected_atoms.cif"))
        gemmi.read_structure(str(job / "XYZOUT.pdb"))
        xml = ET.parse(job / "program.xml")
        rworks = [float(e.text) for e in xml.iter("r_factor")]
        rfrees = [float(e.text) for e in xml.iter("r_free")]
        assert rworks[-1] < 0.23
        assert rfrees[-1] < 0.25


# TODO: Get task working with mmCIF
# def test_8xfm(cif8xfm, mtz8xfm):
#     structure = gemmi.read_structure(cif8xfm)
#     for model in structure:
#         for chain in model:
#             for i, residue in reversed(list(enumerate(chain))):
#                 if residue.name == "A1LU6":
#                     del chain[i]
#     with NamedTemporaryFile(suffix=".cif", delete=False) as temp:
#         xyzin = temp.name
#     structure.make_mmcif_document().write_file(xyzin)
#     args = ["SubstituteLigand"]
#     args += ["--XYZIN", xyzin]
#     args += ["--F_SIGF_IN", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
#     args += ["--FREERFLAG_IN", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FREE]"]
#     args += ["--DICTIN", str(Path(environ["CLIBD_MON"], "a", "A1LU6.cif"))]
#     args += ["--OBSAS", "MERGED"]
#     args += ["--LIGANDAS", "DICT"]
#     args += ["--PIPELINE", "PHASER_RNP"]
#     try:
#         with i2run(args) as job:
#             assert hasLongLigandName(job / "XYZOUT.cif")
#     finally:
#         Path(xyzin).unlink(missing_ok=True)
