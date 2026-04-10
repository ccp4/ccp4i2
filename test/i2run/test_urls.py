import gemmi
import pytest
import rdkit.Chem
from . import urls
from .utils import download


def read_fasta(path):
    with open(path, encoding="utf-8") as f:
        gemmi.read_pir_or_fasta(f.read())


@pytest.mark.parametrize(
    "url_func,code,read_func",
    [
        (urls.pdbe_fasta, "8xfm", read_fasta),
        (urls.pdbe_mmcif, "1gyu", gemmi.read_structure),
        (urls.pdbe_pdb, "6ndn", gemmi.read_structure),
        (urls.pdbe_sfcif, "2ceu", gemmi.cif.read),
        (urls.rcsb_ligand_cif, "A1LU6", gemmi.cif.read),
        (urls.rcsb_ligand_sdf, "A1LU6", rdkit.Chem.SDMolSupplier),
        (urls.rcsb_mmcif, "4dl8", gemmi.read_structure),
        (urls.rcsb_pdb, "8rk1", gemmi.read_structure),
        (urls.redo_cif, "8xfm", gemmi.read_structure),
        (urls.redo_mtz, "8xfm", gemmi.read_mtz_file),
    ],
)
def test(url_func, code, read_func):
    with download(url_func(code)) as path:
        read_func(path)
