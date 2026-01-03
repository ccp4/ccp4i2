import gemmi
import pytest
from . import urls
from .utils import download


def read_fasta(path):
    with open(path, encoding="utf-8") as f:
        gemmi.read_pir_or_fasta(f.read())


def read_sdf(path):
    """Lazy-load RDKit to avoid pickle contamination during pytest collection."""
    import rdkit.Chem
    return rdkit.Chem.SDMolSupplier(path)


@pytest.mark.parametrize(
    "url_func,code,read_func",
    [
        (urls.pdbe_fasta, "8xfm", read_fasta),
        (urls.pdbe_mmcif, "1gyu", gemmi.read_structure),
        (urls.pdbe_pdb, "6ndn", gemmi.read_structure),
        (urls.pdbe_sfcif, "2ceu", gemmi.cif.read),
        (urls.rcsb_ligand_cif, "A1LU6", gemmi.cif.read),
        (urls.rcsb_ligand_sdf, "A1LU6", read_sdf),
        (urls.rcsb_mmcif, "4dl8", gemmi.read_structure),
        (urls.rcsb_pdb, "8rk1", gemmi.read_structure),
        (urls.redo_cif, "8xfm", gemmi.read_structure),
        (urls.redo_mtz, "8xfm", gemmi.read_mtz_file),
    ],
)
def test(url_func, code, read_func):
    with download(url_func(code)) as path:
        read_func(path)
