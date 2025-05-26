from pytest import fixture
from .urls import pdbe_fasta, redo_cif, redo_mtz, rcsb_mmcif, rcsb_sfcif
from .utils import download


@fixture(scope="session")
def cif8xfm():
    with download(redo_cif("8xfm")) as path:
        yield path


@fixture(scope="session")
def mtz8xfm():
    with download(redo_mtz("8xfm")) as path:
        yield path


@fixture(scope="session")
def seq8xfm():
    with download(pdbe_fasta("8xfm")) as path:
        yield path


# For test_servalcat.py
@fixture(scope="session")
def cif8ola():
    with download(rcsb_mmcif("8ola")) as path:
        yield path


@fixture(scope="session")
def mtz8ola():
    with download(redo_mtz("8ola")) as path:
        yield path

