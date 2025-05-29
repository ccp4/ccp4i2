from pytest import fixture
from .urls import pdbe_fasta, pdbe_mmcif, redo_cif, redo_mtz, rcsb_mmcif, rcsb_pdb
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


@fixture(scope="session")
def cif8ola():
    with download(redo_cif("8ola")) as path:
        yield path


@fixture(scope="session")
def mtz8ola():
    with download(pdbe_mmcif("8ola")) as path:
        yield path


@fixture(scope="session")
def cif8olf():
    with download(redo_cif("8olf")) as path:
        yield path


@fixture(scope="session")
def mtz8olf():
    with download(pdbe_mmcif("8olf")) as path:
        yield path


@fixture(scope="session")
def cif8rk1():
    with download(rcsb_mmcif("8rk1")) as path:
        yield path


@fixture(scope="session")
def pdb8rk1():
    with download(rcsb_pdb("8rk1")) as path:
        yield path


@fixture(scope="session")
def mtz8c4y():
    with download(redo_mtz("8c4y")) as path:
        yield path


@fixture(scope="session")
def cif7beq():
    with download(rcsb_mmcif("7beq")) as path:
        yield path


@fixture(scope="session")
def mtz7beq():
    with download(redo_mtz("7beq")) as path:
        yield path