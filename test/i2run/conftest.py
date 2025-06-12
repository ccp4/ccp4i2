from pytest import fixture
from .urls import (
    pdbe_fasta,
    rcsb_ligand_cif,
    rcsb_ligand_sdf,
    rcsb_mmcif,
    redo_cif,
    redo_mtz,
)
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
def sdfA1LU6():
    with download(rcsb_ligand_sdf("A1LU6")) as path:
        yield path


@fixture(scope="session")
def cifA1LU6():
    with download(rcsb_ligand_cif("A1LU6")) as path:
        yield path


@fixture(scope="session")
def sdfOEX():
    with download(rcsb_ligand_sdf("OEX")) as path:
        yield path


@fixture(scope="session")
def cifOEX():
    with download(rcsb_ligand_cif("OEX")) as path:
        yield path


@fixture(scope="session")
def cif4ub6():
    with download(rcsb_mmcif("4ub6")) as path:
        yield path


@fixture(scope="session")
def cifAF3():
    with download(rcsb_ligand_cif("AF3")) as path:
        yield path


@fixture(scope="session")
def cif4dl8():
    with download(rcsb_mmcif("4dl8")) as path:
        yield path
