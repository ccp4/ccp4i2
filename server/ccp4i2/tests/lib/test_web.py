from pytest import mark
from ...lib.web import (
    download,
    pdbe_fasta,
    pdbe_mmcif,
    pdbe_pdb,
    pdbe_sfcif,
    redo_cif,
    redo_mtz,
    redo_pdb,
)


@mark.parametrize(
    "url",
    [
        pdbe_fasta,
        pdbe_mmcif,
        pdbe_pdb,
        pdbe_sfcif,
        redo_cif,
        redo_mtz,
        redo_pdb,
    ],
)
def test_download(url):
    with download(url("1o6a")) as path:
        assert path.exists()
        assert path.stat().st_size > 0
    assert not path.exists()
