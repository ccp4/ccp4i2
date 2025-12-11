"Tests for the parse module"

import gemmi
from pytest import mark
from ...lib import links
from ...lib.parse import parse
from ...lib.web import download


@mark.parametrize(
    ("link", "cls"),
    [
        (links.pdbe_eds_map, gemmi.Ccp4Map),
        (links.pdbe_fasta, list),
        (links.pdbe_mmcif_updated, gemmi.Structure),
        (links.pdbe_mmcif, gemmi.Structure),
        (links.pdbe_pdb, gemmi.Structure),
        (links.pdbe_pdb_gz, gemmi.Structure),
        (links.pdbe_sfcif, gemmi.ReflnBlocks),
        (links.rcsb_fasta, list),
        (links.rcsb_mmcif_gz, gemmi.Structure),
        (links.rcsb_mmcif, gemmi.Structure),
        (links.rcsb_pdb_gz, gemmi.Structure),
        (links.rcsb_pdb, gemmi.Structure),
        (links.rcsb_sfcif_gz, gemmi.ReflnBlocks),
        (links.rcsb_sfcif, gemmi.ReflnBlocks),
        (links.redo_cif, gemmi.Structure),
        (links.redo_mtz, gemmi.Mtz),
        (links.redo_pdb, gemmi.Structure),
    ],
)
def test_parse(link, cls):
    with download(link("1o6a")) as path:
        parsed = parse(path)
    message = f"{path} was not parsed as {cls}."
    if isinstance(parsed, dict):
        for format_, error in parsed.items():
            message += f"\n{format_}: {error}"
    assert isinstance(parsed, cls), message
