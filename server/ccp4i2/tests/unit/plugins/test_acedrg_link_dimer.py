"""AceDRG link mode builds a regularised linked pair; the wrapper must find it.

This is the picture a user can actually look at to judge whether the link
they asked for is chemically sensible -- JLigand drives libcheck and refmac
to produce the same thing. AceDRG does it internally (via servalcat), so all
this wrapper has to do is point at the result.

It had been pointing at "UNL_for_link", a name AceDRG stopped using, so both
outputs resolved to nothing and the job's Moorhen view had nothing to show.
Hence: find by suffix, never by the full name.
"""
import pytest

from ccp4i2.wrappers.AcedrgLink.script.AcedrgLink import _find_linked_dimer


@pytest.fixture
def work(tmp_path):
    (tmp_path / "LYS-PLP_TMP").mkdir()
    return tmp_path


def test_finds_the_pair_where_acedrg_leaves_it(work):
    # Coordinates in the work directory, dictionary one level down: not a
    # tidy arrangement, but AceDRG's, and the wrapper has to cope with it.
    (work / "LIG_for_link.pdb").write_text("END\n")
    (work / "LYS-PLP_TMP" / "LIG_for_link.cif").write_text("data_comp_list\n")
    pdb, cif = _find_linked_dimer(work, "LYS-PLP")
    assert pdb == work / "LIG_for_link.pdb"
    assert cif == work / "LYS-PLP_TMP" / "LIG_for_link.cif"


def test_the_prefix_is_not_hardcoded(work):
    # The old name. AceDRG renamed this once already; that is the whole
    # reason the lookup is by suffix.
    (work / "UNL_for_link.pdb").write_text("END\n")
    (work / "LYS-PLP_TMP" / "UNL_for_link.cif").write_text("data_comp_list\n")
    pdb, cif = _find_linked_dimer(work, "LYS-PLP")
    assert pdb is not None and pdb.name == "UNL_for_link.pdb"
    assert cif is not None and cif.name == "UNL_for_link.cif"


def test_the_link_dictionary_is_not_mistaken_for_the_pair(work):
    # <LINK_ID>_link.cif is CIF_OUT, a different output. It ends in
    # "_link.cif" but not "_for_link.cif", and must not be picked up here.
    (work / "LYS-PLP_link.cif").write_text("data_comp_list\n")
    pdb, cif = _find_linked_dimer(work, "LYS-PLP")
    assert pdb is None
    assert cif is None


def test_nothing_found_is_none_not_a_guess(work):
    assert _find_linked_dimer(work, "LYS-PLP") == (None, None)


def test_a_missing_tmp_directory_is_survivable(tmp_path):
    (tmp_path / "LIG_for_link.pdb").write_text("END\n")
    pdb, cif = _find_linked_dimer(tmp_path, "LYS-PLP")
    assert pdb is not None
    assert cif is None
