"""Which kind of chain a search model covers, by alignment.

Opus, driving a two-lobe kinase through the facade (2026-10-09), split a
566-residue model into lobes of 214 and 276 residues at 100% identity and
was told each "covers none" of the AU contents: the global alignment score
charged the unaligned rest of the target as gaps (-140 and -22 against
+402 for the whole model). The score is now over the aligned span of the
target, with BLOSUM62, so a fragment of the right chain counts and a wrong
pairing still does not.
"""
from pathlib import Path

import gemmi
import pytest

from ccp4i2.lib.utils.formats.model_coverage import MIN_MATCHED, chain_assignments, covered_kinds

DEMO = Path(__file__).resolve().parents[3] / "demo_data" / "beta_blip"


def _seq(name):
    return "".join(l.strip() for l in (DEMO / f"{name}.seq").read_text().splitlines()
                   if not l.startswith(">"))


KINDS = [{"name": "BETA", "sequence": _seq("beta"), "polymerType": "PROTEIN"},
         {"name": "BLIP", "sequence": _seq("blip"), "polymerType": "PROTEIN"}]


def _fragment(path, name, first, last):
    """Residues first..last (1-based, by position) of a demo chain."""
    structure = gemmi.read_structure(str(DEMO / f"{name}.pdb"))
    chain = structure[0][0]
    for index in range(len(chain) - 1, -1, -1):
        if not first - 1 <= index <= last - 1:
            del chain[index]
    structure.setup_entities()
    structure.write_pdb(str(path))
    return path


def _kinds(assignments):
    """(kind, identity) per chain: the demo files' chains have no name."""
    return [(name, identity) for _chain, name, identity in assignments]


def test_whole_chains_are_assigned_to_their_kind():
    assert _kinds(chain_assignments(DEMO / "beta.pdb", KINDS)) == [("BETA", 100.0)]
    assert _kinds(chain_assignments(DEMO / "blip.pdb", KINDS)) == [("BLIP", 100.0)]


@pytest.mark.parametrize("first, last", [(1, 120), (120, 263), (1, 60)])
def test_a_fragment_of_the_right_chain_covers_its_kind(tmp_path, first, last):
    # the lobes: a fragment at 100% identity is not "none of them"
    covered = covered_kinds(_fragment(tmp_path / "lobe.pdb", "beta", first, last), KINDS)
    assert covered == {"BETA": 100.0}


def test_a_wrong_pairing_is_not_coverage(tmp_path):
    # BLIP against beta-lactamase's sequence alone: chance matches at ~25%
    assert _kinds(chain_assignments(DEMO / "blip.pdb", KINDS[:1])) == [(None, 0.0)]
    assert _kinds(chain_assignments(DEMO / "beta.pdb", KINDS[1:])) == [(None, 0.0)]


def test_too_few_matched_residues_is_not_coverage(tmp_path):
    short = _fragment(tmp_path / "short.pdb", "beta", 1, MIN_MATCHED - 5)
    assert covered_kinds(short, KINDS) == {}
