"""A search model covers the AU kinds its chains align to (docs/multi-component-mr.md).

On the beta-lactamase / BLIP demo: each single-chain model covers its own
kind and not the other (the wrong pairing scores negative with chance
matches at ~20% identity), and a file holding both chains covers both.
"""
from pathlib import Path

import gemmi
import pytest

from ccp4i2.lib.utils.formats.model_coverage import (
    chain_assignments, covered_kinds, describe_coverage,
)

DEMO = Path(__file__).resolve().parents[3] / "demo_data" / "beta_blip"


def _seq(name):
    return "".join(l.strip() for l in (DEMO / f"{name}.seq").read_text().splitlines()
                   if not l.startswith(">"))


KINDS = [{"name": "BETA", "sequence": _seq("beta"), "copies": 1, "polymerType": "PROTEIN"},
         {"name": "BLIP", "sequence": _seq("blip"), "copies": 1, "polymerType": "PROTEIN"}]


def _complex(path):
    """beta and BLIP as chains A and B of one file: a complex template."""
    out = gemmi.Structure()
    model = gemmi.Model("1")
    for name, chain_id in (("beta", "A"), ("blip", "B")):
        source = gemmi.read_structure(str(DEMO / f"{name}.pdb"))
        chain = source[0][0]
        chain.name = chain_id
        model.add_chain(chain)
    out.add_model(model)
    out.setup_entities()
    out.write_pdb(str(path))
    return path


def test_a_single_chain_model_covers_its_own_kind_only():
    assert covered_kinds(DEMO / "beta.pdb", KINDS) == {"BETA": 100.0}
    assert covered_kinds(DEMO / "blip.pdb", KINDS) == {"BLIP": 100.0}


def test_a_complex_template_covers_every_kind_it_holds(tmp_path):
    path = _complex(tmp_path / "beta_blip.pdb")
    assert covered_kinds(path, KINDS) == {"BETA": 100.0, "BLIP": 100.0}
    assert [(c, k) for c, k, _ in chain_assignments(path, KINDS)] == [("A", "BETA"), ("B", "BLIP")]


def test_a_chain_matching_nothing_covers_nothing():
    kinds = [{"name": "OTHER", "sequence": "MKVLAAGIVALLLAAGCSSAKEE" * 8,
              "copies": 1, "polymerType": "PROTEIN"}]
    assert covered_kinds(DEMO / "beta.pdb", kinds) == {}


def test_description_reads_as_a_sentence():
    assert describe_coverage("search model 1", {"CDK4": 100.0, "CyclinD1": 92.0}) == \
        "search model 1 covers CDK4 (100%), CyclinD1 (92%)"
    assert describe_coverage("search model 2", {}) == "search model 2 covers none of them"
