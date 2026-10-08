"""Which kinds of chain in the AU contents a search model covers.

Molecular replacement of a complex has three routes (docs/multi-component-mr.md):
one model holding every component, searched as one rigid body; one model per
component, searched together; or the components in turn, with the placed
ones fixed. Whether a job accounts for the whole AU is therefore a question
of what its models COVER, not of how many models it has: CDK4/cyclin D1 was
solved with one search model holding both chains, and a check that counted
models told that correct job it had searched for one of two kinds of chain.

A model chain covers the AU sequence it aligns to best, when that alignment
is worth anything: gemmi's score is positive and at least MIN_MATCHED
residues match. A wrong pairing scores negative (beta-lactamase against
BLIP's sequence: -182, with 52 chance matches at 20% identity), a right one
positive whatever the identity (the demo models: +263 and +165 at 100%); a
remote homologue scores lower but still positive. Identity alone would
misassign, since chance pairings reach 20-28%. The MIN_MATCHED floor keeps
a peptide or a stray short chain from claiming a kind.

gemmi only: no CCP4, no Django.
"""
import tempfile
from typing import Dict, List, Optional, Sequence, Tuple

import gemmi

#: Residues a chain must match for its best alignment to count.
MIN_MATCHED = 30

_POLYMER = {
    "PROTEIN": (gemmi.ResidueKind.AA, gemmi.PolymerType.PeptideL),
    "DNA": (gemmi.ResidueKind.DNA, gemmi.PolymerType.Dna),
    "RNA": (gemmi.ResidueKind.RNA, gemmi.PolymerType.Rna),
}


def asu_kinds(asu_file) -> List[dict]:
    """The AU contents as [{name, sequence, copies, polymerType}], one per
    kind of chain with at least one copy, from a CAsuDataFile."""
    kinds = []
    for seq in asu_file.getFileContent().seqList:
        copies = int(seq.nCopies) if seq.nCopies.isSet() else 1
        if copies <= 0:
            continue
        kinds.append({"name": str(seq.name), "sequence": str(seq.sequence),
                      "copies": copies, "polymerType": str(seq.polymerType or "PROTEIN")})
    return kinds


def chain_assignments(model_path, kinds: Sequence[dict]) -> List[Tuple[str, Optional[str], float]]:
    """[(chain, kind name or None, identity %)] for each polymer chain of the
    model: the kind it aligns to best, or None when no alignment is worth
    anything."""
    structure = gemmi.read_structure(str(model_path), format=gemmi.CoorFormat.Detect)
    structure.setup_entities()
    if len(structure) == 0:
        return []
    expanded = []
    for kind in kinds:
        residue_kind, polymer_type = _POLYMER.get(kind.get("polymerType", "PROTEIN"),
                                                 _POLYMER["PROTEIN"])
        try:
            full = gemmi.expand_one_letter_sequence(kind["sequence"], residue_kind)
        except Exception:
            continue
        expanded.append((kind["name"], full, polymer_type))
    out = []
    for chain in structure[0]:
        polymer = chain.get_polymer()
        if len(polymer) == 0:
            continue
        best = None
        for name, full, polymer_type in expanded:
            try:
                result = gemmi.align_sequence_to_polymer(
                    full, polymer, polymer_type, gemmi.AlignmentScoring())
            except Exception:
                continue
            if best is None or result.score > best[1]:
                best = (name, result.score, result.match_count, result.calculate_identity(2))
        if best and best[1] > 0 and best[2] >= MIN_MATCHED:
            out.append((chain.name, best[0], round(best[3], 1)))
        else:
            out.append((chain.name, None, 0.0))
    return out


def covered_kinds(model_path, kinds: Sequence[dict]) -> Dict[str, float]:
    """{kind name: best identity %} for the kinds the model's chains cover."""
    covered: Dict[str, float] = {}
    for _chain, name, identity in chain_assignments(model_path, kinds):
        if name is not None and identity >= covered.get(name, -1.0):
            covered[name] = identity
    return covered


def search_model_coverage(pdb_file, kinds: Sequence[dict]) -> Dict[str, float]:
    """covered_kinds of a CPdbDataFile, as its atom selection leaves it."""
    if not pdb_file.isSet():
        return {}
    if not pdb_file.isSelectionSet():
        return covered_kinds(pdb_file.getFullPath(), kinds)
    with tempfile.TemporaryDirectory() as work:
        return covered_kinds(pdb_file.getSelectedAtomsFile("selected", work), kinds)


def describe_coverage(label: str, covered: Dict[str, float]) -> str:
    """'model 1 covers CDK4 (100%), CyclinD1 (92%)' or 'model 1 covers none of them'."""
    if not covered:
        return f"{label} covers none of them"
    parts = ", ".join(f"{name} ({identity:.0f}%)" for name, identity in covered.items())
    return f"{label} covers {parts}"
