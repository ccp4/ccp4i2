"""Complex templates among MrParse's hits: one entry matching several targets.

MrParse searches the sequences of a complex merged into one query, and cuts
every hit to one chain. When hits from ONE PDB entry match DIFFERENT targets
(6p8e chains B and A matched CDK4 and cyclin D1), that entry is a template
for the complex: its matching chains, kept in the entry's own frame, make
one search model to place as a single rigid body (docs/multi-component-mr.md,
route A). CDK4/cyclin D1 was solved so, and nothing offered it: every hit
was a chain on its own.

The template is assembled from the per-chain cuts MrParse wrote, which are
the entry's coordinates with residues removed and nothing moved, so the
chains sit exactly as they do in the entry. One chain per target: an entry
holding two copies of the complex gives one copy; the copies are a number
to search for.

gemmi only: no CCP4, no Django.
"""
import os
from typing import Dict, List, Optional, Sequence, Tuple

import gemmi


def read_targets(seqin_path) -> List[Tuple[str, str]]:
    """[(name, sequence)] of a FASTA file as MrParse merges it: in order,
    identical sequences dropped after the first (mr_sequence.py
    merge_multiple_sequences)."""
    targets: List[Tuple[str, str]] = []
    name, lines = None, []

    def flush():
        if name is None:
            return
        seq = "".join(lines).replace(" ", "").upper()
        if seq and seq not in {s for _, s in targets}:
            targets.append((name, seq))

    with open(str(seqin_path)) as stream:
        for line in stream:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                name, lines = line[1:].strip().split()[0] if line[1:].strip() else "", []
            elif name is not None:
                lines.append(line.strip())
    flush()
    return targets


def target_of(hit: dict, targets: Sequence[Tuple[str, str]]) -> Optional[int]:
    """Index of the target a hit lies on, by the midpoint of its query range
    in the merged sequence (1-based, as MrParse reports query_start/stop)."""
    try:
        start, stop = int(hit["query_start"]), int(hit["query_stop"])
    except (KeyError, TypeError, ValueError):
        return None
    midpoint = (start + stop) / 2.0
    offset = 0
    for index, (_name, seq) in enumerate(targets):
        if offset < midpoint <= offset + len(seq):
            return index
        offset += len(seq)
    return None


def find_complexes(hits: Sequence[dict], targets: Sequence[Tuple[str, str]]) -> List[dict]:
    """Entries whose hits cover more than one target, as
    [{entry, chains: [{target, target_index, hit, chain, identity, pdb_file}]}],
    best identity first within an entry, entries in the order first seen.
    One chain per target: the best-identity hit for it."""
    by_entry: Dict[str, Dict[int, dict]] = {}
    order: List[str] = []
    for hit in hits:
        entry = str(hit.get("pdb_id") or "").lower()
        if not entry or not hit.get("pdb_file"):
            continue
        index = target_of(hit, targets)
        if index is None:
            continue
        chain = {"target": targets[index][0], "target_index": index,
                 "hit": hit.get("name"), "chain": hit.get("chain_id"),
                 "identity": float(hit.get("seq_ident") or 0.0),
                 "pdb_file": hit["pdb_file"]}
        if entry not in by_entry:
            order.append(entry)
        held = by_entry.setdefault(entry, {})
        if index not in held or chain["identity"] > held[index]["identity"]:
            held[index] = chain
    found = []
    for entry in order:
        chains = sorted(by_entry[entry].values(), key=lambda c: c["target_index"])
        if len(chains) > 1:
            found.append({"entry": entry, "chains": chains})
    return found


def write_complex(entry: dict, cut_dir, out_dir) -> Optional[str]:
    """One file holding the entry's matching chains, from MrParse's per-chain
    cuts (all in the entry's frame). Returns its path, or None if a cut is
    missing or unreadable."""
    structure = gemmi.Structure()
    model = gemmi.Model("1")
    cell = None
    spacegroup = None
    for chain in entry["chains"]:
        path = os.path.join(str(cut_dir), os.path.basename(chain["pdb_file"]))
        if not os.path.isfile(path):
            return None
        try:
            source = gemmi.read_structure(path, format=gemmi.CoorFormat.Detect)
        except Exception:
            return None
        if len(source) == 0 or len(source[0]) == 0:
            return None
        if cell is None and source.cell.is_crystal():
            cell, spacegroup = source.cell, source.spacegroup_hm
        for source_chain in source[0]:
            if len(source_chain) == 0:
                continue
            copy = gemmi.Chain(source_chain.name or chain["chain"] or "A")
            for residue in source_chain:
                copy.add_residue(residue)
            model.add_chain(copy)
    if len(model) < 2:
        return None
    structure.add_model(model)
    if cell is not None:
        structure.cell = cell
        structure.spacegroup_hm = spacegroup or ""
    structure.setup_entities()
    chains = "".join(c["chain"] or "" for c in entry["chains"])
    out = os.path.join(str(out_dir), f"{entry['entry']}_{chains}_complex.pdb")
    structure.write_pdb(out)
    return out


def describe(entry: dict) -> str:
    """'6p8e chains B (CDK4 92%), A (CyclinD1 100%)'"""
    parts = ", ".join(f"{c['chain']} ({c['target']} {100 * c['identity']:.0f}%)"
                      for c in entry["chains"])
    return f"{entry['entry']} chains {parts}"
