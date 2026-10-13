"""AU-contents rows from a file's digest, and merging them into a table.

The server's copy of what ProvideAsuContents' page did in the browser
(client/renderer/components/task/task-elements/mmcif-sequence-parser.ts:
chainsFromDigest, deduplicateChains, mergeAsuEntries). The page now leaves it
to the server, so a sequence, model or AU-contents file fills the table the
same way whoever sets it: the page, an agent through the API, or i2run.

Entries are plain dicts: name, sequence, polymerType, description, nCopies.
"""
import re


def _nucleic_type(sequence):
    return "RNA" if re.search(r"U", sequence, re.I) else "DNA"


def polymer_type_from_molecule(molecule_type, sequence):
    upper = (molecule_type or "").upper()
    if upper in ("PROTEIN", "DNA", "RNA"):
        return upper
    if upper == "NUCLEIC":
        return _nucleic_type(sequence)
    return "OTHER"


def _classify_chain(chain_id, composition, sequence=""):
    composition = composition or {}
    if chain_id in (composition.get("peptides") or []):
        return "PROTEIN"
    detail = next((d for d in composition.get("chainDetails") or [] if d.get("id") == chain_id), None)
    if detail and detail.get("type") == "protein":
        return "PROTEIN"
    if chain_id in (composition.get("nucleics") or []) or (detail and detail.get("type") == "nucleic"):
        return _nucleic_type(sequence)
    return "OTHER"


def chains_from_digest(digest, label):
    """Every polymer chain a sequence or model digest describes, in file order."""
    if not digest:
        return []
    if digest.get("moleculeType") and digest.get("sequence"):
        sequence = re.sub(r"\s", "", str(digest["sequence"]))
        return [{"chainId": digest.get("name") or label, "sequence": sequence,
                 "polymerType": polymer_type_from_molecule(digest["moleculeType"], sequence),
                 "description": digest.get("description") or label}]
    if digest.get("composition") and digest.get("sequences"):
        composition = digest["composition"]
        ids = list(composition.get("peptides") or []) + list(composition.get("nucleics") or [])
        return [{"chainId": cid, "sequence": str(digest["sequences"][cid]),
                 "polymerType": _classify_chain(cid, composition, str(digest["sequences"][cid])),
                 "description": f"{label} chain {cid}"}
                for cid in ids if digest["sequences"].get(cid)]
    return []


def deduplicate_chains(chains):
    """Identical chains (same polymer type and sequence) as one entry with copies."""
    groups = {}
    for chain in chains:
        groups.setdefault(f"{chain['polymerType']}:{chain['sequence']}", []).append(chain)
    entries = []
    for group in groups.values():
        ids = [c["chainId"] for c in group]
        first = group[0]
        entries.append({"name": "Chain_" + "_".join(ids), "sequence": first["sequence"],
                        "polymerType": first["polymerType"],
                        "description": first.get("description") or ", ".join(f"Chain {i}" for i in ids),
                        "nCopies": len(group)})
    return entries


def entries_from_digest(digest, label):
    """(entries, chains read): an AU-contents digest lists rows already;
    sequence and model digests list chains."""
    if isinstance((digest or {}).get("sequences"), list) and not digest.get("composition"):
        entries = [{"name": s.get("name"), "sequence": s.get("sequence"),
                    "polymerType": s.get("polymerType"), "description": s.get("description"),
                    "nCopies": s.get("nCopies") or 1} for s in digest["sequences"]]
        return entries, len(entries)
    chains = chains_from_digest(digest, label)
    return deduplicate_chains(chains), len(chains)


def _key(entry):
    sequence = re.sub(r"\s", "", str(entry.get("sequence") or ""))
    return f"{str(entry.get('polymerType')).upper()}:{sequence}"


def merge_asu_entries(existing, incoming):
    """``existing`` with ``incoming`` added: a sequence already present gains
    copies rather than a second row."""
    merged = [dict(e) for e in existing]
    by_key = {_key(e): e for e in merged}
    for entry in incoming:
        hit = by_key.get(_key(entry))
        if hit is not None:
            hit["nCopies"] = int(hit.get("nCopies") or 0) + int(entry.get("nCopies") or 0)
        else:
            merged.append(dict(entry))
            by_key[_key(entry)] = merged[-1]
    return merged
