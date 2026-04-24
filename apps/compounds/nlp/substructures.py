"""Curated substructure / scaffold / functional-group catalog.

Maps chemist-typed names to SMARTS patterns. The LLM emits the user's
typed name verbatim ("pyrimidine" / "amide" / "quinoline"); the backend
resolves the name against this catalog and does an on-the-fly RDKit
`HasSubstructMatch` against each compound's SMILES (§19.3 architecture,
slice 13 design). No ingestion-time tagging, no Scaffold table — the
authority is this one file + RDKit.

Curator workflow:
- Edit the ``SCAFFOLDS`` tuple below.
- Each entry: canonical ``name``, ``aliases`` tuple, ``smarts`` string.
- Aliases are matched by normalisation (lowercase, strip non-alphanumeric
  — the same helper used by target resolution). Include common plurals,
  the -yl suffix (*"pyrimidinyl"* → *"pyrimidine"*), and common
  abbreviations when those would be typed in prompts.
- Commit and deploy — no migration, no re-indexing.

If a SMARTS is refined, ALL compounds pick up the new semantics
immediately on the next query; no retag is needed.

**Scale note.** Per-compound ``HasSubstructMatch`` at DDU scale
(a few thousand compounds per target) is sub-second. If the set grows
an order of magnitude and this becomes a hot path, add pattern
fingerprints to MolecularProperties and screen via FP-bit-subset before
the exact substructure match. Not needed yet.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Tuple

from compounds.utils import normalize_ref as normalize


@dataclass(frozen=True)
class Scaffold:
    name: str                          # canonical — this is what lands in the URL / scope sentence
    aliases: Tuple[str, ...]           # additional typed forms that map to this entry
    smarts: str                        # SMARTS pattern used by RDKit HasSubstructMatch


# Seed catalog — ordered for human reading, not algorithmically significant.
# Grouped by family. Curator refines / adds entries as eval surfaces gaps.
SCAFFOLDS: Tuple[Scaffold, ...] = (
    # ────── Aromatic N-heterocycles ──────
    Scaffold("pyridine",
             ("pyridinyl", "pyridines", "pyridyl"),
             "c1ccncc1"),
    Scaffold("pyrimidine",
             ("pyrimidinyl", "pyrimidines"),
             "c1cncnc1"),
    Scaffold("pyrazine",
             ("pyrazinyl", "pyrazines"),
             "c1cnccn1"),
    Scaffold("pyridazine",
             ("pyridazinyl", "pyridazines"),
             "c1ccnnc1"),
    # Five-membered
    Scaffold("pyrazole",
             ("pyrazolyl", "pyrazoles"),
             "c1cc[nH]n1"),
    Scaffold("imidazole",
             ("imidazolyl", "imidazoles"),
             "c1c[nH]cn1"),
    Scaffold("triazole",
             ("triazolyl", "triazoles", "1,2,3-triazole", "1,2,4-triazole"),
             "c1nnc[nH]1"),
    Scaffold("tetrazole",
             ("tetrazolyl", "tetrazoles"),
             "c1nnn[nH]1"),
    Scaffold("oxadiazole",
             ("oxadiazolyl", "oxadiazoles"),
             "c1nnco1"),
    Scaffold("thiazole",
             ("thiazolyl", "thiazoles"),
             "c1cscn1"),
    Scaffold("isoxazole",
             ("isoxazolyl", "isoxazoles"),
             "c1conc1"),
    # Bicyclics
    Scaffold("indole",
             ("indolyl", "indoles"),
             "c1ccc2[nH]ccc2c1"),
    Scaffold("benzimidazole",
             ("benzimidazolyl", "benzimidazoles"),
             "c1ccc2[nH]cnc2c1"),
    Scaffold("quinoline",
             ("quinolinyl", "quinolines"),
             "c1ccc2ncccc2c1"),
    Scaffold("quinazoline",
             ("quinazolinyl", "quinazolines"),
             "c1ccc2ncncc2c1"),
    Scaffold("purine",
             ("purinyl", "purines"),
             "c1ncc2[nH]cnc2n1"),

    # ────── Saturated N-heterocycles ──────
    Scaffold("piperidine",
             ("piperidinyl", "piperidines"),
             "C1CCNCC1"),
    Scaffold("piperazine",
             ("piperazinyl", "piperazines"),
             "C1CNCCN1"),
    Scaffold("morpholine",
             ("morpholinyl", "morpholines"),
             "C1COCCN1"),
    Scaffold("pyrrolidine",
             ("pyrrolidinyl", "pyrrolidines"),
             "C1CCNC1"),

    # ────── Functional groups ──────
    Scaffold("amide",
             ("amides", "carboxamide"),
             "C(=O)N"),
    Scaffold("sulfonamide",
             ("sulfonamides", "sulphonamide"),
             "S(=O)(=O)N"),
    Scaffold("sulfone",
             ("sulfones", "sulphone"),
             "S(=O)(=O)"),
    Scaffold("carboxylic acid",
             ("carboxylic acids", "carboxyl"),
             "C(=O)[OH]"),
    Scaffold("ester",
             ("esters",),
             "C(=O)O[C,c]"),
    Scaffold("nitrile",
             ("nitriles", "cyano"),
             "C#N"),
    Scaffold("primary amine",
             ("primary amines",),
             "[NX3;H2]"),
    Scaffold("secondary amine",
             ("secondary amines",),
             "[NX3;H1;!$(N-C=O)]"),
    Scaffold("trifluoromethyl",
             ("cf3", "trifluoromethyls"),
             "C(F)(F)F"),

    # ────── Carbocycles / generic aryl ──────
    Scaffold("phenyl",
             ("phenyls", "benzene", "benzenes"),
             "c1ccccc1"),
)


def _build_pools() -> Dict[str, Tuple[Scaffold, str]]:
    """Flatten the catalog into normalised-string → (scaffold, matched_via)."""
    pool: Dict[str, Tuple[Scaffold, str]] = {}
    for scaffold in SCAFFOLDS:
        name_norm = normalize(scaffold.name)
        if name_norm:
            pool.setdefault(name_norm, (scaffold, "name"))
        for alias in scaffold.aliases:
            alias_norm = normalize(alias)
            if alias_norm:
                pool.setdefault(alias_norm, (scaffold, "alias"))
    return pool


_POOLS = _build_pools()


def lookup_scaffold_by_typed(as_typed: str) -> Optional[Tuple[Scaffold, str]]:
    """Exact normalised match — returns ``(scaffold, matched_via)`` or None."""
    norm = normalize(as_typed or "")
    if not norm:
        return None
    return _POOLS.get(norm)


def lookup_scaffolds_by_substring(as_typed: str) -> Iterable[Tuple[Scaffold, str]]:
    """Partial match — normalised query is a substring of any pool entry.
    Used when exact lookup misses. Returns an iterable of (scaffold, 'substring')."""
    norm = normalize(as_typed or "")
    if not norm:
        return []
    seen = set()
    hits = []
    for entry_norm, (scaffold, _via) in _POOLS.items():
        if norm in entry_norm and scaffold.name not in seen:
            seen.add(scaffold.name)
            hits.append((scaffold, "substring"))
    return hits


def lookup_scaffold_by_id(scaffold_id: str) -> Optional[Scaffold]:
    """Pinning lookup — `scaffold_id` here is the canonical name (we don't
    assign UUIDs; names are already stable identifiers)."""
    if not scaffold_id:
        return None
    for scaffold in SCAFFOLDS:
        if scaffold.name == scaffold_id:
            return scaffold
    return None


def all_scaffold_names() -> Tuple[str, ...]:
    return tuple(s.name for s in SCAFFOLDS)
