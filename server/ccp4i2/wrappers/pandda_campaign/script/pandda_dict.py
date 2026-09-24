"""Normalise a ligand dictionary so every PanDDA build reads its bond orders.

Design note: docs/pandda-campaign-design.md §4.6.1 and §4.7.

The PanDDA 2 build that CCP4 ships reads a ligand's bond orders from
``_chem_comp_bond.type`` only. A dictionary that spells the same column
``_chem_comp_bond.value_order`` (the PDBx/wwPDB form, with
``pdbx_aromatic_flag``) makes that reader's bond loop run zero times, and the
ligand is built with no bonds at all: RDKit's conformer search then stacks
the atoms on top of each other and PanDDA reports nothing wrong. Upstream
fixed the reader (PR #99, via gemmi ``ChemComp``); the bundle predates it.

This is a *normaliser*, not a converter: whichever spelling the producer
used, the staged ``dict.cif`` carries both, consistently. A ``type``-only
dictionary is returned untouched, because both readers already handle it.
It lives next to ``prepare_mtz_for_pandda`` (lib/pandda_export.py), which
does the same job for the FreeR label, and it is called at staging so the
task works with whichever PanDDA the user has.

Pure gemmi; no CCP4 binary, no Django.
"""
import logging
from pathlib import Path

import gemmi

logger = logging.getLogger(f"ccp4i2:{__name__}")

BOND_PREFIX = "_chem_comp_bond."

#: PDBx ``value_order`` tokens -> the ``type`` tokens the old reader's map
#: accepts (``single, double, triple, aromatic, deloc`` and upper-case
#: variants of the first three). Keys are matched case-insensitively.
VALUE_ORDER_TO_TYPE = {
    "sing": "single",
    "doub": "double",
    "trip": "triple",
    "arom": "aromatic",
    "delo": "deloc",
}

#: What an unrecognised ``value_order`` token (``quad``, ``pi``, ``poly``)
#: becomes. A bond of *some* order keeps the molecule connected, which is the
#: property the old reader loses; the exact order is left to the newer reader,
#: which keeps reading ``value_order``. Logged when it happens.
FALLBACK_TYPE = "single"


def bond_spellings(block) -> set:
    """Which of ``{'type', 'value_order'}`` the block's bond table carries."""
    found = set()
    for spelling in ("type", "value_order"):
        if block.find_values(BOND_PREFIX + spelling):
            found.add(spelling)
    return found


def _bond_blocks(doc):
    return [b for b in doc if b.find_values(BOND_PREFIX + "atom_id_1")]


def type_token(value_order: str, aromatic_flag: str = "") -> str:
    """The ``type`` token equivalent to a PDBx ``value_order`` token.

    ``pdbx_aromatic_flag == 'Y'`` wins over the order token, as the old
    reader treats ``aromatic`` as an order in its own right.
    """
    if aromatic_flag.strip().upper() == "Y":
        return "aromatic"
    token = VALUE_ORDER_TO_TYPE.get(value_order.strip().lower())
    if token is None:
        logger.warning(
            "prepare_dict_for_pandda: unrecognised value_order %r, "
            "writing type=%s", value_order, FALLBACK_TYPE)
        token = FALLBACK_TYPE
    return token


def add_type_column(block) -> bool:
    """Give ``block``'s bond table a ``type`` column derived from
    ``value_order``, and an ``aromatic`` column from ``pdbx_aromatic_flag``
    when it has one. Returns True if anything was added."""
    spellings = bond_spellings(block)
    if "type" in spellings or "value_order" not in spellings:
        return False
    table = block.find(BOND_PREFIX, ["value_order"])
    table.ensure_loop()
    loop = table.loop
    has_flag = bool(block.find_values(BOND_PREFIX + "pdbx_aromatic_flag"))
    has_aromatic = bool(block.find_values(BOND_PREFIX + "aromatic"))
    new_tags = [BOND_PREFIX + "type"]
    if has_flag and not has_aromatic:
        new_tags.append(BOND_PREFIX + "aromatic")
    loop.add_columns(new_tags, value="?")

    wanted = ["value_order", "type"]
    if has_flag:
        wanted.append("pdbx_aromatic_flag")
    if BOND_PREFIX + "aromatic" in new_tags:
        wanted.append("aromatic")
    rows = block.find(BOND_PREFIX, wanted)
    for row in rows:
        flag = row[2] if has_flag else ""
        row[1] = type_token(row[0], flag)
        if BOND_PREFIX + "aromatic" in new_tags:
            row[3] = "y" if flag.strip().upper() == "Y" else "n"
    return True


def needs_normalising(src_path) -> bool:
    """True if any bond table in the file has ``value_order`` but no ``type``."""
    doc = gemmi.cif.read(str(src_path))
    return any(bond_spellings(b) == {"value_order"} for b in _bond_blocks(doc))


def prepare_dict_for_pandda(src_path, staging_dir) -> Path:
    """Return a path to a dictionary every PanDDA build reads the same way.

    If every bond table already carries ``_chem_comp_bond.type`` (or there is
    no bond table to read), ``src_path`` is returned unchanged and nothing is
    written. Otherwise a copy with the ``type`` (and, where derivable,
    ``aromatic``) column added is written to ``staging_dir`` and its path is
    returned. The original columns are kept, so the newer reader sees what it
    saw before.

    Mirrors ``prepare_mtz_for_pandda``: same signature, same contract.
    """
    src_path = Path(src_path)
    doc = gemmi.cif.read(str(src_path))
    changed = False
    for block in _bond_blocks(doc):
        changed = add_type_column(block) or changed
    if not changed:
        return src_path
    staging_dir = Path(staging_dir)
    staging_dir.mkdir(parents=True, exist_ok=True)
    out_path = staging_dir / f"{src_path.stem}_pandda.cif"
    doc.write_file(str(out_path))
    logger.info("prepare_dict_for_pandda: added _chem_comp_bond.type to %s -> %s",
                src_path.name, out_path.name)
    return out_path
