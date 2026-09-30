"""AceDRG's link instruction language, written from a declared description.

AceDRG reads a link instruction as one stream of words (line breaks mean
nothing) and modifies each monomer imperatively:

    DELETE ATOM <name> <monomer>               one item per DELETE
    DELETE BOND <atom> <atom> <monomer>
    CHANGE BOND <atom> <atom> <order> <monomer>   CHANGE runs on until the
    CHANGE CHARGE <monomer> <atom> <charge>        next section keyword
    ADD ATOM <name> <element> <charge> <monomer>   ADD likewise
    ADD BOND <atom> <atom> <order> <monomer>

The monomer number comes last everywhere except CHARGE, where it comes
first. A second item after one DELETE stops AceDRG with "Unknown keyword".
An unrecognised word inside a CHANGE or ADD section is worse: AceDRG's parser
neither consumes it nor reports it, and loops forever (checked against
covLink.py in ccp4-20260702: `CHANGE OE1 DOUBLE 2` never returns).

So MakeLink does not let anyone write this language by hand where it can be
avoided. The edits to each monomer are held as data (MonomerEdits: the atoms
gone, the order each changed bond should have, the charge each changed atom
should carry); edit_problems() says where a description contradicts itself,
and edit_words() writes it in one canonical form. The one place free text
remains, the Advanced "extra instructions" box, is checked word by word
against the grammar above before AceDRG ever sees it.

Pure Python with no CCP4 dependency, so all of it is unit-tested in CI.
"""
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

BOND_ORDERS = ("SINGLE", "DOUBLE", "TRIPLE")

# What AceDRG accepts after CHANGE BOND / ADD BOND. Wider than BOND_ORDERS,
# which is what the task offers: the extra box may use the others.
_ACEDRG_BOND_ORDERS = BOND_ORDERS + ("AROMATIC", "DELOC")


def _bond_key(atom_1: str, atom_2: str) -> Tuple[str, str]:
    """A bond named without regard to which atom was given first."""
    return tuple(sorted((atom_1, atom_2)))


@dataclass
class MonomerEdits:
    """What should differ in one monomer, declared rather than commanded.

    Lists rather than sets or dicts so that a contradiction -- one bond given
    two orders, one atom two charges -- survives to be reported instead of
    being silently resolved by whichever came last.
    """
    deletes: List[str] = field(default_factory=list)
    bond_orders: List[Tuple[str, str, str]] = field(default_factory=list)
    charges: List[Tuple[str, int]] = field(default_factory=list)

    def is_empty(self) -> bool:
        return not (self.deletes or self.bond_orders or self.charges)


def _bad_name(name: str) -> bool:
    # A name is one word to AceDRG: blank or spaced would shift every
    # argument after it into the wrong slot.
    return not name or any(ch.isspace() for ch in name)


def edit_problems(edits: MonomerEdits, link_atom: str = "") -> List[str]:
    """Every way this description of one monomer contradicts itself.

    An empty list means edit_words() can be written and AceDRG will read
    what was meant. Chemistry (valence, whether an atom exists) is AceDRG's
    to judge, and it does so loudly.
    """
    problems = []
    for name in edits.deletes:
        if _bad_name(name):
            problems.append(f'"{name}" is not an atom name')
    for atom_1, atom_2, order in edits.bond_orders:
        if _bad_name(atom_1) or _bad_name(atom_2):
            problems.append(f'"{atom_1}-{atom_2}" is not a bond between two named atoms')
        elif atom_1 == atom_2:
            problems.append(f"A bond needs two different atoms, not {atom_1} twice")
        if order.upper() not in BOND_ORDERS:
            problems.append(
                f"{atom_1}-{atom_2}: bond order must be one of {', '.join(BOND_ORDERS)}, not \"{order}\"")
    for name, _charge in edits.charges:
        if _bad_name(name):
            problems.append(f'"{name}" is not an atom name')
    if problems:
        return problems

    deleted = set(edits.deletes)
    if link_atom and link_atom in deleted:
        problems.append(f"{link_atom} is the linking atom, so it cannot also be deleted")

    orders: Dict[Tuple[str, str], str] = {}
    for atom_1, atom_2, order in edits.bond_orders:
        key = _bond_key(atom_1, atom_2)
        order = order.upper()
        if key in orders and orders[key] != order:
            problems.append(
                f"Bond {atom_1}-{atom_2} is given two orders, {orders[key]} and {order}")
        orders.setdefault(key, order)
        gone = [a for a in key if a in deleted]
        if gone:
            problems.append(
                f"Bond {atom_1}-{atom_2} cannot change order: {gone[0]} is deleted")

    charges: Dict[str, int] = {}
    for name, charge in edits.charges:
        if name in charges and charges[name] != charge:
            problems.append(
                f"{name} is given two charges, {charges[name]:+d} and {charge:+d}")
        charges.setdefault(name, charge)
        if name in deleted:
            problems.append(f"{name} cannot change charge: it is deleted")
    return problems


def edit_words(edits: MonomerEdits, monomer: int) -> List[str]:
    """The AceDRG words for one monomer's edits, in one canonical order.

    Every deletion gets its own DELETE, and the bond and charge changes share
    one CHANGE section, which ends at the next keyword. Repeats are written
    once. Assumes edit_problems() found nothing.
    """
    if monomer not in (1, 2):
        raise ValueError(f"monomer must be 1 or 2, not {monomer!r}")
    serial = str(monomer)
    words: List[str] = []

    seen = set()
    for name in edits.deletes:
        if name not in seen:
            seen.add(name)
            words += ["DELETE", "ATOM", name, serial]

    changes: List[str] = []
    seen_bonds = set()
    for atom_1, atom_2, order in edits.bond_orders:
        key = _bond_key(atom_1, atom_2)
        if key not in seen_bonds:
            seen_bonds.add(key)
            changes += ["BOND", atom_1, atom_2, order.upper(), serial]
    seen_atoms = set()
    for name, charge in edits.charges:
        if name not in seen_atoms:
            seen_atoms.add(name)
            # The one section whose monomer number comes first.
            changes += ["CHARGE", serial, name, str(int(charge))]
    if changes:
        words += ["CHANGE"] + changes
    return words


def _is_int(word: str) -> bool:
    try:
        int(word)
    except ValueError:
        return False
    return True


# For each section item: its keyword and a checker for each argument, in
# AceDRG's order. "monomer" is 1 or 2; "name" any single word.
_MONOMER = ("monomer number (1 or 2)", lambda w: w in ("1", "2"))
_NAME = ("atom name", lambda w: True)
_ORDER = ("bond order (" + ", ".join(_ACEDRG_BOND_ORDERS) + ")",
          lambda w: w.upper() in _ACEDRG_BOND_ORDERS)
_CHARGE = ("whole-number charge", _is_int)
_ELEMENT = ("element symbol", lambda w: w.isalpha())

_ITEMS = {
    "DELETE": {"ATOM": (_NAME, _MONOMER), "BOND": (_NAME, _NAME, _MONOMER)},
    "CHANGE": {"BOND": (_NAME, _NAME, _ORDER, _MONOMER),
               "CHARGE": (_MONOMER, _NAME, _CHARGE)},
    "ADD": {"ATOM": (_NAME, _ELEMENT, _CHARGE, _MONOMER),
            "BOND": (_NAME, _NAME, _ORDER, _MONOMER)},
}


def instruction_words(text: str) -> List[str]:
    """The words of free-text instructions, as AceDRG would read them."""
    words = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        words += [w for w in line.split() if w.upper() != "LINK:"]
    return words


def extra_instruction_problems(text: str) -> List[str]:
    """Why AceDRG would misread these extra instructions, if it would.

    Checked against the grammar in the module docstring. The text must
    begin with a section keyword, so that nothing in it can be read as a
    continuation of the CHANGE section edit_words() may have ended with.
    Returns at most one problem: after the first, positions are guesswork.
    """
    words = instruction_words(text)
    section = None
    single_item = False
    i = 0
    while i < len(words):
        word = words[i].upper()
        if word in _ITEMS:
            section, single_item = word, word == "DELETE"
            i += 1
            if i == len(words):
                return [f'"{words[i - 1]}" must be followed by what to {word.lower()}']
            continue
        if section is None:
            where = "must begin with DELETE, CHANGE or ADD" if i == 0 else "is not a keyword here"
            return [f'Extra instructions: "{words[i]}" {where}']
        arguments = _ITEMS[section].get(word)
        if arguments is None:
            allowed = " or ".join(_ITEMS[section])
            return [f'Extra instructions: after {section}, expected {allowed}, not "{words[i]}"']
        item = words[i + 1:i + 1 + len(arguments)]
        if len(item) < len(arguments):
            return [f"Extra instructions: {section} {word} needs {len(arguments)} values after it"]
        for (what, ok), value in zip(arguments, item):
            if not ok(value):
                return [f'Extra instructions: in {section} {word}, "{value}" is not a {what}']
        i += 1 + len(arguments)
        if single_item:
            # AceDRG takes exactly one item per DELETE.
            section = None
    return []
