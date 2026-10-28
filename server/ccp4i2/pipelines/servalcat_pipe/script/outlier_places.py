"""Servalcat's geometry outliers as places to look at in Moorhen.

The report lists the outliers of the final cycle as tables: atoms, value,
ideal, Z. Looking at one meant finding it by hand in a viewer. Here each
becomes a place: the atoms it names as a Coot selection, and a scene of the
job's refined model and maps centred on them. The report offers the worst as
buttons, and a viewer already open on the project takes the click and moves
there (renderer/lib/viewer-channel.ts), so stepping through them is a camera
move each, not a new window.

Servalcat names atoms gemmi's way, ``G/DC 2012/OP2.A`` (chain, residue name
and number with any insertion code, atom, and alternate location after the
dot), in tags atom1..atom4, atomc (a chiral centre) and atom (a plane).

The same goes for the pipeline's other two judges of the model. MolProbity
(validate_protein's ``Molprobity`` block) names residues, ``chain`` and
``seqnum`` with any insertion code, and clashes as "A 187 CYS O" pairs; the ADP
analysis names atoms gemmi's way, like Servalcat. Every place has a kind and an
index within it, which is all outlierScene is asked for.

A prototype: the selection, the scene and the choice of the worst are here,
and tested; which outliers deserve a button, and how many, is for the people
reading the report to say.
"""
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional

KINDS = ("vdw", "bond", "angle", "torsion", "chir", "plane", "staca", "stacd")
KIND_LABELS = {
    "vdw": "Clash", "bond": "Bond", "angle": "Angle", "torsion": "Torsion",
    "chir": "Chirality", "plane": "Plane", "staca": "Stacking angle",
    "stacd": "Stacking distance",
}

_ATOM_TAG = re.compile(r"^atom(\d+|c)?$")
_ATOM_LABEL = re.compile(
    r"^(?P<chain>[^/]+)/(?P<resname>\S+) (?P<seq>-?\d+)(?P<icode>[A-Za-z]?)"
    r"/(?P<atom>[^.]+?)(?:\.(?P<alt>\S))?$")


def atom_cid(label: str) -> Optional[str]:
    """A Coot CID for one of Servalcat's atom labels, or None if it is not
    one. The alternate location is dropped: the camera centres on the atom
    whichever conformer it is drawn in."""
    match = _ATOM_LABEL.match(label.strip())
    if not match:
        return None
    icode = f".{match['icode']}" if match["icode"] else ""
    return f"//{match['chain']}/{match['seq']}{icode}/{match['atom']}"


def residue_cid(label: str) -> Optional[str]:
    """The CID of the residue one of Servalcat's atom labels belongs to."""
    cid = atom_cid(label)
    return cid.rsplit("/", 1)[0] if cid else None


# The marker, over Moorhen's bonds of radius 0.1 with balls of 0.2: fat enough
# to enclose them (a representation of the same radii loses the depth test to
# the one beneath and does not show), the atoms involved fatter still.
MARKER_COLOUR = "#ff00ff"
RESIDUE_GEOMETRY = {"bondRadius": 0.16, "ballRadius": 0.32}
ATOM_GEOMETRY = {"bondRadius": 0.16, "ballRadius": 0.45}


def marker_representations(atoms: List[str]) -> List[Dict]:
    """How an outlier is marked on the model: its residues drawn again as
    fatter sticks with magenta carbons, element colours kept so the chemistry
    stays readable; then the atoms involved as fatter magenta balls, whatever
    their element (a plain colour reaches carbons only, and an outlier's atoms
    are mostly not carbon: NE2, ZN, O). One representation per CID."""
    return _markers([residue_cid(a) for a in atoms], [atom_cid(a) for a in atoms])


def _markers(residue_cids: List[Optional[str]], atom_cids: List[Optional[str]]) -> List[Dict]:
    reps = []
    for cid in dict.fromkeys(c for c in residue_cids if c):
        reps.append({"style": "CBs", "selection": cid, "colour": MARKER_COLOUR,
                     "geometry": RESIDUE_GEOMETRY})
    for cid in dict.fromkeys(c for c in atom_cids if c):
        reps.append({"style": "CBs", "selection": cid,
                     "colour": {"raw": {"ruleType": "molecule", "args": [cid, MARKER_COLOUR],
                                        "isMultiColourRule": False,
                                        "applyColourToNonCarbonAtoms": True}},
                     "geometry": ATOM_GEOMETRY})
    return reps


def final_block(root: ET.Element) -> ET.Element:
    """The refinement whose outliers describe the final model: the
    after-waters run when there was one, else the main one; the root itself
    for the plain servalcat task's program.xml."""
    for tag in ("SERVALCAT_WATERS", "SERVALCAT_FIRST"):
        block = root.find(tag)
        if block is not None and block.find("cycle") is not None:
            return block
    return root


def outliers(root: ET.Element) -> List[Dict]:
    """Every outlier of the final cycle: kind, its index within that kind
    (what outlierScene is asked for), atom labels, and Z."""
    found = []
    block = final_block(root)
    for kind in KINDS:
        for index, element in enumerate(block.findall(f"./cycle[last()]/geom/outliers/{kind}")):
            atoms = [child.text.strip() for child in element
                     if _ATOM_TAG.match(child.tag) and child.text and child.text.strip()]
            try:
                z = float(element.findtext("z"))
            except (TypeError, ValueError):
                z = math.nan
            found.append({"kind": kind, "index": index, "atoms": atoms, "z": z})
    return found


def worst(root: ET.Element, limit: int = 12) -> List[Dict]:
    """The outliers with the largest |Z|, of every kind together."""
    usable = [o for o in outliers(root) if o["atoms"] and not math.isnan(o["z"])]
    return sorted(usable, key=lambda o: -abs(o["z"]))[:limit]


def describe(outlier: Dict) -> str:
    """'Chirality: D/VAL 173 (CA, N, C, CB) (Z -25.2)'; 'Bond: A/HIS 188/NE2 –
    A/ZN 206/ZN (Z +16.2)'. Atoms of one residue name the residue once.
    MolProbity's and the ADP analysis's places carry their own label."""
    if "label" in outlier:
        return outlier["label"]
    atoms = outlier["atoms"]
    residues = {atom.rsplit("/", 1)[0] for atom in atoms}
    if len(residues) == 1 and len(atoms) > 1:
        names = ", ".join(atom.rsplit("/", 1)[1] for atom in atoms)
        where = f"{residues.pop()} ({names})"
    elif len(atoms) <= 2:
        where = " – ".join(atoms)
    else:
        where = f"{atoms[0]} … {atoms[-1]}"
    return f"{KIND_LABELS.get(outlier['kind'], outlier['kind'])}: {where} (Z {outlier['z']:+.1f})"


# MolProbity's residue outliers: section, kind, label, and whether a larger
# score is worse (a C-beta deviation in Angstrom) or better (the Ramachandran
# and rotamer probabilities). Non-planar omegas and flips have no score.
MOLPROBITY_SECTIONS = (
    ("Ramachandran_outliers", "rama", "Ramachandran", False),
    ("Rotamer_outliers", "rotamer", "Rotamer", False),
    ("CBeta_outliers", "cbeta", "C\u03b2 deviation", True),
    ("Nonplanar_omegas", "omega", "Non-planar \u03c9", False),
    ("Side_chain_flips", "flip", "Side-chain flip", False),
)
_SEQNUM = re.compile(r"^\s*(?P<seq>-?\d+)\s*(?P<icode>[A-Za-z]?)\s*$")
# "A 187 CYS O", "A 100A ASN HD21", "A 52 SER OG B" (an alternate location last)
_CLASH_ATOM = re.compile(
    r"^\s*(?P<chain>\S+)\s+(?P<seq>-?\d+)(?P<icode>[A-Za-z]?)\s+(?P<resname>\S+)\s+(?P<atom>\S+)")


def _cid(chain: str, seqnum: str, atom: Optional[str] = None) -> Optional[str]:
    match = _SEQNUM.match(seqnum or "")
    if not chain or not chain.strip() or not match:
        return None
    icode = f".{match['icode']}" if match["icode"] else ""
    cid = f"//{chain.strip()}/{match['seq']}{icode}"
    return f"{cid}/{atom}" if atom else cid


def _score(text: Optional[str]) -> float:
    try:
        return float(text)
    except (TypeError, ValueError):
        return math.nan


def molprobity_places(root: ET.Element) -> List[Dict]:
    """MolProbity's residue outliers, section by section and worst first
    within each, then its clashes, worst overlap first. Each place names its
    residues (and, for a clash, the two atoms) as CIDs."""
    found = []
    block = root.find(".//Molprobity")
    if block is None:
        return found
    for section, kind, title, larger_is_worse in MOLPROBITY_SECTIONS:
        places = []
        for index, element in enumerate(block.findall(f"./{section}/Outlier")):
            cid = _cid(element.get("chain", ""), element.get("seqnum", ""))
            if cid is None:
                continue
            score = _score(element.get("score"))
            where = f"{element.get('chain')}/{element.get('name', '')} {element.get('seqnum', '').strip()}"
            if math.isnan(score):
                label = f"{title}: {where}"
            elif larger_is_worse:
                label = f"{title}: {where} ({score:.2f} \u00c5)"
            else:  # a probability, which MolProbity's XML rounds to 0.001
                label = f"{title}: {where} (p {'< 0.001' if score < 0.0005 else f'{score:.3f}'})"
            places.append({"kind": kind, "index": index, "label": label, "residues": [cid], "atoms_cids": [],
                           "score": score})
        if any(not math.isnan(p["score"]) for p in places):
            places.sort(key=lambda p: (math.isnan(p["score"]),
                                       -p["score"] if larger_is_worse else p["score"]))
        found.extend(places)
    clashes = []
    for index, element in enumerate(block.findall("./Clashes/Outlier")):
        atoms = [_CLASH_ATOM.match(element.get(end, "")) for end in ("first_atom", "second_atom")]
        if not all(atoms):
            continue
        overlap = _score(element.get("overlap"))
        names = " – ".join(f"{m['chain']}/{m['resname']} {m['seq']}{m['icode']}/{m['atom']}" for m in atoms)
        clashes.append({
            "kind": "clash", "index": index,
            "label": f"Clash: {names}" + ("" if math.isnan(overlap) else f" ({abs(overlap):.2f} \u00c5)"),
            "residues": [_cid(m["chain"], m["seq"] + m["icode"]) for m in atoms],
            "atoms_cids": [_cid(m["chain"], m["seq"] + m["icode"], m["atom"]) for m in atoms],
            "score": overlap})
    # overlaps are negative: the most negative is the worst
    clashes.sort(key=lambda p: (math.isnan(p["score"]), -abs(p["score"]) if not math.isnan(p["score"]) else 0))
    return found + clashes


def adp_places(root: ET.Element) -> List[Dict]:
    """The ADP analysis's atoms outside its limits, furthest outside first,
    low and high together."""
    block = root.find(".//ADP_ANALYSIS/outliers")
    if block is None:
        return []
    found = []
    for side, kind in (("low", "adp_low"), ("high", "adp_high")):
        limit = _score(block.findtext(f"adp_limit_{side}"))
        for index, element in enumerate(block.findall(f"./{side}/data")):
            atom = (element.findtext("atom") or "").strip()
            adp = _score(element.findtext("adp"))
            cid = atom_cid(atom)
            if cid is None or math.isnan(adp):
                continue
            beyond = (limit - adp if side == "low" else adp - limit) if not math.isnan(limit) else 0.0
            label = f"B {side}: {atom} (B {adp:.1f}" + ("" if math.isnan(limit) else f", limit {limit:.1f}") + ")"
            found.append({"kind": kind, "index": index, "label": label,
                          "residues": [residue_cid(atom)], "atoms_cids": [cid], "score": beyond})
    return sorted(found, key=lambda p: -p["score"])


def places(root: ET.Element) -> List[Dict]:
    """Every place to look at: Servalcat's geometry outliers, MolProbity's,
    and the ADP analysis's."""
    return outliers(root) + molprobity_places(root) + adp_places(root)


def find(root: ET.Element, kind: str, index: int) -> Optional[Dict]:
    for outlier in places(root):
        if outlier["kind"] == kind and outlier["index"] == int(index):
            return outlier
    return None


def _place_cids(outlier: Dict):
    """(residue CIDs, atom CIDs) of a place: from Servalcat's atom labels, or
    as MolProbity and the ADP analysis places carry them."""
    if "residues" in outlier:
        return [c for c in outlier["residues"] if c], [c for c in outlier["atoms_cids"] if c]
    return ([c for c in (residue_cid(a) for a in outlier["atoms"]) if c],
            [c for c in (atom_cid(a) for a in outlier["atoms"]) if c])


def outlier_scene(outlier: Dict, file_ids: Dict[str, int], project_uuid: str) -> Dict:
    """A scene of the refined model and its maps, centred and clipped on the
    outlier's atoms. ``file_ids`` maps XYZOUT (required), FPHIOUT and
    DIFFPHIOUT to database file ids; a viewer that already holds those files
    reuses them, so moving between outliers reloads nothing."""
    residues, atoms = _place_cids(outlier)
    # A residue outlier names no atoms: look at the residue
    cids = atoms or residues
    if not cids:
        raise ValueError(f"nothing to look at in {outlier.get('atoms') or outlier.get('label')}")
    selection = "||".join(dict.fromkeys(cids))
    files = [{"name": "model", "fileId": file_ids["XYZOUT"], "projectId": project_uuid}]
    maps = []
    for param, name, difference in (("FPHIOUT", "2mFo-DFc", False), ("DIFFPHIOUT", "mFo-DFc", True)):
        if file_ids.get(param):
            files.append({"name": param.lower(), "kind": "mtz", "fileId": file_ids[param],
                          "projectId": project_uuid})
            entry = {"name": name, "file": param.lower(), "columns": {"F": "F", "PHI": "PHI"},
                     "radius": 12}
            if difference:
                entry["isDifference"] = True
            maps.append(entry)
    scene = {
        "scene": describe(outlier),
        "version": 1,
        "files": files,
        "elements": [{
            "file": "model",
            "representations": [{"style": "CBs"}] + _markers(residues, atoms),
        }],
        "maps": maps,
        # Moorhen's zoom of 1 shows a crystal; 0.3 shows a residue and its neighbours
        "view": {"centre": {"file": "model", "selection": selection},
                 "zoom": 0.3,
                 "slab": {"file": "model", "selection": selection, "pad": 6}},
    }
    if maps:
        # Moorhen's refinement tools act on the active map; without one, its
        # context menu offers them greyed out. The 2mFo-DFc is the one to fit to.
        scene["activeMap"] = maps[0]["name"]
    return scene


def read_program_xml(work_directory) -> Optional[ET.Element]:
    path = Path(work_directory) / "program.xml"
    if not path.is_file():
        return None
    try:
        return ET.parse(path).getroot()
    except ET.ParseError:
        return None
