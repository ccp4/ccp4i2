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
    A/ZN 206/ZN (Z +16.2)'. Atoms of one residue name the residue once."""
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


def find(root: ET.Element, kind: str, index: int) -> Optional[Dict]:
    for outlier in outliers(root):
        if outlier["kind"] == kind and outlier["index"] == int(index):
            return outlier
    return None


def outlier_scene(outlier: Dict, file_ids: Dict[str, int], project_uuid: str) -> Dict:
    """A scene of the refined model and its maps, centred and clipped on the
    outlier's atoms. ``file_ids`` maps XYZOUT (required), FPHIOUT and
    DIFFPHIOUT to database file ids; a viewer that already holds those files
    reuses them, so moving between outliers reloads nothing."""
    cids = [cid for cid in (atom_cid(a) for a in outlier["atoms"]) if cid]
    if not cids:
        raise ValueError(f"no atoms to look at in {outlier['atoms']}")
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
    return {
        "scene": describe(outlier),
        "version": 1,
        "files": files,
        "elements": [{
            "file": "model",
            "representations": [
                {"style": "CBs"},
                # The outlier's own atoms, drawn over the bonds so they stand out
                {"style": "VdwSpheres", "selection": selection, "colour": "#ff00ff",
                 "geometry": {"vdwScale": 0.35}},
            ],
        }],
        "maps": maps,
        "view": {"centre": {"file": "model", "selection": selection},
                 "slab": {"file": "model", "selection": selection, "pad": 6}},
    }


def read_program_xml(work_directory) -> Optional[ET.Element]:
    path = Path(work_directory) / "program.xml"
    if not path.is_file():
        return None
    try:
        return ET.parse(path).getroot()
    except ET.ParseError:
        return None
