"""
Fragment-campaign summary-scene service.

Builds a Moorhen *scene* (see ``client/renderer/types/moorhen-scene.md``)
that overlays every discovered fragment hit on the campaign's parent
reference structure: the parent drawn as a ribbon, each hit dataset's
ligand drawn as sticks, scoped to its own restraint dictionary.

Two jobs this module does that the scene format then exploits:

  * **Hit detection** — a dataset counts as a hit when its refined
    coordinates actually contain the ligand its restraint dictionary
    describes. We read the dictionary's comp_id(s) with gemmi and look
    for those residue names in the coordinates. That single check both
    decides "is this a hit?" and tells us the exact CID to render
    (``//*/(CODE)``). When a member has no dictionary we fall back to
    the conventional placeholder codes (LIG/DRG/UNL).

  * **Dictionary scoping** — the scene lists each hit's dictionary in
    the element's ``dictionaries:`` block, so the resolver loads it and
    re-associates it with *that* molecule's molNo. Two datasets whose
    ligands share a code but differ in chemistry stay correct.

Like ``pandda_export``, this module is deliberately Django-light: it
takes a ``ProjectGroup`` and reads the related ``Job``/``File`` tables,
but does no request/response handling, so it is callable from a
synchronous ViewSet today and a background worker tomorrow.
"""
import logging
from pathlib import Path
from typing import Optional

import gemmi

from ..db import models
from . import pandda_export

logger = logging.getLogger(f"ccp4i2:{__name__}")


# Refinement tasks whose latest finished job carries the ligand-bound
# coordinates. Fragment campaigns here refine with servalcat (via
# servalcat_pipe), so that must come first / be present — otherwise the
# latest finished job falls back to an early apo refmac/dimple model and the
# soaked ligand is missed. order doesn't affect selection (we take the
# highest-id finished job among these), but servalcat_pipe is the canonical
# campaign refinement and its XYZOUT is the final ligand-bound model.
REFINE_TASK_NAMES = (
    "servalcat_pipe",
    "prosmart_refmac",
    "refmac",
    "i2Refmac",
    "i2Dimple",
    "dimple",
)

# Crystallisation additives, cryoprotectants and ions that commonly appear
# as HET residues but are NOT a fragment hit. A refmac LIBOUT can also carry
# standard monomers; those are filtered via gemmi's residue tabulation
# (amino acid / nucleic acid / water) rather than enumerated here.
COMMON_NON_LIGANDS = frozenset({
    "HOH", "WAT", "DOD",
    "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "P6G", "PG0", "MPD", "BME",
    "SO4", "PO4", "ACT", "ACY", "FMT", "CIT", "FLC", "TLA", "MES", "EPE", "TRS",
    "IPA", "DMS", "NH4", "NO3", "CAC",
    "NA", "K", "MG", "CA", "ZN", "MN", "FE", "FE2", "CL", "BR", "IOD",
    "CD", "NI", "CO", "CU", "CU1", "HG",
})

# Parent ribbon colour — a muted grey so the coloured fragment sticks
# read clearly against it.
PARENT_RIBBON_COLOUR = "#b0bec5"

DICT_FILE_TYPE = "application/refmac-dictionary"
COORD_FILE_TYPES = ("chemical/x-pdb", "chemical/x-cif", "chemical/x-mmcif")


# --------------------------------------------------------------------------
# Ligand detection (pure; gemmi only, no DB) -- unit-testable
# --------------------------------------------------------------------------

def _clean_code(value: str) -> str:
    """Normalise a comp_id / residue name: unquote, strip, upper-case."""
    return gemmi.cif.as_string(value).strip().upper()


def dictionary_comp_ids(dict_path) -> set:
    """Return the set of comp_ids declared by a restraint CIF.

    Reads ``_chem_comp.id`` (pair or loop) and ``_chem_comp_atom.comp_id``
    across every block, plus the ``comp_<X>`` block-name convention.
    Mirrors the extraction in ``FileViewSet.molblock`` but collects *all*
    codes (a dict file may declare several monomers in one shot).
    """
    ids: set = set()
    try:
        doc = gemmi.cif.read_file(str(dict_path))
    except Exception as exc:  # noqa: BLE001 - malformed dict shouldn't kill the scene
        logger.warning("Could not read dictionary %s: %s", dict_path, exc)
        return ids
    for block in doc:
        single = block.find_value("_chem_comp.id")
        if single:
            ids.add(_clean_code(single))
        for cid in block.find_loop("_chem_comp.id"):
            ids.add(_clean_code(cid))
        for cid in block.find_loop("_chem_comp_atom.comp_id"):
            ids.add(_clean_code(cid))
        name = block.name or ""
        if name.lower().startswith("comp_"):
            token = name[len("comp_"):]
            if 1 <= len(token) <= 5:
                ids.add(token.upper())
    # Drop empties and the dict's catch-all list block.
    ids.discard("")
    ids.discard("LIST")
    return ids


def _is_fragment_code(code: str) -> bool:
    """True if ``code`` looks like a soaked fragment rather than a standard
    biopolymer residue, water, ion, or crystallisation additive.

    A refmac LIBOUT may declare standard monomers alongside the ligand, so
    "the dict mentions it and the coords contain it" is not sufficient on its
    own — we still exclude anything gemmi recognises as an amino acid /
    nucleic acid / water, plus a curated additive list.
    """
    info = gemmi.find_tabulated_residue(code)
    if info is not None and (
        info.is_amino_acid() or info.is_nucleic_acid() or info.is_water()
    ):
        return False
    return code not in COMMON_NON_LIGANDS


def coordinate_residue_names(coord_path) -> set:
    """Return the set of residue names present in a coordinate file."""
    names: set = set()
    st = gemmi.read_structure(str(coord_path))
    for model in st:
        for chain in model:
            for res in chain:
                names.add(res.name.strip().upper())
    return names


def detect_ligands(coord_path, dict_path: Optional[object] = None) -> list:
    """Return the ligand codes a dataset should be judged a hit on.

    A dataset is a hit if its refined coordinates contain a *fragment-like*
    HET residue: anything that isn't a standard amino acid / nucleic acid /
    water / common crystallisation additive (see ``_is_fragment_code``). This
    catches real three-letter ligand codes (e.g. ``NUT`` for Nutlin), not just
    the conventional placeholders, and works even when no restraint dictionary
    is present — fragment campaigns refined with servalcat often have none.

    When a dictionary *is* available, the dictionary-confirmed subset is
    preferred: it pins the comp_id we render and gives the scoped dictionary
    that lets same-named ligands with different chemistry coexist. If the
    dictionary covers none of the candidates, we still fall back to the
    coordinate-derived candidates.

    Returns a sorted list (possibly empty -> not a hit).
    """
    try:
        coord_names = coordinate_residue_names(coord_path)
    except Exception as exc:  # noqa: BLE001 - unreadable coords -> not a hit
        logger.warning("Could not read coordinates %s: %s", coord_path, exc)
        return []

    candidates = {name for name in coord_names if _is_fragment_code(name)}
    if not candidates:
        return []

    if dict_path is not None:
        confirmed = sorted(candidates & dictionary_comp_ids(dict_path))
        if confirmed:
            return confirmed

    return sorted(candidates)


# --------------------------------------------------------------------------
# Scene assembly (touches the DB)
# --------------------------------------------------------------------------

def _parent_coord_file(group):
    """The parent project's reference XYZOUT coordinate File, or None."""
    parent_membership = group.memberships.filter(
        type=models.ProjectGroupMembership.MembershipType.PARENT
    ).select_related("project").first()
    if not parent_membership:
        return None, None
    parent_project = parent_membership.project
    coord_file = (
        models.File.objects.filter(
            job__project=parent_project,
            job_param_name="XYZOUT",
            type__name="chemical/x-pdb",
        )
        .order_by("-id")
        .first()
    )
    return parent_project, coord_file


def _member_coord_file(job):
    """Pick a refined-coordinate File for a refinement job.

    Prefer the XYZOUT PDB (the canonical refined-model output); fall back
    to any coordinate-typed output of the job.
    """
    xyzout = (
        models.File.objects.filter(job=job, job_param_name="XYZOUT")
        .order_by("-id")
    )
    pdb = xyzout.filter(type__name="chemical/x-pdb").first()
    if pdb:
        return pdb
    any_xyz = xyzout.filter(type__name__in=COORD_FILE_TYPES).first()
    if any_xyz:
        return any_xyz
    return (
        models.File.objects.filter(job=job, type__name__in=COORD_FILE_TYPES)
        .order_by("-id")
        .first()
    )


def _member_dict_file(project, refine_job):
    """Find the restraint dictionary File for a member.

    First the dictionary emitted by the refinement job itself (guaranteed
    to pair with its coordinates); otherwise the dictionary File of the
    member's latest finished acedrg job.
    """
    own = (
        models.File.objects.filter(job=refine_job, type__name=DICT_FILE_TYPE)
        .order_by("-id")
        .first()
    )
    if own:
        return own
    acedrg_job = pandda_export._latest_finished_job(
        project, pandda_export.ACEDRG_TASK_NAMES
    )
    if not acedrg_job:
        return None
    return (
        models.File.objects.filter(job=acedrg_job, type__name=DICT_FILE_TYPE)
        .order_by("-id")
        .first()
    )


def _safe_name(raw: str, used: set) -> str:
    """Slug a project name into a unique, scene-safe file identifier."""
    base = "".join(c if c.isalnum() else "_" for c in (raw or "")).strip("_")
    base = base or "dataset"
    name = base
    n = 2
    while name in used:
        name = f"{base}_{n}"
        n += 1
    used.add(name)
    return name


def build_summary_scene(group) -> dict:
    """Build the fragment-campaign summary scene for ``group``.

    Returns ``{"scene": <MoorhenScene>, "stats": {...}}``. The scene is a
    plain JSON-serialisable dict matching the MoorhenScene TypeScript
    shape; the frontend applies it directly via the scene resolver.
    """
    files: list = []
    elements: list = []
    used_names: set = set()
    stats = {
        "members_total": 0,
        "hits": 0,
        "skipped": [],            # [{project, reason}]
        "parent_present": False,
    }

    # -- Parent reference: ribbon -----------------------------------------
    parent_project, parent_coord = _parent_coord_file(group)
    if parent_coord is not None:
        ref_name = _safe_name("reference", used_names)
        files.append(
            {
                "name": ref_name,
                "kind": "coordinates",
                "fileId": parent_coord.id,
                "projectId": str(parent_project.uuid),
            }
        )
        elements.append(
            {
                "file": ref_name,
                "representations": [
                    {
                        "style": "CRs",
                        "selection": "/*/*/*/*",
                        "colour": PARENT_RIBBON_COLOUR,
                    }
                ],
            }
        )
        stats["parent_present"] = True
    else:
        ref_name = None

    # -- Member hits: ligand sticks, each with its scoped dictionary ------
    member_memberships = group.memberships.filter(
        type=models.ProjectGroupMembership.MembershipType.MEMBER
    ).select_related("project")

    for membership in member_memberships:
        project = membership.project
        stats["members_total"] += 1

        refine_job = pandda_export._latest_finished_job(project, REFINE_TASK_NAMES)
        if not refine_job:
            stats["skipped"].append(
                {"project": project.name, "reason": "no finished refinement job"}
            )
            continue

        coord_file = _member_coord_file(refine_job)
        if coord_file is None or not coord_file.path.exists():
            stats["skipped"].append(
                {"project": project.name, "reason": "no refined coordinate file"}
            )
            continue

        dict_file = _member_dict_file(project, refine_job)
        dict_path = None
        if dict_file is not None and dict_file.path.exists():
            dict_path = dict_file.path
        else:
            # Disk-only acedrg dictionary (no File record): inline its text.
            acedrg_job = pandda_export._latest_finished_job(
                project, pandda_export.ACEDRG_TASK_NAMES
            )
            disk_dict = pandda_export._find_dictionary_cif(acedrg_job)
            if disk_dict is not None:
                dict_path = disk_dict

        codes = detect_ligands(coord_file.path, dict_path)
        if not codes:
            stats["skipped"].append(
                {"project": project.name, "reason": "no ligand in refined model"}
            )
            continue

        ds_name = _safe_name(project.name, used_names)
        files.append(
            {
                "name": ds_name,
                "kind": "coordinates",
                "fileId": coord_file.id,
                "projectId": str(project.uuid),
            }
        )

        dict_ref_name = None
        if dict_file is not None and dict_file.path.exists():
            dict_ref_name = _safe_name(f"{ds_name}_dict", used_names)
            files.append(
                {
                    "name": dict_ref_name,
                    "kind": "dictionary",
                    "fileId": dict_file.id,
                    "projectId": str(project.uuid),
                }
            )
        elif dict_path is not None:
            try:
                cif_text = Path(dict_path).read_text()
                dict_ref_name = _safe_name(f"{ds_name}_dict", used_names)
                files.append(
                    {
                        "name": dict_ref_name,
                        "kind": "dictionary",
                        "cifText": cif_text,
                    }
                )
            except OSError as exc:
                logger.warning("Could not inline dictionary %s: %s", dict_path, exc)

        selection = "||".join(f"//*/({code})" for code in codes)
        element = {
            "file": ds_name,
            "representations": [{"style": "CBs", "selection": selection}],
        }
        if dict_ref_name:
            element["dictionaries"] = [dict_ref_name]
        elements.append(element)
        stats["hits"] += 1

    scene = {
        "scene": f"{group.name} - fragment summary",
        "version": 1,
        "authoredIn": {"projectName": group.name},
        "files": files,
        "elements": elements,
        "resolver": {"onMissingResidues": "clamp-and-log"},
    }

    view = _first_site_view(group)
    if view:
        scene["view"] = view

    return {"scene": scene, "stats": stats}


def _first_site_view(group) -> Optional[dict]:
    """Camera from the campaign's first saved binding site, if any."""
    sites = getattr(group, "sites", None) or []
    if not sites:
        return None
    site = sites[0]
    view: dict = {}
    if site.get("origin"):
        view["origin"] = site["origin"]
    if site.get("quat"):
        view["quat"] = site["quat"]
    if site.get("zoom") is not None:
        view["zoom"] = site["zoom"]
    return view or None
