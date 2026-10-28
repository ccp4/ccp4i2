"""
Fragment-campaign scene service: the whole-campaign summary and one site.

Builds a Moorhen *scene* (see ``client/renderer/types/moorhen-scene.md``)
that overlays fragment hits on the campaign's parent reference structure:
the parent drawn as a ribbon, each hit dataset's ligand drawn as sticks,
scoped to its own restraint dictionary. ``build_summary_scene`` takes every
hit the coordinates reveal; ``build_site_scene`` takes the datasets judged
a hit at one site, and fits them on that site's pocket
(``docs/campaign-site-scene-design.md``).

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
from types import SimpleNamespace
from typing import Optional

import gemmi

from ..db import models
from . import campaign_events, pandda_export, superposition

logger = logging.getLogger(f"ccp4i2:{__name__}")


# Which job carries a member's ligand-bound coordinates is the campaign's one
# definition of a dataset's current model -- campaign_matrix.REFINEMENT_TASKS
# and choose_current_model, shared with the overview's site matrix -- so the
# scene and the overview cannot disagree about which model a dataset has.
# (This used to be its own REFINE_TASK_NAMES list, taking the highest-id
# finished job of refinement and dimple tasks alike at any level.)

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
# The pocket sticks are a near neighbour of the ribbon grey, so the context
# reads as context and the coloured fragments still carry the picture.
POCKET_STICK_COLOUR = "#90a4ae"
# One colour per hit at a site, cycled, so the hits are tellable apart.
HIT_COLOURS = (
    "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e",
    "#9467bd", "#8c564b", "#e377c2", "#17becf",
)
# Unclear verdicts, when asked for, share one muted colour: a site view is a
# claim about what is there, and the doubtful must not look like the confident.
UNCLEAR_COLOUR = "#d7ccc8"

# Residues with any atom this close to the site origin are drawn as the
# pocket. This radius decides only what is shown as sticks -- never which
# datasets are hits -- so a wrong value shows a residue too many or too few,
# visibly, rather than losing a ligand silently.
ENVIRONMENT_RADIUS = 8.0

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


def residue_names(structure: gemmi.Structure) -> set:
    """Return the set of residue names present in a parsed structure."""
    names: set = set()
    for model in structure:
        for chain in model:
            for res in chain:
                names.add(res.name.strip().upper())
    return names


def coordinate_residue_names(coord_path) -> set:
    """Return the set of residue names present in a coordinate file."""
    return residue_names(gemmi.read_structure(str(coord_path)))


def detect_ligands(coord_path, dict_path: Optional[object] = None,
                   structure: Optional[gemmi.Structure] = None) -> list:
    """Return the ligand codes a dataset should be judged a hit on.

    ``structure`` is the already-parsed ``coord_path``, when the caller has
    it. A scene builder reads each member's coordinates for this check and
    again to fit them onto the reference; passing the structure in makes
    that one parse rather than two, which on a campaign of hundreds of
    datasets is the difference between a page and a timeout.

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
        coord_names = (
            residue_names(structure) if structure is not None
            else coordinate_residue_names(coord_path)
        )
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
# Site geometry (pure; gemmi only, no DB) -- unit-testable
# --------------------------------------------------------------------------

def _residue_cid(chain_name: str, res: gemmi.Residue) -> str:
    """One residue as a Coot CID.

    mmdb's residue id puts the insertion code after a dot (``ParseResID``):
    ``//A/45.A``, not ``//A/45A``, which would parse as no residue at all.
    """
    cid = f"//{chain_name}/{res.seqid.num}"
    if res.seqid.icode != " ":
        cid += f".{res.seqid.icode}"
    return cid


def site_position(site) -> tuple:
    """The site's centre as a REAL-SPACE coordinate.

    ``CampaignSite.origin`` is not one. It is Moorhen's view origin, which is
    the *negation* of the point at the centre of the screen: Moorhen's own
    code negates it whenever it hands the value to something that wants real
    coordinates (see ``GetMonomer``, which passes
    ``origin.map(coord => -coord)`` to ``get_monomer_and_position_at``). The
    site save/restore path stores and restores it raw, so the two negations
    cancel and the camera round-trips correctly -- which is why the sign error
    stayed invisible until something treated the stored value as a position.

    Measured on the BAZ2B demo campaign: the stored origin is 40.9 A from the
    nearest CA of the reference with zero CAs inside 30 A, while its negation
    is 3.75 A away with 14 CAs inside 8 A. The negation is the pocket.

    Anything selecting atoms or centring a fit must go through here.
    ``view.origin`` in the scene must NOT -- that is handed back to Moorhen,
    which wants its own convention.
    """
    return (-site.origin_x, -site.origin_y, -site.origin_z)


def pocket_residue_cids(structure: gemmi.Structure, origin,
                        radius: float = ENVIRONMENT_RADIUS) -> list:
    """CIDs of the residues with any atom within ``radius`` of ``origin``.

    One entry per residue however many of its atoms are in range, in
    chain/sequence order, and an origin in empty solvent gives an empty list
    -- which a caller must turn into *no* representation, because an empty
    ``selection:`` draws the whole molecule.

    Waters and fragment-like residues are left out: the pocket is what the
    ligands bind *to*. When the exemplar is itself a hit (no parent), its own
    ligand is already drawn in that hit's colour, and drawing it a second
    time in the pocket grey would fight it for the same atoms.

    A plain distance loop rather than ``gemmi.NeighborSearch``: the search
    would also return symmetry images, and a residue that lines the pocket
    only as a symmetry mate would be drawn where it sits in the asymmetric
    unit, nowhere near the site.
    """
    if len(structure) == 0:
        return []
    centre = gemmi.Position(*origin)
    found: dict = {}
    for chain in structure[0]:
        for res in chain:
            name = res.name.strip().upper()
            info = gemmi.find_tabulated_residue(name)
            if info is not None and info.is_water():
                continue
            if _is_fragment_code(name):
                continue
            if any(atom.pos.dist(centre) <= radius for atom in res):
                key = (chain.name, res.seqid.num, res.seqid.icode)
                found[key] = _residue_cid(chain.name, res)
    return [found[key] for key in sorted(found)]


def fragment_centroids(structure: gemmi.Structure) -> list:
    """``[(code, centroid)]`` for every fragment-like residue in the first model."""
    out: list = []
    if len(structure) == 0:
        return out
    for chain in structure[0]:
        for res in chain:
            code = res.name.strip().upper()
            if len(res) == 0 or not _is_fragment_code(code):
                continue
            n = len(res)
            centroid = gemmi.Position(
                sum(a.pos.x for a in res) / n,
                sum(a.pos.y for a in res) / n,
                sum(a.pos.z for a in res) / n,
            )
            out.append((code, centroid))
    return out


def nearest_fragment_distance(structure: gemmi.Structure, origin,
                              transform: Optional[gemmi.Transform] = None
                              ) -> Optional[float]:
    """Distance from ``origin`` to the closest fragment-like residue's centroid.

    ``None`` -- not ``inf`` -- when the structure has no such residue. With a
    ``transform`` the centroids are moved first, so the distance is measured
    in the frame the scene draws the dataset in: after superposition a
    fragment still far from the site is a real anomaly, not a frame artefact.

    This is a diagnostic, never a filter. A hit verdict whose nearest fragment
    is 40 A away is either recorded against density nobody has modelled yet
    or a numbering mismatch that defeated the fit; either way the user should
    see the number, not lose the ligand.
    """
    centre = gemmi.Position(*origin)
    best: Optional[float] = None
    for _code, centroid in fragment_centroids(structure):
        if transform is not None:
            centroid = gemmi.Position(transform.apply(centroid))
        dist = centroid.dist(centre)
        if best is None or dist < best:
            best = dist
    return best


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
    # Match what the member path accepts. This used to demand
    # job_param_name="XYZOUT" AND chemical/x-pdb, which silently excluded a
    # reference imported as mmCIF or held under any other parameter -- and a
    # missing parent file costs the whole scene its reference structure.
    candidates = models.File.objects.filter(
        job__project=parent_project, type__name__in=COORD_FILE_TYPES
    )
    for queryset in (
        candidates.filter(job_param_name="XYZOUT", type__name="chemical/x-pdb"),
        candidates.filter(job_param_name="XYZOUT"),
        candidates,
    ):
        for coord_file in queryset.order_by("-id"):
            # A File row whose file is gone would put a fileId in the scene
            # that 404s on fetch, which looks exactly like no ribbon at all.
            if coord_file.path.exists():
                return parent_project, coord_file
    return parent_project, None


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


def _ribbon() -> dict:
    return {"style": "CRs", "selection": "/*/*/*/*", "colour": PARENT_RIBBON_COLOUR}


def _resolve_hit(project, used_names: set):
    """Turn one member into the scene pieces that draw its ligand.

    Returns ``(hit, None)`` or ``(None, reason)``. ``hit`` carries the
    ``files[]`` entries (coordinates, and the dictionary by fileId or as
    inlined CIF text), the element with its ``//*/(CODE)`` sticks, and the
    coordinate path later steps fit and measure on. ``sticks`` is the ligand
    representation itself, so a caller can colour it without knowing where
    it sits once a ribbon has been pushed in front of it.

    The summary scene and the site scene both go through here. They differ
    in which members they include, what the exemplar draws and where the
    camera goes -- never in how a hit is resolved, so that a fix to the hit
    rule cannot land in one and not the other.
    """
    refine_job = campaign_events.current_model_job(project)
    if not refine_job:
        return None, "no finished refinement job"

    coord_file = _member_coord_file(refine_job)
    if coord_file is None or not coord_file.path.exists():
        return None, "no refined coordinate file"

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

    # Parsed once here and carried on the hit: the fit onto the reference and
    # the site distances read the same file, and a campaign of hundreds of
    # datasets cannot afford to parse each of them three times.
    try:
        structure = gemmi.read_structure(str(coord_file.path))
    except Exception as exc:  # noqa: BLE001 - unreadable coords -> not a hit
        logger.warning("Could not read coordinates %s: %s", coord_file.path, exc)
        return None, "no ligand in refined model"

    codes = detect_ligands(coord_file.path, dict_path, structure=structure)
    if not codes:
        return None, "no ligand in refined model"

    ds_name = _safe_name(project.name, used_names)
    files = [
        {
            "name": ds_name,
            "kind": "coordinates",
            "fileId": coord_file.id,
            "projectId": str(project.uuid),
        }
    ]

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

    sticks = {
        "style": "CBs",
        "selection": "||".join(f"//*/({code})" for code in codes),
    }
    element = {"file": ds_name, "representations": [sticks]}
    if dict_ref_name:
        element["dictionaries"] = [dict_ref_name]

    return (
        {
            "project": project,
            "name": ds_name,
            "files": files,
            "element": element,
            "sticks": sticks,
            "coord_path": coord_file.path,
            "structure": structure,
            "codes": codes,
        },
        None,
    )


def _parent_reference(group, used_names: set, files: list, elements: list,
                      stats: dict):
    """Draw the parent's coordinates as the ribbon, if it has any.

    Returns ``(name, path, element)`` -- all None when there is no parent
    file, in which case a builder promotes a hit (``_promote_hit``). Claimed
    first so the file is called ``reference`` whatever the members are named.
    """
    parent_project, parent_coord = _parent_coord_file(group)
    if parent_coord is None:
        return None, None, None
    ref_name = _safe_name("reference", used_names)
    files.append(
        {
            "name": ref_name,
            "kind": "coordinates",
            "fileId": parent_coord.id,
            "projectId": str(parent_project.uuid),
        }
    )
    element = {"file": ref_name, "representations": [_ribbon()]}
    elements.append(element)
    stats["parent_present"] = True
    stats["reference"] = {"project": parent_project.name, "is_parent": True}
    return ref_name, parent_coord.path, element


def _promote_hit(hits: list, stats: dict):
    """No parent coordinates: the first hit carries the ribbon.

    A campaign whose parent project was never populated (or whose reference
    model has not been imported yet) otherwise renders as fragments floating
    in space with nothing to place them against, and says nothing about why.
    Every member is a full structure, so one of them can carry the ribbon --
    drawn from the SAME files[] entry that already carries its ligand, so
    this costs no extra download and no second copy in the viewer.

    The frame is then that dataset's rather than the campaign's. Members are
    near-isomorphous, so the picture is right; ``stats["reference"]`` records
    whose frame it is, because a caller that cannot tell cannot explain it.

    Returns ``(name, path, element, movers)``: the promoted hit is the frame,
    so ``movers`` is every other hit, the ones still to be fitted onto it.
    """
    if not hits:
        return None, None, None, hits
    exemplar = hits[0]
    exemplar["element"]["representations"].insert(0, _ribbon())
    stats["reference"] = {"project": exemplar["project"].name, "is_parent": False}
    return exemplar["name"], exemplar["coord_path"], exemplar["element"], hits[1:]


def _hit_structure(hit: dict) -> gemmi.Structure:
    """The hit's parsed coordinates, read on first use and kept on the hit."""
    structure = hit.get("structure")
    if structure is None:
        structure = gemmi.read_structure(str(hit["coord_path"]))
        hit["structure"] = structure
    return structure


def _superpose(ref_name: str, ref_path, movers: list, stats: dict,
               centre=None, min_radius=None, reference=None) -> list:
    """Fit every mover onto the reference; return the ``superpose[]`` entries.

    Each fit is recorded in ``stats["superpose"]`` with the evidence it
    rests on, and kept on the hit as ``fit`` so later steps can measure in
    the frame the scene draws. A dataset that cannot be fitted is still
    drawn, untransformed, and stats say so.

    With a ``centre`` the fit is local to the site. When the count gate
    cannot be met inside the radius cap, a global fit stands in and is
    named as such: past the cap a fit is no longer local in any meaningful
    sense, and drawing the dataset in its own frame would be worse than a
    global fit that says what it is. The scene shows the difference too --
    a global fit's provenance carries no ``radius``.

    ``reference`` is the parsed ``ref_path`` when the caller already has
    it; otherwise it is read here, once, rather than once per mover. Each
    mover's structure comes from the hit (``_hit_structure``), so a member
    is parsed a single time however many steps look at it.
    """
    entries: list = []
    if reference is None:
        try:
            reference = gemmi.read_structure(str(ref_path))
        except Exception as exc:  # noqa: BLE001 - an unreadable reference fails every fit
            reason = f"could not read coordinates: {exc}"
            for hit in movers:
                fit = superposition.FitResult(ok=False, reason=reason)
                hit["fit"] = fit
                stats["superpose"].append({
                    "project": hit["project"].name,
                    "ok": False, "atoms": fit.atoms, "radius": fit.radius,
                    "rmsd": fit.rmsd, "reason": reason,
                })
            return entries

    for hit in movers:
        try:
            moving = _hit_structure(hit)
        except Exception as exc:  # noqa: BLE001 - unreadable file is a failed fit
            fit = superposition.FitResult(
                ok=False, reason=f"could not read coordinates: {exc}"
            )
        else:
            fit = superposition.fit_structures(
                reference, moving, centre=centre, min_radius=min_radius
            )
        record = {
            "project": hit["project"].name,
            "ok": fit.ok,
            "atoms": fit.atoms,
            "radius": fit.radius,
            "rmsd": fit.rmsd,
            "reason": fit.reason,
        }
        if not fit.ok and centre is not None and hit.get("structure") is not None:
            fallback = superposition.fit_structures(reference, hit["structure"])
            if fallback.ok:
                record.update(
                    ok=True,
                    atoms=fallback.atoms,
                    radius=None,
                    rmsd=fallback.rmsd,
                    fallback="global",
                )
                fit = fallback
        stats["superpose"].append(record)
        hit["fit"] = fit
        if fit.ok:
            entries.append(fit.scene_entry(hit["name"], ref_name))
    return entries


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
        # Which structure ended up as the reference ribbon, and whether it was
        # the campaign's parent or a stand-in. A caller that cannot tell those
        # apart cannot explain a scene whose frame is one dataset's.
        "reference": None,        # {"project": name, "is_parent": bool}
        # One record per drawn hit: whether it was fitted onto the reference
        # and on what evidence (CA count, RMSD), or why it could not be. A
        # dataset left unfitted is still drawn, in its own frame; this is
        # where a reader learns that.
        "superpose": [],          # [{project, ok, atoms, radius, rmsd, reason}]
    }

    ref_name, ref_path, _ref_element = _parent_reference(
        group, used_names, files, elements, stats
    )

    # -- Member hits: ligand sticks, each with its scoped dictionary ------
    member_memberships = group.memberships.filter(
        type=models.ProjectGroupMembership.MembershipType.MEMBER
    ).select_related("project")

    hits: list = []
    for membership in member_memberships:
        project = membership.project
        stats["members_total"] += 1
        hit, reason = _resolve_hit(project, used_names)
        if hit is None:
            stats["skipped"].append({"project": project.name, "reason": reason})
            continue
        files.extend(hit["files"])
        elements.append(hit["element"])
        stats["hits"] += 1
        hits.append(hit)

    movers = hits
    if ref_name is None:
        ref_name, ref_path, _ref_element, movers = _promote_hit(hits, stats)

    # -- Put every hit in the reference's frame ----------------------------
    #
    # Members ought to share a frame already (molecular-replaced from one
    # reference), but in practice they do not: PDB imports carry their own
    # origin choice, and even in-campaign refinements drift. Drawn as they
    # come, the ribbons fan out by an Angstrom or so and binding events that
    # are probably equivalent look different. A whole-campaign summary has
    # no single pocket to align on, so the fit is global, on every CA the
    # two structures share; the site scene fits locally.
    superpose: list = []
    if ref_path is not None:
        superpose = _superpose(ref_name, ref_path, movers, stats)

    scene = {
        "scene": f"{group.name} - fragment summary",
        "version": 1,
        "authoredIn": {"projectName": group.name},
        "files": files,
        **({"superpose": superpose} if superpose else {}),
        "elements": elements,
        "resolver": {"onMissingResidues": "clamp-and-log"},
    }

    view = _first_site_view(group)
    if view:
        scene["view"] = view

    annotations, site_stats = _site_markers(group, member_memberships)
    if annotations:
        scene["annotations"] = annotations
    stats["sites"] = site_stats

    return {"scene": scene, "stats": stats}


# Markers over the binding sites of the summary scene: a sphere, which reads
# the same from whatever angle shows the site distribution best. A site with
# hits is warm, a site nobody has found anything at is grey; the banner says
# which.
MARKER_HIT_COLOUR = "#e07a5f"
MARKER_EMPTY_COLOUR = "#9e9e9e"
MARKER_RADIUS = 1.2        # Angstrom


def site_hit_counts(sites, member_ids: set) -> dict:
    """``{site id: number of current members judged a hit there}``.

    A member counts when its cell's effective verdict at the site, across
    reviewers (``lib.site_verdicts.cell_verdicts``), is ``hit``: two reviewers
    disagreeing is "to reconcile", not a hit. Only current members count, for
    the reason ``build_sites_scene`` gives: an evaluation outlives the
    membership it was recorded under.
    """
    from .site_verdicts import cell_verdicts
    counts = {}
    for site in sites:
        n = 0
        for (project_id, _), cell in cell_verdicts(site.evaluations.all()).items():
            if project_id in member_ids and cell["verdict"] == "hit":
                n += 1
        counts[site.id] = n
    return counts


def _marker_text(name: str, hits: int) -> str:
    if hits == 0:
        return f"{name}: no hits"
    return f"{name}: {hits} hit" + ("" if hits == 1 else "s")


def _site_markers(group, member_memberships) -> tuple:
    """A ``sphere`` annotation over every binding site of the campaign
    carrying its name and the number of members judged a hit there, plus
    ``[{id, name, hits}]`` for ``stats``.

    The marker sits at the site's real-space position (``site_position``: the
    stored origin is Moorhen's negated view origin). It is drawn in the
    reference frame the hits are fitted into, which is the frame the sites
    were saved in, since a site is saved from the campaign viewer looking at
    the parent or a stand-in for it.
    """
    from .campaign_site_order import ordered_sites
    sites = ordered_sites(group.site_set.all())
    if not sites:
        return [], []
    member_ids = set(member_memberships.values_list("project_id", flat=True))
    counts = site_hit_counts(sites, member_ids)
    annotations = []
    site_stats = []
    for site in sites:
        hits = counts.get(site.id, 0)
        annotations.append({
            "kind": "sphere",
            "at": [round(c, 3) for c in site_position(site)],
            "radius": MARKER_RADIUS,
            "colour": MARKER_HIT_COLOUR if hits else MARKER_EMPTY_COLOUR,
            "text": _marker_text(site.name, hits),
        })
        site_stats.append({"id": site.id, "name": site.name, "hits": hits})
    return annotations, site_stats


def build_site_scene(group, site, include_unclear: bool = False,
                     superpose: bool = True) -> dict:
    """Build the scene for one binding site of a campaign.

    The exemplar (the parent, else a hit) drawn once as a ribbon with the
    residues around the site origin as sticks, and every dataset judged a
    hit *at this site* drawn on top as its ligand's sticks, each fitted onto
    the exemplar locally, on the CAs around the site.

    Membership is the verdict and nothing else: a dataset is in iff its
    cell's effective verdict at this site (across reviewers,
    lib/site_verdicts.py) is ``hit``, plus ``unclear`` ones when asked for,
    in a muted colour. A cell whose reviewers disagree counts as unclear. ``empty`` and unevaluated members are out, and there
    is deliberately no fall-back to ligand detection when a site has no
    verdicts -- that would show every ligand in the campaign under the name
    of one site. A site nobody has evaluated gets an honest exemplar-only
    scene, and ``stats`` says so.

    Each hit's ligand is drawn as every copy of its code (``//*/(CODE)``),
    as the summary draws it; the camera and slab frame the site, and a copy
    bound in an adjacent subsite is an observation, not clutter (see
    ``docs/campaign-site-scene-design.md``, "Draw every copy").
    """
    return build_sites_scene(group, [site], include_unclear, superpose)


# How strongly a verdict claims a dataset for a site view: a dataset judged a
# hit at one chosen site and empty at another is a hit of the view.
_VERDICT_STRENGTH = {"hit": 3, "unclear": 2, "empty": 1}


def _union_cids(structure: gemmi.Structure, positions: list) -> list:
    """The pocket residues around any of ``positions``, in structure order."""
    wanted = set()
    for position in positions:
        wanted.update(pocket_residue_cids(structure, position))
    ordered = []
    for chain in structure[0]:
        for res in chain:
            cid = _residue_cid(chain.name, res)
            if cid in wanted and cid not in ordered:
                ordered.append(cid)
    return ordered


def build_sites_scene(group, sites: list, include_unclear: bool = False,
                      superpose: bool = True) -> dict:
    """Build one scene for several binding sites of a campaign, taken together.

    For a site and its subsites, or neighbouring sites worth seeing at once
    (docs/campaign-site-scene-design.md, "Several sites in one view"). With
    one site this is the site view of :func:`build_site_scene`, unchanged.
    With more:

      * a dataset is in when it is in for *any* chosen site, drawn once, at
        its strongest verdict across them (hit over unclear over empty), and
        ``stats["drawn"]`` names the sites that claim it;
      * the pocket is every residue around any chosen site, drawn once;
      * every hit is fitted on ONE region -- the CAs around the centroid of
        the chosen sites, in a sphere that starts wide enough to reach the
        farthest of them -- so the view has a single frame;
      * the camera looks at the centroid, oriented as the first site saved
        its view, and the slab covers the whole pocket.
    """
    files: list = []
    elements: list = []
    used_names: set = set()
    single = len(sites) == 1
    stats = {
        "site": {"id": sites[0].id, "name": sites[0].name} if single else None,
        "sites": [{"id": site.id, "name": site.name} for site in sites],
        "hits_claimed": 0,        # datasets judged a hit at a chosen site
        "hits_drawn": 0,
        "empty_verdicts": 0,      # somebody looked and found nothing
        "unclear_verdicts": 0,
        "unclear_drawn": 0,
        "skipped": [],            # [{project, reason, nearest}]
        "parent_present": False,
        "reference": None,        # {"project": name, "is_parent": bool}
        "pocket_residues": 0,
        # Per drawn dataset, the distance from the site origin to its nearest
        # fragment-like residue, measured after superposition. A hit whose
        # fragment is far from the site is worth a look; it is not hidden.
        "drawn": [],              # [{project, verdict, nearest, sites?}]
        "superpose": [],          # [{project, ok, atoms, radius, rmsd, reason, fallback?}]
    }
    if single:
        del stats["sites"]
    else:
        del stats["site"]

    ref_name, ref_path, ref_element = _parent_reference(
        group, used_names, files, elements, stats
    )

    # -- Membership: the verdicts at the chosen sites -----------------------
    #
    # Only current members. An evaluation outlives the membership it was
    # recorded under (it hangs off the project and the site, not the
    # membership row), so a dataset removed from the campaign would otherwise
    # still be drawn in its site views.
    Verdict = models.SiteEvaluation.Verdict
    member_ids = set(
        group.memberships.filter(
            type=models.ProjectGroupMembership.MembershipType.MEMBER
        ).values_list("project_id", flat=True)
    )
    # One entry per dataset, from what its cell says across reviewers
    # (lib/site_verdicts.py): one reviewer's hit beside another's empty is
    # "to reconcile", whose effective verdict is unclear, so it is drawn only
    # when unclears are asked for. Across several sites a dataset keeps its
    # strongest verdict and the names of the sites that gave it. Ordered by
    # project id so the fallback exemplar, when one is needed, is the
    # lowest-id hit: stable, arbitrary, and recorded as such in stats.
    from .site_verdicts import cell_verdicts
    claims: dict = {}             # project_id -> {"verdict", "sites": [site]}
    for site in sites:
        for (project_id, _), cell in cell_verdicts(site.evaluations.all()).items():
            if project_id not in member_ids:
                continue
            verdict = cell["verdict"]
            claim = claims.setdefault(project_id, {"verdict": verdict, "sites": []})
            if _VERDICT_STRENGTH.get(verdict, 0) > _VERDICT_STRENGTH.get(claim["verdict"], 0):
                claim["verdict"] = verdict
                claim["sites"] = []
            if verdict == claim["verdict"]:
                claim["sites"].append(site)
    projects = {p.id: p for p in models.Project.objects.filter(id__in=list(claims))}

    hits: list = []
    for project_id, claim in sorted(claims.items()):
        verdict = claim["verdict"]
        if verdict == Verdict.EMPTY:
            stats["empty_verdicts"] += 1
            continue
        if verdict == Verdict.UNCLEAR:
            stats["unclear_verdicts"] += 1
            if not include_unclear:
                continue
        else:
            stats["hits_claimed"] += 1

        project = projects[project_id]
        hit, reason = _resolve_hit(project, used_names)
        if hit is None:
            stats["skipped"].append(
                {"project": project.name, "reason": reason, "nearest": None}
            )
            continue
        hit["verdict"] = verdict
        hit["sites"] = claim["sites"]
        if verdict == Verdict.UNCLEAR:
            hit["sticks"]["colour"] = UNCLEAR_COLOUR
            stats["unclear_drawn"] += 1
        else:
            hit["sticks"]["colour"] = HIT_COLOURS[stats["hits_drawn"] % len(HIT_COLOURS)]
            stats["hits_drawn"] += 1
        files.extend(hit["files"])
        elements.append(hit["element"])
        hits.append(hit)

    movers = hits
    if ref_name is None:
        ref_name, ref_path, ref_element, movers = _promote_hit(hits, stats)

    positions = [site_position(site) for site in sites]
    centre = tuple(sum(p[i] for p in positions) / len(positions) for i in range(3))
    # How far the farthest chosen site lies from the centre: zero for one.
    spread = max(
        sum((p[i] - centre[i]) ** 2 for i in range(3)) ** 0.5 for p in positions
    )

    # -- The pocket: the exemplar's residues around the chosen sites -------
    #
    # Taken from the exemplar alone: it is the frame, it is drawn once, and
    # the pocket from each hit in turn would pile up six copies of the same
    # side chains. The residue list is spelled out in the scene because a
    # CID cannot express a sphere, and so that a user can edit it.
    pocket_selection = None
    exemplar = None
    if ref_path is not None:
        try:
            exemplar = gemmi.read_structure(str(ref_path))
        except Exception as exc:  # noqa: BLE001 - no pocket beats no scene
            logger.warning("Could not read exemplar %s: %s", ref_path, exc)
            exemplar = None
        if exemplar is not None:
            cids = (
                pocket_residue_cids(exemplar, positions[0]) if single
                else _union_cids(exemplar, positions)
            )
            stats["pocket_residues"] = len(cids)
            if cids:
                pocket_selection = "||".join(cids)
                # Straight after the ribbon; a promoted exemplar keeps its
                # own ligand sticks after both.
                ref_element["representations"].insert(
                    1,
                    {
                        "style": "CBs",
                        "selection": pocket_selection,
                        "colour": POCKET_STICK_COLOUR,
                    },
                )

    # -- Fit every hit onto the exemplar, locally, on the chosen region -----
    #
    # A global fit distributes its residual over the whole molecule and none
    # of it is guaranteed to land anywhere but the pocket; the frame that
    # must coincide here is the site's. The fit sphere is over-determined,
    # so a site origin saved in an older frame (a few Angstroms off, see the
    # design note) trades a few CAs at one edge for a few at the other.
    # Several sites share one fit, on a sphere about their centroid that
    # starts wide enough to take in every one of their pockets.
    superpose_entries: list = []
    if superpose and ref_path is not None:
        superpose_entries = _superpose(
            ref_name, ref_path, movers, stats, centre=centre,
            min_radius=None if single else spread + ENVIRONMENT_RADIUS,
            reference=exemplar,
        )

    for hit in hits:
        fit = hit.get("fit")
        transform = fit.transform() if fit is not None and fit.ok else None
        nearest = None
        try:
            structure = _hit_structure(hit)
            distances = [
                nearest_fragment_distance(structure, site_position(site), transform)
                for site in hit["sites"]
            ]
            distances = [d for d in distances if d is not None]
            nearest = min(distances) if distances else None
        except Exception as exc:  # noqa: BLE001 - a diagnostic, not the picture
            logger.warning("Could not measure %s: %s", hit["coord_path"], exc)
        entry = {"project": hit["project"].name, "verdict": hit["verdict"], "nearest": nearest}
        if not single:
            entry["sites"] = [site.name for site in hit["sites"]]
        stats["drawn"].append(entry)

    # -- Camera: the site's own, with depth clipped to the pocket -----------
    #
    # `origin` sets the camera and `slab` only the clip depth, so both are
    # needed to frame the site (see the grammar). The slab hangs off the
    # pocket residues rather than the origin because it takes a selection.
    # Several sites: the centroid of their saved origins (Moorhen's
    # convention, so averaged as stored), seen the way the first was saved.
    view = _site_view(sites[0])
    if not single:
        view["origin"] = [
            sum(site.origin[i] for site in sites) / len(sites) for i in range(3)
        ]
    if pocket_selection:
        view["slab"] = {
            "file": ref_name,
            "selection": pocket_selection,
            "pad": ENVIRONMENT_RADIUS,
        }

    title = sites[0].name if single else " + ".join(site.name for site in sites)
    scene = {
        "scene": f"{group.name} - site{'' if single else 's'} {title}",
        "version": 1,
        "authoredIn": {"projectName": group.name},
        "files": files,
        **({"superpose": superpose_entries} if superpose_entries else {}),
        "elements": elements,
        "view": view,
        "resolver": {"onMissingResidues": "clamp-and-log"},
    }
    return {"scene": scene, "stats": stats}


def _site_view(site) -> dict:
    """The camera a site row saved.

    origin_x/y/z are non-null columns, so a site always has an origin; quat
    and zoom are navigation extras and may not have been saved.
    """
    view: dict = {"origin": site.origin}
    if site.quat:
        view["quat"] = site.quat
    if site.zoom is not None:
        view["zoom"] = site.zoom
    return view


def _first_site_view(group) -> Optional[dict]:
    """Camera from the campaign's first binding site, if any.

    Reads the ``CampaignSite`` table in display order -- alphabetical by
    name, as every list of sites is shown (``lib.campaign_site_order``) --
    so "first" here is the first site the user sees in the campaign panel.

    This used to read a ``sites`` JSON attribute on the group. Migration 0024
    moved sites into their own table and ``ProjectGroup.sites`` stopped
    existing, but the read was a ``getattr(group, "sites", None)`` whose
    default silently turned a missing attribute into "this campaign has no
    sites" -- so every summary scene came back with no camera at all, and
    nothing failed to say so. Reach for a real attribute, not a defaulted
    getattr, precisely so that the next such move breaks loudly.
    """
    from .campaign_site_order import ordered_sites
    sites = ordered_sites(group.site_set.all())
    if not sites:
        return None
    return _site_view(sites[0])


# ---------------------------------------------------------------------------
# A site of a PanDDA run, across every dataset that has an event there
# ---------------------------------------------------------------------------

#: Members are drawn from their receipts by job number and parameter rather
#: than by file id, the way a receipt's own scenes reference their outputs:
#: it survives a project move, and it does not require a File row to have
#: been gleaned under a name this module would have to guess.
POSE_STYLE = "CBs"


def _run_site_member(member, used_names: set, files: list, elements: list):
    """The scene pieces for one event of one dataset at a run's site.

    Returns ``(entry, None)`` or ``(None, reason)``. The pose is the drawn
    thing; the receipt's apo model is what the fit is computed from, because
    a fit needs a protein and the pose is three dozen atoms of ligand. Both
    come from the same receipt, so they are in the same frame by
    construction.

    The receipt's dictionary travels with the pose. PanDDA writes a pose with
    the component code it built from and a viewer with no dictionary guesses
    the bonding -- which is how a receipt came to be read as a ligand with
    the wrong bonds rather than as the right ligand drawn badly.
    """
    receipt = models.Job.objects.filter(uuid=member["receipt"]["uuid"]).first()
    if receipt is None:
        return None, "receipt job not found"
    if not member.get("has_pose"):
        return None, "no autobuilt pose for this event"

    apo = models.File.objects.filter(job=receipt, job_param_name="XYZIN_APO").first()
    if apo is None or not apo.path.exists():
        return None, "receipt has no apo model to fit on"

    project_uuid = member["project"]["uuid"]
    name = _safe_name(f"{member['project']['name']}_e{member['event_idx']}", used_names)
    position = member["position"]
    files.append({
        "name": name,
        "kind": "coordinates",
        "job": receipt.number,
        "param": f"EVENTS[{position}].POSE",
        "projectId": project_uuid,
    })
    element = {"file": name, "representations": [{"style": POSE_STYLE}]}

    dictionary = models.File.objects.filter(job=receipt, job_param_name="DICT").first()
    if dictionary is not None:
        dict_name = _safe_name(f"{name}_dict", used_names)
        files.append({
            "name": dict_name,
            "kind": "dictionary",
            "job": receipt.number,
            "param": "DICT",
            "projectId": project_uuid,
        })
        element["dictionaries"] = [dict_name]

    elements.append(element)
    return {
        "name": name,
        "element": element,
        "sticks": element["representations"][0],
        "coord_path": apo.path,
        "project": SimpleNamespace(name=member["project"]["name"]),
        "member": member,
    }, None


def _member_site_point(entry) -> Optional[list]:
    """Where this member's event sits in the frame the scene draws it in.

    The event centroid is in its own dataset's frame; the fit is what carries
    it into the exemplar's. Averaging these over the members is how the site
    gets a position in the exemplar frame at all: PanDDA's own site centroid
    is in the frame of whichever dataset it chose as its reference, which is
    not the campaign's parent and is not recorded in a form this can use.
    """
    centroid = entry["member"].get("centroid")
    if not centroid:
        return None
    fit = entry.get("fit")
    if fit is None or not fit.ok:
        return list(centroid)
    moved = fit.transform().apply(gemmi.Position(*centroid))
    return [moved.x, moved.y, moved.z]


def build_run_site_scene(group, site: dict, run: Optional[dict] = None,
                         superpose: bool = True) -> dict:
    """Every ligand PanDDA built at one site of one run, in one frame.

    The panoptic view the curated site scene cannot give before anyone has
    voted: the exemplar drawn once as a ribbon, and the autobuilt pose of
    every event at this site drawn on top, each fitted onto the exemplar.
    Membership is the run's own ``site_idx`` and nothing else -- no verdict is
    required, and none is implied. These are candidates to be judged, which is
    the whole point of looking at them together.

    ``site`` is a site of ``pandda_site_index.build_site_index``; the caller
    chooses the run, so that one index serves the panel and the scene and the
    two cannot disagree about which run they are showing.

    **No maps.** A pose is coordinates and can be moved into the exemplar's
    frame by a matrix; an event map cannot, so a scene of twenty datasets'
    event maps would be twenty maps in twenty frames, all but one of them
    wrong. The receipt's own event scene is where a map belongs, at the
    contour that receipt recorded for it, in the frame it was computed in --
    which is why each member here carries the reference to open it.
    """
    files: list = []
    elements: list = []
    used_names: set = set()
    stats = {
        "site": {"site_idx": site["site_idx"], "centroid": site.get("centroid")},
        "run": run,
        "members_claimed": len(site.get("members", [])),
        "members_drawn": 0,
        "skipped": [],            # [{project, dtag, event_idx, reason}]
        "parent_present": False,
        "reference": None,
        "pocket_residues": 0,
        "drawn": [],              # [{project, dtag, event_idx, score, nearest}]
        "superpose": [],
        "centre": None,           # the site in the exemplar's frame
    }

    ref_name, ref_path, ref_element = _parent_reference(
        group, used_names, files, elements, stats
    )

    entries: list = []
    for member in site.get("members", []):
        entry, reason = _run_site_member(member, used_names, files, elements)
        if entry is None:
            stats["skipped"].append({
                "project": member["project"]["name"], "dtag": member["dtag"],
                "event_idx": member["event_idx"], "reason": reason,
            })
            continue
        entry["sticks"]["colour"] = HIT_COLOURS[stats["members_drawn"] % len(HIT_COLOURS)]
        stats["members_drawn"] += 1
        entries.append(entry)

    movers = entries
    if ref_name is None:
        ref_name, ref_path, ref_element, movers = _promote_hit(entries, stats)

    # -- Two passes, because the site's position is not known in advance ----
    #
    # A global fit first, to learn where this site is in the exemplar's
    # frame; then a fit local to that point, which is the one the scene
    # keeps. The curated site scene skips the first pass because a
    # CampaignSite already carries an origin in the right frame. A run site
    # carries a centroid in PanDDA's own reference frame, which is a
    # different thing and cannot be used as if it were this one.
    exemplar = None
    if ref_path is not None:
        try:
            exemplar = gemmi.read_structure(str(ref_path))
        except Exception as exc:  # noqa: BLE001 - no pocket beats no scene
            logger.warning("Could not read exemplar %s: %s", ref_path, exc)

    superpose_entries: list = []
    if superpose and ref_path is not None and movers:
        _superpose(ref_name, ref_path, movers, stats, reference=exemplar)
        points = [p for p in (_member_site_point(e) for e in entries) if p]
        centre = None
        if points:
            centre = [sum(axis) / len(points) for axis in zip(*points)]
            stats["centre"] = centre
        stats["superpose"] = []
        superpose_entries = _superpose(
            ref_name, ref_path, movers, stats, centre=centre, reference=exemplar
        )
    else:
        points = [p for p in (_member_site_point(e) for e in entries) if p]
        if points:
            stats["centre"] = [sum(axis) / len(points) for axis in zip(*points)]

    centre = stats["centre"]

    # -- The pocket, from the exemplar alone -------------------------------
    pocket_selection = None
    if ref_path is not None and centre is not None:
        if exemplar is not None:
            cids = pocket_residue_cids(exemplar, centre)
            stats["pocket_residues"] = len(cids)
            if cids:
                pocket_selection = "||".join(cids)
                ref_element["representations"].insert(
                    1,
                    {"style": "CBs", "selection": pocket_selection,
                     "colour": POCKET_STICK_COLOUR},
                )

    for entry in entries:
        member = entry["member"]
        stats["drawn"].append({
            "project": member["project"]["name"],
            "dtag": member["dtag"],
            "event_idx": member["event_idx"],
            "score": member.get("score"),
            "hit_probability": member.get("hit_probability"),
            "receipt": member["receipt"],
            "position": member["position"],
        })

    view: dict = {}
    if centre is not None:
        # Moorhen's view origin is the negation of the point to centre on
        # (see site_position for the trap this is the other side of).
        view["origin"] = [-centre[0], -centre[1], -centre[2]]
    if pocket_selection:
        view["slab"] = {"file": ref_name, "selection": pocket_selection,
                        "pad": ENVIRONMENT_RADIUS}

    label = f"site {site['site_idx']}"
    if run and run.get("number"):
        label += f" of run {run.get('project') or ''}/{run['number']}".rstrip("/")
    scene = {
        "scene": f"{group.name} - {label}",
        "version": 1,
        "authoredIn": {"projectName": group.name},
        "files": files,
        **({"superpose": superpose_entries} if superpose_entries else {}),
        "elements": elements,
        **({"view": view} if view else {}),
        "resolver": {"onMissingResidues": "clamp-and-log"},
    }
    return {"scene": scene, "stats": stats}
