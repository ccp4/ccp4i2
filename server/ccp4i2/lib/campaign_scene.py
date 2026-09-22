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
from typing import Optional

import gemmi

from ..db import models
from . import pandda_export, superposition

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
    refine_job = pandda_export._latest_finished_job(project, REFINE_TASK_NAMES)
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

    codes = detect_ligands(coord_file.path, dict_path)
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


def _superpose(ref_name: str, ref_path, movers: list, stats: dict,
               centre=None) -> list:
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
    """
    entries: list = []
    for hit in movers:
        fit = superposition.fit_files(ref_path, hit["coord_path"], centre=centre)
        record = {
            "project": hit["project"].name,
            "ok": fit.ok,
            "atoms": fit.atoms,
            "radius": fit.radius,
            "rmsd": fit.rmsd,
            "reason": fit.reason,
        }
        if not fit.ok and centre is not None:
            fallback = superposition.fit_files(ref_path, hit["coord_path"])
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

    return {"scene": scene, "stats": stats}


def build_site_scene(group, site, include_unclear: bool = False,
                     superpose: bool = True) -> dict:
    """Build the scene for one binding site of a campaign.

    The exemplar (the parent, else a hit) drawn once as a ribbon with the
    residues around the site origin as sticks, and every dataset judged a
    hit *at this site* drawn on top as its ligand's sticks, each fitted onto
    the exemplar locally, on the CAs around the site.

    Membership is the verdict and nothing else: a dataset is in iff it has a
    ``hit`` evaluation at this site, plus ``unclear`` ones when asked for,
    in a muted colour. ``empty`` and unevaluated members are out, and there
    is deliberately no fall-back to ligand detection when a site has no
    verdicts -- that would show every ligand in the campaign under the name
    of one site. A site nobody has evaluated gets an honest exemplar-only
    scene, and ``stats`` says so.

    Each hit's ligand is drawn as every copy of its code (``//*/(CODE)``),
    as the summary draws it; the camera and slab frame the site, and a copy
    bound in an adjacent subsite is an observation, not clutter (see
    ``docs/campaign-site-scene-design.md``, "Draw every copy").
    """
    files: list = []
    elements: list = []
    used_names: set = set()
    stats = {
        "site": {"id": site.id, "name": site.name},
        "hits_claimed": 0,        # hit verdicts found at this site
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
        "drawn": [],              # [{project, verdict, nearest}]
        "superpose": [],          # [{project, ok, atoms, radius, rmsd, reason, fallback?}]
    }

    ref_name, ref_path, ref_element = _parent_reference(
        group, used_names, files, elements, stats
    )

    # -- Membership: the verdicts at this site ------------------------------
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
    # Ordered by project id so the fallback exemplar, when one is needed, is
    # the lowest-id hit: stable, arbitrary, and recorded as such in stats.
    evaluations = site.evaluations.select_related("project").order_by("project_id")

    hits: list = []
    for evaluation in evaluations:
        if evaluation.project_id not in member_ids:
            continue
        verdict = evaluation.verdict
        if verdict == Verdict.EMPTY:
            stats["empty_verdicts"] += 1
            continue
        if verdict == Verdict.UNCLEAR:
            stats["unclear_verdicts"] += 1
            if not include_unclear:
                continue
        else:
            stats["hits_claimed"] += 1

        project = evaluation.project
        hit, reason = _resolve_hit(project, used_names)
        if hit is None:
            stats["skipped"].append(
                {"project": project.name, "reason": reason, "nearest": None}
            )
            continue
        hit["verdict"] = verdict
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

    # -- The pocket: the exemplar's residues around the site origin --------
    #
    # Taken from the exemplar alone: it is the frame, it is drawn once, and
    # the pocket from each hit in turn would pile up six copies of the same
    # side chains. The residue list is spelled out in the scene because a
    # CID cannot express a sphere, and so that a user can edit it.
    pocket_selection = None
    if ref_path is not None:
        try:
            exemplar = gemmi.read_structure(str(ref_path))
        except Exception as exc:  # noqa: BLE001 - no pocket beats no scene
            logger.warning("Could not read exemplar %s: %s", ref_path, exc)
            exemplar = None
        if exemplar is not None:
            cids = pocket_residue_cids(exemplar, site_position(site))
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

    # -- Fit every hit onto the exemplar, locally, on the site --------------
    #
    # A global fit distributes its residual over the whole molecule and none
    # of it is guaranteed to land anywhere but the pocket; the frame that
    # must coincide here is the site's. The fit sphere is over-determined,
    # so a site origin saved in an older frame (a few Angstroms off, see the
    # design note) trades a few CAs at one edge for a few at the other.
    superpose_entries: list = []
    if superpose and ref_path is not None:
        superpose_entries = _superpose(
            ref_name, ref_path, movers, stats, centre=site_position(site)
        )

    for hit in hits:
        fit = hit.get("fit")
        transform = fit.transform() if fit is not None and fit.ok else None
        nearest = None
        try:
            nearest = nearest_fragment_distance(
                gemmi.read_structure(str(hit["coord_path"])), site_position(site), transform
            )
        except Exception as exc:  # noqa: BLE001 - a diagnostic, not the picture
            logger.warning("Could not measure %s: %s", hit["coord_path"], exc)
        stats["drawn"].append(
            {"project": hit["project"].name, "verdict": hit["verdict"], "nearest": nearest}
        )

    # -- Camera: the site's own, with depth clipped to the pocket -----------
    #
    # `origin` sets the camera and `slab` only the clip depth, so both are
    # needed to frame the site (see the grammar). The slab hangs off the
    # pocket residues rather than the origin because it takes a selection.
    view = _site_view(site)
    if pocket_selection:
        view["slab"] = {
            "file": ref_name,
            "selection": pocket_selection,
            "pad": ENVIRONMENT_RADIUS,
        }

    scene = {
        "scene": f"{group.name} - site {site.name}",
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

    Reads the ``CampaignSite`` table, ordered by the site's display ``order``
    (``CampaignSite.Meta``), so "first" here is the first site the user sees
    in the campaign panel.

    This used to read a ``sites`` JSON attribute on the group. Migration 0024
    moved sites into their own table and ``ProjectGroup.sites`` stopped
    existing, but the read was a ``getattr(group, "sites", None)`` whose
    default silently turned a missing attribute into "this campaign has no
    sites" -- so every summary scene came back with no camera at all, and
    nothing failed to say so. Reach for a real attribute, not a defaulted
    getattr, precisely so that the next such move breaks loudly.
    """
    site = group.site_set.first()
    if site is None:
        return None
    return _site_view(site)
