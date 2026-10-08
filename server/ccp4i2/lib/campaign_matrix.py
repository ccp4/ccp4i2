"""The campaign's dataset x site matrix: pure rules, no database.

The campaign overview has one row per dataset and, since this module, one
column per curated ``CampaignSite``. A cell says two independent things:

* whether PanDDA found an event close to that site in that dataset -- a
  machine observation, read from ``CampaignEvent`` rows;
* what a person decided there -- the ``SiteEvaluation`` verdict, or nothing.

Nothing here writes a verdict, and nothing here infers one from an event:
an event near a site is evidence for a person to look at, not a decision
(docs/pandda-campaign-design.md, section 9.1, invariant 1).

Everything in this module is a pure function of plain values, so the rules
can be tested without Django, a database or CCP4. The callers that touch the
database live in :mod:`ccp4i2.lib.campaign_events`.
"""
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

# ---------------------------------------------------------------------------
# Which job holds a dataset's current model
# ---------------------------------------------------------------------------

#: Tasks whose finished, top-level run leaves a refined model (and its maps)
#: as the dataset's model of record. ``SubstituteLigand`` is here because it
#: is the campaign's own route from data to a ligand-bound model: it runs
#: servalcat as a subjob and publishes the refined model, maps and dictionary
#: as its own outputs, so for a campaign member it *is* the refinement.
#: Refinement subjobs of a pipeline are deliberately not counted: the
#: pipeline's outputs are its answer, and a subjob is one step towards it.
REFINEMENT_TASKS = (
    "servalcat_pipe",
    "servalcat",
    "prosmart_refmac",
    "refmac",
    "i2Refmac",      # legacy name; no longer registered, but imported jobs carry it
    "buster",
    "lorestr_i2",
    "pdb_redo_api",
    "SubstituteLigand",
)

#: DIMPLE: rigid-body refinement of the reference model against a dataset.
#: The fallback model, and also what PanDDA takes as input
#: (``pandda_export.DIMPLE_TASK_NAMES`` is this tuple).
DIMPLE_TASKS = ("i2Dimple", "dimple")   # "dimple": legacy name, as i2Refmac

#: ``Job.Status.FINISHED``. Spelled as a number so this module needs no
#: Django import; ``test_campaign_matrix`` pins it to the model's value.
FINISHED = 6


def _field(job, name):
    return job.get(name) if isinstance(job, dict) else getattr(job, name, None)


def choose_current_model(jobs: Iterable) -> Optional[object]:
    """The job holding a dataset's current model, or None.

    ``jobs`` are a project's jobs, as model instances or dicts with ``id``,
    ``task_name``, ``status`` and ``parent_id``. The rule:

    1. the latest (highest id) FINISHED *top-level* job of a refinement task;
    2. else the latest FINISHED dimple run, at any level -- on the merged
       SubstituteLigand route DIMPLE runs as a subjob, and a pipeline that
       failed after it leaves that subjob's model as the best there is;
    3. else None.

    A refinement is preferred to a later DIMPLE run: DIMPLE re-fits the
    reference model, so a newer DIMPLE job is not a better model of the
    dataset than an older refinement of it.
    """
    refined = None
    dimple = None
    for job in jobs:
        if _field(job, "status") != FINISHED:
            continue
        task = _field(job, "task_name")
        job_id = _field(job, "id") or 0
        if task in REFINEMENT_TASKS and _field(job, "parent_id") is None:
            if refined is None or job_id > (_field(refined, "id") or 0):
                refined = job
        elif task in DIMPLE_TASKS:
            if dimple is None or job_id > (_field(dimple, "id") or 0):
                dimple = job
    return refined if refined is not None else dimple


# ---------------------------------------------------------------------------
# Is a dataset in the parent's frame?
# ---------------------------------------------------------------------------

#: Largest relative difference in any cell edge, and largest difference in
#: any cell angle (degrees), at which a dataset's coordinates are still
#: compared directly with a site placed in the parent's frame.
CELL_EDGE_TOLERANCE = 0.02
CELL_ANGLE_TOLERANCE = 2.0
#: How far (A) a cell difference may move a point at the farthest site
#: before the dataset is not compared directly. A cell edge off by a fraction
#: d moves a point at distance r from the origin by about d * r: BAZ2B's
#: 5e9l, 2.5% off the parent in a, moves its pocket ~27 A out by 0.7 A, well
#: inside an 8 A radius, and a 2% edge rule wrongly hid its event.
POSITION_TOLERANCE = 1.5
#: A cell edge off by more than this is another crystal form or setting,
#: whatever the sites' distance from the origin.
CELL_EDGE_LIMIT = 0.10

_EDGES = ("a", "b", "c")
_ANGLES = ("alpha", "beta", "gamma")


def frame_mismatch(dataset_cell: Optional[Sequence[float]],
                   parent_cell: Optional[Sequence[float]],
                   reach: Optional[float] = None) -> Optional[str]:
    """A short reason the dataset is not in the parent's frame, or None.

    Event centroids are stated in each dataset's own frame and site origins
    in the parent's. Comparing them directly is right only because a
    campaign's members are near-isomorphous. What matters is how far a cell
    difference moves a point where the sites are: *reach* is the farthest
    site's distance (A) from the origin, and an edge off by a fraction d
    moves such a point by about d * reach, allowed up to POSITION_TOLERANCE.
    An edge off by more than CELL_EDGE_LIMIT, or an angle by more than
    CELL_ANGLE_TOLERANCE, is another form or setting regardless. Without a
    reach, an edge is held to CELL_EDGE_TOLERANCE.

    When either cell is unknown nothing can be said, and None is returned:
    the caller compares directly, as it would for a matching cell.
    """
    if not dataset_cell or not parent_cell:
        return None
    if len(dataset_cell) != 6 or len(parent_cell) != 6:
        return None
    for i, name in enumerate(_EDGES):
        mine, ref = float(dataset_cell[i]), float(parent_cell[i])
        if ref <= 0:
            return None
        off = abs(mine - ref) / ref
        too_far = (off * reach > POSITION_TOLERANCE if reach is not None
                   else off > CELL_EDGE_TOLERANCE)
        if off > CELL_EDGE_LIMIT or too_far:
            shift = f", ~{off * reach:.1f} A at the sites" if reach is not None else ""
            return (f"cell {name} {mine:.2f} A vs parent {ref:.2f} A "
                    f"({100.0 * off:.1f}% off{shift})")
    for i, name in enumerate(_ANGLES):
        mine, ref = float(dataset_cell[3 + i]), float(parent_cell[3 + i])
        if abs(mine - ref) > CELL_ANGLE_TOLERANCE:
            return f"cell {name} {mine:.1f} vs parent {ref:.1f} degrees"
    return None


def site_reach(sites) -> Optional[float]:
    """The farthest site's distance (A) from the frame origin."""
    distances = [math.sqrt(sum(float(c) ** 2 for c in s["position"]))
                 for s in sites if s.get("position") is not None]
    return max(distances) if distances else None


# ---------------------------------------------------------------------------
# Which event, if any, sits at each site
# ---------------------------------------------------------------------------

#: A site's radius when it has none recorded (A). Matches the model default.
DEFAULT_SITE_RADIUS = 8.0


def _distance(a: Sequence[float], b: Sequence[float]) -> float:
    return math.sqrt(sum((float(a[i]) - float(b[i])) ** 2 for i in range(3)))


def match_events_to_sites(events: Iterable[dict], sites: Iterable[dict],
                          mismatch: Optional[str] = None) -> Dict[str, Optional[dict]]:
    """Each site's cell: the nearest event within the site's radius, or None.

    ``events`` are dicts with ``event_idx``, ``centroid`` (``[x, y, z]`` or
    None), ``hit_probability`` and ``has_pose``. ``sites`` are dicts with
    ``uuid``, ``position`` (the centre in real space, ``[x, y, z]`` -- never
    ``CampaignSite.origin``, which is its negation) and ``radius``. The result has one
    key per site uuid, always, so a caller can render every column.

    * An event counts for a site when its centroid is within the site's
      radius of the site's position, the boundary included.
    * An event may count for more than one site: radii are per site and may
      overlap, and which pocket an event "really" belongs to is a judgement
      this function does not make.
    * Of several events within a site's radius the nearest wins; a tie goes
      to the higher hit probability, then the lower event number, so the
      answer does not depend on input order.
    * An event without a centroid is never matched.
    * With ``mismatch`` set (see :func:`frame_mismatch`) nothing is matched:
      distances across frames are not distances.
    """
    sites = list(sites)
    cells: Dict[str, Optional[dict]] = {str(site["uuid"]): None for site in sites}
    if mismatch:
        return cells
    events = [e for e in events if e.get("centroid") is not None]
    for site in sites:
        position = site["position"]
        radius = site.get("radius")
        radius = DEFAULT_SITE_RADIUS if radius is None else float(radius)
        best = None
        best_key = None
        for event in events:
            d = _distance(event["centroid"], position)
            if d > radius:
                continue
            probability = event.get("hit_probability")
            key = (d, -(probability if probability is not None else -1.0),
                   event.get("event_idx") if event.get("event_idx") is not None else 0)
            if best_key is None or key < best_key:
                best, best_key = event, key
        if best is not None:
            cells[str(site["uuid"])] = {
                "event_idx": best.get("event_idx"),
                "hit_probability": best.get("hit_probability"),
                "distance": round(best_key[0], 2),
                "has_pose": bool(best.get("has_pose")),
            }
    return cells


# ---------------------------------------------------------------------------
# Reading a unit cell cheaply
# ---------------------------------------------------------------------------

_CIF_CELL_KEYS = ("_cell.length_a", "_cell.length_b", "_cell.length_c",
                  "_cell.angle_alpha", "_cell.angle_beta", "_cell.angle_gamma")


def _cif_number(text: str) -> Optional[float]:
    # mmCIF numbers may carry an uncertainty in brackets: 60.123(4)
    match = re.match(r"^\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)", text or "")
    return float(match.group(1)) if match else None


def read_model_cell(path) -> Optional[List[float]]:
    """``[a, b, c, alpha, beta, gamma]`` from a coordinate file's header, or None.

    Reads only as far as the cell: the ``CRYST1`` record of a PDB file, the
    ``_cell.*`` items of an mmCIF one. Standard library only, so it costs a
    few kilobytes of reading rather than a structure parse, and it works on
    the CCP4-free server. A file it cannot read, or one with no cell (or the
    1 A placeholder cell of a model with no crystal), gives None.
    """
    path = Path(path)
    try:
        handle = path.open("r", errors="replace")
    except OSError:
        return None
    cell: Optional[List[float]] = None
    with handle:
        found: Dict[str, float] = {}
        for line in handle:
            if line.startswith("CRYST1"):
                try:
                    cell = [float(line[6:15]), float(line[15:24]), float(line[24:33]),
                            float(line[33:40]), float(line[40:47]), float(line[47:54])]
                except ValueError:
                    cell = None
                break
            if line.startswith(("ATOM", "HETATM")):
                break
            if line.startswith("_cell."):
                parts = line.split(None, 1)
                if len(parts) == 2 and parts[0] in _CIF_CELL_KEYS:
                    value = _cif_number(parts[1])
                    if value is not None:
                        found[parts[0]] = value
                if len(found) == 6:
                    cell = [found[k] for k in _CIF_CELL_KEYS]
                    break
            elif found and line.startswith("_") and not line.startswith("_cell."):
                # The _cell category is contiguous; it has ended incomplete.
                break
            elif line.startswith("loop_") and found:
                break
    if cell is None:
        return None
    if all(abs(v - 1.0) < 1e-6 for v in cell[:3]):
        return None
    return cell
