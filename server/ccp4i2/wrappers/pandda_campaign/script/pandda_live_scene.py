"""Look at a PanDDA run's events while it runs.

Receipts deliver each event to its dataset's project with a scene, but only
after the run ends and fan-out runs. Until then the events are in the run's
own tree, written dataset by dataset, and that is where these read them:

- ``peek_events`` lists the events the finished datasets report, from each
  dataset's ``events.yaml`` (the run-level table comes only at the end).
  A dataset is finished when its Z-map exists; finished datasets do not
  change, so a cache passed in is kept and only new datasets are read.
- ``live_event_scene`` builds the scene for one event, as the receipt would
  (pandda_events/script/pandda_scene.py), but naming the run's files by
  their path in the project, since the database records none of them yet.
  It reads the event map to choose a contour, so it runs when someone asks
  to look, never when the report is drawn.
"""
import os
from pathlib import Path
from typing import Dict, List, Optional

from ccp4i2.wrappers.pandda_events.script.pandda_scene import event_scene
from ccp4i2.wrappers.pandda_events.script.pandda_tree import (
    EVENTS_YAML, PROCESSED_DIR, DatasetNotFound, display_contour, read_dataset, read_events_yaml)

#: A dataset is finished when PanDDA has written its Z-map.
ZMAP_SUFFIX = "-z_map.native.ccp4"


def _finished(dataset_dir: Path) -> bool:
    return (dataset_dir / f"{dataset_dir.name}{ZMAP_SUFFIX}").is_file()


def _float(value) -> Optional[float]:
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


def _dataset_events(dataset_dir: Path) -> List[Dict]:
    path = dataset_dir / EVENTS_YAML
    if not path.is_file():
        return []
    try:
        records = read_events_yaml(path)
    except Exception:       # noqa: BLE001 - one unreadable file is not worth the list
        return []
    events = []
    for idx in sorted(records):
        record = records[idx]
        centroid = record.get("Centroid")
        if not (isinstance(centroid, (list, tuple)) and len(centroid) == 3):
            centroid = None
        events.append({
            "dtag": dataset_dir.name,
            "event_idx": idx,
            "score": _float(record.get("Score")),
            "bdc": _float(record.get("BDC")),
            "centroid": [_float(c) for c in centroid] if centroid else None,
        })
    return events


def peek_events(tree, cache: Optional[Dict[str, List[Dict]]] = None) -> List[Dict]:
    """Every event the finished datasets report, best score first.

    ``cache`` maps a finished dataset to its events and is filled in place;
    pass the same dict on the next call and only datasets finished since are
    read. Unfinished datasets are never cached (their events.yaml may not be
    written yet)."""
    cache = {} if cache is None else cache
    processed = Path(tree) / PROCESSED_DIR
    if not processed.is_dir():
        return []
    for entry in os.scandir(processed):
        if entry.name in cache or not entry.is_dir():
            continue
        dataset_dir = Path(entry.path)
        if _finished(dataset_dir):
            cache[entry.name] = _dataset_events(dataset_dir)
    events = [event for found in cache.values() for event in found]
    events.sort(key=lambda e: (-(e["score"] if e["score"] is not None else float("-inf")),
                               e["dtag"], e["event_idx"]))
    return events


#: What PanDDA names every residue it builds, whatever the dictionary said.
PANDDA_RESIDUE = "LIG"


def dictionary_for_coot(path) -> Optional[str]:
    """The run's ligand dictionary, rewritten so Coot knows ``LIG``.

    PanDDA builds every pose as residue ``LIG``, and staging gives it a
    ``comp_LIG`` block to read (pandda_dict.add_lig_alias). PanDDA opens a
    dictionary by block name; Coot registers compounds by the
    ``_chem_comp_atom.comp_id`` in the rows and the ``comp_list`` header,
    which in that alias still name the true code. So to Coot the pose's LIG
    had no dictionary, and Moorhen fetched the library's ``LIG`` -- a real,
    unrelated CCD compound -- and drew the pose with its bonds.

    Here the ``comp_LIG`` block's rows say ``LIG`` and ``comp_list`` lists
    it, beside the true code, unchanged. The file on disk is not touched.
    None when the file has no restraint block to offer.
    """
    import gemmi

    try:
        doc = gemmi.cif.read(str(path))
    except Exception:       # noqa: BLE001 - no dictionary is better than no scene
        return None
    restraint = [b for b in doc if b.name != "comp_list" and b.find_values("_chem_comp_atom.atom_id")]
    if not restraint:
        return None
    alias = doc.find_block(f"comp_{PANDDA_RESIDUE}")
    if alias is None:
        alias = doc.add_copied_block(restraint[0])
        alias.name = f"comp_{PANDDA_RESIDUE}"
    atom_ids = alias.find_values("_chem_comp_atom.comp_id")
    code = atom_ids[0] if len(atom_ids) else None
    for item in alias:
        if item.loop is not None:
            for tag in item.loop.tags:
                if tag.endswith(".comp_id") or tag == "_chem_comp.id":
                    column = alias.find_values(tag)
                    for i in range(len(column)):
                        column[i] = PANDDA_RESIDUE
        elif item.pair is not None and (item.pair[0].endswith(".comp_id") or item.pair[0] == "_chem_comp.id"):
            alias.set_pair(item.pair[0], PANDDA_RESIDUE)
    header = doc.find_block("comp_list")
    ids = header.find_values("_chem_comp.id") if header is not None else None
    if ids is not None and len(ids) and PANDDA_RESIDUE not in list(ids):
        loop = ids.get_loop()
        if loop is not None:
            width = loop.width()
            rows = [list(loop.values[r * width:(r + 1) * width]) for r in range(loop.length())]
            id_col = loop.tags.index("_chem_comp.id")
            values = next((r for r in rows if r[id_col] == code), rows[0])
            for c, tag in enumerate(loop.tags):
                if tag in ("_chem_comp.id", "_chem_comp.three_letter_code"):
                    values[c] = PANDDA_RESIDUE
            loop.add_row(values)
    return doc.as_string()


def _safe_dtag(dtag: str) -> str:
    dtag = str(dtag)
    if not dtag or "/" in dtag or "\\" in dtag or dtag in (".", ".."):
        raise ValueError(f"not a dataset name: {dtag!r}")
    return dtag


def live_event_scene(tree, staged_datasets, project_dir, project_uuid: str, dtag: str,
                     event_idx: Optional[int] = None) -> Dict:
    """The scene for one event of a run in progress (or finished, before fan-out).

    The apo model, the event's map at its display contour and PanDDA's
    candidate pose, centred on the event; the Z-map is loaded but hidden. The
    ligand's dictionary is attached to the pose and to the apo model both. With no ``event_idx``, the dataset alone: the apo model
    and its Z-map, shown. Files outside the project directory (a link the
    share cannot follow) are left out rather than named.
    """
    dtag = _safe_dtag(dtag)
    staged = Path(staged_datasets) / dtag if staged_datasets else None
    try:
        outputs = read_dataset(tree, dtag, staged_dir=staged)
    except DatasetNotFound:
        raise ValueError(f"PanDDA has not written dataset {dtag} yet")
    project_dir = Path(project_dir).resolve()

    def in_project(path) -> Optional[str]:
        if path is None:
            return None
        resolved = Path(path).resolve()
        if project_dir not in resolved.parents:
            return None
        return resolved.relative_to(project_dir).as_posix()

    event = None
    if event_idx is not None:
        event = next((e for e in outputs.events if e.idx == int(event_idx)), None)
        if event is None:
            raise ValueError(f"{dtag} has no event {event_idx}")

    paths = {
        "XYZIN_APO": in_project(outputs.apo_model),
        "ZMAP": in_project(outputs.zmap),
        "DICT": in_project(outputs.dictionary),
    }
    if event is not None:
        paths["EVENTS[0].EVENT_MAP"] = in_project(event.event_map)
        paths["EVENTS[0].POSE"] = in_project(event.pose)

    # The dictionary travels inline, rewritten so Coot knows the pose's LIG.
    dictionary_text = dictionary_for_coot(outputs.dictionary) if paths["DICT"] else None

    def ref(name: str, role: str, kind: str) -> Dict:
        if role == "DICT" and dictionary_text:
            return {"name": name, "kind": kind, "cifText": dictionary_text}
        return {"name": name, "kind": kind, "projectPath": paths[role], "projectId": project_uuid}

    has_pose = event is not None and paths["EVENTS[0].POSE"] is not None
    if event is None:
        scene = event_scene(dtag, "", project_uuid, {"idx": 0, "position": 0},
                            apo=paths["XYZIN_APO"] is not None, zmap=paths["ZMAP"] is not None,
                            dictionary=paths["DICT"] is not None, ref=ref)
        scene["scene"] = f"{dtag}: Z-map"
        scene.pop("activeMap", None)
        return scene

    has_map = paths["EVENTS[0].EVENT_MAP"] is not None
    record = {
        "idx": event.idx,
        "position": 0,
        "has_map": has_map,
        "has_pose": has_pose,
        "display_contour": display_contour(event.event_map,
                                           event.build.optimal_contour if event.build else None)
        if has_map else None,
        "centroid": event.centroid,
        "ligand_id": outputs.ligand_id,
        "score": event.score,
    }
    # The dictionary whenever there is one, pose or not: it rides on the apo
    # model too, so a ligand built or merged there is read correctly.
    scene = event_scene(dtag, "", project_uuid, record,
                        apo=paths["XYZIN_APO"] is not None, zmap=paths["ZMAP"] is not None,
                        dictionary=paths["DICT"] is not None, ref=ref)
    # The event map is what to look at; the Z-map is there to switch on.
    for entry in scene["maps"]:
        if entry.get("file") == "zmap":
            entry["visible"] = False
    return scene
