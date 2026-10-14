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


def _safe_dtag(dtag: str) -> str:
    dtag = str(dtag)
    if not dtag or "/" in dtag or "\\" in dtag or dtag in (".", ".."):
        raise ValueError(f"not a dataset name: {dtag!r}")
    return dtag


def live_event_scene(tree, staged_datasets, project_dir, project_uuid: str, dtag: str,
                     event_idx: Optional[int] = None) -> Dict:
    """The scene for one event of a run in progress (or finished, before fan-out).

    The apo model, the event's map at its display contour and PanDDA's
    candidate pose with its dictionary, centred on the event; the Z-map is
    loaded but hidden. With no ``event_idx``, the dataset alone: the apo model
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

    def ref(name: str, role: str, kind: str) -> Dict:
        return {"name": name, "kind": kind, "projectPath": paths[role], "projectId": project_uuid}

    has_pose = event is not None and paths["EVENTS[0].POSE"] is not None
    if event is None:
        scene = event_scene(dtag, "", project_uuid, {"idx": 0, "position": 0},
                            apo=paths["XYZIN_APO"] is not None, zmap=paths["ZMAP"] is not None,
                            dictionary=False, ref=ref)
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
    scene = event_scene(dtag, "", project_uuid, record,
                        apo=paths["XYZIN_APO"] is not None, zmap=paths["ZMAP"] is not None,
                        dictionary=has_pose and paths["DICT"] is not None, ref=ref)
    # The event map is what to look at; the Z-map is there to switch on.
    for entry in scene["maps"]:
        if entry.get("file") == "zmap":
            entry["visible"] = False
    return scene
