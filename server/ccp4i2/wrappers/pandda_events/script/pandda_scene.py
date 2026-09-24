"""Scenes for a PanDDA receipt (design note section 11).

A receipt knows what its outputs should look like: the apo model is the model
of record, the Z-map is in z units and reads at z = 3, an event map reads at
the display contour recorded for it (absolute map units, section 7.3), the
pose is a candidate drawn with the dictionary PanDDA built it from, and the
camera belongs on the event centroid. This module says so as scene documents
(the format lives in client/renderer/lib/scene; the grammar is generated,
this authors against it), referencing the receipt's own outputs by job
number and parameter, which survive a project move.
"""
from typing import Dict, List, Optional

Z_CONTOUR = 3.0            # a Z-map is in z units; contour at z = 3
EVENT_MAP_COLOUR = "#3f51b5"
EVENT_ZOOM = 0.35          # a binding pocket, not the whole crystal
MAP_RADIUS = 12.0


def _ref(name: str, param: str, kind: str, job_number: str, project_id: str) -> Dict:
    return {"name": name, "kind": kind, "job": _job(job_number), "param": param, "projectId": project_id}


def _job(job_number: str):
    """The scene's job reference is the job number; an integer where it is
    one (top-level jobs), the dotted string otherwise."""
    text = str(job_number)
    return int(text) if text.isdigit() else text


def _base(name: str, job_number: str, project_id: str, dtag: str, *, apo: bool, zmap: bool,
          dictionary: bool) -> Dict:
    scene: Dict = {"scene": name, "version": 1, "files": [], "elements": [], "maps": []}
    if apo:
        scene["files"].append(_ref("apo", "XYZIN_APO", "coordinates", job_number, project_id))
        scene["elements"].append({"file": "apo", "representations": [{"style": "CRs"}, {"style": "CBs"}]})
    if dictionary:
        scene["files"].append(_ref("dict", "DICT", "dictionary", job_number, project_id))
    if zmap:
        scene["files"].append(_ref("zmap", "ZMAP", "map", job_number, project_id))
        scene["maps"].append({"name": "Z-map", "file": "zmap", "isDifference": True,
                              "contourLevel": Z_CONTOUR, "radius": MAP_RADIUS})
    return scene


def _add_event(scene: Dict, event: Dict, job_number: str, project_id: str, *, dictionary: bool,
               visible: bool) -> None:
    idx = event["idx"]
    if event.get("has_map"):
        scene["files"].append(_ref(f"event{idx}_map", f"EVENTS[{event['position']}].EVENT_MAP", "map",
                                   job_number, project_id))
        entry = {"name": f"Event {idx} map", "file": f"event{idx}_map", "radius": MAP_RADIUS,
                 "colour": EVENT_MAP_COLOUR, "visible": visible}
        if event.get("display_contour") is not None:
            entry["contourLevel"] = float(event["display_contour"])
        scene["maps"].append(entry)
    if event.get("has_pose"):
        scene["files"].append(_ref(f"event{idx}_pose", f"EVENTS[{event['position']}].POSE", "coordinates",
                                   job_number, project_id))
        element = {"file": f"event{idx}_pose", "representations": [{"style": "CBs"}]}
        if dictionary:
            element["dictionaries"] = ["dict"]
        scene["elements"].append(element)


def _look_at(scene: Dict, event: Dict) -> None:
    centroid = event.get("centroid")
    if centroid is not None:
        scene["view"] = {"origin": [float(c) for c in centroid], "zoom": EVENT_ZOOM}
    if event.get("has_map"):
        scene["activeMap"] = f"Event {event['idx']} map"


def event_scene(dtag: str, job_number: str, project_id: str, event: Dict, *, apo: bool, zmap: bool,
                dictionary: bool) -> Dict:
    """One event: the apo model, the Z-map, this event's map at its display
    contour and its pose, centred on the event."""
    ligand = f" ({event['ligand_id']})" if event.get("ligand_id") else ""
    scene = _base(f"{dtag} event {event['idx']}{ligand}", job_number, project_id, dtag,
                  apo=apo, zmap=zmap, dictionary=dictionary)
    _add_event(scene, event, job_number, project_id, dictionary=dictionary, visible=True)
    _look_at(scene, event)
    return scene


def focus_event(events: List[Dict]) -> Optional[Dict]:
    """The event the overview opens on: the highest event score, else the
    first. A score is a machine opinion, not a verdict, but as the place to
    look first it is the best one going."""
    if not events:
        return None
    scored = [e for e in events if e.get("score") is not None]
    return max(scored, key=lambda e: e["score"]) if scored else events[0]


def overview_scene(dtag: str, job_number: str, project_id: str, events: List[Dict], *, apo: bool,
                   zmap: bool, dictionary: bool) -> Dict:
    """The whole receipt: every event's map and pose, the focus event's map
    shown and the rest hidden, centred on the focus event."""
    scene = _base(f"{dtag}: PanDDA events", job_number, project_id, dtag,
                  apo=apo, zmap=zmap, dictionary=dictionary)
    focus = focus_event(events)
    for event in events:
        _add_event(scene, event, job_number, project_id, dictionary=dictionary, visible=(event is focus))
    if focus is not None:
        _look_at(scene, focus)
    return scene


def dump(scene: Dict) -> str:
    import yaml
    return yaml.safe_dump(scene, sort_keys=False, allow_unicode=False)
