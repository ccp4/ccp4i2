"""What PhaserTNG Picard wrote, read into one record a judgement can use.

Picard keeps its results in a database directory of numbered nodes. The
numbers that decide whether the run worked sit in three places:

* ``best.N.dag.cards``: one block per solution, best first, ``====`` between
  blocks. ``phaserdag node <key> <value>`` lines give the whole solution's
  TFZ (``zscore``), LLG, R-factor (``rfactor``, percent), space group,
  whether it packs, twinning and tNCS flags, how many components were
  sought (``seek_components``) and whether all were placed
  (``seek_complete``); one ``phaserdag node pose`` line per placed component
  ends ``tfz <z>``; and ``annotation`` narrates the search, one ``FIND=<tag>
  RF=.../<z>z ... TF=.../<z>z/<packs>`` run per component.
* ``<node>-rfac/solutions.log`` and the run log's ``Best Solution:`` line:
  ``R-factor = <Rwork>/<Rfree> dmin = <d> poses = <n> sg = <sg>``, the only
  place R-free is written.
* ``<node>.result.json`` of the last node: the data (resolution, anomalous,
  Wilson B), the Matthews analysis, twinning and tNCS, and wall time.

The wrapper writes these into program.xml (``program_xml``) and the report
draws from the same parse. Stdlib only.
"""
import json
import os
import re
from pathlib import Path
from typing import Dict, List, Optional
from xml.etree import ElementTree as ET

_BEST = re.compile(r"R-factor\s*=\s*([\d.]+)\s*/\s*([\d.]+)\s+dmin\s*=\s*([\d.]+)\s+poses\s*=\s*(\d+)\s+sg\s*=\s*(\S+)")
_FIND = re.compile(r"FIND=(\S+)(.*?)(?=FIND=|$)")
_RF = re.compile(r"\bRF=\S*?/(\d+(?:\.\d+)?)z")
_TF = re.compile(r"\bTF=\S*?/(\d+(?:\.\d+)?)z/(\S)")
_TFZ = re.compile(r"\btfz\s+(-?[\d.]+)\s*$")


def parse_dag_cards(text: str) -> List[Dict]:
    """Solutions of a best.N.dag.cards file, best first. Each is the first
    value of every ``phaserdag node`` key, plus ``poses`` ([{id, tag, tfz}])
    and the annotation's per-component search record under ``components``."""
    solutions = []
    for block in re.split(r"^={4,}\s*$", text, flags=re.MULTILINE):
        block = block.strip()
        if not block:
            continue
        sol: Dict = {"poses": []}
        for line in block.splitlines():
            line = line.strip()
            if not line.startswith("phaserdag node "):
                continue
            rest = line[len("phaserdag node "):]
            parts = rest.split(None, 1)
            if len(parts) != 2:
                continue
            key, value = parts
            if key == "pose":
                tag = re.search(r'tag\s+"([^"]*)"', value)
                ident = re.search(r"\bid\s+(\d+)", value)
                tfz = _TFZ.search(value)
                sol["poses"].append({"id": ident.group(1) if ident else None,
                                     "tag": tag.group(1) if tag else None,
                                     "tfz": float(tfz.group(1)) if tfz else None})
                continue
            if value.startswith('"') and value.endswith('"'):
                value = value[1:-1]
            sol.setdefault(key, value)
        if len(sol) > 1 or sol["poses"]:
            sol["components"] = annotation_components(sol.get("annotation", ""))
            solutions.append(sol)
    return solutions


def annotation_components(annotation: str) -> List[Dict]:
    """[{tag, rfz, tfz, packs}] from the search narration, one per FIND."""
    found = []
    for match in _FIND.finditer(annotation or ""):
        tag, rest = match.group(1), match.group(2)
        rf = _RF.search(rest)
        tf = _TF.search(rest)
        found.append({"tag": tag,
                      "rfz": float(rf.group(1)) if rf else None,
                      "tfz": float(tf.group(1)) if tf else None,
                      "packs": tf.group(2) if tf else None})
    return found


def parse_best_solution(text: str) -> Optional[Dict]:
    """Rwork, Rfree (percent), dmin, poses and space group from the last
    ``Best Solution:`` (or ``Solution:``) line of a log."""
    lines = [l for l in text.splitlines() if "Best Solution:" in l] or \
            [l for l in text.splitlines() if l.strip().startswith("Solution:")]
    if not lines:
        return None
    match = _BEST.search(lines[-1])
    if not match:
        return None
    return {"rwork": float(match.group(1)), "rfree": float(match.group(2)),
            "dmin": float(match.group(3)), "poses": int(match.group(4)),
            "spacegroup": match.group(5)}


def _flat(obj, prefix=""):
    if isinstance(obj, dict):
        for key, value in obj.items():
            yield from _flat(value, prefix + "/" + key)
    else:
        yield prefix, obj


def read_result_json(path) -> Dict:
    """The run-level facts of a node's result.json: resolution, anomalous,
    Wilson B, the Matthews analysis, twinning, tNCS, wall time."""
    try:
        with open(path) as stream:
            data = json.load(stream)
    except (OSError, ValueError):
        return {}
    flat = {key.replace("/.", "/").rstrip("."): value for key, value in _flat(data)}

    def get(name):
        for key, value in flat.items():
            if key.endswith(name):
                return value
        return None

    matthews = get("voyager.matthews")
    chosen = None
    if isinstance(matthews, list) and matthews:
        chosen = max(matthews, key=lambda m: m.get("probability", 0))
    return {
        "resolution": get("voyager.notifications.hires"),
        "anomalous": get("voyager.data.anomalous"),
        "wilson_b": get("voyager.anisotropy.wilson_bfactor"),
        "matthews_z": chosen.get("z") if chosen else None,
        "matthews_vm": chosen.get("vm") if chosen else None,
        "matthews_probability": chosen.get("probability") if chosen else None,
        "twinning_indicated": get("voyager.twinning.indicated"),
        "tncs_indicated": get("voyager.tncs_result.indicated"),
        "wall_seconds": get("voyager.time.cumulative.wall"),
    }


def last_result_json(db_dir) -> Optional[Path]:
    nodes = sorted(p for p in Path(db_dir).iterdir() if p.is_dir() and p.name[:10].isdigit())
    for node in reversed(nodes):
        found = list(node.glob("*.result.json"))
        if found:
            return found[0]
    return None


def picard_record(db_dir, log_path=None) -> Dict:
    """Everything above, from a Picard database directory and the job log."""
    db_dir = Path(db_dir)
    record: Dict = {"solutions": [], "best": None, "run": {}}
    cards = sorted(db_dir.glob("best.*.dag.cards"))
    if cards:
        record["solutions"] = parse_dag_cards(cards[0].read_text(encoding="utf-8", errors="replace"))
    text = ""
    if log_path and os.path.exists(log_path):
        text = Path(log_path).read_text(encoding="utf-8", errors="replace")
    else:
        logs = sorted(db_dir.glob("*-rfac/solutions.log"))
        if logs:
            text = logs[-1].read_text(encoding="utf-8", errors="replace")
    record["best"] = parse_best_solution(text)
    result = last_result_json(db_dir)
    if result:
        record["run"] = read_result_json(result)
    return record


def _set(element, name, value, places=None):
    if value is None or value == "":
        return
    if places is not None:
        try:
            value = f"{float(value):.{places}f}"
        except (TypeError, ValueError):
            pass
    element.set(name, str(value).lower() if isinstance(value, bool) else str(value))


def program_xml(record: Dict) -> ET.Element:
    """program.xml for a Picard job: run facts, then the solutions."""
    root = ET.Element("phasertng_picard")
    run = record.get("run") or {}
    data = ET.SubElement(root, "Data")
    _set(data, "resolution", run.get("resolution"), 2)
    _set(data, "anomalous", run.get("anomalous"))
    _set(data, "wilson_b", run.get("wilson_b"), 1)
    _set(data, "wall_seconds", run.get("wall_seconds"), 0)
    composition = ET.SubElement(root, "Composition")
    _set(composition, "z", run.get("matthews_z"))
    _set(composition, "vm", run.get("matthews_vm"), 2)
    _set(composition, "probability", run.get("matthews_probability"), 2)
    _set(ET.SubElement(root, "Twinning"), "indicated", run.get("twinning_indicated"))
    _set(ET.SubElement(root, "TNCS"), "indicated", run.get("tncs_indicated"))

    solutions = record.get("solutions") or []
    best = record.get("best") or {}
    node = ET.SubElement(root, "Solutions")
    node.set("count", str(len(solutions)))
    if solutions:
        _set(node, "components_sought", solutions[0].get("seek_components"))
        _set(node, "complete", solutions[0].get("seek_complete"))
    for rank, sol in enumerate(solutions, start=1):
        element = ET.SubElement(node, "Solution")
        element.set("rank", str(rank))
        _set(element, "spacegroup", sol.get("hermann_mauguin"))
        _set(element, "cell", " ".join(sol.get("unitcell", "").split()))
        _set(element, "tfz", sol.get("zscore"), 2)
        _set(element, "llg", sol.get("llg"), 1)
        _set(element, "rfactor", sol.get("rfactor"), 2)
        if rank == 1 and best.get("rfree") is not None:
            _set(element, "rfree", best["rfree"], 2)
            _set(element, "dmin", best.get("dmin"), 2)
        _set(element, "packs", sol.get("packs"))
        _set(element, "twinned", sol.get("twinned"))
        _set(element, "tncs", sol.get("tncs_indicated"))
        _set(element, "scattering_fraction", _full_fs(sol), 3)
        components = ET.SubElement(element, "Components")
        poses = sol.get("poses") or []
        components.set("count", str(len(poses)))
        searched = {c["tag"]: c for c in sol.get("components") or []}
        for pose in poses:
            component = ET.SubElement(components, "Component")
            _set(component, "tag", pose.get("tag"))
            _set(component, "tfz", pose.get("tfz"), 2)
            search = searched.get(pose.get("tag"))
            if search:
                _set(component, "rfz", search.get("rfz"), 1)
                _set(component, "search_tfz", search.get("tfz"), 1)
                _set(component, "packs", search.get("packs"))
        annotation = sol.get("annotation")
        if annotation:
            ET.SubElement(element, "Annotation").text = annotation.strip()
    return root


def _full_fs(sol):
    """The 'fs' of the solution's 'full' line: fraction of scattering."""
    full = sol.get("full")
    if not full:
        return None
    match = re.search(r"\bfs\s+([\d.]+)", full)
    return float(match.group(1)) if match else None
