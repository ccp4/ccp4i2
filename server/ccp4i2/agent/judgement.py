"""A task's judgement file, and the verdict it gives on a job.

The file is ``<task>.agent.yaml`` beside the task's def.xml (format:
docs/agentic-knowledge.md, section 3). ``judge`` reads the results it names
from a job's directory, and returns the first verdict whose condition holds.
It reads files only: no server, no CCP4.
"""
import xml.etree.ElementTree as ET
from pathlib import Path

import yaml

from ..core.tasks import locate_def_xml
from . import condition

DRAFT_NOTE = ("This judgement is a draft that no crystallographer has reviewed; "
              "treat its thresholds as a guide, and check them.")


def judgement_path(task_name):
    def_xml = locate_def_xml(task_name)
    if def_xml is None:
        return None
    return def_xml.with_name(f"{task_name}.agent.yaml")


def load(task_name=None, path=None):
    """The judgement for a task (or from a file), or None if it has none."""
    path = Path(path) if path else judgement_path(task_name)
    if path is None or not path.is_file():
        return None
    with open(path, encoding="utf-8") as stream:
        return yaml.safe_load(stream)


def _number(text):
    """A number as a program wrote it: "4.01" or "4.01%" (a percentage as text)."""
    text = str(text).strip()
    return float(text[:-1] if text.endswith("%") else text)


_TYPES = {"float": _number, "int": lambda s: int(_number(s)), "str": str, "string": str}


def read_result(spec, job_dir, kpis=None):
    """One named result from a job, converted to its type; None if absent."""
    convert = _TYPES.get(spec.get("type", "float"), float)
    source = spec.get("file", "program.xml")
    if source == "kpi":
        raw = (kpis or {}).get(spec.get("kpi"))
    else:
        path = Path(job_dir) / source
        if not path.is_file():
            return None
        # An xpath, or a list of them: the first that finds an element wins
        # (e.g. the R-free after waters were added, else before)
        xpaths = spec.get("xpath")
        xpaths = xpaths if isinstance(xpaths, list) else [xpaths]
        try:
            tree = ET.parse(path)
        except (ET.ParseError, OSError):
            return None
        node = None
        for xpath in xpaths:
            try:
                node = tree.find(xpath)
            except (SyntaxError, KeyError, TypeError):
                node = None
            if node is not None and (spec.get("attribute") or (node.text or "").strip()):
                break
            node = None
        if node is None:
            return None
        attribute = spec.get("attribute")
        raw = node.get(attribute) if attribute else node.text
    if raw is None:
        return None
    try:
        return convert(str(raw).strip()) if not isinstance(raw, (int, float)) else convert(raw)
    except ValueError:
        return None


def read_results(judgement, job_dir, kpis=None):
    return {name: read_result(spec, job_dir, kpis)
            for name, spec in (judgement.get("results") or {}).items()}


def _when(entry):
    """An entry's condition as text (YAML reads a bare true as a boolean)."""
    when = entry["when"]
    if isinstance(when, bool):
        return "true" if when else "false"
    return str(when)


def problems(judgement):
    """What is malformed in a judgement file (empty if nothing is)."""
    found = []
    results = judgement.get("results") or {}
    for name, spec in results.items():
        if spec.get("file", "program.xml") != "kpi" and not spec.get("xpath"):
            found.append(f"result {name}: no xpath")
        elif spec.get("xpath"):
            xpaths = spec["xpath"] if isinstance(spec["xpath"], list) else [spec["xpath"]]
            for xpath in xpaths:
                try:  # a path ElementTree cannot compile would read as missing, always
                    ET.fromstring("<x/>").find(xpath)
                except (SyntaxError, KeyError, TypeError) as err:
                    found.append(f"result {name}: xpath {xpath!r} is not one ElementTree "
                                 f"reads ({err}); keep to tags, /, //, [n], [last()], [@a='v'], [tag='v']")
        if spec.get("type", "float") not in _TYPES:
            found.append(f"result {name}: type {spec['type']!r} is not one of {sorted(_TYPES)}")
    for section in ("verdict", "next"):
        for i, entry in enumerate(judgement.get(section) or []):
            if "when" not in entry:
                found.append(f"{section}[{i}]: no condition")
                continue
            try:
                tree = condition.parse(_when(entry))
            except condition.ConditionError as err:
                found.append(f"{section}[{i}]: {err}")
                continue
            unknown = condition.names(tree) - set(results) - {"outcome"}
            if unknown:
                found.append(f"{section}[{i}]: unknown result {sorted(unknown)}")
            if section == "verdict" and tree != ("lit", True) and not entry.get("basis"):
                found.append(f"verdict[{i}]: a threshold with no basis")
    return found


def judge(task_name, job_dir, kpis=None, judgement=None):
    """The verdict on a job: its results, the outcome, and why."""
    judgement = judgement or load(task_name)
    if judgement is None:
        return {"task": task_name, "outcome": None,
                "note": f"No judgement has been written for {task_name}."}
    values = read_results(judgement, job_dir, kpis)
    optional = {n for n, spec in (judgement.get("results") or {}).items() if spec.get("optional")}
    verdict = {"task": task_name, "results": values, "outcome": None,
               "missing": sorted(n for n, v in values.items() if v is None and n not in optional)}
    for entry in judgement.get("verdict") or []:
        if condition.holds(_when(entry), values):
            verdict.update(outcome=entry.get("outcome"), because=_when(entry),
                           basis=entry.get("basis"))
            break
    with_outcome = dict(values, outcome=verdict["outcome"])
    verdict["next"] = [entry for entry in judgement.get("next") or []
                       if condition.holds(_when(entry), with_outcome)]
    if judgement.get("status") != "reviewed":
        verdict["note"] = DRAFT_NOTE
    return verdict
