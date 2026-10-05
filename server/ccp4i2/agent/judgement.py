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

DRAFT_NOTE = ("A draft no crystallographer has reviewed, written from a handful of runs "
              "on CCP4i2's test projects (Gamma, MDM2, BetaBlip, Thaumatin and a few "
              "more). Each threshold's basis says what it rests on: where that is one or "
              "two jobs, treat the verdict as a guide and weigh the numbers yourself.")


# How an outcome reads at a glance, for the app: good (the job did what it
# was for), caution (usable, but something to check or do next), bad (it did
# not work), unknown (not judged). A verdict entry may set its own "tone".
TONES = {
    "good": {"acceptable", "added", "already_present", "analysed_merged", "built",
             "imported", "improved", "models_found", "plausible", "refined",
             "sites_found", "solved", "trimmed", "usable"},
    "caution": {"ambiguous", "ambiguous_copies", "check_chains", "check_symmetry",
                "free_set_not_kept", "free_set_remade", "low_solvent", "needs_attention",
                "needs_building", "no_better", "no_better_than_input", "not_improving",
                "nothing_to_add", "overpacked", "partial", "phased", "placed", "tncs",
                "too_sparse", "twinning_suspected", "uncertain", "unchecked",
                "unconverged"},
    "bad": {"empty", "failed", "no_free_set", "none_found", "not_fitted", "poor",
            "unusable"},
}


def tone(outcome, entry=None):
    if entry and entry.get("tone"):
        return entry["tone"]
    for name, outcomes in TONES.items():
        if outcome in outcomes:
            return name
    return "unknown"


def pin_references(steps, task_name, number):
    """Next steps with references to this task's latest job ("crank2[-1].X")
    pinned to this job ("[5].X"). An agent working forward wants the latest
    job; the app, showing an older job's judgement, means that job."""
    import copy
    latest = f"{task_name}[-1]"
    pinned = copy.deepcopy(steps)
    for step in pinned:
        for key, value in (step.get("inputs") or {}).items():
            if isinstance(value, str) and value.startswith(latest):
                step["inputs"][key] = f"[{number}]" + value[len(latest):]
    return pinned


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
            if section == "verdict" and entry.get("tone") and entry["tone"] not in (*TONES, "unknown"):
                found.append(f"verdict[{i}]: tone {entry['tone']!r} is not one of "
                             f"{sorted((*TONES, 'unknown'))}")
            if section == "next" and entry.get("rerun") and entry.get("task") not in (
                    None, judgement.get("task")):
                found.append(f"next[{i}]: a rerun is of this task, not {entry['task']!r}")
            if section == "next" and entry.get("rerun") and not entry.get("inputs"):
                found.append(f"next[{i}]: a rerun with nothing changed")
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
    fired = None
    for entry in judgement.get("verdict") or []:
        if condition.holds(_when(entry), values):
            fired = entry
            verdict.update(outcome=entry.get("outcome"), because=_when(entry),
                           basis=entry.get("basis"),
                           clauses=condition.clauses(_when(entry), values))
            break
    verdict["tone"] = tone(verdict["outcome"], fired)
    verdict["meanings"] = {name: " ".join(str(spec.get("meaning") or "").split())
                           for name, spec in (judgement.get("results") or {}).items()}
    with_outcome = dict(values, outcome=verdict["outcome"])
    steps = []
    for entry in judgement.get("next") or []:
        if condition.holds(_when(entry), with_outcome):
            step = dict(entry)
            if step.get("rerun"):
                # the same task again, as a clone with these inputs changed
                step["task"] = task_name
            steps.append(step)
    verdict["next"] = steps
    if judgement.get("status") != "reviewed":
        verdict["note"] = DRAFT_NOTE
    return verdict


CACHE_NAME = "judgement.json"


def judgement_version(task_name):
    """A short hash of the task's judgement file and of the code that
    evaluates it (this module and the condition language), or None without
    a judgement: the verdict on a finished job changes only when one of
    those does."""
    import hashlib
    from pathlib import Path
    path = judgement_path(task_name)
    if path is None or not path.is_file():
        return None
    digest = hashlib.sha256(path.read_bytes())
    for engine in (Path(__file__), Path(__file__).with_name("condition.py")):
        digest.update(engine.read_bytes())
    return digest.hexdigest()[:12]


def judge_finished(task_name, job_dir, kpis=None):
    """The verdict on a finished job, worked out once and kept in the job's
    directory as judgement.json.

    A finished job's files do not change, so its verdict changes only when
    its judgement does; the cache records the judgement's version and is
    used only while that is unchanged (the judgements are drafts, edited
    often). It is also the record of which version of a judgement said what
    about a job. A job still running is judged afresh each time, and not
    kept.
    """
    import json
    from pathlib import Path
    version = judgement_version(task_name)
    cache = Path(job_dir) / CACHE_NAME
    try:
        kept = json.loads(cache.read_text())
        if kept.get("judgement_version") == version and version is not None:
            return kept["verdict"]
    except (OSError, ValueError, KeyError, TypeError):
        pass
    verdict = judge(task_name, job_dir, kpis=kpis)
    verdict["judgement_version"] = version
    if version is not None:
        try:
            cache.write_text(json.dumps({"judgement_version": version, "verdict": verdict},
                                        indent=1, default=str))
        except OSError:
            pass  # a read-only project still gets its verdict
    return verdict
