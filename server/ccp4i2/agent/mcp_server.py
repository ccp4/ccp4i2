"""CCP4i2 as MCP tools: the facade an agent drives the app through.

    CCP4I2_URL=http://127.0.0.1:<port>/api/ccp4i2 CCP4I2_TOKEN=<token> \\
        python -m ccp4i2.agent.mcp_server

The port and token are the ones Help > About shows (the token is needed
only where the server asks for one). Speaks MCP over stdio. It is a client
of the REST API and nothing else — no Django, no CCP4 — so the server's
validation decides everything, here as in the app, and it can run on any
Python with the ``mcp`` package (``pip install ccp4i2[agent]``).

Deliberately absent: deleting anything, changing a job's status, and any
file outside a project. See docs/agentic-knowledge.md.
"""
import json
import mimetypes
import os
import urllib.error
import urllib.parse
import urllib.request
import uuid
import xml.etree.ElementTree as ET
from contextvars import ContextVar
from pathlib import Path

from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError

STATUS = {0: "Unknown", 1: "Pending", 2: "Queued", 3: "Running", 4: "Interrupted",
          5: "Failed", 6: "Finished", 7: "Running remotely", 8: "File holder",
          9: "To delete", 10: "Unsatisfactory"}

INSTRUCTIONS = """\
CCP4i2 runs crystallographic tasks as jobs in projects. Work like this:

0. Reflection data from outside: import them first by running the TASK
   import_merged (create_job, upload_file, run_job; for unmerged data the
   task aimless_pipe), which checks them and makes or keeps the free-R set
   every later job needs. Models and sequences can be uploaded where used.
   Jobs are named by project and job NUMBER (as project_jobs lists them and
   [n].PARAM references use), in every tool.
1. describe_task(task) before using a task: its judgement says when to use
   it, which inputs need thought, how its result is judged, its traps and
   what comes next. A task with no judgement is one nobody has written up:
   be more careful with it. Prefer a pipeline to running its stages one by
   one (e.g. shelx or crank2, not ShelxCD then the rest): the pipeline holds
   the decisions between stages. Go stage by stage only as a fallback a
   judgement points you to.
2. create_job(project, task). A new job is filled in from the project's
   latest results, as in the app: read job_parameters(job, only_set=True)
   and check those are the files you mean.
3. Set inputs: set_file(job, path, "[n].PARAM[i]") for a file another job
   made or imported (n is a job number; -1 the latest), upload_file for a
   file from disk, set_parameter for anything else. Paths are as
   job_parameters lists them, e.g. inputData.XYZIN.
4. validate(job) and fix every error before run_job(job). Warnings are
   advice; say why you go on if you do.
5. wait_for_job(job) until it is Finished, Failed or Unsatisfactory; it
   returns after max_seconds with waited_out true, so if your own tool calls
   time out, give a max_seconds below that limit and call it again (jobs
   can run for an hour). Then judge_job(job). "Finished" means it ran, not
   that it worked: the judgement decides, and a draft judgement is a
   guide, not a rule. Do not argue a verdict away with the numbers it
   already weighed (a low R-free does not excuse geometry the verdict
   faulted); act on its advice, or report that you stopped short and why.
6. To change a finished job, clone_job it and edit the clone.
7. When a job fails, or its judgement cannot read its numbers, job_errors
   says what it recorded and shows its log. file_summary("[n].XYZOUT")
   says what a model holds (chains, residues built, UNK, ligands, heavy
   atoms) without downloading it.

Report what you ran, the numbers that decided each step, and anything you
were unsure of. Never delete or overwrite a user's work.
"""

from .request_state import request_state_security  # noqa: E402

server = MCPServer("ccp4i2", instructions=INSTRUCTIONS,
                   request_state_security=request_state_security())

# Served over HTTP by the app itself (agent/http.py), each call goes back to
# the same server's REST API with the caller's own Authorization; over stdio
# the environment says where the server is and what token to use.
_caller = ContextVar("ccp4i2_mcp_caller", default=None)


class ApiError(ToolError):
    """The server refused or failed: its reason, for the agent to read."""


def _base():
    caller = _caller.get()
    if caller is not None:
        return caller["base"].rstrip("/")
    return os.environ.get("CCP4I2_URL", "http://127.0.0.1:3421/api/ccp4i2").rstrip("/")


def _auth_headers():
    """Who the REST call is for. Over HTTP: the caller's own Authorization,
    or none (the API then refuses), never the server's token; and the
    caller's address as X-Forwarded-For, so the API does not see every agent
    as 127.0.0.1. Over stdio: CCP4I2_TOKEN, the command line's own."""
    caller = _caller.get()
    headers = {}
    if caller is not None:
        if caller.get("authorization"):
            headers["Authorization"] = caller["authorization"]
        if caller.get("forwarded_for"):
            headers["X-Forwarded-For"] = caller["forwarded_for"]
    elif os.environ.get("CCP4I2_TOKEN"):
        headers["Authorization"] = f"Bearer {os.environ['CCP4I2_TOKEN']}"
    return headers


def _request(method, path, body=None, query=None, content_type="application/json"):
    url = f"{_base()}/{path.strip('/')}/"
    if query:
        url += "?" + urllib.parse.urlencode(query, doseq=True)
    data = None
    headers = {}
    if body is not None:
        data = json.dumps(body).encode() if content_type == "application/json" else body
        headers["Content-Type"] = content_type
    headers.update(_auth_headers())
    req = urllib.request.Request(url, data=data, method=method, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=300) as response:
            payload = json.loads(response.read() or b"null")
    except urllib.error.HTTPError as err:
        try:
            payload = json.loads(err.read())
        except ValueError:
            raise ApiError(f"{method} {path}: HTTP {err.code}") from None
        raise ApiError(payload.get("error") or payload.get("detail") or str(payload)) from None
    if isinstance(payload, dict) and payload.get("success") is False:
        raise ApiError(payload.get("error") or str(payload))
    if isinstance(payload, dict) and "data" in payload and "success" in payload:
        return payload["data"]
    return payload


def _text(path):
    """A file the server serves as is (a log), as text; "" if there is none."""
    url = f"{_base()}/{path.lstrip('/')}"
    headers = _auth_headers()
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=60) as r:
            return r.read().decode("utf-8", "replace")
    except urllib.error.HTTPError:
        return ""


def _get(path, **query):
    return _request("GET", path, query={k: v for k, v in query.items() if v not in (None, "", False)})


def _post(path, body=None):
    return _request("POST", path, body=body if body is not None else {})


def _jid(project_id, job):
    """The server's id for job number ``job`` (e.g. "3", or a sub-job "3.1")
    of this project: a job outside the project cannot be reached by mistake."""
    number = str(job).strip()

    def find(jobs):
        for entry in jobs:
            if str(entry["number"]) == number:
                return entry["id"]
            found = find(entry.get("children") or [])
            if found is not None:
                return found
        return None

    found = find(_get(f"projects/{project_id}/job_tree")["job_tree"])
    if found is None:
        raise ApiError(f"No job {number} in project {project_id}: give the job NUMBER "
                       f"as project_jobs lists it")
    return found


def _job(job_id):
    return _get(f"jobs/{job_id}")


def _errors(xml_text):
    """An error report list as plain dicts (empty if there are none)."""
    if not xml_text:
        return []
    root = ET.fromstring(xml_text)
    out = []
    for report in root.iter("errorReport"):
        entry = {child.tag: (child.text or "").strip() for child in report}
        out.append({k: v for k, v in entry.items() if v})
    return out


# -- projects and tasks ------------------------------------------------------

@server.tool()
def list_projects() -> dict:
    """The projects, with their ids."""
    return {"projects": [{"id": p["id"], "name": p["name"]} for p in _get("projects")]}


@server.tool()
def create_project(name: str) -> dict:
    """Make a new project (its directory is chosen by the server)."""
    project = _request("POST", "projects", body={"name": name})
    return {"id": project["id"], "name": project["name"]}


@server.tool()
def list_tasks(query: str = "") -> dict:
    """The tasks that can be run, optionally only those whose name, title or
    description contains ``query``. Superseded and interactive (in-app
    window) tasks are left out. ``judged`` says whether the task has a
    judgement (describe_task gives it)."""
    query = query.lower()
    out = []
    for name, task in _get("task_lookup").items():
        if task.get("supersededBy") or task.get("interactive"):
            continue
        text = " ".join(str(task.get(k) or "") for k in ("TASKTITLE", "DESCRIPTION", "shortTitle"))
        if query and query not in f"{name} {text}".lower():
            continue
        out.append({"task": name, "title": task.get("TASKTITLE"), "description": task.get("DESCRIPTION"),
                    "judged": bool(task.get("hasJudgement"))})
    return {"tasks": out}


@server.tool()
def describe_task(task: str, full: bool = False) -> dict:
    """What to know before running a task: its title and its judgement —
    when to use it and when not, the inputs needing thought, the outcomes it
    can be judged to, its traps and next steps. Read this before create_job.
    A summary by default; full=True gives everything, including where each
    number is read and the basis of every threshold."""
    data = _get(f"agent/tasks/{task}")
    judgement = data.get("judgement")
    if full or not judgement:
        return data

    def short(text, n=400):
        text = " ".join(str(text or "").split())
        return text if len(text) <= n else text[:n].rsplit(" ", 1)[0] + " ..."

    data["judgement"] = {
        "status": judgement.get("status"),
        "purpose": short(judgement.get("purpose"), 600),
        "use_when": [short(u) for u in judgement.get("use_when") or []],
        "not_when": [{"text": short(n.get("text") if isinstance(n, dict) else n),
                      "instead": n.get("instead") if isinstance(n, dict) else None}
                     for n in judgement.get("not_when") or []],
        "inputs": [{"param": i.get("param"), "advice": short(i.get("advice"))}
                   for i in judgement.get("inputs") or []],
        "outcomes": sorted({str(v.get("outcome")) for v in judgement.get("verdict") or []
                            if v.get("outcome")}),
        "traps": [short(t, 300) for t in judgement.get("traps") or []],
        "next": [{"when": n.get("when"), "task": n.get("task"), "advice": short(n.get("advice"), 300)}
                 for n in judgement.get("next") or []],
        "note": "Summary: describe_task(task, full=True) for where each number is read and why.",
    }
    return data


@server.tool()
def project_jobs(project_id: int) -> dict:
    """A project's jobs, newest first: job (its number, what every job tool
    and every [n].PARAM reference takes), task, status, key numbers, and the
    files each made or imported, each with a file_id that set_file and file
    fields accept."""
    tree = _get(f"projects/{project_id}/job_tree")
    out = []
    for job in tree["job_tree"]:
        files = [{"param": f.get("job_param_name"), "name": f.get("name"),
                  "type": f.get("type"), "annotation": f.get("annotation") or None,
                  "file_id": f.get("uuid")}
                 for f in job.get("files") or []]
        out.append({"job": job["number"], "task": job["task_name"],
                    "title": job.get("title"), "status": STATUS.get(job["status"], job["status"]),
                    "kpis": {**(job.get("kpis") or {}).get("float_values", {}),
                             **(job.get("kpis") or {}).get("char_values", {})},
                    "files": files})
    return {"jobs": out}


# -- making and running a job ------------------------------------------------

@server.tool()
def create_job(project_id: int, task: str) -> dict:
    """Create a job of ``task`` in a project. Its inputs are filled in from
    the project's latest results, as in the app: check them with
    job_parameters(only_set=True)."""
    job = _post(f"projects/{project_id}/create_task", {"task_name": task})["new_job"]
    return {"project_id": job["project"], "job": job["number"], "task": job["task_name"]}


@server.tool()
def clone_job(project_id: int, job: str | int) -> dict:
    """Copy a job, with its inputs, as a new pending job to edit and run."""
    job_id = _jid(project_id, job)
    job = _post(f"jobs/{job_id}/clone")
    return {"project_id": job["project"], "job": job["number"], "task": job["task_name"]}


@server.tool()
def job_parameters(project_id: int, job: str | int, section: str = "", only_set: bool = False, query: str = "") -> dict:
    """A job's parameters, one entry each: path, label, value, whether set
    and required, choices, default; for a list, item_fields: what each item
    holds. ``section`` limits to inputData, controlParameters or outputData;
    ``query`` to paths or labels containing it. Use the path with
    set_parameter / set_file."""
    job_id = _jid(project_id, job)
    data = _get(f"jobs/{job_id}/parameters", section=section or None,
                only_set="1" if only_set else None, query=query or None)
    return {"parameters": data["parameters"]}


def _value_at(job_id, path):
    """A parameter's value now, for a path that may reach into a list item
    (inputData.ENSEMBLES[1].pdbItemList[0].structure)."""
    import re
    section, _, rest = path.partition(".")
    top = re.split(r"[.\[]", rest, maxsplit=1)[0]
    entries = _get(f"jobs/{job_id}/parameters", query=top)["parameters"]
    entry = next((e for e in entries if e["path"] == f"{section}.{top}"), None)
    if entry is None:
        return None, False
    value = entry.get("value")
    for step in re.findall(r"\[(\d+)\]|\.(\w+)", rest[len(top):]):
        index, key = step
        try:
            value = value[int(index)] if index else value.get(key)
        except (IndexError, KeyError, TypeError, AttributeError):
            return None, False
    return value, value is not None


def _resolve_files(value, project_id):
    """Replace each {"file": <reference or file id>} (or a value echoed from
    job_parameters, which carries fileId) with the file as the server takes it."""
    if isinstance(value, list):
        return [_resolve_files(v, project_id) for v in value]
    if not isinstance(value, dict):
        return value
    target = value.get("fileId") or value.get("file_id") or (
        value.get("file") if set(value) <= {"file", "annotation", "fileId", "file_id"} else None)
    if target is not None and set(value) <= {"file", "annotation", "fileId", "file_id"}:
        resolved = _get(f"projects/{project_id}/resolve_fileuse", fileuse=str(target))
        resolved.pop("fullPath", None)  # a path the client is not to hand back
        return resolved
    return {k: _resolve_files(v, project_id) for k, v in value.items()}


@server.tool()
def set_parameter(project_id: int, job: str | int, path: str,
                  value: str | int | float | bool | list | dict | None) -> dict:
    """Set a parameter of a pending job (path as job_parameters gives it).
    A list is set whole, as a JSON array of items shaped as its item_fields
    say; a file anywhere in a value is {"file": "<reference like [3].XYZOUT[0],
    or a file id>"}. For a single file parameter set_file is simpler."""
    job_id = _jid(project_id, job)
    project_id = _job(job_id)["project"]
    _post(f"jobs/{job_id}/set_parameter",
          {"object_path": path, "value": _resolve_files(value, project_id)})
    now, is_set = _value_at(job_id, path)
    return {"path": path, "now": now, "set": is_set}


@server.tool()
def set_file(project_id: int, job: str | int, path: str, reference: str) -> dict:
    """Set a file input of a pending job to a file already in the project:
    ``reference`` is a file_id from project_jobs, or "[n].PARAM" /
    "[n].PARAM[i]" (n a job NUMBER, or -1 the latest job, -2 the one before;
    PARAM the producing job's parameter), or "task[-1].PARAM" (the latest
    job of a task)."""
    job_id = _jid(project_id, job)
    job = _job(job_id)
    resolved = _get(f"projects/{job['project']}/resolve_fileuse", fileuse=reference)
    resolved.pop("fullPath", None)  # a path the client is not to hand back
    _post(f"jobs/{job_id}/set_parameter", {"object_path": path, "value": resolved})
    return {"path": path, "file": resolved.get("baseName"), "annotation": resolved.get("annotation")}


@server.tool()
def upload_file(project_id: int, job: str | int, path: str, local_path: str, column_labels: str = "") -> dict:
    """Import a file from this computer into a pending job's file input
    (an MTZ needs ``column_labels`` when it holds more than one data set,
    e.g. "/*/*/[FP,SIGFP]"). For a file inside a list item, set the list
    first with set_parameter (items without the file), then upload to the
    item's path, e.g. inputData.ENSEMBLES[0].pdbItemList[0].structure.
    Reflection data from outside the project: import them with import_merged
    instead, which checks them and makes the free-R set."""
    job_id = _jid(project_id, job)
    source = Path(local_path).expanduser()
    if not source.is_file():
        raise ApiError(f"no file {source}")
    boundary = uuid.uuid4().hex
    fields = {"objectPath": path}
    if column_labels:
        fields["column_selector"] = column_labels
    parts = []
    for key, value in fields.items():
        parts.append(f'--{boundary}\r\nContent-Disposition: form-data; name="{key}"\r\n\r\n{value}\r\n'.encode())
    kind = mimetypes.guess_type(source.name)[0] or "application/octet-stream"
    parts.append((f'--{boundary}\r\nContent-Disposition: form-data; name="file"; '
                  f'filename="{source.name}"\r\nContent-Type: {kind}\r\n\r\n').encode()
                 + source.read_bytes() + b"\r\n")
    parts.append(f"--{boundary}--\r\n".encode())
    _request("POST", f"jobs/{job_id}/upload_file_param", body=b"".join(parts),
             content_type=f"multipart/form-data; boundary={boundary}")
    now, _ = _value_at(job_id, path)
    return {"path": path, "now": now}


@server.tool()
def validate(project_id: int, job: str | int) -> dict:
    """Check a pending job before running it: the errors (which block the
    run) and warnings (advice). Includes the slower checks run at
    submission."""
    job_id = _jid(project_id, job)
    quick = _errors(_get(f"jobs/{job_id}/validation").get("xml"))
    try:
        slow = _errors(_get(f"jobs/{job_id}/run_time_validation").get("xml"))
    except ApiError as err:
        slow = [{"description": f"run-time validation failed: {err}", "severity": "4"}]
    seen, reports = set(), []
    for report in quick + slow:
        key = json.dumps(report, sort_keys=True)
        if key not in seen:
            seen.add(key)
            reports.append(report)

    def severity(report):
        # The server writes the name ("ERROR", "WARNING"); older reports a number
        text = str(report.get("severity", "0")).strip().upper()
        names = {"ERROR": 4, "WARNING": 2, "INFO": 1, "OK": 0}
        if text in names:
            return names[text]
        try:
            return int(text)
        except ValueError:
            return 4  # a severity nobody recognises is not safe to call advice
    return {"errors": [r for r in reports if severity(r) >= 4],
            "warnings": [r for r in reports if 2 <= severity(r) < 4]}


@server.tool()
def run_job(project_id: int, job: str | int) -> dict:
    """Start a pending job, then wait_for_job. As in the app, a job whose
    validation has errors is not started: the errors come back instead."""
    job_id = _jid(project_id, job)
    errors = validate(project_id, job)["errors"]
    if errors:
        raise ApiError("Not started: validation has errors. " + "; ".join(
            f"{e.get('objectPath') or e.get('name') or ''}: {e.get('description') or e.get('details') or ''}".strip(": ")
            for e in errors))
    _post(f"jobs/{job_id}/run")
    return job_status(project_id, job)


@server.tool()
def job_status(project_id: int, job: str | int) -> dict:
    """A job's status and key numbers."""
    job_id = _jid(project_id, job)
    job = _job(job_id)
    return {"project_id": job["project"], "job": job["number"], "task": job["task_name"],
            "status": STATUS.get(job["status"], job["status"]),
            "kpis": {**(job.get("float_values") or {}), **(job.get("char_values") or {})}}


TERMINAL = {"Finished", "Failed", "Unsatisfactory", "Interrupted"}


@server.tool()
def wait_for_job(project_id: int, job: str | int, max_seconds: int = 600) -> dict:
    """Wait until a job ends (Finished, Failed, Unsatisfactory, Interrupted)
    or max_seconds pass (at most 1800), then give its status. Call again if
    it is still running."""
    job_id = _jid(project_id, job)
    import time
    deadline = time.monotonic() + max(0, min(int(max_seconds), 1800))
    while True:
        status = job_status(project_id, job)
        if status["status"] in TERMINAL or time.monotonic() >= deadline:
            status["waited_out"] = status["status"] not in TERMINAL
            return status
        time.sleep(10)


@server.tool()
def job_errors(project_id: int, job: str | int, log_lines: int = 40) -> dict:
    """Why a job failed, or what it warned about: the errors and warnings it
    recorded, the end of its log, and the same for any of its steps
    (sub-jobs) that failed. Read this before trying again differently."""
    job_id = _jid(project_id, job)
    tree = _get(f"projects/{project_id}/job_tree")["job_tree"]

    def node(entries, number):
        for entry in entries:
            if str(entry["number"]) == number:
                return entry
            found = node(entry.get("children") or [], number)
            if found:
                return found
        return None

    def report(jid, number):
        reports = []
        for r in _errors(_get(f"jobs/{jid}/diagnostic_xml").get("xml")):
            keep = {k: r[k] for k in ("code", "description", "details", "severityName", "class")
                    if r.get(k)}
            if keep and keep not in reports:
                reports.append(keep)
        directory = "/".join(f"job_{n}" for n in str(number).split("."))
        log = _text(f"projects/{project_id}/files_by_path/CCP4_JOBS/{directory}/log.txt")
        tail = [line for line in log.splitlines() if line.strip()][-max(0, min(int(log_lines), 200)):]
        return {"job": number, "errors": reports, "log_tail": "\n".join(tail)}

    this = node(tree, str(job).strip()) or {}
    out = report(job_id, str(job).strip())
    out["status"] = STATUS.get(this.get("status"), this.get("status"))
    failed = [c for c in this.get("children") or [] if c.get("status") in (4, 5, 10)]
    out["failed_steps"] = [dict(report(c["id"], c["number"]), task=c["task_name"]) for c in failed]
    return out


@server.tool()
def judge_job(project_id: int, job: str | int) -> dict:
    """Did a finished job work? Reads the numbers the task's judgement names
    from the job's files, and gives the outcome, the reason, and the next
    steps the judgement suggests. ``missing`` lists numbers that could not be
    read; an outcome resting on missing numbers is not to be trusted."""
    job_id = _jid(project_id, job)
    return _get(f"jobs/{job_id}/judgement")


@server.tool()
def what_next(project_id: int, job: str | int) -> dict:
    """Next steps after a job: the judgement's (if the task has one) and the
    app's usual follow-on tasks."""
    job_id = _jid(project_id, job)
    verdict = judge_job(project_id, job)
    usual = _get(f"jobs/{job_id}/what_next").get("result", [])
    return {"from_judgement": verdict.get("next") or [],
            "usual_next_tasks": [t.get("taskName") for t in usual]}


@server.tool()
def file_summary(project_id: int, file: str) -> dict:
    """What a file holds, read by the server: for a model its chains, residue
    ranges and counts, sequence per chain (X where a residue is UNK),
    ligands, waters and elements (so: which residues were left unbuilt, and
    whether heavy atoms or a ligand are in it); for reflection data its
    cell, space group, resolution and columns; for a sequence file its
    sequences. ``file`` is a reference "[n].PARAM" (n a job number, -1 the
    latest; e.g. "[7].XYZOUT") or a file id from job_parameters."""
    resolved = _get(f"projects/{project_id}/resolve_fileuse", fileuse=str(file))
    uuid = resolved.get("dbFileId")
    if not uuid:
        raise ApiError(f"{file}: no file found")
    return {"file": {"name": resolved.get("baseName"), "annotation": resolved.get("annotation"),
                     "job_directory": resolved.get("relPath")},
            "summary": _get(f"files_by_uuid/{uuid}/digest")}


def main():
    server.run("stdio")


if __name__ == "__main__":
    main()
