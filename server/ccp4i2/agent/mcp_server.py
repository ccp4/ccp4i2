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

0. Reflection data from outside: import them first with import_merged,
   which checks them and makes (or keeps) the free-R set every later job
   needs; models and sequences can be uploaded where they are used.
1. describe_task(task) before using a task: its judgement says when to use
   it, which inputs need thought, how its result is judged, its traps and
   what comes next. A task with no judgement is one nobody has written up:
   be more careful with it.
2. create_job(project, task). A new job is filled in from the project's
   latest results, as in the app: read job_parameters(job, only_set=True)
   and check those are the files you mean.
3. Set inputs: set_file(job, path, "[n].PARAM[i]") for a file another job
   made or imported (n is a job number; -1 the latest), upload_file for a
   file from disk, set_parameter for anything else. Paths are as
   job_parameters lists them, e.g. inputData.XYZIN.
4. validate(job) and fix every error before run_job(job). Warnings are
   advice; say why you go on if you do.
5. wait_for_job(job) until it is Finished, Failed or Unsatisfactory.
   Then judge_job(job). "Finished" means it ran, not that it worked: the
   judgement decides, and a draft judgement is a guide, not a rule.
6. To change a finished job, clone_job it and edit the clone.

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


def _request(method, path, body=None, query=None, content_type="application/json"):
    url = f"{_base()}/{path.strip('/')}/"
    if query:
        url += "?" + urllib.parse.urlencode(query, doseq=True)
    data = None
    headers = {}
    if body is not None:
        data = json.dumps(body).encode() if content_type == "application/json" else body
        headers["Content-Type"] = content_type
    caller = _caller.get()
    if caller is not None:
        if caller.get("authorization"):
            headers["Authorization"] = caller["authorization"]
    elif os.environ.get("CCP4I2_TOKEN"):
        headers["Authorization"] = f"Bearer {os.environ['CCP4I2_TOKEN']}"
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


def _get(path, **query):
    return _request("GET", path, query={k: v for k, v in query.items() if v not in (None, "", False)})


def _post(path, body=None):
    return _request("POST", path, body=body if body is not None else {})


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
    window) tasks are left out."""
    query = query.lower()
    out = []
    for name, task in _get("task_lookup").items():
        if task.get("supersededBy") or task.get("interactive"):
            continue
        text = " ".join(str(task.get(k) or "") for k in ("TASKTITLE", "DESCRIPTION", "shortTitle"))
        if query and query not in f"{name} {text}".lower():
            continue
        out.append({"task": name, "title": task.get("TASKTITLE"), "description": task.get("DESCRIPTION")})
    return {"tasks": out}


@server.tool()
def describe_task(task: str) -> dict:
    """What to know before running a task: its title, and its judgement —
    when to use it and when not, the inputs needing thought, how its results
    are judged, its traps and next steps. Read this before create_job."""
    return _get(f"agent/tasks/{task}")


@server.tool()
def project_jobs(project_id: int) -> dict:
    """A project's jobs, newest first: id (what job tools take), number (what
    references use), task, status, key numbers, and the files each made or
    imported, each with a file_id that set_file and file fields accept."""
    tree = _get(f"projects/{project_id}/job_tree")
    out = []
    for job in tree["job_tree"]:
        files = [{"param": f.get("job_param_name"), "name": f.get("name"),
                  "type": f.get("type"), "annotation": f.get("annotation") or None,
                  "file_id": f.get("uuid")}
                 for f in job.get("files") or []]
        out.append({"number": job["number"], "id": job["id"], "task": job["task_name"],
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
    return {"id": job["id"], "number": job["number"], "task": job["task_name"]}


@server.tool()
def clone_job(job_id: int) -> dict:
    """Copy a job, with its inputs, as a new pending job to edit and run."""
    job = _post(f"jobs/{job_id}/clone")
    return {"id": job["id"], "number": job["number"], "task": job["task_name"]}


@server.tool()
def job_parameters(job_id: int, section: str = "", only_set: bool = False, query: str = "") -> dict:
    """A job's parameters, one entry each: path, label, value, whether set
    and required, choices, default; for a list, item_fields: what each item
    holds. ``section`` limits to inputData, controlParameters or outputData;
    ``query`` to paths or labels containing it. Use the path with
    set_parameter / set_file."""
    data = _get(f"jobs/{job_id}/parameters", section=section or None,
                only_set="1" if only_set else None, query=query or None)
    return {"parameters": data["parameters"]}


def _value_at(job_id, path):
    """A parameter's value now, for a path that may reach into a list item
    (inputData.ENSEMBLES[1].pdbItemList[0].structure)."""
    import re
    section, _, rest = path.partition(".")
    top = re.split(r"[.\[]", rest, 1)[0]
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
def set_parameter(job_id: int, path: str,
                  value: str | int | float | bool | list | dict | None) -> dict:
    """Set a parameter of a pending job (path as job_parameters gives it).
    A list is set whole, as a JSON array of items shaped as its item_fields
    say; a file anywhere in a value is {"file": "<reference like [3].XYZOUT[0],
    or a file id>"}. For a single file parameter set_file is simpler."""
    project_id = _job(job_id)["project"]
    _post(f"jobs/{job_id}/set_parameter",
          {"object_path": path, "value": _resolve_files(value, project_id)})
    now, is_set = _value_at(job_id, path)
    return {"path": path, "now": now, "set": is_set}


@server.tool()
def set_file(job_id: int, path: str, reference: str) -> dict:
    """Set a file input of a pending job to a file already in the project:
    ``reference`` is a file_id from project_jobs, or "[n].PARAM" /
    "[n].PARAM[i]" (n a job NUMBER, or -1 the latest job, -2 the one before;
    PARAM the producing job's parameter), or "task[-1].PARAM" (the latest
    job of a task)."""
    job = _job(job_id)
    resolved = _get(f"projects/{job['project']}/resolve_fileuse", fileuse=reference)
    resolved.pop("fullPath", None)  # a path the client is not to hand back
    _post(f"jobs/{job_id}/set_parameter", {"object_path": path, "value": resolved})
    return {"path": path, "file": resolved.get("baseName"), "annotation": resolved.get("annotation")}


@server.tool()
def upload_file(job_id: int, path: str, local_path: str, column_labels: str = "") -> dict:
    """Import a file from this computer into a pending job's file input
    (an MTZ needs ``column_labels`` when it holds more than one data set,
    e.g. "/*/*/[FP,SIGFP]"). For a file inside a list item, set the list
    first with set_parameter (items without the file), then upload to the
    item's path, e.g. inputData.ENSEMBLES[0].pdbItemList[0].structure.
    Reflection data from outside the project: import them with import_merged
    instead, which checks them and makes the free-R set."""
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
def validate(job_id: int) -> dict:
    """Check a pending job before running it: the errors (which block the
    run) and warnings (advice). Includes the slower checks run at
    submission."""
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
def run_job(job_id: int) -> dict:
    """Start a pending job. Poll job_status until it ends."""
    _post(f"jobs/{job_id}/run")
    return job_status(job_id)


@server.tool()
def job_status(job_id: int) -> dict:
    """A job's status and key numbers."""
    job = _job(job_id)
    return {"id": job["id"], "number": job["number"], "task": job["task_name"],
            "status": STATUS.get(job["status"], job["status"]),
            "kpis": {**(job.get("float_values") or {}), **(job.get("char_values") or {})}}


TERMINAL = {"Finished", "Failed", "Unsatisfactory", "Interrupted"}


@server.tool()
def wait_for_job(job_id: int, max_seconds: int = 600) -> dict:
    """Wait until a job ends (Finished, Failed, Unsatisfactory, Interrupted)
    or max_seconds pass (at most 1800), then give its status. Call again if
    it is still running."""
    import time
    deadline = time.monotonic() + max(0, min(int(max_seconds), 1800))
    while True:
        status = job_status(job_id)
        if status["status"] in TERMINAL or time.monotonic() >= deadline:
            status["waited_out"] = status["status"] not in TERMINAL
            return status
        time.sleep(10)


@server.tool()
def judge_job(job_id: int) -> dict:
    """Did a finished job work? Reads the numbers the task's judgement names
    from the job's files, and gives the outcome, the reason, and the next
    steps the judgement suggests. ``missing`` lists numbers that could not be
    read; an outcome resting on missing numbers is not to be trusted."""
    return _get(f"jobs/{job_id}/judgement")


@server.tool()
def what_next(job_id: int) -> dict:
    """Next steps after a job: the judgement's (if the task has one) and the
    app's usual follow-on tasks."""
    verdict = judge_job(job_id)
    usual = _get(f"jobs/{job_id}/what_next").get("result", [])
    return {"from_judgement": verdict.get("next") or [],
            "usual_next_tasks": [t.get("taskName") for t in usual]}


def main():
    server.run("stdio")


if __name__ == "__main__":
    main()
