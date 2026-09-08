# -*- coding: utf-8 -*-
"""CCP4i2 data-interaction layer for code running inside Coot.

This module is the shared half of the Coot bridge: it knows how to talk
to the CCP4i2 REST API, how to turn a job's input_params.xml into a
"load plan", how to shape the project hierarchy into a toolkit-neutral
browse model, and how to hand out harvestable save paths in the job's
COOT_FILE_DROP directory. It contains NO GUI code and NO Coot calls:
those live in the per-generation renderers (coot1_gui.py for Coot 1.x /
GTK4, coot09_loader.py for Coot 0.9).

Hard constraints (do not break these):

* Runs on Python 2.7 (Coot 0.9's embedded interpreter) and Python 3.x
  (Coot 1.x). No f-strings, no pathlib, no annotations.
* Standard library only. Inside a CCP4-shipped Coot the ccp4-python
  site-packages are available, but a standalone Coot build only
  guarantees stdlib.
* No imports from the rest of ccp4i2. Coot 0.9 stubs load this file
  directly by path (imp.load_source), bypassing the package.

The connection contract - everything this module is given - is a set of
environment variables exported by the task wrapper when it launches
Coot (credentials/contact details plus the job identity, nothing else):

    CCP4I2_API_URL            e.g. http://127.0.0.1:43434 (no trailing /)
    CCP4I2_ACCESS_TOKEN       bearer token (Electron per-launch secret);
                              absent in open localhost dev mode
    CCP4I2_JOB_ID             database integer id of the launching job
    CCP4I2_JOB_UUID           uuid of the launching job
    CCP4I2_PROJECT_ID         database integer id of the project
    CCP4I2_PROJECT_UUID       uuid of the project
    CCP4I2_PROJECT_NAME       project name (labels only)
    CCP4I2_PROJECT_DIRECTORY  project root on this machine (local mode)
    CCP4I2_JOB_DIRECTORY      the job's work directory
    CCP4I2_DROP_DIR           where saves go to be harvested
"""

from __future__ import absolute_import, print_function

import glob
import json
import os
import re
import xml.etree.ElementTree as ET

try:  # Python 3
    from urllib.request import Request, urlopen
    from urllib.parse import urlencode
except ImportError:  # Python 2
    from urllib2 import Request, urlopen  # type: ignore
    from urllib import urlencode  # type: ignore

DEFAULT_TIMEOUT = 20.0
API_PREFIX = "/api/ccp4i2/"

# ---------------------------------------------------------------------------
# Configuration (the environment contract)
# ---------------------------------------------------------------------------


class BridgeConfig(object):
    """Reads the connection/identity contract from the environment."""

    def __init__(self, environ=None):
        env = environ if environ is not None else os.environ
        self.api_url = self._resolve_api_url(env)
        self.token = env.get("CCP4I2_ACCESS_TOKEN") or env.get(
            "CCP4I2_LOCAL_SESSION_TOKEN"
        )
        self.job_id = env.get("CCP4I2_JOB_ID")
        self.job_uuid = env.get("CCP4I2_JOB_UUID")
        self.project_id = env.get("CCP4I2_PROJECT_ID")
        self.project_uuid = env.get("CCP4I2_PROJECT_UUID")
        self.project_name = env.get("CCP4I2_PROJECT_NAME")
        self.project_directory = env.get("CCP4I2_PROJECT_DIRECTORY")
        self.job_directory = env.get("CCP4I2_JOB_DIRECTORY")
        self.drop_dir = env.get("CCP4I2_DROP_DIR")
        if not self.drop_dir and self.job_directory:
            self.drop_dir = os.path.join(self.job_directory, "COOT_FILE_DROP")

    @staticmethod
    def _resolve_api_url(env):
        url = env.get("CCP4I2_API_URL")
        if url:
            return url.rstrip("/")
        port = env.get("UVICORN_PORT")
        if port:
            return "http://127.0.0.1:{0}".format(port)
        # Plain `manage.py runserver` development default.
        return "http://127.0.0.1:8000"

    def has_job_context(self):
        return bool(self.job_id or self.job_uuid)

    def cache_dir(self):
        """Directory for files fetched over the API (no shared filesystem)."""
        base = self.job_directory or self.drop_dir or os.getcwd()
        return os.path.join(base, "COOT_DOWNLOADS")


# ---------------------------------------------------------------------------
# HTTP client
# ---------------------------------------------------------------------------


class BridgeError(Exception):
    pass


class CootBridgeClient(object):
    """Thin JSON client for the CCP4i2 REST API.

    ``opener`` is injectable for tests: any callable with the signature
    of ``urllib.request.urlopen`` (request, timeout=...).
    """

    def __init__(self, config=None, opener=None, timeout=DEFAULT_TIMEOUT):
        self.config = config or BridgeConfig()
        self._opener = opener or urlopen
        self._timeout = timeout

    # -- transport ----------------------------------------------------------

    def _url(self, path, params=None):
        url = self.config.api_url + API_PREFIX + path.lstrip("/")
        if params:
            url += "?" + urlencode(params)
        return url

    def _headers(self):
        headers = {"Accept": "application/json"}
        if self.config.token:
            headers["Authorization"] = "Bearer " + self.config.token
        return headers

    def get_raw(self, path, params=None):
        request = Request(self._url(path, params), headers=self._headers())
        response = self._opener(request, timeout=self._timeout)
        try:
            return response.read()
        finally:
            try:
                response.close()
            except Exception:
                pass

    def get_json(self, path, params=None):
        payload = json.loads(self.get_raw(path, params).decode("utf-8"))
        return unwrap_envelope(payload)

    # -- browse -------------------------------------------------------------

    def projects(self):
        return self.get_json("projects/")

    def project_by_uuid(self, project_uuid):
        wanted = str(project_uuid).replace("-", "").lower()
        for project in self.projects():
            have = str(project.get("uuid", "")).replace("-", "").lower()
            if have == wanted:
                return project
        return None

    def job_tree(self, project_id):
        """The project's job tree as a list of root job nodes.

        The endpoint answers {"job_tree": [...], "total_jobs": ...};
        tolerate a bare list too.
        """
        data = self.get_json("projects/{0}/job_tree/".format(project_id))
        if isinstance(data, dict):
            return data.get("job_tree") or []
        return data or []

    def files(self, job_id=None, file_type=None):
        params = {}
        if job_id is not None:
            params["job"] = job_id
        if file_type is not None:
            params["type"] = file_type
        return self.get_json("files/", params or None)

    # -- single files -------------------------------------------------------

    def file_local_path(self, file_id):
        """Absolute path of a file on this machine (desktop mode), or None."""
        try:
            data = self.get_json("files/{0}/file_path/".format(file_id))
            path = data.get("path")
            if path and os.path.exists(path):
                return path
        except Exception:
            pass
        return None

    def file_by_uuid(self, file_uuid):
        return self.get_json("files_by_uuid/{0}/".format(file_uuid))

    def download_file(self, file_id, destination):
        content = self.get_raw("files/{0}/download/".format(file_id))
        _ensure_dir(os.path.dirname(destination))
        handle = open(destination, "wb")
        try:
            handle.write(content)
        finally:
            handle.close()
        return destination

    # -- job data -----------------------------------------------------------

    def params_xml(self, job_id):
        data = self.get_json("jobs/{0}/params_xml/".format(job_id))
        return data.get("xml") if isinstance(data, dict) else data


def unwrap_envelope(payload):
    """CCP4i2 endpoints answer either bare JSON or the api_success
    envelope {"success": true, "data": ...}; accept both."""
    if isinstance(payload, dict) and "success" in payload:
        if not payload.get("success", False):
            raise BridgeError(str(payload.get("error", "API error")))
        return payload.get("data")
    return payload


# ---------------------------------------------------------------------------
# The load plan: job input_params.xml -> what Coot should load
# ---------------------------------------------------------------------------

# Parameter names understood by the Coot tasks, in load order, with the
# toolkit-neutral kind each maps to. Renderers dispatch on the kind, so
# adding an entry here lights it up in every Coot generation at once.
INPUT_PARAM_KINDS = (
    ("DICT", "dictionary"),          # dictionaries first: models may need them
    ("XYZIN_LIST", "coordinates"),
    ("XYZIN", "coordinates"),
    ("FPHIIN_LIST", "map_2fofc"),
    ("FPHIIN", "map_2fofc"),
    ("DELFPHIIN_LIST", "map_fofc"),
    ("DELFPHIIN", "map_fofc"),
    ("DELFPHIINANOM_LIST", "map_anom"),
    ("DELFPHIINANOM", "map_anom"),
)

_FILE_CHILD_TAGS = ("baseName", "relPath", "project", "dbFileId",
                    "subType", "annotation", "contentFlag")


def parse_input_params(xml_text):
    """Extract file entries from a params/input_params XML document.

    Returns a list of dicts:
        {"param", "kind", "base_name", "rel_path", "project",
         "db_file_id", "sub_type", "annotation"}
    tolerating both single-file parameters (<XYZIN><baseName>..) and
    lists (<XYZIN_LIST><CPdbDataFile><baseName>..).
    """
    root = ET.fromstring(xml_text)
    body = root.find(".//ccp4i2_body")
    if body is None:
        body = root
    input_data = body.find("inputData")
    if input_data is None:
        return []
    entries = []
    for param_name, kind in INPUT_PARAM_KINDS:
        element = input_data.find(param_name)
        if element is None:
            continue
        if element.find("baseName") is not None:
            items = [element]
        else:
            items = [child for child in list(element)
                     if child.find("baseName") is not None]
        for item in items:
            entry = _read_file_element(item)
            if entry is None:
                continue
            entry["param"] = param_name
            entry["kind"] = kind
            entries.append(entry)
    return entries


def _read_file_element(element):
    values = {}
    for tag in _FILE_CHILD_TAGS:
        child = element.find(tag)
        values[tag] = (child.text or "").strip() if child is not None else ""
    if not values["baseName"]:
        return None
    return {
        "base_name": values["baseName"],
        "rel_path": values["relPath"],
        "project": values["project"],
        "db_file_id": values["dbFileId"],
        "sub_type": values["subType"],
        "annotation": values["annotation"],
    }


def resolve_entry_path(entry, config, client=None):
    """Turn a parsed file entry into a local path, or None.

    Resolution ladder:
      1. Same project as the launching job (or no project recorded):
         join the local project directory with relPath/baseName.
      2. dbFileId known: ask the API for the file's local path
         (desktop, shared filesystem), else download it to the cache.
    """
    project = (entry.get("project") or "").replace("-", "").lower()
    own = (config.project_uuid or "").replace("-", "").lower()
    if config.project_directory and (not project or project == own):
        path = os.path.join(
            config.project_directory,
            *_split_rel(entry.get("rel_path", "")) + [entry["base_name"]]
        )
        if os.path.exists(path):
            return path
    db_file_id = entry.get("db_file_id")
    if db_file_id and client is not None:
        try:
            meta = client.file_by_uuid(db_file_id)
            file_id = meta.get("id") if isinstance(meta, dict) else None
            if file_id is not None:
                path = client.file_local_path(file_id)
                if path:
                    return path
                destination = os.path.join(
                    config.cache_dir(), entry["base_name"])
                return client.download_file(file_id, destination)
        except Exception:
            pass
    return None


def _split_rel(rel_path):
    if not rel_path:
        return []
    return [part for part in re.split(r"[\\/]+", rel_path) if part]


def load_plan(config, client=None):
    """Build the ordered list of {kind, path, label} for the launching job.

    Prefers the API (jobs/{id}/params_xml/); falls back to the
    input_params.xml sitting in the job directory when the API is not
    reachable, so a network-less session still loads its data.
    """
    xml_text = None
    if client is not None and config.job_id:
        try:
            xml_text = client.params_xml(config.job_id)
        except Exception:
            xml_text = None
    if xml_text is None and config.job_directory:
        for name in ("input_params.xml", "params.xml"):
            candidate = os.path.join(config.job_directory, name)
            if os.path.exists(candidate):
                handle = open(candidate, "r")
                try:
                    xml_text = handle.read()
                finally:
                    handle.close()
                break
    if not xml_text:
        return []
    plan = []
    for entry in parse_input_params(xml_text):
        path = resolve_entry_path(entry, config, client)
        if not path:
            continue
        label = entry.get("annotation") or entry["base_name"]
        plan.append({"kind": entry["kind"], "path": path, "label": label})
    return plan


# ---------------------------------------------------------------------------
# The browse model: project hierarchy -> toolkit-neutral menu structure
# ---------------------------------------------------------------------------

#: File types a Coot session can load, mapped to kinds. MTZ map
#: coefficients dispatch further on sub_type (1=2Fo-Fc, 2=Fo-Fc, 3=anom).
LOADABLE_TYPES = {
    "chemical/x-pdb": "coordinates",
    "chemical/x-cif": "coordinates",
    "chemical/x-mmcif": "coordinates",
    "application/CCP4-mtz-map": "map_coeffs",
    "application/CCP4-map": "map",
    "application/refmac-dictionary": "dictionary",
}

_MAP_SUBTYPE_KINDS = {1: "map_2fofc", 2: "map_fofc", 3: "map_anom"}

#: Order in which kinds must be loaded. Dictionaries first so that
#: coordinate parsing understands ligand geometry, then coordinates,
#: then maps - the ordering Moorhen's job loader uses (moorhen-wrapper
#: fetchJobFiles). A file whose kind is absent here sorts last.
_LOAD_PRIORITY = {
    "dictionary": 0,
    "coordinates": 1,
    "map_2fofc": 2,
    "map_fofc": 2,
    "map_anom": 2,
    "map": 2,
}


def classify_file(file_record):
    """Kind for a File API record, or None if Coot cannot load it."""
    kind = LOADABLE_TYPES.get(file_record.get("type"))
    if kind == "map_coeffs":
        sub_type = file_record.get("sub_type")
        try:
            sub_type = int(sub_type)
        except (TypeError, ValueError):
            sub_type = 1
        return _MAP_SUBTYPE_KINDS.get(sub_type, "map_2fofc")
    return kind


def order_for_load(files):
    """Stable-sort loadable-file dicts into safe load order (dicts first).

    ``files`` are {"kind", ...} records as produced by browse_model.
    Loading a coordinate file before its ligand dictionary leaves the
    ligand without geometry, so dictionaries must precede coordinates.
    """
    return sorted(
        files, key=lambda record: _LOAD_PRIORITY.get(record.get("kind"), 99))


def filter_rows(rows, query, key):
    """Rows whose ``key(row)`` contains ``query`` (case-insensitive).

    Empty/blank query returns all rows. Shared by both browsers so the
    filter behaves identically on Coot 1.x and 0.9.
    """
    needle = (query or "").strip().lower()
    if not needle:
        return list(rows)
    return [row for row in rows if needle in (key(row) or "").lower()]


def display_label(label, source_project_name, own_project_name):
    """Molecule label for a loaded file, prefixed with the project name
    when it comes from a project other than the launching job's - the
    legacy patchMoleculeName cross-project cue - so two jobs' outputs
    stay distinguishable in Coot's display manager.
    """
    if source_project_name and source_project_name != own_project_name:
        return "{0}: {1}".format(source_project_name, label)
    return label


def browse_model(job_tree, include_empty=False):
    """Flatten a projects/{id}/job_tree/ response into browser rows.

    Returns a list of job dicts, depth-first (a sub-job follows its
    parent), each:
        {"job_id", "job_uuid", "label", "status", "depth",
         "files": [{"file_id", "file_uuid", "label", "kind"}]}
    Jobs with no loadable files are skipped unless include_empty.
    """
    rows = []

    def visit(node, depth):
        files = []
        for record in node.get("files", []):
            kind = classify_file(record)
            if kind is None:
                continue
            files.append({
                "file_id": record.get("id"),
                "file_uuid": record.get("uuid"),
                "label": record.get("annotation") or record.get("name"),
                "kind": kind,
            })
        files = order_for_load(files)
        if files or include_empty:
            rows.append({
                "job_id": node.get("id"),
                "job_uuid": node.get("uuid"),
                "label": "{0}: {1}".format(
                    node.get("number"), node.get("title") or
                    node.get("task_name")),
                "status": node.get("status"),
                "depth": depth,
                "files": files,
            })
        for child in node.get("children", []):
            visit(child, depth + 1)

    for node in job_tree:
        visit(node, 0)
    return rows


# ---------------------------------------------------------------------------
# Saving: the COOT_FILE_DROP harvest contract
# ---------------------------------------------------------------------------

_OUTPUT_RE = re.compile(r"output(\d+)\.(pdb|cif)$")


def next_output_number(drop_dir):
    """1 + the highest output<N>.pdb|cif already in the drop directory."""
    highest = 0
    for path in glob.glob(os.path.join(drop_dir, "output*.pdb")) + \
            glob.glob(os.path.join(drop_dir, "output*.cif")):
        match = _OUTPUT_RE.search(os.path.basename(path))
        if match:
            highest = max(highest, int(match.group(1)))
    return highest + 1


def output_path(drop_dir, number, extension="pdb"):
    _ensure_dir(drop_dir)
    return os.path.join(
        drop_dir, "output{0}.{1}".format(number, extension.lstrip(".")))


def harvestable_outputs(drop_dir):
    """(number, path) pairs in the drop directory, in save order."""
    found = []
    for path in glob.glob(os.path.join(drop_dir, "output*.pdb")) + \
            glob.glob(os.path.join(drop_dir, "output*.cif")):
        match = _OUTPUT_RE.search(os.path.basename(path))
        if match:
            found.append((int(match.group(1)), path))
    found.sort()
    return found


def _ensure_dir(path):
    if path and not os.path.isdir(path):
        try:
            os.makedirs(path)
        except OSError:
            pass
