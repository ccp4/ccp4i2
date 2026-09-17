"""Fetch a file from a public repository straight into a project.

The repository download modes the file selector offers (PDBe, RCSB, UniProt,
PDB-REDO) are client-driven: the browser fetches the bytes through the Next
proxy and uploads them. That is fine for a 2 MB mmCIF and wrong for a
cryo-EM map, which is hundreds of MB, arrives gzipped, and on a served
deployment would then have to be staged past the 100 MB body caps. So maps
are fetched here, on the server: download from the archive into the
project's import directory, gunzip on the way, set the subtype, and set the
parameter through the same import path an upload uses.

Only EMDB is implemented. The shape is repository-agnostic on purpose: a
``pdbe`` repository could carry large structure-factor files the same way.

Stdlib only (``urllib``, ``gzip``): this runs on the slim server too.
See docs/emdb-map-fetch-plan.md.
"""

import gzip
import json
import logging
import pathlib
import re
import shutil
import tempfile
import urllib.request
from dataclasses import dataclass, field
from typing import Optional

logger = logging.getLogger(f"ccp4i2:{__name__}")

EMDB_API = "https://www.ebi.ac.uk/emdb/api/entry"
EMDB_ARCHIVE = "https://ftp.ebi.ac.uk/pub/databases/emdb/structures"
USER_AGENT = "ccp4i2 (repository fetch)"

#: Kinds of map an EMDB entry can hold, each with the archive directory it
#: lives in and the CMapDataFile subtype it becomes.
EMDB_KINDS = {
    "map": {"directory": "map", "sub_type": 1, "label": "Main map"},
    "half_map": {"directory": "other", "sub_type": 5, "label": "Half map"},
    "mask": {"directory": "masks", "sub_type": 4, "label": "Mask"},
}

_FILE_NAME_RE = re.compile(r"^emd_\d+[A-Za-z0-9_.-]*\.(map|mrc|ccp4)(\.gz)?$")


class RepositoryError(Exception):
    """A request the fetch cannot honour, with the HTTP status it earns."""

    def __init__(self, status, message):
        super().__init__(message)
        self.status = status


# ---------------------------------------------------------------------------
# EMDB
# ---------------------------------------------------------------------------

def normalise_emdb_entry(text) -> str:
    """``11638``, ``EMD-11638``, ``emd_11638`` or ``EMD11638`` -> ``EMD-11638``."""
    match = re.fullmatch(r"\s*(?:emd[-_ ]?)?(\d{4,5})\s*", str(text or ""), re.IGNORECASE)
    if not match:
        raise RepositoryError(400, f"Not an EMDB entry: {text!r} (expected e.g. EMD-11638)")
    return f"EMD-{match.group(1)}"


def fetch_emdb_entry(entry: str, opener=None) -> dict:
    """The entry's JSON from the EMDB API. 404 becomes a RepositoryError."""
    url = f"{EMDB_API}/{entry}"
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    opener = opener or urllib.request.urlopen
    try:
        with opener(request, timeout=30) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as err:
        if err.code == 404:
            raise RepositoryError(404, f"{entry} is not an EMDB entry")
        raise RepositoryError(502, f"EMDB answered {err.code} for {entry}")
    except (urllib.error.URLError, OSError) as err:
        raise RepositoryError(502, f"Could not reach EMDB: {err}")


def _value(node, *keys, default=None):
    """Walk ``node[k1][k2]...`` tolerating missing keys and EMDB's
    ``{"valueOf_": x, "units": u}`` leaves."""
    for key in keys:
        if isinstance(node, dict) and key in node:
            node = node[key]
        elif isinstance(node, list) and isinstance(key, int) and len(node) > key:
            node = node[key]
        else:
            return default
    if isinstance(node, dict) and "valueOf_" in node:
        return node["valueOf_"]
    return node


def emdb_resolution(entry_json: dict):
    """The reported resolution in Å, or None."""
    value = _value(entry_json, "structure_determination_list", "structure_determination", 0,
                   "image_processing", 0, "final_reconstruction", "resolution")
    try:
        return float(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def _describe_map(node: dict) -> dict:
    """Pixel spacing, dimensions and the author contour of one map node."""
    spacing = _value(node, "pixel_spacing", "x")
    dims = [_value(node, "dimensions", axis) for axis in ("col", "row", "sec")]
    contour = None
    for item in _value(node, "contour_list", "contour", default=[]) or []:
        if item.get("primary", False) or contour is None:
            contour = item.get("level")
    return {
        "pixel_spacing": float(spacing) if spacing is not None else None,
        "dimensions": dims if all(d is not None for d in dims) else None,
        "size_kbytes": node.get("size_kbytes"),
        "contour_level": contour,
        "format": node.get("format"),
    }


def emdb_files(entry_json: dict) -> list:
    """The maps an entry actually has, as the client and the fetch see them.

    Each: ``{kind, file, sub_type, label, index, directory, pixel_spacing,
    dimensions, size_kbytes, contour_level}``. Main map first, then half maps
    in order, then masks. Entries without half maps or masks simply list
    fewer items; the caller must not assume any kind is present.
    """
    files = []
    main = entry_json.get("map") or {}
    if main.get("file"):
        files.append({"kind": "map", "file": main["file"], "index": 1, **_describe_map(main)})
    halves = _value(entry_json, "interpretation", "half_map_list", "half_map", default=[]) or []
    for i, node in enumerate(halves, start=1):
        if node.get("file"):
            files.append({"kind": "half_map", "file": node["file"], "index": i, **_describe_map(node)})
    masks = _value(entry_json, "interpretation", "segmentation_list", "segmentation", default=[]) or []
    for i, node in enumerate(masks, start=1):
        if node.get("file"):
            files.append({"kind": "mask", "file": node["file"], "index": i, **_describe_map(node)})
    for item in files:
        kind = EMDB_KINDS[item["kind"]]
        item["sub_type"] = kind["sub_type"]
        item["directory"] = kind["directory"]
        item["label"] = kind["label"] + (f" {item['index']}" if item["kind"] != "map" else "")
    return files


def emdb_entry_summary(entry_json: dict) -> dict:
    """What the fetch dialog shows: title, resolution, files, fitted models."""
    pdb_ids = [ref.get("pdb_id") for ref in
               (_value(entry_json, "crossreferences", "pdb_list", "pdb_reference", default=[]) or [])
               if ref.get("pdb_id")]
    return {
        "repository": "emdb",
        "entry": entry_json.get("emdb_id"),
        "title": _value(entry_json, "admin", "title"),
        "resolution": emdb_resolution(entry_json),
        "files": emdb_files(entry_json),
        "pdb_ids": pdb_ids,
    }


def emdb_archive_url(entry: str, item: dict) -> str:
    """The archive URL for one listed file. Built from the kind's directory,
    never from anything a client supplied beyond the entry and file name,
    and the file name must look like an EMDB map."""
    if not _FILE_NAME_RE.match(item["file"]):
        raise RepositoryError(400, f"Not an EMDB map file name: {item['file']}")
    return f"{EMDB_ARCHIVE}/{entry}/{item['directory']}/{item['file']}"


def emdb_annotation(entry: str, item: dict, entry_json: dict) -> str:
    """``EMD-11638 half map 1, 0.53 Å/px, 256³, 1.22 Å``."""
    parts = [f"{entry} {item['label'].lower()}"]
    if item.get("pixel_spacing"):
        parts.append(f"{item['pixel_spacing']:.2f} A/px")
    dims = item.get("dimensions")
    if dims:
        parts.append("x".join(str(d) for d in dims) if len(set(dims)) > 1 else f"{dims[0]}^3")
    resolution = emdb_resolution(entry_json)
    if resolution:
        parts.append(f"{resolution:g} A")
    return ", ".join(parts)


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def download_to(url: str, destination: pathlib.Path, gunzip: bool, opener=None,
                chunk_bytes: int = 4 * 1024 * 1024) -> pathlib.Path:
    """Stream ``url`` to ``destination``, decompressing when ``gunzip``.
    Never holds the file in memory."""
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    opener = opener or urllib.request.urlopen
    try:
        with opener(request, timeout=120) as response:
            source = gzip.GzipFile(fileobj=response) if gunzip else response
            with open(destination, "wb") as handle:
                shutil.copyfileobj(source, handle, chunk_bytes)
    except urllib.error.HTTPError as err:
        raise RepositoryError(502 if err.code != 404 else 404,
                              f"Archive answered {err.code} for {url}")
    except (urllib.error.URLError, OSError, EOFError) as err:
        raise RepositoryError(502, f"Download of {url} failed: {err}")
    return destination


# ---------------------------------------------------------------------------
# The fetch
# ---------------------------------------------------------------------------

@dataclass
class RepositoryFetch:
    repository: str
    entry: str
    file: str
    object_path: str
    sub_type: Optional[int] = None
    description: str = ""
    extra: dict = field(default_factory=dict)


def fetch_repository_file(job, spec: RepositoryFetch, *, fetch_entry=None, download=None) -> dict:
    """Fetch one listed file of a repository entry into the job's project and
    set ``spec.object_path`` to it. Returns what ``upload_file_param`` returns,
    plus ``source_url`` and ``annotation``.

    ``fetch_entry`` and ``download`` exist for tests; the defaults talk to EBI.
    """
    from .upload_param import ImportSpec, _LocalPathUpload, import_file_for_param

    if spec.repository != "emdb":
        raise RepositoryError(400, f"Unknown repository {spec.repository!r}; only 'emdb' is supported")
    entry = normalise_emdb_entry(spec.entry)
    entry_json = (fetch_entry or fetch_emdb_entry)(entry)
    listed = {item["file"]: item for item in emdb_files(entry_json)}
    item = listed.get(spec.file)
    if item is None:
        raise RepositoryError(400, f"{entry} does not list a file called {spec.file!r}")
    sub_type = int(spec.sub_type) if spec.sub_type is not None else item["sub_type"]
    if sub_type != item["sub_type"]:
        raise RepositoryError(
            400, f"{spec.file} is a {item['label'].lower()} (subtype {item['sub_type']}), not subtype {sub_type}")

    url = emdb_archive_url(entry, item)
    gunzip = spec.file.endswith(".gz")
    local_name = spec.file[:-3] if gunzip else spec.file

    # A scratch directory inside the project's import directory: the import
    # copies from there into place on the same filesystem, and the scratch
    # is removed whatever happens.
    import_dir = pathlib.Path(job.project.directory) / "CCP4_IMPORTED_FILES"
    import_dir.mkdir(parents=True, exist_ok=True)
    scratch = pathlib.Path(tempfile.mkdtemp(prefix=".fetch-", dir=import_dir))
    try:
        target = scratch / local_name
        logger.info("Fetching %s -> %s (gunzip=%s)", url, target, gunzip)
        (download or download_to)(url, target, gunzip)
        annotation = emdb_annotation(entry, item, entry_json)
        description = spec.description.strip() or f"Fetched from {url}"
        result = import_file_for_param(job, ImportSpec(
            object_path=spec.object_path,
            files=[_LocalPathUpload(target)],
            provenance_description=description,
            sub_type=sub_type,
            annotation=annotation,
        ))
        result.update({"source_url": url, "annotation": annotation, "entry": entry, "kind": item["kind"]})
        return result
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
