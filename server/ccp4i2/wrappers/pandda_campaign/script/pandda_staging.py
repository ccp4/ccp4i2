"""Stage a declared list of datasets into the PanDDA invocation contract's
input tree, and write the manifest that lets fan-out find its way back.

Design note: docs/pandda-campaign-design.md §3.1-§3.3, §4.7, §8.4, §12.

The tree is the one the Materia invocation contract specifies, unchanged::

    <staging>/
    ├── datasets/
    │   ├── xtal-0000/{final.pdb, final.mtz, dict.cif?}
    │   └── xtal-0001/...
    ├── Projects.csv          # "Dataset, Project", one row per xtal-NNNN
    └── manifest.json         # uuid-keyed record of what went where

Names are clean ``xtal-NNNN`` because PanDDA takes the crystal number from the
last run of digits in the directory name, and a project-derived name with a
date in it parses as ~20 million and is silently dropped by the range filter.
``xtal-NNNN`` *is* the dtag PanDDA reports; the user's label survives in
``Projects.csv`` and the manifest.

Bytes move through ``link_or_copy`` only: a hardlink when the two paths share
a filesystem, a copy otherwise, never a symlink (Windows privilege, and the
project export/move/restore code would each have to decide whether to follow
one). The two normalisers, ``prepare_mtz_for_pandda`` for the FreeR label and
``prepare_dict_for_pandda`` for bond orders, run on the way in, so a
relabelled MTZ or a re-spelled dictionary is a real file by construction.

``stage_datasets`` is pure path work: give it ``DatasetSpec`` records and a
directory. ``campaign_dataset_specs`` is the one function here that reads the
database; it builds those records for a campaign from the same
``collect_pandda_datasets`` the export endpoint uses, so the two agree on
which datasets a campaign has.
"""
import csv
import hashlib
import json
import logging
import os
import shutil
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from ccp4i2.lib.pandda_export import prepare_mtz_for_pandda
from .pandda_dict import prepare_dict_for_pandda

logger = logging.getLogger(f"ccp4i2:{__name__}")

MANIFEST_VERSION = 1
DATASETS_DIR = "datasets"
PROJECTS_CSV = "Projects.csv"
MANIFEST_JSON = "manifest.json"
MODEL_NAME = "final.pdb"
REFLECTIONS_NAME = "final.mtz"
DICT_NAME = "dict.cif"


@dataclass(frozen=True)
class DatasetSpec:
    """One dataset as the orchestrator's ``DATASETS`` list declares it.

    ``label`` is the user's DTAG: what ``Projects.csv`` records against the
    ``xtal-NNNN`` name and what a person recognises. Identity for fan-out is
    ``project_uuid``, never the label, because project names are renameable.
    """
    label: str
    xyzin: Path
    hklin: Path
    dictionary: Optional[Path] = None
    project_uuid: Optional[str] = None
    source_job_uuid: Optional[str] = None
    #: ``{'xyzin'|'hklin'|'dict': File uuid}`` where known.
    source_file_uuids: Dict[str, str] = field(default_factory=dict)


def xtal_name(index: int) -> str:
    return f"xtal-{index:04d}"


def link_or_copy(src, dst) -> str:
    """Put the bytes of ``src`` at ``dst``; returns ``'link'`` or ``'copy'``.

    Hardlinks the *resolved* source (a symlinked source would otherwise give
    a link to the link). Falls back to a copy whenever the link fails: a
    different filesystem, a filesystem without hardlinks, or Windows refusing.
    """
    src = Path(src).resolve()
    dst = Path(dst)
    if dst.exists() or dst.is_symlink():
        # os.link would fail on this and the copy fallback would then
        # overwrite silently; staging never means "replace".
        raise FileExistsError(f"refusing to overwrite {dst}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(src, dst)
        return "link"
    except OSError:
        shutil.copy2(src, dst)
        return "copy"


def sha256_of(path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _check_specs(specs: List[DatasetSpec]) -> None:
    if not specs:
        raise ValueError("no datasets to stage")
    seen = set()
    for spec in specs:
        if not spec.label or not str(spec.label).strip():
            raise ValueError("a dataset has an empty label")
        if spec.label in seen:
            raise ValueError(f"duplicate dataset label {spec.label!r}")
        seen.add(spec.label)
        for role, path in (("model", spec.xyzin), ("reflections", spec.hklin),
                           ("dictionary", spec.dictionary)):
            if path is None:
                continue
            if not Path(path).is_file():
                raise FileNotFoundError(
                    f"dataset {spec.label!r}: {role} file not found: {path}")


def _stage_file(src: Path, dst: Path, *, source_uuid: Optional[str],
                prepare=None) -> dict:
    """Stage one file, through ``prepare`` (a normaliser with the
    ``prepare_*_for_pandda`` signature) when given. Returns its manifest entry."""
    src = Path(src)
    prepared = Path(prepare(src, dst.parent)) if prepare else src
    rewritten = prepared != src
    if rewritten:
        prepared.replace(dst)
        how = "rewritten"
    else:
        how = link_or_copy(src, dst)
    return {
        "source": str(src),
        "source_file_uuid": source_uuid,
        "sha256": sha256_of(dst),
        "how": how,
    }


def write_projects_csv(rows: Iterable, path) -> None:
    """``Projects.csv`` exactly as the export endpoint writes it (and as the
    Reinspect ingest reads it): header ``Dataset, Project``, one row each."""
    with open(path, "w", newline="") as handle:
        handle.write("Dataset, Project\n")
        for xtal, label in rows:
            handle.write(f"{xtal}, {label}\n")


def read_projects_csv(path) -> List[tuple]:
    with open(path, newline="") as handle:
        reader = csv.reader(handle, skipinitialspace=True)
        header = next(reader)
        if [h.strip() for h in header] != ["Dataset", "Project"]:
            raise ValueError(f"{path}: unexpected header {header!r}")
        return [(row[0].strip(), row[1].strip()) for row in reader if row]


def write_manifest(manifest: dict, path) -> None:
    with open(path, "w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=False)
        handle.write("\n")


def read_manifest(path) -> dict:
    with open(path) as handle:
        manifest = json.load(handle)
    version = manifest.get("manifest_version")
    if version != MANIFEST_VERSION:
        raise ValueError(f"{path}: manifest_version {version!r}, "
                         f"expected {MANIFEST_VERSION}")
    return manifest


def stage_datasets(specs: Iterable[DatasetSpec], staging_root,
                   provenance: Optional[dict] = None) -> dict:
    """Build the contract's input tree under ``staging_root``.

    ``xtal-NNNN`` numbers follow the order of ``specs``, so the list as
    submitted is the record of the assignment. Refuses to stage over an
    existing ``datasets/`` directory: what to do with a previous attempt is
    the caller's decision, not something to guess here.

    Returns the manifest (also written to ``manifest.json``). ``provenance``
    is stored as given; the orchestrator fills it with what produced the
    run (§4.4) once that is known.
    """
    specs = list(specs)
    _check_specs(specs)
    staging_root = Path(staging_root)
    datasets_dir = staging_root / DATASETS_DIR
    if datasets_dir.exists():
        raise FileExistsError(f"staging tree already exists: {datasets_dir}")
    datasets_dir.mkdir(parents=True)

    entries = []
    for index, spec in enumerate(specs):
        xtal = xtal_name(index)
        dataset_dir = datasets_dir / xtal
        dataset_dir.mkdir()
        uuids = spec.source_file_uuids or {}
        files = {
            MODEL_NAME: _stage_file(
                spec.xyzin, dataset_dir / MODEL_NAME,
                source_uuid=uuids.get("xyzin")),
            REFLECTIONS_NAME: _stage_file(
                spec.hklin, dataset_dir / REFLECTIONS_NAME,
                source_uuid=uuids.get("hklin"),
                prepare=prepare_mtz_for_pandda),
        }
        if spec.dictionary is not None:
            files[DICT_NAME] = _stage_file(
                spec.dictionary, dataset_dir / DICT_NAME,
                source_uuid=uuids.get("dict"),
                prepare=prepare_dict_for_pandda)
        entries.append({
            "xtal": xtal,
            "label": spec.label,
            "project_uuid": spec.project_uuid,
            "source_job_uuid": spec.source_job_uuid,
            "files": files,
        })
        logger.info("staged %s as %s (%s)", spec.label, xtal,
                    ", ".join(f"{name}:{info['how']}" for name, info in files.items()))

    write_projects_csv(((e["xtal"], e["label"]) for e in entries),
                       staging_root / PROJECTS_CSV)
    manifest = {
        "manifest_version": MANIFEST_VERSION,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "datasets_dir": DATASETS_DIR,
        "projects_csv": PROJECTS_CSV,
        "provenance": dict(provenance or {}),
        "datasets": entries,
    }
    write_manifest(manifest, staging_root / MANIFEST_JSON)
    return manifest


def campaign_dataset_specs(group, skipped: Optional[list] = None) -> List[DatasetSpec]:
    """``DatasetSpec`` records for a campaign's member projects.

    The one database-reading function in this module. Uses the same
    ``collect_pandda_datasets`` as the export endpoint, so the datasets a
    campaign prefills into the orchestrator's ``DATASETS`` list are the
    datasets its ZIP export would contain. Members without a finished dimple,
    or whose dimple outputs are missing on disk, are skipped with a warning,
    exactly as the export skips them. ``skipped``, if given, collects
    ``(project_name, reason)`` for each, so a caller can say why.
    """
    from ccp4i2.db import models
    from ccp4i2.lib.pandda_export import (
        _find_dictionary_cif, collect_pandda_datasets)

    specs = []
    for project, dimple_job, acedrg_job in collect_pandda_datasets(group):
        if not dimple_job:
            logger.warning("campaign %s: %s has no finished dimple, skipped",
                           group.name, project.name)
            if skipped is not None:
                skipped.append((project.name, "no finished dimple job"))
            continue
        dimple_dir = Path(dimple_job.directory)
        xyzin = dimple_dir / MODEL_NAME
        hklin = dimple_dir / REFLECTIONS_NAME
        if not xyzin.is_file() or not hklin.is_file():
            logger.warning("campaign %s: dimple outputs missing for %s, skipped",
                           group.name, project.name)
            if skipped is not None:
                skipped.append((project.name, f"dimple job {dimple_job.number} has no final.pdb/final.mtz on disk"))
            continue
        uuids = {}
        for role, path in (("xyzin", xyzin), ("hklin", hklin)):
            row = models.File.objects.filter(job=dimple_job, name=path.name).first()
            if row is not None:
                uuids[role] = str(row.uuid)
        dict_cif = _find_dictionary_cif(acedrg_job)
        if dict_cif is not None:
            row = models.File.objects.filter(job=acedrg_job, name=dict_cif.name).first()
            if row is not None:
                uuids["dict"] = str(row.uuid)
        specs.append(DatasetSpec(
            label=project.name,
            xyzin=xyzin,
            hklin=hklin,
            dictionary=dict_cif,
            project_uuid=str(project.uuid),
            source_job_uuid=str(dimple_job.uuid),
            source_file_uuids=uuids,
        ))
    return specs
