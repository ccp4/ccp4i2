"""
Job / project bibliography.

Compiles the bibliographic references relevant to a job (its task + subjobs,
optionally input progenitors) or a whole project, from the MedLine reference
files that already power task reports.

Design (see docs/BIBLIOGRAPHY_PLAN.md):

- Source of truth is ``I2_TOP/references/{key}.medline.txt`` — the same MedLine
  files ``report.metadata.ReferenceGroup.loadFromMedLine`` consumes. We parse
  them directly here (no dependency on the heavier report layer).
- A task maps to one or more *citation keys* via ``TASK_CITES``. Most tasks cite
  themselves (identity). Pipelines that drive programs *internally* (not as i2
  subjobs) — e.g. crank2 runs shelxc/d/e, buccaneer, parrot… by direct program
  calls, so subjob-walking never sees them — declare the programs they cite.
  Thin variant tasks (phaser_MR_FRF, crank2_faest…) alias to their program's
  canonical key so we keep one physical file per program.
- References are deduped by PMID (falling back to a normalised title).
"""

import logging
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

from ccp4i2 import I2_TOP
from ccp4i2.core.CCP4Utils import findReferenceFile as find_reference_file
from ccp4i2.db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")

REFERENCES_DIR = I2_TOP / "references"

# ---------------------------------------------------------------------------
# Task -> citation-key expansion
# ---------------------------------------------------------------------------
# A task expands to the set of citation keys whose {key}.medline.txt files make
# up its bibliography. Absent from this map => the task cites itself (identity),
# which is correct for a plain wrapper. Present => use exactly these keys.
#
# Two uses:
#   1. MONOLITHIC pipelines that invoke programs INTERNALLY (as direct program
#      calls, not i2 subjobs) must enumerate the programs they run so those
#      citations surface — subjob-walking can't see them. These are crank2,
#      modelcraft, mrbump_basic, mrparse, slicendice, arcimboldo. (Toolchains
#      grounded in each pipeline's own source, not assumed.)
#   2. Thin variant tasks alias to the canonical program key (phaser_MR_FRF ->
#      phaser), so we don't duplicate the MedLine file per variant.
#
# NOTE: SUBJOB-COMPOSING pipelines (SubstituteLigand, import_merged, phaser_*
# _pipeline, dr_mr_modelbuild, aimless_pipe, servalcat_pipe, prosmart_refmac…)
# whose steps are each real i2 plugins are DELIBERATELY absent: task_names_for_job
# already burrows Job.children and unions their tasks, so their constituents
# surface honestly from what actually ran. Adding them here would double-count
# and risk drifting from the real execution graph.
# Canonical home is core.citations - Django-free, so report/ shares the very
# same map instead of keeping a second, always-poorer one of its own.
from ccp4i2.core.citations import TASK_CITES  # noqa: E402,F401

# Tasks that are pure CCP4i2 plumbing (import/provide/select/format shims) with
# no citable upstream program. Listed for documentation; they simply resolve to
# an empty bibliography (no {task}.medline.txt, not in TASK_CITES).
# Canonical home is core.citations - Django-free, so report/ can use it too.
from ccp4i2.core.citations import NON_CITABLE as _NON_CITABLE

# ---------------------------------------------------------------------------
# MedLine parsing
# ---------------------------------------------------------------------------


def _parse_medline(text: str) -> List[dict]:
    """Parse a MedLine reference file into a list of reference dicts.

    Mirrors report.metadata.ReferenceGroup.loadFromMedLine: records are split on
    'PMID- ', and each yields title (TI), source (SO), authors (AU, repeatable)
    and a link (URL/UR/doi). A record is kept only if it has a source or title.
    """
    refs: List[dict] = []
    for chunk in re.split(r"\nPMID- ", "\n" + text):
        pmid_m = re.match(r"\s*(\d{4,})", chunk)
        pmid = pmid_m.group(1) if pmid_m else None

        ti = re.search(r"TI  -(.*(?:\n      .*)*)", chunk)
        so = re.search(r"SO  -(.*)", chunk)
        authors = [a.strip() for a in re.findall(r"AU  -(.*)", chunk)]
        link_m = (
            re.search(r"URL -\s*(\S+)", chunk)
            or re.search(r"UR  -\s*(\S+)", chunk)
        )
        link = link_m.group(1).strip() if link_m else None
        if link is None:
            doi = re.search(r"AID - (10\.\S+) \[doi\]", chunk) or re.search(
                r"LID - (10\.\S+) \[doi\]", chunk
            )
            if doi:
                link = f"https://doi.org/{doi.group(1)}"

        title = None
        if ti:
            # collapse the MedLine hanging-indent continuation lines
            title = re.sub(r"\s+", " ", ti.group(1)).strip()
        source = so.group(1).strip() if so else None

        if not title and not source:
            continue
        refs.append(
            {
                "pmid": pmid,
                "title": title,
                "authors": authors,
                "source": source,
                "link": link,
            }
        )
    return refs


def _citation_keys_for_task(task_name: str) -> List[str]:
    """Expand a task name to its citation keys (identity if unmapped)."""
    if task_name in TASK_CITES:
        return TASK_CITES[task_name]
    return [task_name]


def _dedup_key(ref: dict) -> str:
    if ref.get("pmid"):
        return f"pmid:{ref['pmid']}"
    title = (ref.get("title") or "").lower()
    return "title:" + re.sub(r"[^a-z0-9]+", "", title)


def references_for_tasks(task_names: Iterable[str]) -> List[dict]:
    """Return deduped references for a set of task names.

    Each task expands via TASK_CITES to citation keys; each key is read from
    ``references/{key}.medline.txt`` (missing files are silently skipped — that
    just means "no reference for this program"). Returns plain dicts:
    ``{pmid, title, authors[], source, link}``.
    """
    keys: Set[str] = set()
    for task in task_names:
        keys.update(_citation_keys_for_task(task))

    seen: Set[str] = set()
    out: List[dict] = []
    for key in sorted(keys):
        # Case-insensitive: a report class names the same file by TASKNAME,
        # which differs in case from the citation key for AceDRG. Only Linux
        # notices, and it notices by silently finding nothing.
        path = find_reference_file(key)
        if path is None:
            if key not in _NON_CITABLE:
                logger.debug("No reference file for citation key %s", key)
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError as err:
            logger.warning("Failed to read reference file %s: %s", path, err)
            continue
        for ref in _parse_medline(text):
            dk = _dedup_key(ref)
            if dk in seen:
                continue
            seen.add(dk)
            ref = dict(ref)
            ref["cited_by"] = key
            out.append(ref)
    return out


# ---------------------------------------------------------------------------
# Task-name collectors (job / project)
# ---------------------------------------------------------------------------


def task_names_for_job(
    job: models.Job,
    include_subjobs: bool = True,
    include_progenitors: bool = False,
) -> Set[str]:
    """Collect the task names whose references make up this job's bibliography.

    - the job's own task
    - all descendant subjobs (Job.children, recursive) when include_subjobs
    - tasks of input-side FileUse progenitors when include_progenitors (the
      reverse of JobViewSet.dependent_jobs' consumer walk)
    """
    names: Set[str] = {job.task_name}

    if include_subjobs:
        frontier = [job]
        while frontier:
            current = frontier.pop()
            for child in current.children.all():
                if child.task_name:
                    names.add(child.task_name)
                frontier.append(child)

    if include_progenitors:
        _add_progenitor_task_names(job, names)

    return names


def _add_progenitor_task_names(job: models.Job, names: Set[str]) -> None:
    """Walk the input-file provenance graph backwards, adding producer tasks."""
    visited: Set[int] = {job.id}
    frontier = [job]
    while frontier:
        current = frontier.pop()
        input_files = models.File.objects.filter(
            file_uses__job=current, file_uses__role=models.FileUse.Role.IN
        ).distinct()
        for f in input_files:
            producer = f.job  # the job that owns (produced) the file
            if producer is None or producer.id in visited:
                continue
            visited.add(producer.id)
            if producer.task_name:
                names.add(producer.task_name)
            frontier.append(producer)


def task_names_for_project(project: models.Project) -> Set[str]:
    """All distinct task names of jobs in a project (includes subjobs, which are
    themselves Job rows)."""
    return {
        t
        for t in models.Job.objects.filter(project=project)
        .values_list("task_name", flat=True)
        .distinct()
        if t
    }


def bibliography_for_job(
    job: models.Job,
    include_subjobs: bool = True,
    include_progenitors: bool = False,
) -> List[dict]:
    """Compiled, deduped bibliography for a job."""
    return references_for_tasks(
        task_names_for_job(job, include_subjobs, include_progenitors)
    )


def bibliography_for_project(project: models.Project) -> List[dict]:
    """Compiled, deduped bibliography for a whole project."""
    return references_for_tasks(task_names_for_project(project))
