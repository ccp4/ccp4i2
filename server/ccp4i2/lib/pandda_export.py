"""
PANDDA export service.

Builds a PANDDA-ready ZIP from a campaign's member projects, collecting
``final.pdb`` + ``final.mtz`` from each member's most recent finished
i2Dimple job and (optionally) the dictionary CIF from its most recent
finished LidiaAcedrgNew job. MTZs whose FreeR column carries a label
PANDDA2 doesn't recognise are rewritten with an accepted label.

This module is deliberately Django-light: it accepts a
``ProjectGroup`` instance and reads the related ``Job`` table, but
contains no request/response handling. That makes it callable from a
synchronous ViewSet today and from a background worker tomorrow
without further refactoring.
"""
import logging
import os
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import gemmi

from ..db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")


PANDDA_FREER_LABELS = ("FREE", "FreeR_flag", "R-free-flags")
COMMON_FREER_LABELS = (
    "FREER",
    "FREE",
    "FreeR_flag",
    "FreeRflag",
    "FreeR",
    "R-free-flags",
    "free",
)

DIMPLE_TASK_NAMES = ("i2Dimple", "dimple")
ACEDRG_TASK_NAMES = ("LidiaAcedrgNew", "acedrg")
DICT_CIF_CANDIDATES = ("LIG.cif", "DRG.cif")


@dataclass
class PanddaExportResult:
    zip_path: Path
    temp_dir: Path
    included_count: int


def prepare_mtz_for_pandda(src_path: Path, staging_dir: Path) -> Path:
    """Return a path to an MTZ whose FreeR column carries a PANDDA-accepted label.

    PANDDA2 only recognises 'FREE', 'FreeR_flag', or 'R-free-flags'. If the
    source MTZ already has one of those, the original path is returned (no
    rewrite). Otherwise the first column matching a common FreeR convention
    (e.g. CCP4's 'FREER') is renamed to 'FreeR_flag' and the rewritten file
    is placed in ``staging_dir``. If no FreeR-like column is found at all,
    the source path is returned unchanged.
    """
    mtz = gemmi.read_mtz_file(str(src_path))
    labels = [col.label for col in mtz.columns]
    if any(label in PANDDA_FREER_LABELS for label in labels):
        return src_path
    for col in mtz.columns:
        if col.label in COMMON_FREER_LABELS:
            col.label = "FreeR_flag"
            out_path = Path(staging_dir) / f"{src_path.stem}_pandda.mtz"
            mtz.write_to_file(str(out_path))
            return out_path
    return src_path


def _latest_finished_job(project, task_names):
    return (
        models.Job.objects.filter(
            project=project,
            task_name__in=list(task_names),
            status=models.Job.Status.FINISHED,
        )
        .order_by("-id")
        .first()
    )


def _find_dictionary_cif(acedrg_job) -> Optional[Path]:
    if not acedrg_job:
        return None
    acedrg_dir = Path(acedrg_job.directory)
    for name in DICT_CIF_CANDIDATES:
        candidate = acedrg_dir / name
        if candidate.exists():
            return candidate
    return None


def collect_pandda_datasets(group):
    """Yield ``(project, dimple_job, acedrg_job)`` for each member project.

    Used by both the export endpoint (to build the ZIP) and the
    ``pandda_data`` summary endpoint (to report counts to the UI), so
    the two stay consistent.
    """
    member_memberships = group.memberships.filter(
        type=models.ProjectGroupMembership.MembershipType.MEMBER
    ).select_related("project")
    for membership in member_memberships:
        project = membership.project
        dimple_job = _latest_finished_job(project, DIMPLE_TASK_NAMES)
        acedrg_job = _latest_finished_job(project, ACEDRG_TASK_NAMES)
        yield project, dimple_job, acedrg_job


def build_pandda_zip(group, temp_dir: Optional[Path] = None) -> PanddaExportResult:
    """Build a PANDDA-ready ZIP for ``group`` and return its location.

    Creates ``temp_dir`` if not supplied (using ``tempfile.mkdtemp``).
    The caller owns cleanup of ``temp_dir`` once the ZIP has been
    served / persisted. Projects without dimple outputs are skipped
    with a warning; the returned ``included_count`` may be zero, in
    which case the caller is expected to surface that to the user.
    """
    if temp_dir is None:
        temp_dir = Path(tempfile.mkdtemp(prefix="pandda_export_"))
    else:
        temp_dir = Path(temp_dir)
        temp_dir.mkdir(parents=True, exist_ok=True)

    zip_path = temp_dir / f"pandda_{group.name}.zip"
    included_count = 0
    # Maps the synthetic ``xtal-NNNN`` directory name (which is what
    # PANDDA2 reads as the crystal number) back to the original
    # project name, written out as ``Projects.csv`` at the ZIP root.
    # The same CSV is consumed by MoorhenPanddaInspect to surface the
    # original provenance as a subtitle alongside each xtal id.
    manifest_rows = []

    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for project, dimple_job, acedrg_job in collect_pandda_datasets(group):
            if not dimple_job:
                continue

            dimple_dir = Path(dimple_job.directory)
            final_pdb = dimple_dir / "final.pdb"
            final_mtz = dimple_dir / "final.mtz"

            if not final_pdb.exists() or not final_mtz.exists():
                logger.warning(
                    "Dimple outputs not found for project %s", project.name
                )
                continue

            xtal_name = f"xtal-{included_count:04d}"
            dataset_path = f"datasets/{xtal_name}"
            try:
                mtz_to_write = prepare_mtz_for_pandda(final_mtz, temp_dir)
            except Exception as relabel_err:
                logger.warning(
                    "FreeR relabel failed for %s, using original MTZ: %s",
                    project.name, relabel_err,
                )
                mtz_to_write = final_mtz

            zf.write(str(final_pdb), f"{dataset_path}/final.pdb")
            zf.write(str(mtz_to_write), f"{dataset_path}/final.mtz")
            manifest_rows.append((xtal_name, project.name))
            included_count += 1

            dict_cif = _find_dictionary_cif(acedrg_job)
            if dict_cif is not None:
                zf.write(str(dict_cif), f"{dataset_path}/dict.cif")

        if manifest_rows:
            csv_lines = ["Dataset, Project"]
            csv_lines.extend(f"{xtal}, {name}" for xtal, name in manifest_rows)
            zf.writestr("Projects.csv", "\n".join(csv_lines) + "\n")

    return PanddaExportResult(
        zip_path=zip_path,
        temp_dir=temp_dir,
        included_count=included_count,
    )
