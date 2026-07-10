"""
Job export utilities.

Provides functions to export jobs as ZIP archives, and to serve a job's
exportable data files (MTZ) via the per-plugin export contract or a
generic fallback. See docs/EXPORT_MTZ_PLAN.md for the design.
"""

import logging
import os
from pathlib import Path
from django.http import FileResponse
from ccp4i2.db import models
from ccp4i2.lib.response import Result, api_error
from ccp4i2.db.export_project import export_project_to_zip

logger = logging.getLogger(f"ccp4i2:{__name__}")

# FileType.name prefix shared by every MTZ variant (full, mini, observed,
# freerflag, phases, map). Used by the generic export fallback.
MTZ_TYPE_PREFIX = "application/CCP4-mtz"

# Tasks whose exportJobFile reconstructs via the (now-absent) runCad. Until a
# gemmi join utility lands (EXPORT_MTZ_PLAN Phase 2) these would crash, so we
# skip their plugin export and let the generic fallback handle them instead.
_RECONSTRUCTOR_TASKS = {"aimless_pipe", "import_merged"}


def export_job(job: models.Job, output_path: Path) -> Result[Path]:
    """
    Export a single job as a ZIP archive.

    Creates a downloadable ZIP archive containing the job data,
    associated files, and dependency information.

    Args:
        job: Job model instance to export
        output_path: Path where ZIP file should be created

    Returns:
        Result containing path to created ZIP file

    Example:
        >>> result = export_job(job, Path("/tmp/job_export.zip"))
        >>> if result.success:
        ...     print(f"Exported to: {result.data}")
    """
    try:
        project = job.project
        job_selection = {str(job.number)}

        # Use existing export functionality
        result_path = export_project_to_zip(
            project=project,
            output_path=output_path,
            job_selection=job_selection
        )

        logger.info(
            "Exported job %s (number: %s) from project %s to %s",
            job.uuid, job.number, project.name, result_path
        )

        return Result.ok(result_path)

    except Exception as err:
        logger.exception("Failed to export job %s", job.uuid, exc_info=err)
        return Result.fail(
            f"Failed to export job: {str(err)}",
            details={
                "job_id": str(job.uuid),
                "job_number": job.number,
                "output_path": str(output_path),
                "error_type": type(err).__name__
            }
        )


# ---------------------------------------------------------------------------
# Data-file (MTZ) export — per-plugin contract + generic fallback
# ---------------------------------------------------------------------------

# Prefix for generic-fallback modes, to distinguish them from plugin modes.
_FALLBACK_PREFIX = "file:"


def _output_mtz_files(job: models.Job):
    """Output MTZ File rows produced by this job, whose file exists on disk."""
    files = (
        models.File.objects.filter(
            job=job, type__name__startswith=MTZ_TYPE_PREFIX
        )
        .select_related("type", "job")
        .order_by("job_param_name", "name")
    )
    result = []
    for f in files:
        try:
            if f.path and Path(f.path).exists():
                result.append(f)
        except Exception:
            continue
    return result


def _generic_menu(job: models.Job):
    """Fallback export menu: one entry per tracked output MTZ file on disk."""
    menu = []
    for f in _output_mtz_files(job):
        label = f.annotation or f.job_param_name or f.name
        menu.append([f"{_FALLBACK_PREFIX}{f.uuid}", f"MTZ file: {label}", "application/CCP4-mtz"])
    return menu


def export_job_file_menu(job: models.Job):
    """Return the list of exportable-file menu items for a job.

    Prefers the plugin's module-level ``exportJobFileMenu(jobId)`` contract;
    falls back to the job's tracked output MTZ file(s). Each item is
    ``[mode, label, mimeType]``. Empty list => nothing to export
    (frontend hides the button).
    """
    from ccp4i2.core.tasks import get_plugin_module

    task_name = job.task_name
    module = get_plugin_module(task_name)
    if (
        module is not None
        and task_name not in _RECONSTRUCTOR_TASKS
        and hasattr(module, "exportJobFileMenu")
        and hasattr(module, "exportJobFile")
    ):
        try:
            declared = module.exportJobFileMenu(jobId=str(job.uuid)) or []
        except Exception as err:
            logger.warning(
                "Plugin exportJobFileMenu failed for task %s (job %s): %s",
                task_name, job.uuid, err,
            )
            declared = []
        # The plugin menu is optimistic — it advertises a mode regardless of
        # whether the file exists. For these (non-reconstructor) tasks
        # exportJobFile is a cheap path resolve, so validate each mode by
        # confirming it yields a file on disk before offering it.
        menu = []
        for item in declared:
            mode = item[0] if item else None
            if mode and _plugin_export_path(module, job, mode):
                menu.append(item)
        if menu:
            return menu
    # Generic fallback (also covers reconstructor tasks and any task whose
    # plugin export produced nothing).
    return _generic_menu(job)


def _plugin_export_path(module, job: models.Job, mode):
    """Resolve a plugin export mode to an existing path, or None."""
    try:
        path = module.exportJobFile(jobId=str(job.uuid), mode=mode)
    except Exception as err:
        logger.warning(
            "Plugin exportJobFile(mode=%s) failed for task %s (job %s): %s",
            mode, job.task_name, job.uuid, err,
        )
        return None
    if path and os.path.exists(str(path)):
        return Path(path)
    return None


def export_job_file(pk, mode):
    """Resolve and stream a job's exportable file for the given mode.

    Returns ``(file_response, error_response)`` — exactly one is non-None,
    matching the JobViewSet.export_job_file action contract.
    """
    if not mode:
        return None, api_error("Missing required 'mode' parameter", status=400)

    try:
        job = models.Job.objects.get(id=pk)
    except models.Job.DoesNotExist:
        return None, api_error(f"Job {pk} not found", status=404)

    export_path = None

    # Generic-fallback mode: serve a specific tracked output File as-is.
    if str(mode).startswith(_FALLBACK_PREFIX):
        file_uuid = str(mode)[len(_FALLBACK_PREFIX):]
        try:
            f = models.File.objects.get(uuid=file_uuid, job=job)
        except models.File.DoesNotExist:
            return None, api_error(f"Export file {file_uuid} not found for job", status=404)
        export_path = f.path
    else:
        # Plugin mode: call the module-level exportJobFile contract.
        if job.task_name in _RECONSTRUCTOR_TASKS:
            return None, api_error(
                "MTZ reconstruction for this task is not yet available "
                "(pending a gemmi CAD replacement).",
                status=501,
            )
        from ccp4i2.core.tasks import get_plugin_module

        module = get_plugin_module(job.task_name)
        if module is None or not hasattr(module, "exportJobFile"):
            return None, api_error(
                f"Task {job.task_name} does not support file export", status=400
            )
        export_path = _plugin_export_path(module, job, mode)

    if not export_path or not os.path.exists(str(export_path)):
        return None, api_error("No exportable file was produced", status=404)

    export_path = Path(export_path)
    download_name = f"{job.project.name}_job_{job.number}_{export_path.name}"
    response = FileResponse(
        open(export_path, "rb"),
        as_attachment=True,
        filename=download_name,
        content_type="application/octet-stream",
    )
    return response, None
