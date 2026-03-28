# Copyright (C) 2025-2026 University of York
# Copyright (C) 2025-2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
import logging
import os

from rest_framework.decorators import api_view
from django.http import JsonResponse
from django.db import connection
from ..core.tasks import get_task_tree
from ..db import models
import psutil

logger = logging.getLogger(__name__)


@api_view(["GET"])
def task_tree(request):
    """
    Returns the task tree structure for displaying available plugins.

    Response format:
    {
        "success": true,
        "data": {
            "task_tree": {
                "tree": [[module_name, title, [task_names...]], ...],
                "lookup": {taskName: {metadata...}, ...},
                "iconLookup": {module_name: icon_path, ...}
            }
        }
    }
    """
    task_tree_data = get_task_tree()
    return JsonResponse({"success": True, "data": {"task_tree": task_tree_data}})


@api_view(["GET"])
def active_jobs(request):
    """
    Returns a list of all running/queued jobs in the ccp4i2 job queue.

    Includes jobs that are RUNNING locally, RUNNING_REMOTELY on Azure workers,
    or QUEUED waiting for a worker to pick them up.
    """
    active_statuses = [
        models.Job.Status.RUNNING,
        models.Job.Status.RUNNING_REMOTELY,
        models.Job.Status.QUEUED,
    ]
    running_jobs = models.Job.objects.filter(status__in=active_statuses)
    active_jobs_list = []
    for job in running_jobs:
        pid = job.process_id
        cpu_percent = None
        mem_usage = None

        # Try to get process metrics if we have a local PID
        if pid:
            try:
                proc = psutil.Process(pid)

                def get_total_cpu_percent(process):
                    try:
                        process.cpu_percent(interval=None)
                        children = process.children(recursive=True)
                        for child in children:
                            child.cpu_percent(interval=None)
                        total_cpu = process.cpu_percent(interval=0.5)
                        for child in children:
                            total_cpu += get_total_cpu_percent(child)
                        return total_cpu
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        return 0.0

                cpu_percent = get_total_cpu_percent(proc)
                mem_info = proc.memory_info()
                mem_usage = mem_info.rss
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass

        active_jobs_list.append(
            {
                "project": job.project.name,
                "job_id": job.pk,
                "job_task_name": job.task_name,
                "job_uuid": job.uuid,
                "job_number": job.number,
                "status": job.get_status_display(),
                "pid": pid,
                "cpu_percent": cpu_percent,
                "memory_usage_bytes": mem_usage,
            }
        )
    return JsonResponse({"success": True, "data": {"active_jobs": active_jobs_list}})


def health_check(request):
    """
    Simple health check endpoint for deployment monitoring.
    Returns 200 OK if the service is healthy.
    """
    try:
        # Test database connection
        with connection.cursor() as cursor:
            cursor.execute("SELECT 1")

        return JsonResponse(
            {
                "status": "healthy",
                "service": "ccp4i2-django-api",
                "database": "connected",
            }
        )
    except Exception as e:
        return JsonResponse(
            {
                "status": "unhealthy",
                "service": "ccp4i2-django-api",
                "database": "disconnected",
                "error": str(e),
            },
            status=503,
        )


def version_info(request):
    """
    Returns version information for the server deployment.
    Build timestamp and git commit are set via environment variables during Docker build.
    """
    return JsonResponse(
        {
            "buildTimestamp": os.environ.get("BUILD_TIMESTAMP", "dev"),
            "gitCommit": os.environ.get("GIT_COMMIT", "unknown"),
        }
    )


@api_view(["GET"])
def monomer_info(request, code):
    """
    Return atom and bond information for a CCP4 monomer library entry.

    GET /api/ccp4i2/monomer-info/<code>/

    Response:
    {
        "success": true,
        "data": {
            "code": "LYS",
            "atoms": ["N", "CA", "C", "O", "CB", ...],
            "bonds": [
                {"atom1": "N", "atom2": "CA", "type": "single"},
                ...
            ]
        }
    }
    """
    from ..lib.utils.formats.cif_ligand import extract_monomer_atoms_bonds

    ccp4 = os.environ.get("CCP4", "")
    if not ccp4:
        return JsonResponse(
            {"success": False, "error": "CCP4 environment not configured"},
            status=500,
        )

    # Monomer library: $CCP4/lib/data/monomers/{first_letter}/{CODE}.cif
    code_upper = code.upper()
    first_letter = code_upper[0].lower()
    cif_path = os.path.join(ccp4, "lib", "data", "monomers", first_letter, f"{code_upper}.cif")

    if not os.path.isfile(cif_path):
        return JsonResponse(
            {"success": False, "error": f"Monomer '{code_upper}' not found in library"},
            status=404,
        )

    try:
        result = extract_monomer_atoms_bonds(cif_path)
        return JsonResponse({
            "success": True,
            "data": {
                "code": code_upper,
                "atoms": result["atoms"],
                "bonds": result["bonds"],
            },
        })
    except Exception as e:
        logger.exception("Error parsing monomer CIF for %s", code_upper)
        return JsonResponse(
            {"success": False, "error": f"Error parsing monomer: {str(e)}"},
            status=500,
        )
