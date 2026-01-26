from rest_framework.decorators import api_view
from django.http import JsonResponse
from django.db import connection
from ..db import models
from ..lib.utils.navigation.task_tree import get_task_tree
import psutil


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
    Returns a list of all running jobs in the ccp4i2 job queue.
    """
    running_jobs = models.Job.objects.filter(status=models.Job.Status.RUNNING)
    active_jobs_list = []
    for job in running_jobs:
        pid = job.process_id
        try:
            proc = psutil.Process(pid)

            def get_total_cpu_percent(process):
                try:
                    # Initialize measurement
                    process.cpu_percent(interval=None)
                    children = process.children(recursive=True)
                    for child in children:
                        child.cpu_percent(interval=None)
                    # Wait for interval and measure again
                    total_cpu = process.cpu_percent(interval=0.5)
                    for child in children:
                        total_cpu += get_total_cpu_percent(child)
                    return total_cpu
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    return 0.0

            cpu_percent = get_total_cpu_percent(proc)

            mem_info = proc.memory_info()
            mem_usage = mem_info.rss  # Resident Set Size in bytes
            active_jobs_list.append(
                {
                    "project": job.project.name,
                    "job_id": job.pk,
                    "job_task_name": job.task_name,
                    "job_uuid": job.uuid,
                    "job_number": job.number,
                    "pid": pid,
                    "cpu_percent": cpu_percent,
                    "memory_usage_bytes": mem_usage,
                }
            )
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            active_jobs_list.append(
                {"pid": pid, "error": "Process not found or access denied"}
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
