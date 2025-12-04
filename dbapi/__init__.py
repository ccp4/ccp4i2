"""
CCP4i2 Database API - Environment-aware compatibility layer.

This module provides database access that works in both:
1. Django mode: Uses Django ORM (server/ccp4x/db/models.py)
2. Qt mode: Uses legacy QtSql-based CCP4DbApi

Usage:
    from dbapi import CDbApi, JOB_STATUS_PENDING, ...

The appropriate implementation is selected based on the baselayer
environment detection.
"""

# Job status constants (shared between both backends)
JOB_STATUS_UNKNOWN = 0
JOB_STATUS_PENDING = 1
JOB_STATUS_QUEUED = 2
JOB_STATUS_RUNNING = 3
JOB_STATUS_INTERRUPTED = 4
JOB_STATUS_FAILED = 5
JOB_STATUS_FINISHED = 6
JOB_STATUS_REMOTE = 7
JOB_STATUS_FILE_HOLDER = 8
JOB_STATUS_TO_DELETE = 9
JOB_STATUS_UNSATISFACTORY = 10

# File role constants
FILE_ROLE_OUT = 0
FILE_ROLE_IN = 1
FILE_ROLE_IMPORT = 2

# Job evaluation constants
JOB_EVALUATION_NOT_SET = 0
JOB_EVALUATION_POOR = 1
JOB_EVALUATION_OK = 2
JOB_EVALUATION_GOOD = 3
JOB_EVALUATION_EXCELLENT = 4


def _is_django_mode():
    """Check if we're in Django mode."""
    try:
        from baselayer import DJANGO
        return DJANGO()
    except ImportError:
        return False


if _is_django_mode():
    # Django mode: Use Django ORM models
    try:
        from ccp4x.db.models import (
            Project,
            Job,
            File,
            FileUse,
            FileType,
            ProjectGroup,
            ProjectGroupMembership,
            ProjectTag,
            ProjectExport,
            ProjectImport,
            ServerJob,
            JobValueKey,
            JobFloatValue,
            JobCharValue,
            FileExport,
            FileImport,
            XData,
        )

        # Compatibility class for Django mode
        class CDbApi:
            """
            Django-based database API.

            This provides a compatibility interface for code that expects
            the legacy CDbApi class while using Django ORM under the hood.
            """

            def __init__(self, dbFileName=None, projectId=None):
                self.dbFileName = dbFileName
                self.projectId = projectId

            def getProjects(self):
                """Return list of all projects."""
                return list(Project.objects.all().values())

            def getProject(self, projectId):
                """Return project by ID."""
                try:
                    return Project.objects.get(pk=projectId)
                except Project.DoesNotExist:
                    return None

            def getJob(self, jobId):
                """Return job by ID."""
                try:
                    return Job.objects.get(pk=jobId)
                except Job.DoesNotExist:
                    return None

            def getJobsInProject(self, projectId):
                """Return all jobs in a project."""
                return list(Job.objects.filter(project_id=projectId))

            # Add more compatibility methods as needed

    except ImportError as e:
        # Django models not available - provide stub
        import warnings
        warnings.warn(f"Django models not available: {e}. Using stubs.")

        class CDbApi:
            """Stub CDbApi when Django models unavailable."""
            def __init__(self, *args, **kwargs):
                raise RuntimeError("Django database not configured")

else:
    # Qt mode: Use legacy QtSql-based implementation
    from .CCP4DbApi import *
    from .CCP4DbUtils import *
