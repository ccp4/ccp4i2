"""
CCP4ProjectsManager.py - Modern Django-based Projects Manager

This module provides a modern replacement for the legacy CCP4ProjectsManager
that uses Django ORM for database access instead of the old XML database.

The CProjectsManager class provides methods for:
- Managing project and job information
- Generating file paths for job outputs
- Accessing job directories
- Database operations through Django ORM

Backward Compatibility:
- CPurgeProject is imported from ccp4i2.core.CPurgeProject for legacy plugin code
"""

import os
import logging
from typing import Optional, Dict, Any

# Import CPurgeProject for backward compatibility with legacy plugins
from ccp4i2.core.CPurgeProject import CPurgeProject

logger = logging.getLogger(f"ccp4x:{__name__}")


# Ensure Django is configured early to prevent app registry errors in async/subprocess contexts
def _ensure_django_configured():
    """
    Ensure Django is configured before any model imports.

    This is critical for async workers and subprocesses that may not inherit
    the Django configuration from the main process.

    Safe to call multiple times (idempotent).
    """
    import os
    import sys

    # Check if Django is available
    try:
        import django
    except ImportError:
        # Django not available - this is OK for some contexts
        return

    # Check if already configured
    # NOTE: We must check if apps attribute exists before accessing it
    # because it's not available until django.setup() is called
    try:
        if hasattr(django, 'apps') and django.apps.apps.ready:
            return
    except Exception:
        # If checking ready state fails, assume not configured
        pass

    # Ensure DJANGO_SETTINGS_MODULE is set
    if 'DJANGO_SETTINGS_MODULE' not in os.environ:
        os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'ccp4x.settings')

    # Configure Django
    try:
        django.setup()
        logger.debug("Django configured successfully in CCP4ProjectsManager")
    except Exception as e:
        logger.warning(f"Failed to configure Django in CCP4ProjectsManager: {e}")


# Configure Django on module import
_ensure_django_configured()


class CProjectsManager:
    """
    Modern projects manager that uses Django database access.

    This class provides the same interface as the legacy CProjectsManager
    but uses Django ORM instead of XML database files.
    """

    insts = None  # Singleton instance

    def __init__(self, database=None):
        """
        Initialize the projects manager.

        Args:
            database: Optional database instance (for compatibility).
                     In the modern version, this is managed by Django.
        """
        self._db = None
        if database is not None:
            self.setDatabase(database)

    def setDatabase(self, database):
        """
        Set the database instance.

        Args:
            database: Database instance to use
        """
        self._db = database

    def db(self, label=None):
        """
        Get the database API instance.

        Args:
            label: Optional label (unused, for compatibility)

        Returns:
            CCP4i2DjangoDbApi instance for database operations
        """
        if self._db is None:
            # Lazy load the Django database API
            try:
                # Import Django database API
                # NOTE: Django MUST be configured before this is called (via django.setup() or pytest-django)
                # DO NOT call django.setup() here - it must be done in the main thread before any async workers
                # Import without 'server.' prefix since server/ is in sys.path (added by conftest or main)
                from ccp4x.db.ccp4i2_django_dbapi import CCP4i2DjangoDbApi
                self._db = CCP4i2DjangoDbApi()
            except Exception as e:
                logger.error(f"Failed to initialize Django database: {e}")
                raise
        return self._db

    def makeFileName(self, jobId: Optional[str] = None, mode: str = 'PARAMS',
                    jobInfo: Optional[Dict] = None, qualifier: Optional[str] = None) -> str:
        """
        Generate suitable name for job output file.

        This replicates the behavior of the legacy makeFileName() method,
        generating standard filenames for various job output types.

        Args:
            jobId: Job UUID or identifier
            mode: Type of file to generate name for. Options include:
                  'ROOT', 'PARAMS', 'JOB_INPUT', 'PROGRAMXML', 'LOG',
                  'STDOUT', 'STDERR', 'INTERRUPT', 'DIAGNOSTIC', 'REPORT',
                  'DIAGNOSTIC_REPORT', 'TABLE_RTF', 'TABLES_DIR',
                  'XML_TABLES_DIR', 'COM', 'MGPICDEF', 'PIC', 'RVAPIXML'
            jobInfo: Optional job info dict (unused, for compatibility)
            qualifier: Optional qualifier to add to filename (e.g., for numbered files)

        Returns:
            Full path to the generated filename
        """
        # Get job directory
        myDir = self.jobDirectory(jobId=jobId)

        # Standard filename definitions
        defNames = {
            'ROOT': '',
            'PARAMS': 'params.xml',
            'JOB_INPUT': 'input_params.xml',
            'PROGRAMXML': 'program.xml',
            'LOG': 'log.txt',
            'STDOUT': 'stdout.txt',
            'STDERR': 'stderr.txt',
            'INTERRUPT': 'interrupt_status.xml',
            'DIAGNOSTIC': 'diagnostic.xml',
            'REPORT': 'report.html',
            'DIAGNOSTIC_REPORT': 'diagnostic_report.html',
            'TABLE_RTF': 'tables.rtf',
            'TABLES_DIR': 'tables_as_csv_files',
            'XML_TABLES_DIR': 'tables_as_xml_files',
            'COM': 'com.txt',
            'MGPICDEF': 'report.mgpic.py',
            'PIC': 'report.png',
            'RVAPIXML': 'i2.xml'
        }

        fileName = defNames.get(mode, 'unknown.unk')

        # Add qualifier if specified (e.g., for numbered files)
        if qualifier is not None:
            if '.' in fileName:
                base, ext = fileName.rsplit('.', 1)
                fileName = f"{base}_{qualifier}.{ext}"
            else:
                fileName = f"{fileName}_{qualifier}"

        # Special case: ROOT returns just the directory
        if mode == 'ROOT':
            return myDir

        return os.path.join(myDir, fileName)

    def jobDirectory(self, jobId: Optional[str] = None, jobNumber: Optional[str] = None,
                    projectDirectory: Optional[str] = None, projectId: Optional[str] = None,
                    create: bool = True, subDir: Optional[str] = None) -> Optional[str]:
        """
        Get the directory path for a job.

        Args:
            jobId: Job UUID or identifier
            jobNumber: Job number (e.g., "1.2.3")
            projectDirectory: Project directory path
            projectId: Project UUID
            create: Whether to create the directory if it doesn't exist
            subDir: Optional subdirectory to append

        Returns:
            Full path to the job directory, or None if directory doesn't exist and create=False

        Raises:
            Exception: If directory cannot be created or accessed
        """
        try:
            # Get directory from database
            # We pass jobId and jobNumber; the implementation will look up project info
            directory = self.db().jobDirectory(
                projectId=projectId,
                projectDirectory=projectDirectory,
                jobId=jobId,
                jobNumber=jobNumber,
                projectName=None,  # Will be looked up from jobId if needed
                create=create
            )
        except Exception as e:
            logger.error(
                f"Error getting job directory: jobId={jobId}, jobNumber={jobNumber}, "
                f"projectDirectory={projectDirectory}, projectId={projectId}"
            )
            raise

        # Create directory if it doesn't exist and create=True
        if not os.path.exists(directory):
            if create:
                try:
                    os.makedirs(directory, exist_ok=True)
                    logger.debug(f"Created job directory: {directory}")
                except Exception as e:
                    logger.error(f"Failed to create directory {directory}: {e}")
                    raise
            else:
                logger.warning(f"Job directory does not exist: {directory}")
                return None

        # Handle subdirectory if specified
        if subDir is not None:
            directory = os.path.join(directory, subDir)
            if not os.path.exists(directory):
                try:
                    os.makedirs(directory, exist_ok=True)
                    logger.debug(f"Created subdirectory: {directory}")
                except Exception as e:
                    logger.error(f"Failed to create subdirectory {directory}: {e}")
                    raise

        return directory


def PROJECTSMANAGER() -> CProjectsManager:
    """
    Get the singleton instance of CProjectsManager.

    This function provides the same interface as the legacy PROJECTSMANAGER()
    function, but returns a modern Django-based manager.

    Returns:
        Singleton CProjectsManager instance
    """
    if CProjectsManager.insts is None:
        CProjectsManager.insts = CProjectsManager()
    return CProjectsManager.insts
