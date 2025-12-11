"""
CPurgeProject.py - Job Directory Cleanup Utility

This module provides functionality to clean up non-essential files from job
execution directories. It implements the CPurgeProject class which allows
selective deletion of intermediate, diagnostic, and scratch files based on
configurable purge categories and contexts.

The purge system is popular with users who want to keep their project
directories tidy by removing large intermediate files while preserving
essential outputs and reports.

Usage:
    from ccp4i2.core.CPurgeProject import CPurgeProject

    # Create purge instance for a project
    cleanup = CPurgeProject(projectId)

    # Purge intermediate files from a job
    cleanup.purgeJob(jobId, context="extended_intermediate", reportMode="skip")

    # Purge with custom categories
    cleanup.purgeJob(jobId, purgeCategories=[1, 2, 5])
"""

import os
import glob
import logging
from typing import Optional, List, Tuple, Dict
from pathlib import Path

logger = logging.getLogger(f"ccp4i2:{__name__}")


class CPurgeProject:
    """
    Manages cleanup of non-essential files from job directories.

    Files are assigned to purge categories based on their purpose:
    - Category 1: Scratch files (temporary, always safe to delete)
    - Category 2: Scratch with potential diagnostic value
    - Category 3: Diagnostic files (logs, debug output)
    - Category 4: Files redundant after report is created
    - Category 5: Intermediate data files
    - Category 6: Redundant on project completion
    - Category 7: Large files that could be deleted
    - Category 0: Retained files (overrides default deletion)

    Different contexts define which categories to delete:
    - 'script_finish': [1, 2] - Cleanup at end of script execution
    - 'onfinish': [1, 2, 4] - Cleanup when job finishes normally
    - 'temporary': [1, 2, 4] - Remove temporary files
    - 'intermediate': [1, 2, 4, 5] - Remove intermediate data
    - 'extended_intermediate': [1, 2, 5, 7] - Remove large intermediate files
    - 'project_complete': [1, 2, 3, 4, 5, 6] - Full cleanup on project completion
    """

    # Purge category definitions
    PURGECODES = {
        1: 'Scratch file',
        2: 'Scratch file potentially useful diagnostic',
        3: 'Diagnostic file',
        4: 'File redundant after report created',
        5: 'Intermediate data file',
        6: 'Redundant on project completion',
        7: 'Large file that could be deleted',
        0: 'Retained file - used to override default'
    }

    # Context definitions - which categories to purge in each context
    CONTEXTLOOKUP = {
        'script_finish': [1, 2],
        'onfinish': [1, 2, 4],
        'temporary': [1, 2, 4],
        'intermediate': [1, 2, 4, 5],
        'extended_intermediate': [1, 2, 5, 7],
        'project_complete': [1, 2, 3, 4, 5, 6]
    }

    # Default search list - files common to many tasks
    SEARCHLIST = [
        ['*-coordinates_*.*', 1],
        ['report.previous_*.html', 1],
        ['report_tmp.previous_*.html', 1],
        ['params.previous_*.xml', 1],
        ['diagnostic.previous_*.xml', 1],
        ['hklin.mtz', 1],
        ['hklout.mtz', 1],
        ['*mtzsplit.log', 2],
        ['log_mtz*.txt', 2],
        ['sftools', 2],
        ['ctruncate', 2],
        ['scratch', 2],
        ['stderr.txt', 3],
        ['stdout.txt', 3],
        ['diagnostic.xml', 3],
        ['report_tmp.html', 4],
        ['program.xml', 4],
        ['XMLOUT.xml', 4],
        ['log.txt', 4],
        ['*.scene.xml', 6],
        ['*observed_data_as*', 6],
        ['*phases_as*', 6],
        ['*%*/report.html', 6],
        ['*%*/tables_as_csv_files', 6],
        ['*%*/tables_as_xml_files', 6],
        ['*%*/params.xml', 6],
        ['com.txt', 6]
    ]

    def __init__(self, projectId: Optional[str] = None):
        """
        Initialize purge manager for a project.

        Args:
            projectId: UUID of the project to manage cleanup for
        """
        self.projectId = projectId
        self._db = None

    def db(self):
        """
        Get database API instance (lazy initialization).

        Returns:
            Database API instance for accessing job information
        """
        if self._db is None:
            try:
                # Import Django database API
                # NOTE: Django MUST be configured before this is called (via django.setup() or pytest-django)
                # DO NOT call django.setup() here - it must be done in the main thread before any async workers
                # Import without 'server.' prefix since server/ is in sys.path (added by conftest or main)
                from ccp4i2.db.ccp4i2_django_dbapi import CCP4i2DjangoDbApi
                self._db = CCP4i2DjangoDbApi()
            except Exception as e:
                logger.error(f"Failed to initialize Django database: {e}")
                raise
        return self._db

    def purgeJob(self, jobId: str, context: Optional[str] = None,
                 purgeCategories: Optional[List[int]] = None,
                 reportMode: str = "report") -> Dict[str, int]:
        """
        Purge files from a job directory based on context or categories.

        Args:
            jobId: Job UUID to purge files from
            context: Purge context name (e.g., 'extended_intermediate', 'project_complete').
                    If provided, uses CONTEXTLOOKUP to determine categories.
            purgeCategories: Explicit list of category numbers to purge (overrides context)
            reportMode: How to report deletions:
                       - "report": Log all deletions
                       - "skip": Silent operation
                       - "verbose": Detailed logging

        Returns:
            Dictionary with purge statistics:
            - 'files_deleted': Number of files deleted
            - 'bytes_freed': Total bytes freed
            - 'errors': Number of errors encountered
        """
        # Determine which categories to purge
        if purgeCategories is None:
            if context is None:
                context = 'intermediate'
            purgeCategories = self.CONTEXTLOOKUP.get(context, [1, 2, 4])

        # Get job directory
        try:
            from ccp4i2.core.CCP4ProjectsManager import PROJECTSMANAGER
            job_dir = PROJECTSMANAGER().jobDirectory(jobId=jobId, create=False)
            if job_dir is None or not os.path.exists(job_dir):
                logger.warning(f"Job directory not found for job {jobId}")
                return {'files_deleted': 0, 'bytes_freed': 0, 'errors': 1}
        except Exception as e:
            logger.error(f"Error getting job directory for {jobId}: {e}")
            return {'files_deleted': 0, 'bytes_freed': 0, 'errors': 1}

        # Get task-specific purge list if available
        task_purge_list = self._getTaskPurgeList(jobId)

        # Build combined search list
        search_list = self._buildSearchList(task_purge_list)

        # Track statistics
        stats = {
            'files_deleted': 0,
            'bytes_freed': 0,
            'errors': 0
        }

        # Process each search pattern
        for entry in search_list:
            pattern = entry[0]
            category = entry[1]

            # Skip if category not in purge list
            if category not in purgeCategories:
                continue

            # Handle sub-job patterns (e.g., 'refmac%*/ABCDOUT.mtz')
            if '%' in pattern:
                self._purgeSubJobFiles(job_dir, pattern, category, reportMode, stats)
            else:
                self._purgeFiles(job_dir, pattern, category, reportMode, stats)

        # Log summary
        if reportMode != "skip":
            logger.info(f"Purge complete for job {jobId}: "
                       f"{stats['files_deleted']} files deleted, "
                       f"{stats['bytes_freed']} bytes freed, "
                       f"{stats['errors']} errors")

        return stats

    def _getTaskPurgeList(self, jobId: str) -> List[List]:
        """
        Get task-specific PURGESEARCHLIST from the job's plugin class.

        Args:
            jobId: Job UUID

        Returns:
            Task-specific purge search list, or empty list if not available
        """
        try:
            # Try to get job info from database
            db = self.db()
            job_info = db.getJobInfo(jobId)
            if not job_info:
                return []

            # Get task name and try to import the plugin
            task_name = job_info.get('taskname')
            if not task_name:
                return []

            # Try to get plugin class
            from ccp4i2.core.task_manager.plugin_registry import get_plugin_class
            plugin_class = get_plugin_class(task_name)

            if plugin_class and hasattr(plugin_class, 'PURGESEARCHLIST'):
                return plugin_class.PURGESEARCHLIST
        except Exception as e:
            logger.debug(f"Could not get task purge list for {jobId}: {e}")

        return []

    def _buildSearchList(self, task_purge_list: List[List]) -> List[List]:
        """
        Build combined search list from default and task-specific lists.

        Task-specific entries can override defaults by using category 0.

        Args:
            task_purge_list: Task-specific purge search list

        Returns:
            Combined search list with overrides applied
        """
        # Start with default list
        search_list = list(self.SEARCHLIST)

        # Add task-specific entries
        if task_purge_list:
            # Check for overrides (category 0)
            overridden_patterns = {entry[0] for entry in task_purge_list if entry[1] == 0}

            # Remove overridden patterns from default list
            if overridden_patterns:
                search_list = [entry for entry in search_list
                             if entry[0] not in overridden_patterns]

            # Add non-override task entries
            search_list.extend([entry for entry in task_purge_list if entry[1] != 0])

        return search_list

    def _purgeFiles(self, base_dir: str, pattern: str, category: int,
                   reportMode: str, stats: Dict[str, int]) -> None:
        """
        Purge files matching a pattern in the base directory.

        Args:
            base_dir: Directory to search in
            pattern: Glob pattern for files to delete
            category: Purge category number
            reportMode: Reporting mode
            stats: Statistics dictionary to update
        """
        try:
            # Use glob to find matching files
            full_pattern = os.path.join(base_dir, pattern)
            matches = glob.glob(full_pattern)

            for file_path in matches:
                self._deleteFile(file_path, category, reportMode, stats)
        except Exception as e:
            logger.error(f"Error purging files with pattern {pattern}: {e}")
            stats['errors'] += 1

    def _purgeSubJobFiles(self, base_dir: str, pattern: str, category: int,
                         reportMode: str, stats: Dict[str, int]) -> None:
        """
        Purge files in sub-job directories matching pattern.

        Handles patterns like 'refmac%*/ABCDOUT.mtz' where % indicates
        a sub-job directory.

        Args:
            base_dir: Base job directory
            pattern: Pattern with % wildcards for sub-jobs
            category: Purge category number
            reportMode: Reporting mode
            stats: Statistics dictionary to update
        """
        try:
            # Replace % with * for glob matching
            glob_pattern = pattern.replace('%', '')
            full_pattern = os.path.join(base_dir, glob_pattern)
            matches = glob.glob(full_pattern, recursive=True)

            for file_path in matches:
                self._deleteFile(file_path, category, reportMode, stats)
        except Exception as e:
            logger.error(f"Error purging sub-job files with pattern {pattern}: {e}")
            stats['errors'] += 1

    def _deleteFile(self, file_path: str, category: int, reportMode: str,
                   stats: Dict[str, int]) -> None:
        """
        Delete a single file or directory and update statistics.

        Args:
            file_path: Path to file or directory to delete
            category: Purge category number
            reportMode: Reporting mode
            stats: Statistics dictionary to update
        """
        try:
            if not os.path.exists(file_path):
                return

            # Get size before deletion
            if os.path.isfile(file_path):
                size = os.path.getsize(file_path)
            else:
                # For directories, sum all file sizes
                size = sum(f.stat().st_size for f in Path(file_path).rglob('*') if f.is_file())

            # Delete file or directory
            if os.path.isfile(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                import shutil
                shutil.rmtree(file_path)
            else:
                return  # Skip special files

            # Update statistics
            stats['files_deleted'] += 1
            stats['bytes_freed'] += size

            # Report if requested
            if reportMode == "verbose":
                category_name = self.PURGECODES.get(category, 'Unknown')
                logger.info(f"Deleted ({category_name}): {file_path} ({size} bytes)")
            elif reportMode == "report":
                logger.debug(f"Deleted: {file_path}")

        except Exception as e:
            logger.warning(f"Failed to delete {file_path}: {e}")
            stats['errors'] += 1
