"""
CCP4i2RunnerBase - Refactored version using clean component architecture.

External API remains unchanged for backward compatibility.
Internal implementation delegates to KeywordExtractor, ArgumentBuilder, and PluginPopulator.
"""

import pathlib
import unittest
import argparse
import re
import os
import sys
import logging
import shlex
import numpy
import traceback

from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.core.CCP4ErrorHandling import CException
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4Container
from ccp4i2.core.CCP4TaskManager import TASKMANAGER
from ccp4i2.core import CCP4Data
from ccp4i2.core import CCP4XtalData

from ccp4i2.lib.utils.parameters.set_parameter import set_parameter_container
from ccp4i2.lib.utils.containers.find_objects import find_object_by_path
from ccp4i2.lib.utils.formats.gemmi_split_mtz import gemmi_split_mtz

# Import refactored components
from .i2run_components import KeywordExtractor, ArgumentBuilder, PluginPopulator

logger = logging.getLogger("root")


# ============================================================================
# Legacy helper functions - kept for backward compatibility
# These delegate to the new components but maintain the old function signatures
# ============================================================================

def get_leaf_paths(container):
    """
    Legacy function - delegates to KeywordExtractor.

    Given a container, return a list of unique leaf element dictionaries.
    """
    # Use new component but match old return format
    keywords = KeywordExtractor._get_leaf_paths(container)

    # Old format included className, new format doesn't - add it back
    for kw in keywords:
        kw["className"] = type(kw["object"]).__name__
        kw["minimumPath"] = kw["path"]  # Will be computed later

    return keywords


def compute_minimum_paths(keywords):
    """
    Legacy function - delegates to KeywordExtractor.

    Compute minimum unique paths for keywords.
    """
    return KeywordExtractor._compute_minimum_paths(keywords)


# ============================================================================
# CCP4i2RunnerBase - Refactored to use component architecture
# ============================================================================

class CCP4i2RunnerBase(object):
    """
    Base class for running CCP4i2 tasks from command line.

    Refactored to use clean component architecture:
    - KeywordExtractor: Extract parameters from plugin
    - ArgumentBuilder: Build argparse arguments
    - PluginPopulator: Populate plugin with parsed args

    External API remains unchanged for backward compatibility.
    """

    def __init__(self, the_args=None, command_line=None, parser=None, parent=None):
        """
        Initialize runner with command-line arguments.

        Args:
            the_args: List of command-line arguments
            command_line: String command line (alternative to the_args)
            parser: argparse.ArgumentParser instance (optional)
            parent: Parent object (legacy, unused)
        """
        # parent parameter is legacy/unused - runner classes don't inherit from HierarchicalObject

        assert (
            the_args is not None or command_line is not None
        ), "Need to provide one of args or command_line"

        if the_args is None and command_line is not None:
            the_args = [os.path.expandvars(a) for a in shlex.split(command_line)]

        self.args = the_args

        if parser is None:
            self.parser = argparse.ArgumentParser(usage="i2run")
        else:
            self.parser = parser

        self.task_name = self.args[0]
        self.parser.add_argument("task_name")
        self.parser.add_argument("--project_name", default=None)
        self.parser.add_argument("--project_path", default=None)
        self.parser.add_argument("--delay", action="store_true")
        self.parser.add_argument("--batch", action="store_true")

        self.parsed_args = None
        self.job_id = None
        self.project_id = None
        self.list_map = {}
        self._plugin = None
        self._keywords = None  # Cache keywords

    def keywordsOfTask(self):
        """Get keywords for the current task."""
        return self.keywordsOfTaskName(self.task_name)

    def parseArgs(self, arguments_parsed=False):
        """
        Parse command-line arguments for the task.

        Uses ArgumentBuilder to add task-specific arguments.

        Args:
            arguments_parsed: If True, skip adding task arguments

        Returns:
            Parsed arguments namespace
        """
        if not arguments_parsed:
            CCP4i2RunnerBase.addTaskArguments(
                self.parser, self.task_name, parent=None
            )

        if self.parsed_args is None:
            self.parsed_args = self.parser.parse_args(self.args)

        return self.parsed_args

    def getPlugin(self, jobId=None, arguments_parsed=False):
        """
        Get plugin instance populated with command-line arguments.

        Caches the plugin instance to prevent multiple instantiations.

        Args:
            jobId: Optional job ID
            arguments_parsed: If True, skip parsing arguments

        Returns:
            Configured plugin instance (cached)
        """
        # Return cached plugin if already created
        if self._plugin is not None:
            return self._plugin

        parsed_args = self.parseArgs(arguments_parsed=arguments_parsed)

        # Extract project info if provided
        if parsed_args.project_name is not None:
            self.projectId = self.projectWithName(
                parsed_args.project_name, parsed_args.project_path
            )
            self.jobId = self.projectJobWithTask(
                projectId=self.projectId, task_name=self.task_name
            )
            jobId = self.jobId

        # Create and cache plugin instance
        self._plugin = self.pluginWithArgs(parsed_args=parsed_args, jobId=jobId)
        return self._plugin

    def projectWithName(self, project_name):
        """Abstract method - implemented in subclasses."""
        raise Exception(
            "CCP4i2RunnerBase does not provide projectWithName - implement in subclasses"
        )

    def projectJobWithTask(self, projectId, task_name=None):
        """Abstract method - implemented in subclasses."""
        raise Exception(
            "CCP4i2RunnerBase does not provide projectJobWithTask - implement in subclasses"
        )

    def execute(self):
        """Abstract method - implemented in subclasses."""
        raise Exception(
            "CCP4i2RunnerBase does not provide execute - implement in subclasses"
        )

    # ========================================================================
    # Static methods for keyword extraction and argument building
    # ========================================================================

    @staticmethod
    def keywordsOfContainer(container: CContainer, growingList=None):
        """Legacy method - delegates to KeywordExtractor."""
        return get_leaf_paths(container)

    @staticmethod
    def getCandidatePath(currentPath):
        """Legacy method - kept for backward compatibility."""
        pathElements = currentPath.split(".")
        return ".".join(pathElements[1:])

    @staticmethod
    def minimisePaths(allKeywords):
        """Legacy method - delegates to compute_minimum_paths."""
        compute_minimum_paths(allKeywords)
        return allKeywords

    @staticmethod
    def keywordsOfTaskName(task_name, parent=None):
        """
        Get all keywords for a task by name.

        Uses KeywordExtractor to extract parameters from plugin definition.

        Args:
            task_name: Name of the task/plugin
            parent: Unused (kept for backward compatibility)

        Returns:
            List of keyword dictionaries with metadata
        """
        # Use new component
        keywords = KeywordExtractor.extract_from_task_name(task_name)

        # Compute minimum paths
        keywords = CCP4i2RunnerBase.minimisePaths(keywords)

        return keywords

    @staticmethod
    def addTaskArguments(theParser, task_name, parent=None):
        """
        Add command-line arguments to parser.

        Uses ArgumentBuilder to add arguments with backward-compatible aliases.

        Args:
            theParser: argparse.ArgumentParser instance
            task_name: Name of the task/plugin
            parent: Unused (kept for backward compatibility)
        """
        keywords = CCP4i2RunnerBase.keywordsOfTaskName(task_name, parent=None)
        logger.debug(str(keywords))

        # Use new component to build arguments
        ArgumentBuilder.build_arguments(theParser, keywords)

        return theParser

    # ========================================================================
    # Utility methods
    # ========================================================================

    @staticmethod
    def slugify_variable(variable):
        """Convert variable name to slug format."""
        # Remove any non-alphanumeric characters except for commas and hyphens
        variable = re.sub(r"[/*\[\] ,]+", "-", variable)
        variable = variable.lower()
        variable = variable.strip("-")
        return variable

    # ========================================================================
    # Plugin population methods
    # ========================================================================

    def handleItem(self, thePlugin: CPluginScript, objectPath, value):
        """
        Set a value on a plugin parameter.

        Delegates to PluginPopulator for actual work, but kept as instance method
        for backward compatibility.

        Args:
            thePlugin: CPluginScript instance
            objectPath: Dot-separated path to parameter
            value: Value(s) to set
        """
        # Delegate to new component
        PluginPopulator._handle_item(thePlugin, objectPath, value)

    def fileForFileUse(
        self,
        thePlugin,
        jobNumber,
        jobParamName,
        paramIndex,
        task_name=None,
        listMap=None,
    ):
        """
        Abstract method for file use resolution.

        Implemented in Django subclass.
        """
        raise NotImplementedError("fileForFileUse must be implemented in subclasses")

    def handleItemOrList(self, thePlugin, objectPath, value):
        """
        Handle setting a value that may be a list.

        Args:
            thePlugin: CPluginScript instance
            objectPath: Path to parameter
            value: Single value or list of values
        """
        # Navigate to the target object
        target = find_object_by_path(thePlugin, objectPath)

        if target is None:
            logger.warning(f"Could not find object at path: {objectPath}")
            return

        # Delegate to PluginPopulator
        if isinstance(value, list):
            path_parts = objectPath.split(".")
            parent_path = ".".join(path_parts[:-1])
            parent = find_object_by_path(thePlugin, parent_path) if parent_path else thePlugin
            PluginPopulator._handle_item_or_list(thePlugin, parent, target, value)
        else:
            PluginPopulator._handle_single_value(target, value)

    def pluginWithArgs(self, parsed_args, workDirectory=None, jobId=None):
        """
        Create plugin instance and populate with command-line arguments.

        Abstract method - implemented in subclasses.

        Args:
            parsed_args: Parsed command-line arguments
            workDirectory: Working directory for the plugin
            jobId: Optional job ID for database integration

        Returns:
            Configured plugin instance
        """
        raise NotImplementedError("pluginWithArgs must be implemented in subclasses")

    @staticmethod
    def parseFileUse(fileUse):
        """
        Parse fileUse string into components.

        Supports formats like:
        - task_name[jobIndex].jobParamName[paramIndex]
        - jobIndex.jobParamName[paramIndex]
        - etc.

        Args:
            fileUse: FileUse string to parse

        Returns:
            Dict with task_name, jobIndex, jobParamName, paramIndex
        """
        fileUseParser_tN_jI_jP_pI = re.compile(
            r"^(?P<task_name>.+)\[(?P<jobIndex>.+)\]\.(?P<jobParamName>.+)\[(?P<paramIndex>.+)\].*$"
        )
        fileUseParser_jI_jP_pI = re.compile(
            r"^(?P<jobIndex>.+)\.(?P<jobParamName>.*)\[(?P<paramIndex>.+)\].*$"
        )
        fileUseParser_tN_jI_jP = re.compile(
            r"^(?P<task_name>.+)\[(?P<jobIndex>.+)\]\.(?P<jobParamName>.+).*$"
        )
        fileUseParser_jI_jP = re.compile(r"^(?P<jobIndex>.+)\.(?P<jobParamName>.+).*$")

        result = {
            "task_name": None,
            "jobIndex": None,
            "jobParamName": None,
            "paramIndex": -1,
        }

        logger.debug(f"Trying to match [{fileUse}]")

        try:
            matches = fileUseParser_tN_jI_jP_pI.match(fileUse)
            result.update(matches.groupdict())
            result["jobIndex"] = int(result["jobIndex"])
            result["paramIndex"] = int(result["paramIndex"])
            return result
        except AttributeError:
            try:
                matches = fileUseParser_jI_jP_pI.match(fileUse)
                result.update(matches.groupdict())
                result["jobIndex"] = int(result["jobIndex"])
                result["paramIndex"] = int(result["paramIndex"])
                return result
            except AttributeError:
                try:
                    matches = fileUseParser_tN_jI_jP.match(fileUse)
                    result.update(matches.groupdict())
                    result["jobIndex"] = int(result["jobIndex"])
                    result["paramIndex"] = int(result["paramIndex"])
                    return result
                except AttributeError:
                    matches = fileUseParser_jI_jP.match(fileUse)
                    result.update(matches.groupdict())
                    result["jobIndex"] = int(result["jobIndex"])
                    result["paramIndex"] = int(result["paramIndex"])
                    return result


# ============================================================================
# Unit tests - kept for backward compatibility
# ============================================================================

class TestParse(unittest.TestCase):
    """Legacy unit tests for fileUse parsing."""

    modelPdb = os.path.join(
        os.path.expandvars("$CCP4I2_TOP/demo_data/gamma/gamma_model.pdb")
    )
    gammaMtz = os.path.join(
        os.path.expandvars("$CCP4I2_TOP/demo_data/gamma/merged_intensities_Xe.mtz")
    )
    betaPdb = os.path.join(
        os.path.expandvars("$CCP4I2_TOP/demo_data/beta_blip/beta.pdb")
    )
    blipPdb = os.path.join(
        os.path.expandvars("$CCP4I2_TOP/demo_data/beta_blip/blip.pdb")
    )

    def setUp(self):
        pass

    def test1(self):
        """Test basic fileUse parsing."""
        parsed = CCP4i2RunnerBase.parseFileUse("task_name[1].XYZIN[0]")
        self.assertEqual(parsed["task_name"], "task_name")
        self.assertEqual(parsed["jobIndex"], 1)
        self.assertEqual(parsed["jobParamName"], "XYZIN")
        self.assertEqual(parsed["paramIndex"], 0)

    def test2(self):
        """Test fileUse without task name."""
        parsed = CCP4i2RunnerBase.parseFileUse("1.XYZIN[0]")
        self.assertEqual(parsed["jobIndex"], 1)
        self.assertEqual(parsed["jobParamName"], "XYZIN")
        self.assertEqual(parsed["paramIndex"], 0)

    def test3(self):
        """Test complex fileUse patterns."""
        parsed = CCP4i2RunnerBase.parseFileUse("-1.HKLOUT[0]")
        self.assertEqual(parsed["jobIndex"], -1)
        self.assertEqual(parsed["jobParamName"], "HKLOUT")
        self.assertEqual(parsed["paramIndex"], 0)

    def test4(self):
        """Test another fileUse pattern."""
        parsed = CCP4i2RunnerBase.parseFileUse("prosmart_refmac[-1].XYZOUT[0]")
        self.assertEqual(parsed["task_name"], "prosmart_refmac")
        self.assertEqual(parsed["jobIndex"], -1)
        self.assertEqual(parsed["jobParamName"], "XYZOUT")
        self.assertEqual(parsed["paramIndex"], 0)

    def test5(self):
        """Test file existence patterns."""
        # This test is incomplete in the original
        pass

    def test6(self):
        """Test MTZ splitting."""
        # This test is incomplete in the original
        pass
