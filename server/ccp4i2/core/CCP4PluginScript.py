"""
CPluginScript - Base class for CCP4i2 task wrappers and pipelines.

This is a modernized version that uses the new CData infrastructure
while maintaining backward compatibility with the existing API used
by all .def.xml files and plugin implementations.
"""

from typing import Optional
from pathlib import Path
import xml.etree.ElementTree as ET
import logging
import os
import shutil

from ccp4i2.core.base_object.base_classes import CData, CContainer
from ccp4i2.core.base_object.error_reporting import CErrorReport, SEVERITY_ERROR, SEVERITY_WARNING
from ccp4i2.core.task_manager.def_xml_handler import DefXmlParser
from ccp4i2.core.task_manager.params_xml_handler import ParamsXmlHandler
from ccp4i2.core.CCP4TaskManager import locate_def_xml
from ccp4i2.core.task_manager.plugin_registry import get_plugin_class
from ccp4i2.core.base_object.class_metadata import cdata_class

# Module-level logger
logger = logging.getLogger(__name__)


@cdata_class(
    error_codes={
        '200': {'description': 'Error merging MTZ files'},
        '201': {'description': 'Invalid miniMtzsIn specification'},
        '202': {'description': 'File object has no CONTENT_SIGNATURE_LIST'},
        '203': {'description': 'File object has no path set'},
        '204': {'description': 'MTZ file not found'},
        '205': {'description': 'File object not found in inputData or outputData'},
        '206': {'description': 'Invalid contentFlag for file type'},
        '207': {'description': 'Invalid item format in miniMtzsIn'},
        '208': {'description': 'Conversion method not found on file object'},
        '209': {'description': 'Conversion to requested format not yet implemented'},
        '210': {'description': 'Error during MTZ file conversion'},
        '300': {'description': 'Input MTZ file not found for split operation'},
        '301': {'description': 'miniMtzsOut and programColumnNames length mismatch'},
        '302': {'description': 'Output object not found in container.outputData'},
        '303': {'description': 'Output object has no path set'},
        '304': {'description': 'Error splitting MTZ file'},
    }
)
class CPluginScript(CData):
    """
    Base class for CCP4i2 wrappers and pipelines.

    A CPluginScript wraps a crystallographic program or script, providing:
    - Parameter management through containers
    - Input/output file handling
    - Command generation from templates
    - Process execution and monitoring
    - Error reporting and database integration

    Subclasses should define:
        TASKTITLE: Display title for GUI
        TASKNAME: Unique task identifier
        TASKCOMMAND: Executable name
    """

    # Class attributes to be defined in subclasses
    TASKTITLE = None
    TASKNAME = None
    TASKCOMMAND = None
    ASYNCHRONOUS = False  # Set to True for async execution

    # Status codes
    SUCCEEDED = 0
    FAILED = 1
    RUNNING = 2
    UNSATISFACTORY = 3  # Job completed but with warnings/issues

    def __init__(self,
                 parent=None,
                 name: Optional[str] = None,
                 xmlFile: Optional[str] = None,
                 workDirectory: Optional[str | Path] = None,
                 dummy: bool = False,
                 **kwargs):
        """
        Initialize CPluginScript.

        Args:
            parent: Parent object (usually None for top-level scripts)
            name: Script instance name
            xmlFile: Path to input_params.xml file to load
            workDirectory: Working directory for the script (str or Path, optional)
            dummy: If True, skip def.xml loading and create minimal container (legacy API)
            **kwargs: Additional arguments
        """
        # Initialize CData base class (provides hierarchy and event system)
        super().__init__(parent=parent, name=name or self.TASKNAME, **kwargs)

        # Store dummy flag for later reference
        self._dummy = dummy

        # Create signals (inherits SignalManager from HierarchicalObject)
        from ccp4i2.core.base_object.signal_system import Signal
        # Emitted when the plugin completes execution
        self.finished = self._signal_manager.create_signal("finished", dict)
        # Emitted during streaming execution when progress is available (e.g., new refinement cycle)
        self.progressUpdated = self._signal_manager.create_signal("progressUpdated", dict)

        # Initialize infrastructure components
        self._def_parser = DefXmlParser()
        self._params_handler = ParamsXmlHandler()

        # Create main container
        # CPluginScript is now the parent of the container
        self.container = CContainer(parent=self, name="container")

        # Error report for tracking issues during execution
        self.errorReport = CErrorReport()

        # Process management
        self._process = None
        self._status = None

        # Asynchronous execution control (legacy API)
        # Set to True to run plugin asynchronously (non-blocking)
        # Set to False for synchronous (blocking) execution
        self.doAsync = False

        # Child job counter for sub-plugins (follows legacy convention)
        self._childJobCounter = 0

        # Database integration attributes (for database-backed environments)
        # These are set by the database handler when running in CCP4i2 GUI
        # Extract dbHandler from kwargs if provided (e.g., from get_job_plugin())
        self._dbHandler = kwargs.pop('dbHandler', None)
        self._dbProjectId = None       # Project identifier in database
        self._dbProjectName = None     # Project name
        self._dbJobId = None           # Job identifier in database
        self._dbJobNumber = None       # Job number (e.g., "1.2.3" for nested jobs)

        # Command line for external program
        self.commandLine = []

        # Command script (list of lines to write to stdin or script file)
        self.commandScript = []

        # Working directory and file paths
        if workDirectory is not None:
            self.workDirectory = Path(workDirectory)
        else:
            self.workDirectory = Path.cwd()
        self.defFile = None
        self.paramsFile = None

        # Load DEF file if available (defines container structure)
        # This will create inputData, outputData, controlParameters, guiAdmin
        # as children of self.container
        # Skip def.xml loading if dummy=True (creates minimal container only)
        #
        # IMPORTANT: Handle TWO inheritance patterns:
        #
        # Pattern 1 - XML-based inheritance (phaser, prosmart_refmac):
        #   - Child class inherits from parent (pythonic)
        #   - Child has .def.xml with <file>parent.def.xml</file> tag
        #   - load_nested_xml() automatically merges parent and child .def.xml
        #   - Example: phaser_simple.def.xml contains <file>phaser_pipeline.def.xml</file>
        #
        # Pattern 2 - Pythonic-only inheritance (crank2 sub-wrappers):
        #   - Child class inherits from parent (pythonic)
        #   - Child has NO .def.xml file
        #   - Must load parent's .def.xml using parent's TASKNAME
        #   - Example: crank2_refatompick has no .def.xml, inherits from crank2
        #
        # Strategy: Try to load child's .def.xml. If not found AND parent class exists,
        # try loading parent's .def.xml.
        if self.TASKNAME and not dummy:
            def_path = locate_def_xml(self.TASKNAME)

            if def_path and def_path.exists():
                # Pattern 1: Child has .def.xml (may contain <file> tag for parent)
                logger.debug(f"[__init__] Loading .def.xml for {self.TASKNAME}")
                self._loadDefFile()
            else:
                # Pattern 2: No .def.xml for child - try parent's TASKNAME
                parent_classes = [c for c in self.__class__.__mro__[1:]
                                if c.__name__ != 'CPluginScript'
                                and issubclass(c, CPluginScript)
                                and hasattr(c, 'TASKNAME')]

                if parent_classes and parent_classes[0].TASKNAME:
                    parent_taskname = parent_classes[0].TASKNAME
                    logger.debug(f"[__init__] No .def.xml for {self.TASKNAME}, trying parent {parent_taskname}")

                    # Temporarily swap TASKNAME to load parent's def file
                    original_taskname = self.TASKNAME
                    self.TASKNAME = parent_taskname
                    self._loadDefFile()
                    self.TASKNAME = original_taskname
                else:
                    # No .def.xml and no parent - will use default containers
                    logger.debug(f"[__init__] No .def.xml found for {self.TASKNAME} and no parent class")

        # Create default empty sub-containers ONLY if they don't exist after .def.xml loading
        # This ensures backward compatibility for plugins without .def.xml files
        # Also used for dummy plugins which skip def.xml loading entirely
        self._ensure_standard_containers()

        # Set standard container ordering for serialization
        # Standard containers come first, then any custom containers from .def.xml
        self._set_container_order()

        # Load PARAMS file if provided (actual parameter values)
        if xmlFile:
            self.loadDataFromXml(xmlFile)

    # Getter/setter methods for database-related attributes
    def get_status(self) -> Optional[int]:
        """Get the current plugin status."""
        # Handle case where _status might be a dict (from postProcessWrapper)
        if isinstance(self._status, dict):
            return self._status.get('finishStatus', self.FAILED)
        return self._status

    def get_db_job_id(self) -> Optional[str]:
        """Get the database job UUID."""
        return self._dbJobId

    def set_db_job_id(self, job_id: str) -> None:
        """Set the database job UUID."""
        self._dbJobId = job_id

    def get_db_job_number(self) -> Optional[str]:
        """Get the database job number (e.g., '1.2.3')."""
        return self._dbJobNumber

    def set_db_job_number(self, job_number: str) -> None:
        """Set the database job number."""
        self._dbJobNumber = job_number

    def setDbData(self, handler=None, projectName=None, projectId=None,
                  jobNumber=None, jobId=None) -> None:
        """
        Set database context for this plugin (legacy CCP4i2 API).

        This method provides compatibility with legacy code that calls
        setDbData() to configure the plugin's database context.

        Args:
            handler: Database handler instance (CCP4i2DjangoDbHandler)
            projectName: Project name string
            projectId: Project UUID (may have hyphens removed)
            jobNumber: Job number string (e.g., '1.2.3')
            jobId: Job UUID (may have hyphens removed)
        """
        if handler is not None:
            self._dbHandler = handler

        if projectName is not None:
            self._dbProjectName = projectName

        if projectId is not None:
            # Store as-is (may or may not have hyphens)
            self._dbProjectId = projectId

        if jobNumber is not None:
            self.set_db_job_number(jobNumber)

        if jobId is not None:
            # Restore hyphens if they were removed
            if '-' not in jobId and len(jobId) == 32:
                # UUID without hyphens - restore them
                jobId = f"{jobId[0:8]}-{jobId[8:12]}-{jobId[12:16]}-{jobId[16:20]}-{jobId[20:]}"
            self.set_db_job_id(jobId)
            logger.debug(f"[DEBUG setDbData] After set_db_job_id: self._dbJobId = {self._dbJobId}, get_db_job_id() = {self.get_db_job_id()}")

    def _ensure_standard_containers(self):
        """
        Ensure standard sub-containers exist.

        Creates inputData, outputData, controlParameters, and guiAdmin
        as children of self.container only if they don't already exist.
        This ensures backward compatibility for plugins without .def.xml files.
        """
        from ccp4i2.core.base_object.fundamental_types import CString, CInt

        standard_containers = ['inputData', 'outputData', 'controlParameters', 'guiAdmin']

        for container_name in standard_containers:
            # Check if container exists as a child
            try:
                # Try to access via __getattr__ (which searches children)
                container = getattr(self.container, container_name)
                logger.debug(f"[_ensure_standard_containers] {container_name} already exists with {len(container.children())} children")
            except AttributeError:
                # Container doesn't exist - create it
                logger.debug(f"[_ensure_standard_containers] Creating new {container_name} container")
                container = CContainer(
                    parent=self.container,
                    name=container_name
                )
                # Use setattr to properly register in _data_order for serialization
                setattr(self.container, container_name, container)

            # Add standard fields to guiAdmin
            if container_name == 'guiAdmin':
                # Ensure jobTitle exists (used by arcimboldo, baverage, etc.)
                if not hasattr(container, 'jobTitle'):
                    job_title = CString(
                        parent=container,
                        name='jobTitle'
                    )
                    # Set default value from job database if available
                    if hasattr(self, '_dbJobId') and self._dbJobId:
                        try:
                            from ccp4i2.db.models import Job
                            job = Job.objects.get(uuid=self._dbJobId)
                            if job.name:
                                job_title.value = job.name
                        except Exception:
                            # If we can't get job name from DB, leave it unset
                            pass
                    # Use setattr to properly register in _data_order for serialization
                    setattr(container, 'jobTitle', job_title)

                # Ensure jobStatus exists (stores job completion status)
                # Values: 0=Pending, 1=Running, 2=Finished, 3=Failed, etc.
                if not hasattr(container, 'jobStatus'):
                    job_status = CInt(
                        parent=container,
                        name='jobStatus'
                    )
                    job_status.value = 0  # Default: Pending
                    # Use setattr to properly register in _data_order for serialization
                    setattr(container, 'jobStatus', job_status)

    def _set_container_order(self):
        """
        Set the standard ordering for container children.

        Ensures that standard containers (inputData, outputData, controlParameters,
        guiAdmin) appear first in serialization order, followed by any custom
        containers from .def.xml (e.g., prosmartProtein, prosmartNucleicAcid).

        This order is used by dataOrder() which is called by CCP4i2JsonEncoder
        for JSON serialization.
        """
        # Standard container order (these should appear first)
        standard_order = ['inputData', 'outputData', 'controlParameters', 'guiAdmin']

        # Get all current child names
        all_children = [child.objectName() for child in self.container.children()
                       if hasattr(child, 'objectName') and child.objectName()]

        # Build ordered list: standard containers first, then others
        ordered = []
        for name in standard_order:
            if name in all_children:
                ordered.append(name)

        # Add any remaining children (custom containers from .def.xml)
        for name in all_children:
            if name not in ordered:
                ordered.append(name)

        # Set CONTENT_ORDER on the container
        self.container.CONTENT_ORDER = ordered
        logger.debug(f"[DEBUG _set_container_order] Set CONTENT_ORDER: {ordered}")

    def _loadDefFile(self):
        """
        Load the .def.xml file for this task.

        Locate the .def.xml file, then uses
        DefXmlParser to load the structure into containers.
        """
        logger.debug(f"[_loadDefFile] Looking for .def.xml for task: {self.TASKNAME}")
        def_path = locate_def_xml(self.TASKNAME)

        if def_path and def_path.exists():
            logger.debug(f"[_loadDefFile] Found .def.xml at: {def_path}")
            # Load using DefXmlParser
            error = self.loadContentsFromXml(str(def_path))
            if error:
                logger.error(f"[FATAL] Failed to load .def.xml for task '{self.TASKNAME}': {error}")
                raise RuntimeError(f"Failed to load .def.xml for task '{self.TASKNAME}': {error}")
            else:
                logger.debug(f"[_loadDefFile] Successfully loaded .def.xml")
        else:
            # WARNING: .def.xml not found - plugin will use default containers
            # This is normal for some legacy plugins (e.g., crank2 sub-wrappers)
            # that inherit container structure from parent classes
            warning_msg = (
                f"[WARNING] No .def.xml file found for task '{self.TASKNAME}'\n"
                f"  Looked at: {def_path}\n"
                f"  Default containers (inputData, outputData, etc.) will be created.\n"
                f"  This is normal for some legacy plugins (e.g., crank2 sub-wrappers)."
            )
            logger.warning(warning_msg)

    def loadContentsFromXml(self, fileName: str) -> CErrorReport:
        """
        Load container structure from a DEF file using DefXmlParser.

        Args:
            fileName: Path to .def.xml file

        Returns:
            CErrorReport indicating success or failure
        """
        error = CErrorReport()
        try:
            logger.debug(f"[loadContentsFromXml] Parsing .def.xml file: {fileName}")
            # Use DefXmlParser to parse the .def.xml file
            parsed_container = self._def_parser.parse_def_xml(fileName)

            # The parsed container has the full hierarchy
            # We need to extract ALL sub-containers and attach them to our
            # container (which is parented to this CPluginScript instance)
            # This includes standard containers (inputData, outputData, controlParameters, guiAdmin)
            # AND any custom containers (e.g., prosmartProtein, prosmartNucleicAcid for pipelines)

            # Iterate through all children of the parsed container
            logger.debug(f"[loadContentsFromXml] Attaching children to self.container:")
            for child in parsed_container.children():
                child_name = child.objectName()
                logger.debug(f"[loadContentsFromXml]   - {child_name} (type={type(child).__name__})")
                # Attach the child to our container using __setattr__
                setattr(self.container, child_name, child)
                # Update parent to be our container (which is parented to self)
                child.set_parent(self.container)

            self.defFile = fileName
            logger.debug(f"[loadContentsFromXml] Successfully loaded .def.xml and attached all children")

        except Exception as e:
            logger.error(f"[DEBUG loadContentsFromXml] Exception loading .def.xml: {e}")
            import traceback
            logger.error(f"[DEBUG loadContentsFromXml] Traceback:\n{traceback.format_exc()}")
            error.append(
                klass=self.__class__.__name__,
                code=100,
                details=f"Failed to load DEF file {fileName}: {e}",
                name=str(fileName)
            )

        return error

    def loadContentsFromEtree(self, element: ET.Element) -> CErrorReport:
        """
        Load container structure from an eTree element.

        Args:
            element: eTree element containing container definitions

        Returns:
            CErrorReport indicating success or failure
        """
        error = CErrorReport()
        # This method is kept for API compatibility but delegates to
        # DefXmlParser internally if needed in the future
        error.append(
            klass=self.__class__.__name__,
            code=101,
            details="loadContentsFromEtree not yet implemented",
            name=self.objectName() or ""
        )
        return error

    def loadDataFromXml(self, fileName: str) -> CErrorReport:
        """
        Load parameter values from a PARAMS file using ParamsXmlHandler.

        Args:
            fileName: Path to .params.xml or input_params.xml file

        Returns:
            CErrorReport indicating success or failure
        """
        error = CErrorReport()
        try:
            # Use ParamsXmlHandler to import params
            success = self._params_handler.import_params_xml(
                self.container, fileName)

            if not success:
                error.append(
                    klass=self.__class__.__name__,
                    code=102,
                    details=f"Failed to load data from {fileName}",
                    name=str(fileName)
                )
            else:
                self.paramsFile = fileName

        except Exception as e:
            error.append(
                klass=self.__class__.__name__,
                code=103,
                details=f"Failed to load data from {fileName}: {e}",
                name=str(fileName)
            )

        return error

    def loadDataFromEtree(self, element: ET.Element) -> CErrorReport:
        """
        Load parameter values from an eTree element.

        Args:
            element: eTree element containing parameter values

        Returns:
            CErrorReport indicating success or failure
        """
        error = CErrorReport()
        try:
            self.container.loadDataFromEtree(element)
        except Exception as e:
            error.append(
                klass=self.__class__.__name__,
                code=104,
                details=f"Failed to load data from etree: {e}",
                name=self.objectName() or ""
            )
        return error

    def saveContentsToXml(self, fileName: str) -> CErrorReport:
        """
        Save container structure to a DEF file.

        Args:
            fileName: Path to output .def.xml file

        Returns:
            CErrorReport indicating success or failure
        """
        error = CErrorReport()
        # DEF file saving is typically not done by tasks
        # Kept for API compatibility
        error.append(
            klass=self.__class__.__name__,
            code=105,
            details="saveContentsToXml not yet implemented",
            name=self.objectName() or ""
        )
        return error

    def saveDataToXml(self, fileName: str, exclude_unset: bool = True) -> CErrorReport:
        """
        Save parameter values to a PARAMS file using ParamsXmlHandler.

        Args:
            fileName: Path to output .params.xml file
            exclude_unset: If True, only save parameters that have been explicitly set

        Returns:
            CErrorReport indicating success or failure
        """
        import logging
        logger = logging.getLogger(f"ccp4i2:{__name__}")
        logger.info(f"saveDataToXml called with fileName: {fileName}, exclude_unset: {exclude_unset}")

        error = CErrorReport()
        try:
            # Use ParamsXmlHandler to export params
            logger.info(f"Calling _params_handler.export_params_xml...")
            success = self._params_handler.export_params_xml(
                self.container, fileName, exclude_unset=exclude_unset)
            logger.info(f"export_params_xml returned: {success}")

            if not success:
                error.append(
                    klass=self.__class__.__name__,
                    code=106,
                    details=f"Failed to save data to {fileName}",
                    name=str(fileName)
                )

        except Exception as e:
            logger.exception(f"Exception in saveDataToXml: {e}")
            error.append(
                klass=self.__class__.__name__,
                code=107,
                details=f"Failed to save data to {fileName}: {e}",
                name=str(fileName)
            )

        return error

    def saveParams(self, fileName: Optional[str] = None, exclude_unset: bool = True) -> CErrorReport:
        """
        Save current parameters to PARAMS file using ParamsXmlHandler.

        This is typically called at the end of execution to save the
        final state including any output data that was generated.

        Args:
            fileName: Path to output file (defaults to auto-generated name)
            exclude_unset: If True, only save parameters that have been explicitly set

        Returns:
            CErrorReport indicating success or failure
        """
        if fileName is None:
            # Use TASKNAME (not self.name which may contain spaces/special chars from job title)
            fileName = str(self.workDirectory / f"{self.TASKNAME}.params.xml")
        return self.saveDataToXml(fileName, exclude_unset=exclude_unset)

    # =========================================================================
    # Process workflow methods
    # =========================================================================

    def process(self) -> int:
        """
        Main processing method - orchestrates the entire workflow.

        This method calls the following steps in order:
        1. validity() - validate plugin (allows qualifier adjustments before file checks)
        2. checkOutputData() - set output file names
        3. processInputFiles() - pre-process input files
        4. makeCommandAndScript() - generate command line/file
        5. startProcess() - execute the program

        Returns:
            Status code (SUCCEEDED, FAILED, or RUNNING)
        """
        # Validate input data using validity() which allows plugins to adjust
        # qualifiers for embedded wrappers before checkInputData() runs
        # (e.g., servalcat_pipe sets allowUndefined on metalCoordWrapper.inputData.XYZIN)
        try:
            error = self.validity()
            if error:
                self.errorReport.extend(error)
                if error.maxSeverity() >= 4:  # ERROR level
                    # Emit finished signal so pipelines are notified of failure
                    self.reportStatus(self.FAILED)
                    return self.FAILED
        except Exception as e:
            # Catch Python exceptions in validity() and add to error report
            import traceback
            tb_str = ''.join(traceback.format_exception(type(e), e, e.__traceback__))
            self.errorReport.append(
                klass=self.__class__.__name__,
                code=998,
                details=f'Python exception in validity(): {type(e).__name__}: {str(e)}\n\n{tb_str}',
                name='validity',
                severity=4  # ERROR
            )
            self.reportStatus(self.FAILED)
            return self.FAILED

        # Set up output data
        try:
            error = self.checkOutputData()
            if error:
                self.errorReport.extend(error)
                # Don't fail - checkOutputData should fix issues
        except Exception as e:
            # Catch Python exceptions in checkOutputData() and add to error report
            import traceback
            tb_str = ''.join(traceback.format_exception(type(e), e, e.__traceback__))
            self.errorReport.append(
                klass=self.__class__.__name__,
                code=997,
                details=f'Python exception in checkOutputData(): {type(e).__name__}: {str(e)}\n\n{tb_str}',
                name='checkOutputData',
                severity=4  # ERROR
            )
            self.reportStatus(self.FAILED)
            return self.FAILED

        # Save params.xml after setting output file attributes
        # This ensures output files have project/relPath/baseName set before execution
        # NOTE: We save to params.xml (not input_params.xml) to preserve the original inputs
        # NOTE: We do NOT reload from input_params.xml here because:
        # 1. The async runner already loaded params from input_params.xml (async_run_job.py line 172)
        # 2. The async runner imported files and updated the container (async_run_job.py line 73)
        # 3. The async runner saved params.xml (async_run_job.py line 108)
        # 4. Reloading here would corrupt the file paths that were set during import
        #
        # The file metadata should already be correct from the import step.
        has_method = hasattr(self, 'get_db_job_id')
        job_id = self.get_db_job_id() if has_method else None
        pass  # DEBUG: print(f"[DEBUG process] Checking if should save params: has_method={has_method}, job_id={job_id}")
        pass  # DEBUG: print(f"[DEBUG process] self.workDirectory = {self.workDirectory}")
        if has_method and job_id:
            try:
                # Just save params.xml with current state (which already has correct file paths from import)
                params_path = os.path.join(self.workDirectory, "params.xml")
                pass  # DEBUG: print(f"[DEBUG process] About to call saveDataToXml({params_path})")
                save_error = self.saveDataToXml(params_path)
                if save_error and hasattr(save_error, 'hasError') and save_error.hasError():
                    pass  # DEBUG: print(f"[DEBUG process] Warning: Failed to save params.xml after checkOutputData: {save_error}")
                else:
                    pass  # DEBUG: print(f"[DEBUG process] Saved params.xml after checkOutputData with output file attributes")
            except Exception as e:
                pass  # DEBUG: print(f"[DEBUG process] Warning: Exception saving params.xml: {e}")
                import traceback
                traceback.print_exc()
        else:
            pass  # DEBUG: print(f"[DEBUG process] Skipping params save (no database context)")

        # Pre-process input files if needed
        try:
            result = self.processInputFiles()
            # Handle both modern API (CErrorReport) and legacy API (int)
            if isinstance(result, int):
                # Legacy API: returns SUCCEEDED (0) or FAILED (1)
                if result != self.SUCCEEDED:
                    # Emit finished signal so pipelines are notified of failure
                    self.reportStatus(result)
                    return result
            elif result:
                # Modern API: returns CErrorReport (truthy if has errors)
                self.errorReport.extend(result)
                # Emit finished signal so pipelines are notified of failure
                self.reportStatus(self.FAILED)
                return self.FAILED
        except Exception as e:
            # Catch Python exceptions in processInputFiles() and add to error report
            import traceback
            tb_str = ''.join(traceback.format_exception(type(e), e, e.__traceback__))
            self.errorReport.append(
                klass=self.__class__.__name__,
                code=996,
                details=f'Python exception in processInputFiles(): {type(e).__name__}: {str(e)}\n\n{tb_str}',
                name='processInputFiles',
                severity=4  # ERROR
            )
            self.reportStatus(self.FAILED)
            return self.FAILED

        # Generate command and script
        try:
            error = self.makeCommandAndScript()
            if error:
                self.errorReport.extend(error)
                # Emit finished signal so pipelines are notified of failure
                self.reportStatus(self.FAILED)
                return self.FAILED
        except Exception as e:
            # Catch Python exceptions in makeCommandAndScript() and add to error report
            import traceback
            tb_str = ''.join(traceback.format_exception(type(e), e, e.__traceback__))
            self.errorReport.append(
                klass=self.__class__.__name__,
                code=995,
                details=f'Python exception in makeCommandAndScript(): {type(e).__name__}: {str(e)}\n\n{tb_str}',
                name='makeCommandAndScript',
                severity=4  # ERROR
            )
            self.reportStatus(self.FAILED)
            return self.FAILED

        # Start the process
        # Legacy compatibility: plugins have various startProcess signatures
        # Inspect the signature and call with appropriate arguments
        try:
            result = self.startProcess()

            # Handle both modern API (CErrorReport) and legacy API (int)
            if isinstance(result, int):
                # Legacy API: returns SUCCEEDED (0) or FAILED (1)
                if result != self.SUCCEEDED:
                    # Emit finished signal so pipelines are notified of failure
                    self.reportStatus(result)
                    return result
            elif result:
                # Modern API: returns CErrorReport (truthy if has errors)
                self.errorReport.extend(result)
                # Emit finished signal so pipelines are notified of failure
                self.reportStatus(self.FAILED)
                return self.FAILED
        except Exception as e:
            # Catch Python exceptions in startProcess() and add to error report
            import traceback
            tb_str = ''.join(traceback.format_exception(type(e), e, e.__traceback__))
            self.errorReport.append(
                klass=self.__class__.__name__,
                code=994,
                details=f'Python exception in startProcess(): {type(e).__name__}: {str(e)}\n\n{tb_str}',
                name='startProcess',
                severity=4  # ERROR
            )
            self.reportStatus(self.FAILED)
            return self.FAILED

        # For synchronous execution (subprocess.run), process is complete when startProcess returns
        # Check if the process succeeded by examining exit code
        status, exit_status, exit_code = self.postProcessCheck()
        import logging
        logger = logging.getLogger(__name__)
        logger.debug(f"[process] postProcessCheck returned status: {status}, exitStatus: {exit_status}, exitCode: {exit_code} (SUCCEEDED={self.SUCCEEDED}, FAILED={self.FAILED})")

        # Only proceed with output processing if the process succeeded
        if status == self.SUCCEEDED:
            # Call processOutputFiles to extract output data
            try:
                error = self.processOutputFiles()
                if error:
                    self.errorReport.extend(error)
                    # Don't fail the job for processOutputFiles errors if the process succeeded
                    # Legacy wrappers often have non-fatal issues in processOutputFiles
                    print(f"Warning: processOutputFiles() returned errors but process completed successfully")
                    # status = self.FAILED  # Commented out - don't fail for postprocessing errors
            except Exception as e:
                # Legacy wrappers may not have processOutputFiles implemented
                print(f"Warning: processOutputFiles() exception: {type(e).__name__}: {str(e)}")
                import traceback
                traceback.print_exc()
                # Don't change status - process succeeded even if postprocessing failed

            # Glean output files to database if in database-connected mode
            # This is essential for subjobs created via makePluginObject() which don't go
            # through the async track_job context manager
            self._glean_output_files_sync()

        # Emit finished signal so pipelines can continue
        # This is essential for sub-plugins in pipelines
        logger.debug(f"[process] Calling reportStatus with status: {status}")
        self.reportStatus(status)
        logger.debug(f"[process] Returning status: {status}")

        return status

    def _find_datafile_descendants(self, container) -> list:
        """
        Recursively find all CDataFile descendants in a container hierarchy.

        Handles both CContainer children and CList elements.

        Args:
            container: Container to search

        Returns:
            List of (name, file_obj) tuples for all CDataFile descendants
        """
        from ccp4i2.core.base_object.base_classes import CDataFile
        from ccp4i2.core.base_object.fundamental_types import CList

        results = []

        # Check all children of this container
        for child in container.children():
            # If it's a CDataFile, add it
            if isinstance(child, CDataFile):
                results.append((child.objectName(), child))
            # If it's a CList, check its elements
            elif isinstance(child, CList):
                for i, item in enumerate(child):
                    # If the list item is a CDataFile, add it
                    if isinstance(item, CDataFile):
                        # Use list element name if available, otherwise use index
                        item_name = item.objectName() if item.objectName() else f"{child.objectName()}[{i}]"
                        results.append((item_name, item))
                    # If the list item is a container, recurse into it
                    elif hasattr(item, 'children'):
                        results.extend(self._find_datafile_descendants(item))
            # If it's a container, recurse into it
            elif hasattr(child, 'children'):
                results.extend(self._find_datafile_descendants(child))

        return results

    def _glean_output_files_sync(self):
        """
        Glean output files to database (synchronous wrapper).

        This method is called at the end of process() to register output files
        in the database. It's essential for subjobs created via makePluginObject()
        which don't go through the async track_job context manager.

        Only acts if:
        - Plugin is NOT tracked by track_job (which handles its own gleaning)
        - _dbHandler is set
        - _dbJobId is set
        - container.outputData exists
        """
        # Skip if being tracked by track_job context manager
        # track_job handles gleaning for top-level jobs to avoid double-gleaning
        if getattr(self, '_tracked_by_track_job', False):
            logger.debug(f"[_glean_output_files_sync] Skipping - tracked by track_job context manager")
            return

        # Check if we're in database mode
        if not hasattr(self, '_dbHandler') or self._dbHandler is None:
            return
        if not hasattr(self, '_dbJobId') or self._dbJobId is None:
            return
        if not hasattr(self.container, 'outputData') or self.container.outputData is None:
            return

        logger.debug(f"[_glean_output_files_sync] Gleaning output files for subjob {self._dbJobId}")

        try:
            # Use async_to_sync to call the async glean method
            from asgiref.sync import async_to_sync
            import uuid as uuid_module

            # Normalize job UUID
            job_uuid = self._dbJobId
            if isinstance(job_uuid, str):
                if '-' not in job_uuid and len(job_uuid) == 32:
                    job_uuid = uuid_module.UUID(job_uuid)
                else:
                    job_uuid = uuid_module.UUID(job_uuid)

            # Call glean_job_files
            files_gleaned = async_to_sync(self._dbHandler.glean_job_files)(
                job_uuid,
                self.container.outputData,
                plugin=self
            )
            logger.debug(f"[_glean_output_files_sync] Gleaned {len(files_gleaned)} output files for {self.__class__.__name__}")

            # Also glean KPIs if available
            kpis_gleaned = async_to_sync(self._dbHandler.glean_performance_indicators)(
                job_uuid,
                self.container.outputData
            )
            logger.debug(f"[_glean_output_files_sync] Gleaned {kpis_gleaned} performance indicators")

        except Exception as e:
            # Don't fail the job if gleaning fails
            logger.warning(f"[_glean_output_files_sync] Failed to glean output files: {e}")
            import traceback
            traceback.print_exc()

    def checkInputData(self) -> CErrorReport:
        """
        Check that required input data is set.

        This method recursively checks all CDataFile descendants in the
        entire container hierarchy and verifies that files marked with
        mustExist actually exist.

        Only returns actual errors (SEVERITY_ERROR), not warnings about
        optional files. Warnings are for GUI display, not for blocking
        job execution.

        Returns:
            CErrorReport with any validation errors (no warnings)
        """
        from ccp4i2.core.base_object.error_reporting import SEVERITY_ERROR

        error = CErrorReport()

        # Find all CDataFile descendants recursively in the entire container
        file_items = self._find_datafile_descendants(self.container)

        # Validate each file
        for name, obj in file_items:
            obj_error = obj.validity()
            if obj_error:
                # Only include actual errors, not warnings
                # Warnings are for GUI display, not for blocking execution
                for err_item in obj_error._errors:
                    if err_item.get('severity', 0) >= SEVERITY_ERROR:
                        # Directly append to _errors list since err_item dict uses
                        # 'class' as key but append() takes 'klass' as parameter
                        error._errors.append(err_item)

        return error

    def validity(self) -> CErrorReport:
        """
        Validate the plugin's container and return an error report.

        This method provides comprehensive validation of the plugin's container,
        combining container validation with input file checking. It is the
        canonical way to validate a CPluginScript instance.

        The validation includes:
        - Container validity (all child objects' validity() methods)
        - Input file existence checks (mustExist files)

        Returns:
            CErrorReport containing all validation errors/warnings

        Example:
            >>> plugin = get_job_plugin(job)
            >>> errors = plugin.validity()
            >>> if errors.maxSeverity() >= SEVERITY_ERROR:
            ...     print(f"Validation failed: {errors.report()}")
        """
        error = CErrorReport()

        # Container validation recursively validates all children,
        # including CDataFile objects (allowUndefined, mustExist, contentFlag).
        # No separate checkInputData() call needed - that would duplicate errors.
        if hasattr(self, 'container') and self.container is not None:
            container_errors = self.container.validity()
            if container_errors:
                error.extend(container_errors)

        return error

    def validity_as_xml(self):
        """
        Validate the plugin's container and return an XML Element.

        This is a convenience method that calls validity() and converts
        the CErrorReport to an XML Element tree, suitable for serialization
        or API responses.

        Returns:
            xml.etree.ElementTree.Element: Root element 'errorReportList'
            containing 'errorReport' children

        Example:
            >>> plugin = get_job_plugin(job)
            >>> error_xml = plugin.validity_as_xml()
            >>> xml_str = ET.tostring(error_xml, encoding='unicode')
        """
        from xml.etree import ElementTree as ET
        from ccp4i2.core import CCP4ErrorHandling

        error_report = self.validity()

        # Convert CErrorReport to XML using the same logic as validate_container
        element = ET.Element("errorReportList")

        # Mapping from severity codes to text
        SEVERITY_TEXT = {
            CCP4ErrorHandling.SEVERITY_OK: "OK",
            CCP4ErrorHandling.SEVERITY_UNDEFINED: "UNDEFINED",
            CCP4ErrorHandling.SEVERITY_WARNING: "WARNING",
            CCP4ErrorHandling.SEVERITY_UNDEFINED_ERROR: "UNDEFINED_ERROR",
            CCP4ErrorHandling.SEVERITY_ERROR: "ERROR"
        }

        for item in error_report.getErrors():
            try:
                ele = ET.Element("errorReport")

                e = ET.Element("className")
                class_name = item["class"] if isinstance(item["class"], str) else item["class"].__name__
                e.text = class_name
                ele.append(e)

                e = ET.Element("code")
                e.text = str(item["code"])
                ele.append(e)

                e = ET.Element("description")
                e.text = item["details"]
                ele.append(e)

                e = ET.Element("severity")
                severity = item["severity"]
                e.text = SEVERITY_TEXT.get(severity, f"UNKNOWN({severity})")
                ele.append(e)

                # Add objectPath from 'name' field (which contains the object path)
                if item.get("name"):
                    e = ET.Element("objectPath")
                    e.text = item["name"]
                    ele.append(e)

                element.append(ele)
            except Exception as e:
                logger.exception("Error converting error to XML", exc_info=e)

        ET.indent(element, " ")
        return element

    def checkOutputData(self) -> CErrorReport:
        """
        Set output file names if not already set.

        This method generates appropriate file names for any output files
        that don't have names yet, based on their objectName() and workDirectory.

        For lists of CDataFiles, pre-populates up to maxListElements items.

        Returns:
            CErrorReport with any issues (should fix rather than fail)
        """
        import os
        import re
        from ccp4i2.core.base_object.base_classes import CDataFile
        from ccp4i2.core.base_object.fundamental_types import CList

        logger.debug(f"[DEBUG checkOutputData] Called for task: {self.TASKNAME if hasattr(self, 'TASKNAME') else 'unknown'}")
        logger.debug(f"[DEBUG checkOutputData] _dbProjectId = {getattr(self, '_dbProjectId', 'NOT SET')}")
        logger.debug(f"[DEBUG checkOutputData] _dbJobNumber = {getattr(self, '_dbJobNumber', 'NOT SET')}")
        logger.debug(f"[DEBUG checkOutputData] workDirectory = {self.workDirectory}")

        error = CErrorReport()

        if not hasattr(self.container, 'outputData'):
            return error

        # Get maxListElements setting (default 50)
        max_list_elements = getattr(self, 'maxListElements', 50)

        def slugify(name: str) -> str:
            """Convert objectName to valid filename, removing special chars like braces."""
            # Remove or replace special characters
            name = re.sub(r'[\[\]{}()<>]', '', name)  # Remove braces and brackets
            name = re.sub(r'[^\w\s\-\.]', '_', name)  # Replace other special chars with underscore
            name = re.sub(r'[-\s]+', '_', name)  # Replace spaces and hyphens with underscore
            return name.strip('_')

        def populate_list_outputs(obj, parent_path: str = ""):
            """Pre-populate CList of CDataFiles with proper file paths."""
            if not isinstance(obj, CList):
                return

            # If list is empty or has fewer than max elements, populate it
            current_len = len(obj) if hasattr(obj, '__len__') else 0
            target_len = max_list_elements

            # Get the item type if specified
            item_class = None

            # If we couldn't determine item type, check if there are existing items
            if current_len > 0:
                item_class = type(obj[0])

            # Only populate if we know it's a CDataFile list
            if item_class and issubclass(item_class, CDataFile):
                for i in range(current_len, target_len):
                    # Create new instance
                    new_item = item_class()
                    # Generate file path
                    base_name = slugify(obj.objectName() or obj.name or "output")
                    file_name = f"{base_name}_{i}.mtz"  # Default to .mtz extension
                    file_path = os.path.join(self.workDirectory, file_name)
                    new_item.setFullPath(file_path)
                    # Add to list
                    obj.append(new_item)

        def process_container(container, parent_path: str = ""):
            """Recursively process container to set output file paths."""
            logger.debug(f'[DEBUG checkOutputData] Processing container: {container.objectName() or "unknown"}')
            children = list(container.children())
            logger.debug(f'[DEBUG checkOutputData] Found {len(children)} children')
            for child in children:
                logger.debug(f'[DEBUG checkOutputData] Processing child: {child.objectName() or "unknown"} (type: {type(child).__name__})')
                # Handle CDataFile
                if isinstance(child, CDataFile):
                    # Only set path if baseName has not been set by the user
                    # Check baseName rather than fullPath as it's more fundamental
                    # Also check that baseName is non-empty (DEF files may set it to empty string)
                    basename_is_set = False
                    if hasattr(child, 'baseName') and hasattr(child.baseName, 'isSet'):
                        basename_is_set = child.baseName.isSet('value') and bool(str(child.baseName).strip())

                    obj_name = child.objectName() if hasattr(child, 'objectName') else (child.name if hasattr(child, 'name') else 'unknown')
                    logger.debug(f'[DEBUG checkOutputData]   {obj_name}: basename_is_set={basename_is_set}, baseName={str(child.baseName) if hasattr(child, "baseName") else "N/A"}')
                    if hasattr(child, 'baseName'):
                        if hasattr(child.baseName, '_value_states'):
                            logger.debug(f'[DEBUG checkOutputData]     baseName._value_states={child.baseName._value_states}')
                        if hasattr(child.baseName, 'value'):
                            logger.debug(f'[DEBUG checkOutputData]     baseName.value={child.baseName.value}')

                    if not basename_is_set:
                        # Use setOutputPath if available (database-aware method)
                        # Otherwise fall back to setFullPath
                        if hasattr(child, 'setOutputPath') and self._dbProjectId and self._dbJobNumber:
                            from pathlib import Path
                            # Calculate relPath from job number (e.g., "CCP4_JOBS/job_36" for job 36)
                            relPath = Path("CCP4_JOBS").joinpath(
                                *[f"job_{num}" for num in str(self._dbJobNumber).split(".")]
                            )
                            logger.debug(f'[DEBUG checkOutputData] Calling setOutputPath with projectId={self._dbProjectId}, relPath={relPath}')
                            child.setOutputPath(
                                jobName="",  # No prefix
                                projectId=str(self._dbProjectId),
                                relPath=str(relPath)
                            )
                            retrieved_path = child.getFullPath()
                            logger.debug(f'[DEBUG checkOutputData] Retrieved path after setOutputPath: {retrieved_path}')
                        else:
                            # Fallback: generate simple local path
                            obj_name = child.objectName()
                            if not obj_name:
                                obj_name = 'output'

                            file_name = slugify(obj_name)
                            if not any(file_name.endswith(ext) for ext in ['.mtz', '.pdb', '.cif', '.log', '.xml']):
                                # Get default extension - prefer calling method over qualifier
                                default_ext = '.mtz'  # Fallback default

                                # First try calling fileExtensions() method if it exists (e.g., CPdbDataFile)
                                # This allows classes to determine extension dynamically based on contentFlag
                                if hasattr(child, 'fileExtensions') and callable(child.fileExtensions):
                                    try:
                                        file_exts = child.fileExtensions()
                                        if file_exts and isinstance(file_exts, list) and len(file_exts) > 0:
                                            default_ext = '.' + file_exts[0].lstrip('.')
                                            logger.debug(f'[DEBUG checkOutputData] Using extension {default_ext} from fileExtensions() method')
                                    except Exception as e:
                                        logger.debug(f'[DEBUG checkOutputData] Failed to call fileExtensions() method: {e}')
                                else:
                                    # Fallback: use qualifier from class metadata
                                    from ccp4i2.core.base_object.class_metadata import get_class_metadata_by_type
                                    try:
                                        meta = get_class_metadata_by_type(type(child))
                                        if meta and meta.qualifiers and 'fileExtensions' in meta.qualifiers:
                                            file_exts = meta.qualifiers['fileExtensions']
                                            if file_exts and isinstance(file_exts, list) and len(file_exts) > 0:
                                                default_ext = '.' + file_exts[0].lstrip('.')
                                                logger.debug(f'[DEBUG checkOutputData] Using default extension {default_ext} from fileExtensions qualifier={file_exts}')
                                    except Exception as e:
                                        logger.debug(f'[DEBUG checkOutputData] Failed to get class metadata: {e}')
                                file_name += default_ext

                            file_path = os.path.join(self.workDirectory, file_name)
                            logger.debug(f'[DEBUG checkOutputData] Setting path for {obj_name}: {file_path}')
                            child.setFullPath(file_path)
                            retrieved_path = child.getFullPath()
                            logger.debug(f'[DEBUG checkOutputData] Retrieved path: {retrieved_path}')

                # Handle CList - pre-populate if it contains CDataFiles
                elif isinstance(child, CList):
                    populate_list_outputs(child, parent_path)
                    # Also process any items already in the list
                    for item in child:
                        if isinstance(item, CDataFile):
                            # Check if baseName is already set
                            basename_is_set = False
                            if hasattr(item, 'baseName') and hasattr(item.baseName, 'isSet'):
                                basename_is_set = item.baseName.isSet('value')

                            if not basename_is_set:
                                obj_name = item.objectName() or item.name or 'output'
                                file_name = slugify(obj_name)
                                if not any(file_name.endswith(ext) for ext in ['.mtz', '.pdb', '.cif', '.log', '.xml']):
                                    # Get default extension from class metadata fileExtensions
                                    default_ext = '.mtz'  # Fallback default
                                    from ccp4i2.core.base_object.class_metadata import get_class_metadata_by_type
                                    try:
                                        meta = get_class_metadata_by_type(type(item))
                                        if meta and meta.qualifiers and 'fileExtensions' in meta.qualifiers:
                                            file_exts = meta.qualifiers['fileExtensions']
                                            if file_exts and isinstance(file_exts, list) and len(file_exts) > 0:
                                                default_ext = '.' + file_exts[0].lstrip('.')
                                    except Exception:
                                        pass
                                    file_name += default_ext
                                file_path = os.path.join(self.workDirectory, file_name)
                                item.setFullPath(file_path)
                        elif hasattr(item, 'children'):
                            process_container(item, parent_path)

                # Handle nested containers
                elif hasattr(child, 'children'):
                    process_container(child, parent_path)

        # Process outputData container
        try:
            process_container(self.container.outputData)
        except Exception as e:
            error.append(
                klass=self.__class__.__name__,
                code=50,
                details=f"Error setting output file paths: {str(e)}"
            )

        # Save params.xml after setting output file attributes
        # Do this in checkOutputData() rather than process() because some plugins
        # override process() and bypass the base implementation
        # NOTE: We save to params.xml (not input_params.xml) to preserve the original inputs
        if self.get_db_job_id():
            try:
                params_path = os.path.join(self.workDirectory, "params.xml")
                logger.debug(f"[DEBUG checkOutputData] Saving params.xml to {params_path}")
                save_error = self.saveDataToXml(params_path)
                if save_error and hasattr(save_error, 'hasError') and save_error.hasError():
                    logger.debug(f"[DEBUG checkOutputData] Warning: Failed to save params.xml: {save_error}")
                else:
                    logger.debug(f"[DEBUG checkOutputData] ✅ Saved params.xml with output file attributes")
            except Exception as e:
                logger.debug(f"[DEBUG checkOutputData] Warning: Exception saving params.xml: {e}")

        return error

    def processInputFiles(self) -> CErrorReport:
        """
        Pre-process input files before running main program.

        This is a hook for subclasses to perform any manipulations
        on input data or files before calling the main program.

        Returns:
            CErrorReport with any errors
        """
        return CErrorReport()

    def makeCommandAndScript(self) -> CErrorReport:
        """
        Generate command line and command file for the program.

        Returns:
            CErrorReport with any errors
        """
        return CErrorReport()

    def appendCommandLine(self, wordList=[], clear=False) -> CErrorReport:
        """
        Add text strings or list of strings to the command line.

        Args:
            wordList: String or list of strings to add to command line
            clear: If True, clear command line before appending

        Returns:
            CErrorReport with any errors
        """
        import re
        from ccp4i2.core.base_object.fundamental_types import CList

        error = CErrorReport()

        if clear:
            self.clearCommandLine()

        if not isinstance(wordList, list):
            wordList = [wordList]

        for item in wordList:
            if isinstance(item, (list, CList)):
                for subItem in item:
                    try:
                        myText = str(subItem)
                        self.commandLine.append(myText)
                    except Exception:
                        error.append(
                            klass=self.__class__.__name__,
                            code=2,
                            details="Error converting command line item to string"
                        )
            else:
                try:
                    myText = str(item)
                    # Remove newlines from command line arguments
                    myText = re.sub(r'\n', ' ', myText)
                    self.commandLine.append(myText)
                except Exception:
                    error.append(
                        klass=self.__class__.__name__,
                        code=2,
                        details="Error converting command line item to string"
                    )

        if error:
            self.errorReport.extend(error)

        return error

    def clearCommandLine(self):
        """Clear the command line list."""
        self.commandLine = []

    def makeFileName(self, format='COM', ext='', qualifier=None):
        """
        Generate consistent names for output files.

        Args:
            format: File type format (e.g., 'PROGRAMXML', 'LOG', 'REPORT')
            ext: File extension (unused, kept for compatibility)
            qualifier: Optional qualifier to modify the filename

        Returns:
            Full path to the file in the work directory
        """
        import os

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
            'COM': 'com.txt',
            'MGPICDEF': 'report.mgpic.py',
            'PIC': 'report.png',
            'RVAPIXML': 'i2.xml'
        }

        fileName = defNames.get(format, 'unknown.unk')
        if qualifier is not None:
            base, ext = fileName.split('.', 1)
            fileName = base + '_' + str(qualifier) + '.' + ext
        return os.path.join(self.workDirectory, fileName)

    def renameFile(self, src, dst):
        """
        Atomically rename a file from src to dst.

        This is used for atomic file writes - write to a temp file, then rename it.
        Renaming is atomic on most filesystems, preventing partial reads.

        Args:
            src: Source file path
            dst: Destination file path
        """
        try:
            os.rename(src, dst)
        except OSError as e:
            # If rename fails (e.g., cross-device link), fall back to copy+delete
            shutil.move(src, dst)

    def logFileText(self) -> str:
        """
        Read and return the contents of the log file.

        This is a legacy API method used by plugins to parse their output.

        Returns:
            Contents of the log file as a string, or empty string if file doesn't exist
        """
        log_path = self.makeFileName('LOG')
        try:
            with open(log_path, 'r') as f:
                return f.read()
        except FileNotFoundError:
            return ""
        except Exception as e:
            print(f"Warning: Error reading log file {log_path}: {e}")
            return ""

    def jobNumberString(self) -> str:
        """
        Return a string representation of the job number for annotations.

        In a full CCP4i2 environment with database integration, this would
        return something like "Job #123". For standalone testing without
        database, returns the task name.

        Returns:
            String describing the job (e.g., "parrot" or "Job #123")
        """
        # When running standalone (no database integration), use task name
        # In full CCP4i2, this would query the database for the job number
        if hasattr(self, 'jobId') and self.jobId:
            return f"Job #{self.jobId}"
        return self.TASKNAME or "Job"

    def startProcess(self) -> CErrorReport:
        """
        Start the external program process.

        Runs the program specified by TASKCOMMAND with the command line
        arguments built by makeCommandAndScript(), capturing output to log files.

        If ASYNCHRONOUS is True, uses AsyncProcessManager for non-blocking execution.
        Otherwise uses subprocess.run() for synchronous execution.

        If both commandLine and commandScript are defined, the commandLine is used
        to start the process and commandScript is fed as stdin.

        Returns:
            CErrorReport with any errors
        """

        error = CErrorReport()

        if not self.TASKCOMMAND:
            error.append(
                klass=self.__class__.__name__,
                code=100,
                details="No TASKCOMMAND specified for this plugin"
            )
            return error

        # Check if async execution is requested
        # Priority order for execution mode:
        # 1. doAsync instance variable (explicit override - highest priority)
        #    - doAsync=False OVERRIDES ASYNCHRONOUS=True
        #    - This allows i2run tests to force synchronous execution
        # 2. ASYNCHRONOUS class variable (set in plugin definition)
        if self.doAsync is False:
            # Explicit synchronous override
            return self._startProcessSync()
        if self.ASYNCHRONOUS or self.doAsync:
            return self._startProcessAsync()
        return self._startProcessSync()

    def _prepareProcessExecution(self):
        """
        Common setup for both sync and async process execution.

        Returns:
            dict with keys:
                - command: List of command parts [exe_path, *args]
                - command_script_file: Path to script file or None
                - stdout_path: Path to LOG file
                - stderr_path: Path to STDERR file
                - env: Copy of os.environ
                - stdin_input: String to send to stdin or None
        """
        import os
        from pathlib import Path

        # Ensure working directory exists
        work_dir_path = Path(self.workDirectory)
        if not work_dir_path.exists():
            work_dir_path.mkdir(parents=True, exist_ok=True)
            print(f"Created working directory: {self.workDirectory}")

        # Write command script to file if present
        command_script_file = None
        if self.commandScript:
            command_script_file = self.writeCommandFile()

        # Prepare log file paths
        # Use LOG for program output (not STDOUT) to match CCP4 conventions
        stdout_path = self.makeFileName('LOG')
        stderr_path = self.makeFileName('STDERR')

        # Find full path to executable to ensure subprocess can find it
        exe_path = shutil.which(self.TASKCOMMAND)
        if exe_path:
            # Use full path if found
            command = [exe_path] + self.commandLine
        else:
            # Fall back to command name
            command = [self.TASKCOMMAND] + self.commandLine

        # Copy environment to ensure subprocess inherits all variables
        env = os.environ.copy()

        # Prepare stdin input
        stdin_input = None
        if self.commandScript:
            # Join command script lines into single string
            stdin_input = ''.join(self.commandScript)

        print(f"\n{'='*60}")
        print(f"Command: {' '.join(command)}")
        print(f"Working directory: {self.workDirectory}")
        print(f"Log file: {stdout_path}")
        print(f"Stderr file: {stderr_path}")
        if command_script_file:
            print(f"Command script: {command_script_file}")
        print(f"Environment CBIN: {env.get('CBIN', 'NOT SET')}")
        print(f"Environment CCP4: {env.get('CCP4', 'NOT SET')}")
        print(f"Environment CLIB: {env.get('CLIB', 'NOT SET')}")
        print(f"{'='*60}\n")

        return {
            'command': command,
            'command_script_file': command_script_file,
            'stdout_path': stdout_path,
            'stderr_path': stderr_path,
            'env': env,
            'stdin_input': stdin_input
        }

    def _startProcessSync(self) -> CErrorReport:
        """Synchronous process execution using subprocess.run()."""
        import subprocess
        import os

        error = CErrorReport()

        # Register with PROCESSMANAGER so it can return our exit code
        from ccp4i2.core.CCP4Modules import PROCESSMANAGER
        PROCESSMANAGER().register(self)

        # Prepare execution environment (common setup)
        prep = self._prepareProcessExecution()
        command = prep['command']
        stdout_path = prep['stdout_path']
        stderr_path = prep['stderr_path']
        stdin_input = prep['stdin_input']
        env = prep['env']

        try:

            with open(stdout_path, 'w') as stdout_file, open(stderr_path, 'w') as stderr_file:
                # Write formatted header to stdout
                stdout_file.write("="*70 + "\n")
                stdout_file.write(f"CCP4i2 Task: {self.TASKNAME}\n")
                stdout_file.write("="*70 + "\n\n")

                # Write command line
                stdout_file.write("Command Line:\n")
                stdout_file.write("-" * 70 + "\n")
                stdout_file.write(f"{' '.join(command)}\n\n")

                # Write command script if present
                if self.commandScript:
                    stdout_file.write("Command Script (stdin):\n")
                    stdout_file.write("-" * 70 + "\n")
                    for line in self.commandScript:
                        stdout_file.write(line)
                    stdout_file.write("\n")

                stdout_file.write("="*70 + "\n")
                stdout_file.write("Program Output:\n")
                stdout_file.write("="*70 + "\n\n")
                stdout_file.flush()

                # Run the process
                # Note: When input= is provided, stdin is automatically set to PIPE and closed after writing
                # When input= is None, stdin defaults to inheriting from parent (but files are closed)
                result = subprocess.run(
                    command,
                    cwd=self.workDirectory,
                    input=stdin_input,
                    stdout=stdout_file,
                    stderr=stderr_file,
                    text=True,
                    timeout=300,  # 5 minute timeout
                    env=env
                )

            # Store exit code for PROCESSMANAGER queries
            self._exitCode = result.returncode
            self._exitStatus = 0 if result.returncode == 0 else 1

            # Check return code
            if result.returncode != 0:
                error.append(
                    klass=self.__class__.__name__,
                    code=101,
                    details=f"Process {self.TASKCOMMAND} exited with code {result.returncode}",
                    severity=4  # ERROR
                )

                # Read stderr for error details
                if os.path.exists(stderr_path):
                    with open(stderr_path, 'r') as f:
                        stderr_content = f.read()
                        if stderr_content:
                            print(f"STDERR:\n{stderr_content}")

            else:
                print(f"✅ Process completed successfully (exit code 0)")

        except FileNotFoundError as e:
            error.append(
                klass=self.__class__.__name__,
                code=102,
                details=f"File not found error: {str(e)}. "
                        f"Command: {command[0]}. "
                        f"Working directory: {self.workDirectory}. "
                        f"Make sure CCP4 is set up (source ccp4.setup-sh)"
            )
        except subprocess.TimeoutExpired:
            error.append(
                klass=self.__class__.__name__,
                code=103,
                details=f"Process {self.TASKCOMMAND} timed out after 300 seconds"
            )
        except Exception as e:
            error.append(
                klass=self.__class__.__name__,
                code=104,
                details=f"Error running {self.TASKCOMMAND}: {str(e)}"
            )

        return error

    def _startProcessWithStreaming(self, lineHandler=None) -> CErrorReport:
        """
        Process execution with streaming stdout capture.

        Uses subprocess.Popen() to capture stdout line-by-line, enabling
        real-time processing of program output without file-watching race conditions.

        Args:
            lineHandler: Optional callback function(line: str) called for each line.
                        If None, uses self.onProcessOutput if defined.

        The output is:
        1. Passed to lineHandler for real-time processing (e.g., log scraping)
        2. Written to the LOG file for archival
        3. Available for post-processing after completion

        Example usage in a plugin:
            def onProcessOutput(self, line):
                self.logScraper.processLine(line)

            def startProcess(self):
                return self._startProcessWithStreaming(lineHandler=self.onProcessOutput)
        """
        import subprocess
        import os

        error = CErrorReport()

        # Register with PROCESSMANAGER
        from ccp4i2.core.CCP4Modules import PROCESSMANAGER
        PROCESSMANAGER().register(self)

        # Prepare execution environment
        prep = self._prepareProcessExecution()
        command = prep['command']
        stdout_path = prep['stdout_path']
        stderr_path = prep['stderr_path']
        stdin_input = prep['stdin_input']
        env = prep['env']

        # Determine handler - use provided callback or instance method
        handler = lineHandler
        if handler is None and hasattr(self, 'onProcessOutput'):
            handler = self.onProcessOutput

        try:
            with open(stdout_path, 'w') as stdout_file, open(stderr_path, 'w') as stderr_file:
                # Write formatted header
                stdout_file.write("="*70 + "\n")
                stdout_file.write(f"CCP4i2 Task: {self.TASKNAME}\n")
                stdout_file.write("="*70 + "\n\n")
                stdout_file.write("Command Line:\n")
                stdout_file.write("-" * 70 + "\n")
                stdout_file.write(f"{' '.join(command)}\n\n")

                if self.commandScript:
                    stdout_file.write("Command Script (stdin):\n")
                    stdout_file.write("-" * 70 + "\n")
                    for line in self.commandScript:
                        stdout_file.write(line)
                    stdout_file.write("\n")

                stdout_file.write("="*70 + "\n")
                stdout_file.write("Program Output:\n")
                stdout_file.write("="*70 + "\n\n")
                stdout_file.flush()

                # Start process with stdout captured via PIPE
                process = subprocess.Popen(
                    command,
                    cwd=self.workDirectory,
                    stdin=subprocess.PIPE if stdin_input else None,
                    stdout=subprocess.PIPE,
                    stderr=stderr_file,
                    text=True,
                    bufsize=1,  # Line buffered
                    env=env
                )

                # Write stdin if needed, then close
                if stdin_input:
                    process.stdin.write(stdin_input)
                    process.stdin.close()

                # Read stdout line by line
                for line in process.stdout:
                    # Write to log file immediately
                    stdout_file.write(line)
                    stdout_file.flush()

                    # Call handler for real-time processing
                    if handler:
                        try:
                            handler(line.rstrip('\n'))
                        except Exception as e:
                            logger.warning(f"Error in lineHandler: {e}")

                # Wait for process to complete
                process.wait()

                # Store exit code
                self._exitCode = process.returncode
                self._exitStatus = 0 if process.returncode == 0 else 1

                if process.returncode != 0:
                    error.append(
                        klass=self.__class__.__name__,
                        code=101,
                        details=f"Process {self.TASKCOMMAND} exited with code {process.returncode}",
                        severity=4  # ERROR
                    )

        except subprocess.TimeoutExpired:
            error.append(
                klass=self.__class__.__name__,
                code=102,
                details=f"Process {self.TASKCOMMAND} timed out"
            )
        except FileNotFoundError:
            error.append(
                klass=self.__class__.__name__,
                code=103,
                details=f"Executable {self.TASKCOMMAND} not found"
            )
        except Exception as e:
            logger.exception(f"Error running {self.TASKCOMMAND}: {e}")
            error.append(
                klass=self.__class__.__name__,
                code=104,
                details=f"Error running {self.TASKCOMMAND}: {str(e)}"
            )

        return error

    def _startProcessAsync(self) -> CErrorReport:
        """
        Asynchronous process execution using AsyncProcessManager.

        This method:
        1. Starts the subprocess in non-blocking mode
        2. Returns immediately (RUNNING status)
        3. When subprocess finishes, calls _onProcessFinished handler
        4. Handler calls postProcess() which emits finished signal
        """
        import os
        from ccp4i2.core.async_process_manager import ASYNC_PROCESSMANAGER

        error = CErrorReport()

        # Prepare execution environment (common setup)
        prep = self._prepareProcessExecution()
        command = prep['command']
        command_script_file = prep['command_script_file']
        stdout_path = prep['stdout_path']
        stderr_path = prep['stderr_path']
        env = prep['env']

        try:
            # Get async process manager
            pm = ASYNC_PROCESSMANAGER()

            # Prepare stdin input
            input_file = command_script_file if command_script_file else None

            # Create handler callback
            handler = [self._onProcessFinished, {}]

            # Start async process
            # Note: command[0] is the executable path, command[1:] are the arguments
            self._runningProcessId = pm.startProcess(
                command=command[0],
                args=command[1:],
                inputFile=input_file,
                logFile=stdout_path,
                cwd=str(self.workDirectory),
                env=env,
                handler=handler,
                timeout=-1,  # No timeout
                ifAsync=True
            )

            print(f"✅ Process started asynchronously (PID: {self._runningProcessId})")

        except Exception as e:
            error.append(
                klass=self.__class__.__name__,
                code=110,
                details=f"Error starting async process: {e}"
            )

        return error

    def _onProcessFinished(self, pid: int):
        """
        Called when async subprocess completes.

        This handler:
        1. Checks process exit status
        2. Calls postProcess() to handle completion
        3. postProcess() eventually calls reportStatus() which emits finished signal
        """
        print(f"\n{'='*60}")
        print(f"Process {pid} finished")
        print(f"{'='*60}\n")

        # Call postProcess to handle completion
        # This will call processOutputFiles() and reportStatus()
        # reportStatus() emits the finished signal
        status = self.postProcess()

        return status

    def postProcess(self) -> int:
        """
        Post-processing after program completes.

        This method calls:
        1. postProcessCheck() - check if program succeeded
        2. processOutputFiles() - extract output data
        3. reportStatus() - save params and report to database

        Returns:
            Status code (SUCCEEDED or FAILED)
        """
        # Check if process succeeded
        status, exit_status, exit_code = self.postProcessCheck()

        if status == self.SUCCEEDED:
            # Extract output data
            # Wrap in try/except to handle legacy wrappers that may have incomplete implementations
            try:
                error = self.processOutputFiles()
                if error:
                    self.errorReport.extend(error)
                    status = self.FAILED
            except Exception as e:
                # Legacy wrappers may depend on methods we haven't implemented yet
                # For now, just log a warning and continue
                print(f"Warning: processOutputFiles() exception: {type(e).__name__}: {str(e)}")
                import traceback
                traceback.print_exc()
                # Don't fail the job for this - the main output file should still be created

        # Report status and save params
        self.reportStatus(status)

        return status

    def postProcessCheck(self, processId=None):
        """
        Check if the program process completed successfully.

        Args:
            processId: Optional process ID (for legacy compatibility, not currently used
                      as we check self._runningProcessId internally)

        Checks:
        1. Process exit code (for both sync and async processes)
        2. Presence of expected output files
        3. Log file contents for errors

        Returns:
            tuple: (status, exitStatus, exitCode) where:
                   status: SUCCEEDED or FAILED
                   exitStatus: 0 for success, 1 for failure
                   exitCode: The actual process exit code
        """
        import os
        import logging
        logger = logging.getLogger(__name__)

        exit_code = None
        exit_status = None
        status = self.SUCCEEDED

        # Check exit code from synchronous execution (subprocess.run)
        if hasattr(self, '_exitCode') and self._exitCode is not None:
            logger.debug(f"[postProcessCheck] Checking sync exit code: {self._exitCode}")
            exit_code = self._exitCode
            exit_status = self._exitStatus if hasattr(self, '_exitStatus') else (0 if exit_code == 0 else 1)
            if exit_code != 0:
                # Process failed - error already added by _startProcessSync
                logger.debug(f"[postProcessCheck] Returning FAILED due to non-zero exit code: {exit_code}")
                status = self.FAILED

        # For async processes, check the exit code from the process manager
        elif hasattr(self, '_runningProcessId') and self._runningProcessId is not None:
            from ccp4i2.core.async_process_manager import ASYNC_PROCESSMANAGER
            pm = ASYNC_PROCESSMANAGER()

            exit_code = pm.getJobData(self._runningProcessId, 'exitCode')
            exit_status = pm.getJobData(self._runningProcessId, 'exitStatus')

            if exit_code is not None and exit_code != 0:
                # Process failed - add to error report
                self.errorReport.append(
                    klass=self.__class__.__name__,
                    code=993,
                    details=f'Program {self.TASKCOMMAND} exited with error code {exit_code}. '
                           f'Check {self.makeFileName("LOG")} and {self.makeFileName("STDERR")} for details.',
                    name='postProcessCheck',
                    severity=4  # ERROR
                )
                status = self.FAILED

            elif exit_status is not None and exit_status != 0:
                # Process marked as failed by process manager
                self.errorReport.append(
                    klass=self.__class__.__name__,
                    code=992,
                    details=f'Program {self.TASKCOMMAND} failed during execution. '
                           f'Check {self.makeFileName("LOG")} and {self.makeFileName("STDERR")} for details.',
                    name='postProcessCheck',
                    severity=4  # ERROR
                )
                status = self.FAILED

        # Always return tuple for consistency
        if exit_code is None:
            exit_code = 0
        if exit_status is None:
            exit_status = 0 if status == self.SUCCEEDED else 1

        return status, exit_status, exit_code

    def processOutputFiles(self) -> CErrorReport:
        """
        Extract data from output files and populate outputData container.

        This is a hook for subclasses to read output files, extract
        relevant data, and update the outputData container.

        Returns:
            CErrorReport with any errors
        """
        return CErrorReport()

    def reportStatus(self, status: int):
        """
        Report job completion and save parameters.

        This method:
        1. Saves parameters to PARAMS file
        2. Reports completion to database
        3. Emits finished signal

        Args:
            status: Final status (SUCCEEDED or FAILED)
        """
        import logging
        logger = logging.getLogger(__name__)

        logger.debug(f"[reportStatus] Called with status: {status} (SUCCEEDED={self.SUCCEEDED}, FAILED={self.FAILED})")

        # Save params
        self.saveParams()

        # Report final status to database
        if hasattr(self, '_dbHandler') and self._dbHandler is not None:
            logger.debug(f"[reportStatus] _dbHandler exists: {self._dbHandler}")
            if hasattr(self, '_dbJobId') and self._dbJobId is not None:
                logger.debug(f"[reportStatus] _dbJobId exists: {self._dbJobId}")
                try:
                    logger.debug(f"[reportStatus] Calling updateJobStatus with jobId={self._dbJobId}, finishStatus={status}")
                    self._dbHandler.updateJobStatus(
                        jobId=str(self._dbJobId),
                        finishStatus=status,
                        container=self.container
                    )
                    logger.debug(f"[reportStatus] Successfully updated job status in database")
                except Exception as e:
                    logger.warning(f"Failed to update job status in database: {e}")
                    import traceback
                    traceback.print_exc()
            else:
                logger.debug(f"[reportStatus] _dbJobId not set (hasattr={hasattr(self, '_dbJobId')}, value={getattr(self, '_dbJobId', None)})")
        else:
            logger.debug(f"[reportStatus] _dbHandler not set (hasattr={hasattr(self, '_dbHandler')}, value={getattr(self, '_dbHandler', None)})")

        # Emit finished signal with status dict (modern API)
        # Legacy plugins may expect just int, handled by connectSignal() wrapper
        status_dict = {
            'finishStatus': status,
            'jobId': getattr(self, 'jobId', None)
        }
        self.finished.emit(status_dict)

        # Set the internal status
        self._status = status

    def postProcessWrapper(self, finishStatus):
        """
        Wrapper method for propagating finish status from sub-plugins.

        This is called by pipelines as a callback when a sub-plugin finishes.
        It simply forwards the status to reportStatus() to emit the finished signal.

        Args:
            finishStatus: Can be either int (legacy) or dict (modern)
        """
        # Handle both int and dict status formats
        if isinstance(finishStatus, dict):
            status = finishStatus.get('finishStatus', self.FAILED)
        else:
            status = finishStatus

        self.reportStatus(status)

    # =========================================================================
    # File watching methods (legacy API compatibility)
    # =========================================================================

    def watchedFiles(self):
        """Get dictionary of watched files (legacy API).

        Returns:
            Dictionary mapping file paths to watch metadata
        """
        if not hasattr(self, "_watchedFiles"):
            self._watchedFiles = {}
        return self._watchedFiles

    def watchedDirectories(self):
        """Get dictionary of watched directories (legacy API).

        Returns:
            Dictionary mapping directory paths to watch metadata
        """
        if not hasattr(self, "_watchedDirectories"):
            self._watchedDirectories = {}
        return self._watchedDirectories

    def watchFile(self, fileName, handler, minDeltaSize=0, unwatchWhileHandling=False):
        """Watch a file for changes using background thread polling.

        In the legacy Qt-based implementation, this used QFileSystemWatcher.
        This implementation uses a background thread to poll the file for changes.

        Args:
            fileName: Path to file to watch
            handler: Callback function to call when file changes
            minDeltaSize: Minimum file size change to trigger handler
            unwatchWhileHandling: Whether to temporarily unwatch during handling
        """
        import os
        import threading
        import time

        parentDirectoryPath, fileRoot = os.path.split(fileName)

        # Store watch metadata
        watch_info = {
            "parentDirectoryPath": parentDirectoryPath,
            "handler": handler,
            "maxSizeYet": 0,
            "minDeltaSize": minDeltaSize,
            "unwatchWhileHandling": unwatchWhileHandling,
            "active": True,
            "handling": False
        }
        self.watchedFiles()[fileName] = watch_info

        def poll_file():
            """Background thread to poll file for changes."""
            poll_interval = 2.0  # Check every 2 seconds (reduced from 1s to lower CPU usage)
            last_size = 0

            while watch_info["active"]:
                try:
                    if os.path.exists(fileName):
                        current_size = os.path.getsize(fileName)
                        size_delta = current_size - last_size

                        # Check if file has grown enough to trigger handler
                        if size_delta >= minDeltaSize and current_size > watch_info["maxSizeYet"]:
                            # Skip if we're already handling and unwatchWhileHandling is True
                            if unwatchWhileHandling and watch_info["handling"]:
                                time.sleep(poll_interval)
                                continue

                            watch_info["handling"] = True
                            watch_info["maxSizeYet"] = current_size

                            try:
                                logger.debug(f"File changed: {fileName} (size: {current_size})")
                                handler(fileName)
                            except Exception as e:
                                logger.error(f"Error in file watch handler for {fileName}: {e}")
                            finally:
                                watch_info["handling"] = False

                        last_size = current_size

                except Exception as e:
                    logger.debug(f"Error polling file {fileName}: {e}")

                time.sleep(poll_interval)

        # Start background polling thread
        watcher_thread = threading.Thread(target=poll_file, daemon=True, name=f"FileWatcher-{fileRoot}")
        watcher_thread.start()
        logger.debug(f"Started file watcher for {fileName}")

    def unwatchFile(self, fileName):
        """Stop watching a file for changes.

        Args:
            fileName: Path to file to stop watching
        """
        if fileName in self.watchedFiles():
            self.watchedFiles()[fileName]["active"] = False
            del self.watchedFiles()[fileName]
            logger.debug(f"Stopped file watcher for {fileName}")

    def stopAllFileWatchers(self):
        """Stop all file watchers.

        Called during cleanup to ensure all background polling threads terminate.
        """
        for fileName, watch_info in list(self.watchedFiles().items()):
            watch_info["active"] = False
        self._watchedFiles = {}
        logger.debug("Stopped all file watchers")

    def watchDirectory(self, directoryName, handler):
        """Watch a directory for changes (legacy API compatibility).

        This is a no-op stub for backward compatibility.

        Args:
            directoryName: Path to directory to watch
            handler: Callback function to call when directory changes
        """
        self.watchedDirectories()[directoryName] = {"handler": handler}
        logger.debug(f"watchDirectory called for {directoryName} (compatibility mode - no actual watching)")

    # =========================================================================
    # Dictionary merging for crystallographic refinement
    # =========================================================================

    def joinDicts(self, outfile, infiles):
        """Merge multiple crystallographic dictionary files into one.

        This method combines multiple CIF dictionary files (geometric restraint
        definitions for monomers) into a single dictionary file. It first tries
        to use dictionaryAccumulator for proper merging, then falls back to
        libcheck if that fails.

        Args:
            outfile: CDataFile object for the output merged dictionary
            infiles: List of CDataFile objects for input dictionaries

        Returns:
            tuple: (status_code, error_code) where status_code is SUCCEEDED/FAILED
        """
        import os

        # Handle empty input
        if len(infiles) == 0:
            return self.SUCCEEDED, None

        # Handle single input - just copy reference
        if len(infiles) == 1:
            outfile.set(infiles[0])
            return self.SUCCEEDED, None

        # Multiple inputs - need to merge
        outfile_string = str(outfile.fullPath) if hasattr(outfile, 'fullPath') else str(outfile.getFullPath())

        try:
            # Build list of input file paths
            input_cif_list = []
            for dict_file in infiles:
                path = str(dict_file.fullPath) if hasattr(dict_file, 'fullPath') else str(dict_file.getFullPath())
                input_cif_list.append(path)

            # Try to use dictionaryAccumulator (preferred method)
            try:
                from ccp4i2.utils import dictionaryAccumulator
                dictionaryAccumulator.accumulate(input_cif_list, outfile_string)

                if os.path.exists(outfile_string):
                    outfile.setFullPath(outfile_string)
                    logger.info(f"Merged {len(input_cif_list)} dictionaries using dictionaryAccumulator")
                    return self.SUCCEEDED, None
            except Exception as e:
                logger.warning(f"dictionaryAccumulator failed: {e}, trying libcheck...")

            # Fall back to libcheck
            logger.warning("Using libcheck to merge dictionaries (output may not be self-consistent)")

            # Find libcheck executable
            exe = shutil.which('libcheck')

            # Try to find libcheck in CCP4 installation
            # For now, just use subprocess to call libcheck if it's in PATH
            import subprocess

            # Merge dictionaries pairwise
            outfile_string = ""
            for idx, dict_file in enumerate(infiles):
                dict_path = str(dict_file.fullPath) if hasattr(dict_file, 'fullPath') else str(dict_file.getFullPath())

                if len(outfile_string) == 0:
                    # First file - just use it as starting point
                    outfile_string = dict_path
                else:
                    # Merge with previous result
                    outfile_old = outfile_string
                    outfile_string = os.path.join(self.workDirectory, f'merged_dictionary_tmp{idx}')

                    logger.info(f"Merging dictionaries:")
                    logger.info(f"  DictIn1: {outfile_old}")
                    logger.info(f"  DictIn2: {dict_path}")
                    logger.info(f"  DictOut: {outfile_string}")

                    # Build libcheck command file
                    comFileText = f"_N\n_FILE_L {outfile_old}\n_FILE_L2 {dict_path}\n_FILE_O {outfile_string}\n_END\n"

                    # Run libcheck
                    try:
                        result = subprocess.run(
                            [exe],
                            input=comFileText,
                            capture_output=True,
                            text=True,
                            cwd=str(self.workDirectory)
                        )

                        # libcheck creates .lib extension
                        outfile_string = os.path.join(self.workDirectory, f'merged_dictionary_tmp{idx}.lib')

                        if result.returncode not in [0] or not os.path.exists(outfile_string):
                            logger.error(f"libcheck failed with return code {result.returncode}")
                            return self.FAILED, 32

                    except FileNotFoundError:
                        logger.error("libcheck executable not found")
                        return self.FAILED, 32

            # Copy final merged dictionary to output location
            outfile_final = os.path.join(self.workDirectory, 'merged_dictionary.lib')
            if os.path.exists(outfile_string):
                shutil.copyfile(outfile_string, outfile_final)
            else:
                return self.FAILED, 32

            if os.path.exists(outfile_final):
                outfile.setFullPath(outfile_final)
                logger.info(f"Final merged dictionary: {outfile_final}")
                return self.SUCCEEDED, None
            else:
                return self.FAILED, 32

        except Exception as e:
            logger.error(f"Failed to merge dictionaries: {e}")
            return self.FAILED, 28

    # =========================================================================
    # Utility methods for backward compatibility with old API
    # =========================================================================

    def makePluginObject(self, taskName: str = None,
                         reportToDatabase: bool = True, **kwargs) -> Optional['CPluginScript']:
        """
        Create a sub-plugin (sub-job) instance.

        This is used when a pipeline calls other wrappers as sub-jobs.

        Following legacy CCP4i2 convention, each sub-plugin is assigned:
        - A working directory: parent_workdir/job_N (where N = 1, 2, 3, ...)
        - A unique name: parent_name_N

        The working directory is created if it doesn't exist.

        Args:
            taskName: Name of the task to instantiate (or use legacy pluginName= kwarg)
            reportToDatabase: Whether to report this job to the database (default True).
                            In database-backed environments (CCP4i2 GUI), this controls
                            whether the sub-job is registered in the project database.
                            In standalone mode, this parameter is ignored.
            **kwargs: Additional arguments passed to the plugin constructor
                     (workDirectory will be overridden based on convention)
                     Legacy: also accepts 'pluginName' as alias for taskName

        Returns:
            New CPluginScript instance, or None if plugin not found
        """
        import os
        from pathlib import Path

        # Legacy compatibility: accept pluginName= as alias for taskName=
        if taskName is None and 'pluginName' in kwargs:
            taskName = kwargs.pop('pluginName')
            logger.debug(f"[makePluginObject] Legacy pluginName parameter used: {taskName}")

        if taskName is None:
            raise ValueError("makePluginObject requires taskName parameter (or legacy pluginName)")

        # Increment child job counter
        self._childJobCounter += 1

        # Create subdirectory following convention: job_1, job_2, etc.
        child_work_dir = Path(self.workDirectory) / f"job_{self._childJobCounter}"
        try:
            if not child_work_dir.exists():
                child_work_dir.mkdir(parents=True, exist_ok=True)
                logger.debug(f"[DEBUG makePluginObject] Created subdirectory: {child_work_dir}")
        except Exception as e:
            self.errorReport.append(
                klass=self.__class__.__name__,
                code=24,
                details=f"Failed to create working directory '{child_work_dir}': {e}",
                name=str(child_work_dir)
            )
            # Fall back to parent's work directory
            child_work_dir = Path(self.workDirectory)

        # Create name following convention: parent_name_N
        child_name = f"{self.objectName()}_{self._childJobCounter}" if self.objectName() else f"job_{self._childJobCounter}"

        # Check if dummy mode requested (extract from kwargs)
        dummy_mode = kwargs.get('dummy', False)

        # Get the plugin class
        plugin_class = get_plugin_class(taskName)

        if plugin_class is None:
            # If dummy=True, create a generic CPluginScript instead of failing
            if dummy_mode:
                logger.debug(f"[makePluginObject] Plugin '{taskName}' not found, creating dummy CPluginScript")
                plugin_class = CPluginScript
            else:
                # Raise CException - this makes the error immediately visible
                # rather than returning None and causing NoneType errors downstream
                from ccp4i2.core.base_object.error_reporting import CException
                error_msg = f"Plugin '{taskName}' not found in registry. Check that the plugin exists and all its dependencies can be imported."
                self.errorReport.append(
                    klass=self.__class__.__name__,
                    code=108,
                    details=error_msg,
                    name=taskName
                )
                raise CException(error_msg)

        # Instantiate the plugin with computed workDirectory and name
        # Use automatic workDirectory/name unless explicitly provided in kwargs
        try:
            plugin_kwargs = kwargs.copy()

            # Extract pluginTitle (handled separately, not passed to __init__)
            plugin_title = plugin_kwargs.pop('pluginTitle', None)

            # Use automatic workDirectory unless explicitly overridden
            if 'workDirectory' not in plugin_kwargs:
                plugin_kwargs['workDirectory'] = str(child_work_dir)

            # Use automatic name unless explicitly overridden
            if 'name' not in plugin_kwargs:
                plugin_kwargs['name'] = child_name

            # Always set parent relationship
            plugin_kwargs['parent'] = self

            actual_name = plugin_kwargs['name']
            actual_workdir = plugin_kwargs['workDirectory']

            # Save params.xml BEFORE creating the nested plugin
            # This ensures the nested job can load the parent's current parameters
            if self.get_db_job_id():
                try:
                    params_path = os.path.join(self.workDirectory, "params.xml")
                    logger.debug(f"[DEBUG makePluginObject] Saving parent params.xml to {params_path}")
                    save_error = self.saveDataToXml(params_path)
                    if save_error and hasattr(save_error, 'hasError') and save_error.hasError():
                        logger.debug(f"[DEBUG makePluginObject] Warning: Failed to save params.xml: {save_error}")
                    else:
                        logger.debug(f"[DEBUG makePluginObject] ✅ Saved params.xml for nested job to load")
                except Exception as e:
                    logger.debug(f"[DEBUG makePluginObject] Warning: Exception saving params.xml: {e}")

            plugin_instance = plugin_class(**plugin_kwargs)
            # Ensure child job counter is reset for new instance (prevents state pollution)
            plugin_instance._childJobCounter = 0

            # Set pluginTitle on container header if provided (legacy CCP4i2 API)
            if plugin_title is not None and hasattr(plugin_instance, 'container') and plugin_instance.container is not None:
                if hasattr(plugin_instance.container, 'header'):
                    plugin_instance.container.header.pluginTitle = plugin_title
                    logger.debug(f"[DEBUG makePluginObject] Set pluginTitle to '{plugin_title}'")

            # Propagate database context to nested plugin so it can resolve file paths
            # Nested plugins need the dbHandler to lookup files via dbFileId
            if hasattr(self, '_dbHandler') and self._dbHandler is not None:
                plugin_instance._dbHandler = self._dbHandler
                logger.debug(f"Propagated dbHandler to nested plugin")
            if hasattr(self, '_dbProjectId') and self._dbProjectId is not None:
                plugin_instance._dbProjectId = self._dbProjectId
                logger.debug(f"[DEBUG makePluginObject] Propagated dbProjectId to nested plugin")

            # Handle database job creation for sub-job
            # In database-backed mode with reportToDatabase=True, delegate to dbHandler
            if reportToDatabase and hasattr(self, '_dbHandler') and self._dbHandler is not None and hasattr(self._dbHandler, 'createSubJob'):
                try:
                    parent_job_id = self._dbJobId if hasattr(self, '_dbJobId') else None
                    if parent_job_id:
                        job_number = str(self._childJobCounter)
                        # Delegate to dbHandler for clean separation of concerns
                        new_job_id = self._dbHandler.createSubJob(
                            taskName=taskName,
                            parentJobId=parent_job_id,
                            jobNumber=job_number
                        )
                        plugin_instance._dbJobId = new_job_id
                        logger.debug(f"[DEBUG makePluginObject] Assigned new job ID to sub-job: {new_job_id}")
                    else:
                        logger.debug(f"[DEBUG makePluginObject] No parent job ID, skipping sub-job creation")
                except Exception as e:
                    logger.warning(f"[WARNING makePluginObject] Failed to create database job for sub-job: {e}")
                    # Fall back to propagating parent's job ID (old behavior)
                    if hasattr(self, '_dbJobId') and self._dbJobId is not None:
                        plugin_instance._dbJobId = self._dbJobId
                        logger.debug(f"[DEBUG makePluginObject] Falling back to parent's job ID due to error")
            else:
                # Not in database mode or reportToDatabase is False
                # Still propagate parent's job ID for backward compatibility
                if hasattr(self, '_dbJobId') and self._dbJobId is not None:
                    plugin_instance._dbJobId = self._dbJobId
                    logger.debug(f"[DEBUG makePluginObject] Propagated parent's dbJobId to nested plugin (no database creation)")
            logger.debug(f"[DEBUG makePluginObject] Created sub-plugin '{taskName}' as '{actual_name}' in '{actual_workdir}'")
            return plugin_instance
        except Exception as e:
            self.errorReport.append(
                klass=self.__class__.__name__,
                code=109,
                details=f"Failed to instantiate plugin '{taskName}': {e}",
                name=taskName
            )
            import traceback
            traceback.print_exc()
            return None

    def getErrorReport(self) -> CErrorReport:
        """Get the accumulated error report."""
        return self.errorReport

    def getContainer(self) -> CContainer:
        """Get the main container."""
        return self.container

    # =========================================================================
    # Database integration methods (for database-backed environments)
    # =========================================================================

    def getJobId(self):
        """Get the database job ID.

        Returns:
            Job ID in database, or None if not running in database-backed mode
        """
        return self._dbJobId

    @property
    def jobId(self):
        """Legacy property for backward compatibility.

        Returns:
            Job ID in database (same as getJobId())
        """
        return self._dbJobId

    @jobId.setter
    def jobId(self, value):
        """Legacy setter for backward compatibility.

        Args:
            value: Job ID to set
        """
        self._dbJobId = value

    @property
    def jobTitle(self):
        """Get job title from guiAdmin container.

        Returns:
            Job title string, or empty string if not set
        """
        try:
            return str(self.container.guiAdmin.jobTitle) if self.container.guiAdmin.jobTitle.isSet() else ''
        except AttributeError:
            return ''

    @jobTitle.setter
    def jobTitle(self, value):
        """Set job title in guiAdmin container.

        Args:
            value: Job title string to set
        """
        try:
            self.container.guiAdmin.jobTitle.value = str(value) if value else ''
        except AttributeError:
            logger.warning("Cannot set jobTitle: guiAdmin container not available")

    @property
    def jobStatus(self):
        """Get job status from guiAdmin container.

        Returns:
            Job status integer (0=Pending, 1=Running, 2=Finished, 3=Failed, etc.)
        """
        try:
            return self.container.guiAdmin.jobStatus.value if self.container.guiAdmin.jobStatus.isSet() else 0
        except AttributeError:
            return 0

    @jobStatus.setter
    def jobStatus(self, value):
        """Set job status in guiAdmin container.

        Args:
            value: Job status integer to set
        """
        try:
            self.container.guiAdmin.jobStatus.value = int(value)
        except AttributeError:
            logger.warning("Cannot set jobStatus: guiAdmin container not available")

    def getJobNumber(self):
        """Get the database job number.

        Returns:
            Job number (e.g., "1.2.3"), or None if not running in database-backed mode
        """
        return self._dbJobNumber

    def getProjectId(self):
        """Get the database project ID.

        Returns:
            Project ID in database, or None if not running in database-backed mode
        """
        return self._dbProjectId

    def projectId(self):
        """Legacy alias for getProjectId() (compatibility with old CCP4i2 code).

        Returns:
            Project ID in database, or None if not running in database-backed mode
        """
        return self._dbProjectId

    def relPath(self, jobNumber=None):
        """Get relative path to job directory from project root (legacy CCP4i2 API).

        Args:
            jobNumber: Job number string (e.g., "1.2.3"). If None, uses self._dbJobNumber.

        Returns:
            Relative path string like "CCP4_JOBS/job_1/job_2/job_3"
        """
        import os
        if jobNumber is None:
            jobNumber = self._dbJobNumber
        if jobNumber is None:
            # Fallback: extract from workDirectory if it contains CCP4_JOBS
            if 'CCP4_JOBS' in str(self.workDirectory):
                work_dir = str(self.workDirectory)
                return work_dir[work_dir.index('CCP4_JOBS'):]
            return None
        numList = str(jobNumber).split('.')
        path = os.path.join('CCP4_JOBS', 'job_' + numList[0])
        for num in numList[1:]:
            path = os.path.join(path, 'job_' + num)
        return path

    def getChildJobNumber(self):
        """Get the current child job counter.

        This returns how many sub-plugins have been created via makePluginObject().

        Returns:
            Child job counter (starts at 0, increments with each makePluginObject call)
        """
        return self._childJobCounter

    def getWorkDirectory(self):
        "Get the absolute working directory path as string."
        return str(self.workDirectory)

    def testForInterrupt(self) -> bool:
        """Test if user has requested pipeline interruption.

        Checks for existence of 'INTERRUPT' file in the working directory.
        This is used by pipelines (e.g., crank2) to detect user cancellation.

        Returns:
            True if INTERRUPT file exists, False otherwise
        """
        return (self.workDirectory / "INTERRUPT").exists()

    # connectSignal() is now inherited from HierarchicalObject base class
    # with automatic signature adaptation for legacy int handlers

    # =========================================================================
    # MTZ File Merging Methods (makeHklin family)
    # =========================================================================

    def makeHklinGemmi(
        self,
        file_objects: list,
        output_name: str = 'hklin',
        merge_strategy: str = 'first'
    ) -> Path:
        """
        Merge normalized mini-MTZ files into a single HKLIN file (new Pythonic API).

        This is the modern, Pythonic replacement for makeHklin. It works with
        container attribute names and uses gemmi for MTZ operations.

        Args:
            file_objects: List of file specifications. Each can be either:
                - str: Attribute name in inputData/outputData (e.g., 'HKLIN1')
                       Uses the file's contentFlag to determine columns automatically.
                - dict: Explicit specification with keys:
                    {
                        'name': str,                    # Attribute name (required)
                        'target_contentFlag': int,      # Optional: Convert to this contentFlag if needed
                        'rename': Dict[str, str],       # Optional: Column renaming
                        'display_name': str             # Optional: Name for column prefixes (default: name)
                    }

            output_name: Base name for output file (default: 'hklin')
                        Output written to: self.workDirectory / f"{output_name}.mtz"

            merge_strategy: How to handle column conflicts (default: 'first')
                - 'first': Keep column from first file
                - 'last': Keep column from last file
                - 'error': Raise error on conflicts
                - 'rename': Auto-rename conflicts (F, F_1, F_2, ...)

        Returns:
            Path: Full path to created HKLIN file

        Raises:
            AttributeError: If file object not found in inputData/outputData
            ValueError: If contentFlag unknown or file has no path set
            FileNotFoundError: If MTZ file doesn't exist at specified path
            NotImplementedError: If conversion method not available

        Example:
            >>> # Simple merge using contentFlags
            >>> hklin = self.makeHklinGemmi(['HKLIN1', 'FREERFLAG'])
            >>> # Result: merged.mtz with columns from both files

            >>> # With explicit column renaming
            >>> hklin = self.makeHklinGemmi([
            ...     'HKLIN1',  # Uses contentFlag automatically
            ...     {
            ...         'name': 'HKLIN2',
            ...         'rename': {'F': 'F_deriv', 'SIGF': 'SIGF_deriv'}
            ...     }
            ... ])

            >>> # With automatic conversion
            >>> hklin = self.makeHklinGemmi([
            ...     'HKLIN1',
            ...     {
            ...         'name': 'HKLIN2',
            ...         'target_contentFlag': 4  # Convert to FMEAN if not already
            ...     }
            ... ])
        """
        from ccp4i2.core.CCP4Utils import merge_mtz_files
        from ccp4i2.core.base_object.fundamental_types import CInt

        input_specs = []
        converted_files = []  # Track temporary file objects for cleanup

        for file_spec_idx, file_spec in enumerate(file_objects):
            # Parse spec to get name, display_name, target_contentFlag, and optional rename
            if isinstance(file_spec, str):
                name = file_spec
                display_name = file_spec
                target_content_flag = None
                rename_map = {}
            elif isinstance(file_spec, dict):
                name = file_spec['name']
                display_name = file_spec.get('display_name', name)
                target_content_flag = file_spec.get('target_contentFlag', None)
                rename_map = file_spec.get('rename', {})
            else:
                raise ValueError(
                    f"Invalid file_spec type: {type(file_spec)}. "
                    f"Expected str or dict."
                )

            # Lookup file object in containers
            file_obj = None
            if hasattr(self.container.inputData, name):
                file_obj = getattr(self.container.inputData, name)
            elif hasattr(self.container.outputData, name):
                file_obj = getattr(self.container.outputData, name)
            else:
                raise AttributeError(
                    f"File object '{name}' not found in inputData or outputData"
                )

            # Debug logging
            pass  # DEBUG: print(f"[DEBUG makeHklinGemmi] Processing file '{name}':")
            print(f"  - Type: {type(file_obj).__name__}")
            if hasattr(file_obj, 'baseName'):
                val = file_obj.baseName.value if hasattr(file_obj.baseName, 'value') else file_obj.baseName
                print(f"  - baseName.value: '{val}'")
                print(f"  - baseName.isSet(): {file_obj.baseName.isSet() if hasattr(file_obj.baseName, 'isSet') else 'N/A'}")
            if hasattr(file_obj, 'relPath'):
                val = file_obj.relPath.value if hasattr(file_obj.relPath, 'value') else file_obj.relPath
                print(f"  - relPath.value: '{val}'")
            if hasattr(file_obj, 'dbFileId'):
                val = file_obj.dbFileId.value if hasattr(file_obj.dbFileId, 'value') else file_obj.dbFileId
                print(f"  - dbFileId.value: '{val}'")
            if hasattr(file_obj, 'project'):
                val = file_obj.project.value if hasattr(file_obj.project, 'value') else file_obj.project
                print(f"  - project.value: '{val}'")
            print(f"  - getFullPath(): '{file_obj.getFullPath() if hasattr(file_obj, 'getFullPath') else 'N/A'}'")

            pass  # DEBUG: print(f"[DEBUG makeHklinGemmi] About to check if setContentFlag needed for '{name}'")
            pass  # DEBUG: print(f"[DEBUG makeHklinGemmi]   hasattr(file_obj, 'setContentFlag'): {hasattr(file_obj, 'setContentFlag')}")
            pass  # DEBUG: print(f"[DEBUG makeHklinGemmi]   contentFlag: {file_obj.contentFlag if hasattr(file_obj, 'contentFlag') else 'N/A'}")
            pass  # DEBUG: print(f"[DEBUG makeHklinGemmi]   int(contentFlag): {int(file_obj.contentFlag) if hasattr(file_obj, 'contentFlag') else 'N/A'}")

            # Auto-detect contentFlag from file content to ensure accuracy
            # Check if contentFlag is NOT_SET (isSet() returns False) or has value 0
            if hasattr(file_obj, 'setContentFlag'):
                cf = file_obj.contentFlag if hasattr(file_obj, 'contentFlag') else None
                needs_detection = False
                if cf is None:
                    needs_detection = True
                elif hasattr(cf, 'isSet') and not cf.isSet():
                    needs_detection = True
                elif int(file_obj.contentFlag) == 0:
                    needs_detection = True

                if needs_detection:
                    logger.debug(f"Auto-detecting contentFlag for '{name}'")
                    file_obj.setContentFlag()

            # Check if conversion is needed
            current_content_flag = int(file_obj.contentFlag)
            pass  # DEBUG: print(f"[DEBUG makeHklinGemmi] current_content_flag={current_content_flag}, target_content_flag={target_content_flag}")

            if target_content_flag is not None and current_content_flag != target_content_flag:
                # CONVERSION NEEDED!
                logger.debug(f"[DEBUG makeHklinGemmi] Conversion needed: {name} from contentFlag={current_content_flag} to {target_content_flag}")

                # Validate target contentFlag
                if not hasattr(file_obj, 'CONTENT_SIGNATURE_LIST'):
                    raise ValueError(
                        f"File object '{name}' (class {file_obj.__class__.__name__}) "
                        f"has no CONTENT_SIGNATURE_LIST. Cannot convert."
                    )

                if target_content_flag < 1 or target_content_flag > len(file_obj.CONTENT_SIGNATURE_LIST):
                    raise ValueError(
                        f"Invalid target_contentFlag {target_content_flag} for '{name}'. "
                        f"Valid range: 1-{len(file_obj.CONTENT_SIGNATURE_LIST)}"
                    )

                # Get target content flag name (e.g., 'IPAIR', 'FMEAN')
                target_name = self._get_content_flag_name(file_obj, target_content_flag)

                # Call conversion method (e.g., as_IPAIR(), as_FMEAN())
                method_name = f'as_{target_name}'
                if not hasattr(file_obj, method_name):
                    raise NotImplementedError(
                        f"Conversion method '{method_name}' not found on {file_obj.__class__.__name__}. "
                        f"Cannot convert from contentFlag={current_content_flag} to {target_content_flag}."
                    )

                conversion_method = getattr(file_obj, method_name)
                converted_path = conversion_method(self.workDirectory)
                logger.debug(f"[DEBUG makeHklinGemmi] Converted {name} to {converted_path}")

                # Create a temporary file object pointing to converted file
                temp_name = f"_converted_{name}_{file_spec_idx}"
                temp_file_obj = file_obj.__class__(parent=self.container.inputData, name=temp_name)
                temp_file_obj.setFullPath(str(converted_path))
                temp_file_obj.contentFlag = CInt(target_content_flag)

                # Add to inputData temporarily
                setattr(self.container.inputData, temp_name, temp_file_obj)
                converted_files.append(temp_name)

                # Use temp object for merging
                file_obj = temp_file_obj
                logger.debug(f"[DEBUG makeHklinGemmi] Using temp file object '{temp_name}' with contentFlag={target_content_flag}")

            # Get filesystem path
            pass  # DEBUG: print(f"[DEBUG makeHklinGemmi] Getting path for '{name}'...")
            path = file_obj.getFullPath()
            pass  # DEBUG: print(f"[DEBUG makeHklinGemmi] Got path: {path}")
            logger.debug(f"[DEBUG makeHklinGemmi] Processing '{name}' -> path: {path}")
            if not path:
                raise ValueError(f"File object '{name}' has no path set")

            pass  # DEBUG: print(f"[DEBUG makeHklinGemmi] Getting columns for '{name}'...")
            # Get columns from CONTENT_SIGNATURE_LIST using contentFlag
            # contentFlag is 1-indexed, CONTENT_SIGNATURE_LIST is 0-indexed
            pass  # DEBUG: print(f"[DEBUG makeHklinGemmi] Converting contentFlag to int...")
            content_flag = int(file_obj.contentFlag)
            pass  # DEBUG: print(f"[DEBUG makeHklinGemmi] content_flag={content_flag}")
            pass  # DEBUG: print(f"[DEBUG makeHklinGemmi] Checking for CONTENT_SIGNATURE_LIST...")
            if not hasattr(file_obj, 'CONTENT_SIGNATURE_LIST'):
                raise ValueError(
                    f"File object '{name}' (class {file_obj.__class__.__name__}) "
                    f"has no CONTENT_SIGNATURE_LIST. Is it a CMiniMtzDataFile?"
                )
            pass  # DEBUG: print(f"[DEBUG makeHklinGemmi] CONTENT_SIGNATURE_LIST exists, length={len(file_obj.CONTENT_SIGNATURE_LIST)}")

            if content_flag < 1 or content_flag > len(file_obj.CONTENT_SIGNATURE_LIST):
                raise ValueError(
                    f"Invalid contentFlag {content_flag} for '{name}'. "
                    f"Valid range: 1-{len(file_obj.CONTENT_SIGNATURE_LIST)}"
                )

            pass  # DEBUG: print(f"[DEBUG makeHklinGemmi] Getting columns from CONTENT_SIGNATURE_LIST[{content_flag - 1}]...")
            columns = file_obj.CONTENT_SIGNATURE_LIST[content_flag - 1]
            pass  # DEBUG: print(f"[DEBUG makeHklinGemmi] Got columns: {columns}")

            # Build column_mapping (input_label -> output_label)
            # By default, prepend display_name to column (e.g., HKLIN1_F)
            # unless explicit rename is provided or identity mapping is requested
            column_mapping = {}

            # Check if identity mapping is requested (for legacy makeHklin compatibility)
            use_identity = (rename_map == 'identity')

            for col in columns:
                if use_identity:
                    # Identity mapping: F -> F, SIGF -> SIGF (no prefixing)
                    output_col = col
                elif isinstance(rename_map, dict) and col in rename_map:
                    # Explicit rename provided
                    output_col = rename_map[col]
                else:
                    # Default: prepend display_name with underscore
                    output_col = f"{display_name}_{col}"
                column_mapping[col] = output_col

            # Build spec for merge_mtz_files
            input_specs.append({
                'path': path,
                'column_mapping': column_mapping
            })

        # Call low-level gemmi utility
        output_path = self.workDirectory / f"{output_name}.mtz"
        result = merge_mtz_files(
            input_specs=input_specs,
            output_path=output_path,
            merge_strategy=merge_strategy
        )

        return result

    def _get_content_flag_name(self, file_obj, content_flag: int) -> str:
        """
        Get the name of a content flag from its integer value.

        Args:
            file_obj: File object with CONTENT_FLAG_* class constants
            content_flag: Integer content flag value

        Returns:
            Name of the content flag (e.g., 'IPAIR', 'FMEAN')

        Raises:
            ValueError: If content flag not found
        """
        # Search class constants for matching content flag
        for attr_name in dir(file_obj.__class__):
            if attr_name.startswith('CONTENT_FLAG_'):
                flag_value = getattr(file_obj.__class__, attr_name)
                if flag_value == content_flag:
                    # Extract name after CONTENT_FLAG_
                    return attr_name.replace('CONTENT_FLAG_', '')

        raise ValueError(f"No content flag name found for value {content_flag}")

    def makeHklin(self, miniMtzsIn: list, hklin: str = 'hklin', ignoreErrorCodes: list = []) -> tuple:
        """
        Merge mini-MTZ files into HKLIN (backward-compatible legacy API).

        This is a lightweight wrapper around makeHklinGemmi() that translates
        the old API format to the new API format. All conversion logic is
        handled by makeHklinGemmi().

        Args:
            miniMtzsIn: List of file specifications. Each can be either:
                - str: Attribute name in inputData (e.g., 'HKLIN1')
                       Uses the file object's own contentFlag
                - [str, int]: [attribute_name, target_contentFlag]
                       If file's contentFlag != target_contentFlag,
                       converts file to target format first (handled by makeHklinGemmi)

            hklin: Base name for output file (default: 'hklin')
            ignoreErrorCodes: Error codes to ignore (for compatibility, not used)

        Returns:
            tuple: (hklin_filename, CErrorReport) where:
                - hklin_filename: Path to created HKLIN file (None if error)
                - CErrorReport: Error report (empty if successful)

        Example (old API):
            >>> # Simple merge (uses objects' contentFlags)
            >>> error = self.makeHklin(['HKLIN1', 'FREERFLAG'])

            >>> # Request HKLIN2 in FPAIR format (converts if needed)
            >>> error = self.makeHklin([
            ...     'HKLIN1',
            ...     ['HKLIN2', CObsDataFile.CONTENT_FLAG_FPAIR]
            ... ])
        """
        error = CErrorReport()
        hklin_filename = None

        try:
            # Translate old API to new API
            # Legacy makeHklin() should NOT prefix column names - it uses identity mapping
            # (unlike makeHklin0 which DOES prefix with mtzName_columnName)
            file_objects = []

            for item in miniMtzsIn:
                if isinstance(item, str):
                    # Simple name - use identity mapping (no prefixing)
                    file_objects.append({
                        'name': item,
                        'display_name': item,
                        'rename': 'identity'  # Special value: use identity mapping
                    })

                elif isinstance(item, (list, tuple)) and len(item) == 2:
                    # [name, target_contentFlag] - use identity mapping
                    name, target_flag = item
                    file_objects.append({
                        'name': name,
                        'display_name': name,
                        'target_contentFlag': target_flag,
                        'rename': 'identity'  # Special value: use identity mapping
                    })

                else:
                    error.append(
                        klass=self.__class__.__name__,
                        code=207,
                        details=f"Invalid miniMtzsIn item: {item}. Expected str or [str, int]",
                        name=str(item)
                    )
                    self.errorReport.extend(error)
                    return (None, error)

            # Call new API - it handles all conversions
            output_path = self.makeHklinGemmi(
                file_objects=file_objects,
                output_name=hklin,
                merge_strategy='first'
            )

            # Store the output filename for legacy API compatibility
            hklin_filename = str(output_path)

        except (AttributeError, ValueError, NotImplementedError, FileNotFoundError) as e:
            hklin_filename = None
            # Map specific exceptions to error codes
            error.append(
                klass=self.__class__.__name__,
                code=200,
                details=f"Error in makeHklin: {e}",
                name=hklin
            )
            self.errorReport.extend(error)

        except Exception as e:
            # Catch-all for unexpected errors
            error.append(
                klass=self.__class__.__name__,
                code=200,
                details=f"Unexpected error in makeHklin: {e}",
                name=hklin
            )
            self.errorReport.extend(error)

        # Add errors to plugin's error report
        self.errorReport.extend(error)

        # Print ERROR-level messages to terminal
        if error.maxSeverity() >= SEVERITY_ERROR:
            print(f"\n{'='*60}")
            print(f"ERROR in {self.__class__.__name__}.makeHklin():")
            print(f"{'='*60}")
            print(error.report())
            print(f"{'='*60}\n")

        return (hklin_filename, error)

    def makeHklInput(
        self,
        miniMtzsIn: list = [],
        hklin: str = 'hklin',
        ignoreErrorCodes: list = [],
        extendOutputColnames: bool = True,
        useInputColnames: bool = False
    ) -> tuple:
        """
        Legacy API for makeHklin - returns (outfile, colnames, error).

        This method provides backward compatibility with the old CCP4i2 API.
        It wraps the modern makeHklin() method and returns the expected tuple.

        Args:
            miniMtzsIn: List of file names or [name, contentFlag] pairs
            hklin: Output filename (without extension)
            ignoreErrorCodes: Error codes to ignore (not currently used)
            extendOutputColnames: Whether to extend column names with parameter prefix
            useInputColnames: Whether to use original input column names (identity mapping)
                              When True, preserves standard column names like F, SIGF, FreeR_flag

        Returns:
            Tuple of (outfile_path, column_names, error_report)
        """
        # Choose which makeHklin variant to call based on parameters
        #
        # Legacy ccp4i2 behavior from _buildInputVector:
        # - (False, False): [infile, colout]              - standard output names
        # - (True, True):   [infile, colin, ext_outputCol] - input→prefixed mapping (old makeHklin0)
        # - (True, False):  [infile, ext_outputCol]       - prefixed output names
        # - (False, True):  [infile, colin, colout]       - input→standard mapping
        #
        # The key insight from legacy: when extendOutputColnames=True, output has PREFIXED names
        # Both shelxeMR and phaser_singleMR use (True, True) and expect prefixed column names.
        # shelxeMR uses MTZ_parse which finds columns by TYPE, so works with any naming.
        # phaser_singleMR hardcodes prefixed names (F_SIGF_F, etc.)

        if extendOutputColnames:
            # extendOutputColnames=True means output MTZ has prefixed column names
            # This matches legacy makeHklin0 behavior
            # useInputColnames affects the LABIN (input column names) but output is still prefixed
            outfile, colnames, error = self.makeHklin0(miniMtzsIn, hklin, ignoreErrorCodes)
            return outfile, colnames, error

        elif useInputColnames:
            # Only useInputColnames (without extendOutputColnames): identity mapping
            outfile, error = self.makeHklin(miniMtzsIn, hklin)
            # Get column names from the merged MTZ
            colnames = ""
            try:
                import gemmi
                if outfile:
                    mtz = gemmi.read_mtz_file(str(outfile))
                    column_names = [col.label for col in mtz.columns
                                    if col.label not in ['H', 'K', 'L', 'M/ISYM']]
                    colnames = ','.join(column_names)
            except Exception:
                pass
            return outfile, colnames, error

        else:
            # Neither parameter: identity mapping (default, old makeHklin behavior)
            outfile, error = self.makeHklin(miniMtzsIn, hklin)

        # outfile might be None if there was an error
        if outfile is None:
            outfile = str(self.workDirectory / f"{hklin}.mtz")

        # For now, return empty column names string
        # (Full column introspection would require reading the merged MTZ)
        colnames = ""

        return outfile, colnames, error

    def makeHklin0(
        self,
        miniMtzsIn: list = [],
        hklin: str = 'hklin',
        ignoreErrorCodes: list = []
    ) -> tuple:
        """
        Legacy API for makeHklin that returns prefixed column names.

        This is similar to makeHklInput but follows the original makeHklin0 behavior
        from legacy ccp4i2, which prefixes output column names with the parameter name
        (e.g., "F_SIGF_F,F_SIGF_SIGF,ABCD_HLA,ABCD_HLB,...").

        Args:
            miniMtzsIn: List of file names or [name, contentFlag] pairs
            hklin: Output filename (without extension)
            ignoreErrorCodes: Error codes to ignore (for compatibility)

        Returns:
            Tuple of (outfile_path, column_names_string, error_report)
            where column_names_string is comma-separated prefixed column names

        Example:
            >>> # Returns: ('hklin.mtz', 'F_SIGF_F,F_SIGF_SIGF,ABCD_HLA,ABCD_HLB,ABCD_HLC,ABCD_HLD', error)
            >>> self.makeHklin0(['F_SIGF', 'ABCD'])
        """
        import gemmi

        # Initialize error report
        error = self.errorReport.__class__()

        # Build file_objects list with prefixing (no identity mapping)
        # Unlike makeHklin(), we want column names prefixed with parameter name
        file_objects = []

        for item in miniMtzsIn:
            if isinstance(item, str):
                # Simple name - prefix columns with parameter name
                file_objects.append({
                    'name': item,
                    'display_name': item,
                    # Do NOT use 'identity' - let makeHklinGemmi do default prefixing
                })

            elif isinstance(item, (list, tuple)) and len(item) == 2:
                # [name, target_contentFlag] - prefix columns with parameter name
                name, target_flag = item
                file_objects.append({
                    'name': name,
                    'display_name': name,
                    'target_contentFlag': target_flag,
                    # Do NOT use 'identity' - let makeHklinGemmi do default prefixing
                })

            else:
                error.append(
                    klass=self.__class__.__name__,
                    code=207,
                    details=f"Invalid miniMtzsIn item: {item}. Expected str or [str, int]",
                    name=str(item)
                )
                self.errorReport.extend(error)
                return (None, "", error)

        # Call new API with prefixing enabled
        try:
            output_path = self.makeHklinGemmi(
                file_objects=file_objects,
                output_name=hklin,
                merge_strategy='first'
            )
            outfile = str(output_path)

        except (AttributeError, ValueError, NotImplementedError, FileNotFoundError) as e:
            outfile = str(self.workDirectory / f"{hklin}.mtz")
            error.append(
                klass=self.__class__.__name__,
                code=200,
                details=f"Error in makeHklin0: {e}",
                name=hklin
            )
            self.errorReport.extend(error)
            return outfile, "", error

        # Read the merged MTZ to get actual column names
        try:
            mtz = gemmi.read_mtz_file(str(outfile))
            # Get all data columns (exclude H, K, L, M_ISYM which are special)
            column_names = []
            for col in mtz.columns:
                if col.label not in ['H', 'K', 'L', 'M/ISYM']:
                    column_names.append(col.label)

            # Join with commas
            allColout = ','.join(column_names)

            pass  # DEBUG: print(f"[DEBUG makeHklin0] Created {outfile} with columns: {allColout}")

        except Exception as e:
            print(f"[WARNING makeHklin0] Could not read column names from {outfile}: {e}")
            allColout = ""

        return outfile, allColout, error

    def splitMtz(self, infile: str, outfiles: list, logFile: str = None) -> int:
        """
        Split an MTZ file into multiple mini-MTZ files with selected columns.

        This is a thin CData wrapper around split_mtz_file() from CCP4Utils.
        It handles the legacy outfiles list format and converts it to the
        simple column_mapping dict format.

        Args:
            infile: Path to input MTZ file
            outfiles: List of output specifications, each is a list of:
                     [output_path, input_columns] or
                     [output_path, input_columns, output_columns]
                     where input_columns and output_columns are comma-separated column names
            logFile: Optional path to log file (not used in gemmi implementation)

        Returns:
            SUCCEEDED or FAILED status code

        Example:
            >>> self.splitMtz(
            ...     '/path/to/input.mtz',
            ...     [['/path/to/output.mtz', 'FMEAN,SIGFMEAN', 'F,SIGF']],
            ...     '/path/to/log'
            ... )
        """
        from ccp4i2.core.CCP4Utils import split_mtz_file, MtzSplitError

        logger.debug(f'[DEBUG splitMtz] Splitting {infile} using gemmi')
        logger.debug(f'[DEBUG splitMtz] Output specs: {outfiles}')

        try:
            # Process each output file
            for outfile_spec in outfiles:
                # Parse output specification
                if len(outfile_spec) == 2:
                    output_path, input_cols = outfile_spec
                    output_cols = input_cols  # Use same names for output
                elif len(outfile_spec) == 3:
                    output_path, input_cols, output_cols = outfile_spec
                else:
                    print(f'[ERROR] Invalid outfile spec: {outfile_spec}')
                    return self.FAILED

                # Parse column names
                input_col_names = [c.strip() for c in input_cols.split(',')]
                output_col_names = [c.strip() for c in output_cols.split(',')]

                if len(input_col_names) != len(output_col_names):
                    print(f'[ERROR] Input and output column counts must match')
                    return self.FAILED

                # Build column mapping dict for utility function
                column_mapping = dict(zip(input_col_names, output_col_names))

                logger.debug(f'[DEBUG splitMtz] Creating {output_path}')
                logger.debug(f'[DEBUG splitMtz]   Column mapping: {column_mapping}')

                # Call CData-agnostic utility function
                result_path = split_mtz_file(
                    input_path=infile,
                    output_path=output_path,
                    column_mapping=column_mapping
                )

                import os
                file_size = os.path.getsize(result_path)
                logger.debug(f'[DEBUG splitMtz] Created: {result_path} ({file_size} bytes)')

            return self.SUCCEEDED

        except (FileNotFoundError, ValueError, MtzSplitError) as e:
            print(f'[ERROR] splitMtz failed: {e}')
            return self.FAILED
        except Exception as e:
            print(f'[ERROR] splitMtz unexpected error: {e}')
            import traceback
            traceback.print_exc()
            return self.FAILED

    def splitHklout(
        self,
        miniMtzsOut: list,
        programColumnNames: list,
        outputColumnNames: list = None,
        inFile: str = None,
        logFile: str = None,
        **kwargs
    ) -> 'CErrorReport':
        """
        Split an HKLOUT file into multiple mini-MTZ files (CData-aware API).

        This is the standard CCP4i2 API method that works with object names
        in container.outputData. It wraps the lower-level splitMtz() method.

        **Column Name Standardization (Auto-Inference)**:
        This method automatically standardizes column names when outputColumnNames is not provided.
        It inspects each output file object's CONTENT_SIGNATURE_LIST and uses the first signature
        as the standard column names. This ensures generated files match expected column naming
        conventions (e.g., 'FreeR_flag' -> 'FREER', allowing contentFlag introspection to work).

        Args:
            miniMtzsOut: List of output object names in container.outputData
            programColumnNames: List of comma-separated INPUT column name strings
                              e.g., ['F,SIGF', 'HLA,HLB,HLC,HLD', 'FreeR_flag']
            outputColumnNames: Optional list of comma-separated OUTPUT column name strings
                             for explicit relabeling. If provided, must match length of programColumnNames.
                             e.g., ['F,SIGF', 'A,B,C,D', 'FREER']
                             If None (default), column names are AUTO-INFERRED from each output
                             file object's CONTENT_SIGNATURE_LIST[0], enabling automatic standardization.
            inFile: Input MTZ file path (default: workDirectory/hklout.mtz)
            logFile: Log file path (default: workDirectory/splitmtz.log)
            **kwargs: Legacy compatibility - accepts 'infile' (lowercase) as alias for 'inFile'

        Returns:
            CErrorReport with any errors

        Example:
            >>> # Split with AUTOMATIC column standardization (inferred from CONTENT_SIGNATURE_LIST)
            >>> error = self.splitHklout(
            ...     ['FREEROUT'],
            ...     ['FreeR_flag']  # Will be auto-relabeled to 'FREER' based on CFreeRDataFile.CONTENT_SIGNATURE_LIST
            ... )
            >>>
            >>> # Split with EXPLICIT relabeling
            >>> error = self.splitHklout(
            ...     ['ABCDOUT'],
            ...     ['HLA,HLB,HLC,HLD'],
            ...     ['A,B,C,D']  # Explicitly rename columns
            ... )
            >>>
            >>> # Split without relabeling (no CONTENT_SIGNATURE_LIST available)
            >>> error = self.splitHklout(
            ...     ['FPHIOUT'],
            ...     ['F,phi']  # No relabeling if CONTENT_SIGNATURE_LIST not found
            ... )
        """
        error = CErrorReport()

        # Legacy compatibility: accept 'infile' (lowercase) as alias for 'inFile'
        if 'infile' in kwargs and inFile is None:
            inFile = kwargs['infile']

        # IMPORTANT: Handle legacy wrapper bug where MTZ path is passed as outputColumnNames
        # Some locked legacy wrappers (e.g., servalcat) call splitHklout with 3 positional args:
        #   splitHklout(miniMtzsOut, programColumnNames, mtz_path_string)
        # This incorrectly passes the MTZ path as outputColumnNames instead of inFile.
        #
        # Detection logic:
        # - If outputColumnNames is a string (not a list)
        # - AND it looks like a file path (contains '/' or ends with .mtz)
        # - AND inFile is still None
        # Then treat outputColumnNames as inFile and set outputColumnNames to None for auto-inference
        if (outputColumnNames is not None and
            isinstance(outputColumnNames, str) and
            ('/' in outputColumnNames or outputColumnNames.endswith('.mtz')) and
            inFile is None):
            logger.debug(f'[DEBUG splitHklout] Legacy wrapper bug detected: '
                        f'MTZ path passed as outputColumnNames: {outputColumnNames}')
            inFile = outputColumnNames
            outputColumnNames = None
            logger.debug(f'[DEBUG splitHklout] Corrected: inFile={inFile}, outputColumnNames=None (auto-inference)')

        # Default input file
        if inFile is None:
            inFile = str(self.workDirectory / 'hklout.mtz')

        # Validate input file exists
        from pathlib import Path
        if not Path(inFile).exists():
            error.append(
                klass=self.__class__.__name__,
                code=300,
                details=f"Input file not found: {inFile}"
            )
            return error

        # Validate arguments
        if len(miniMtzsOut) != len(programColumnNames):
            error.append(
                klass=self.__class__.__name__,
                code=301,
                details=f"miniMtzsOut and programColumnNames must have same length. "
                        f"Got {len(miniMtzsOut)} vs {len(programColumnNames)}"
            )
            return error

        # Validate outputColumnNames if provided
        if outputColumnNames is not None and len(outputColumnNames) != len(programColumnNames):
            error.append(
                klass=self.__class__.__name__,
                code=305,
                details=f"outputColumnNames must have same length as programColumnNames. "
                        f"Got {len(outputColumnNames)} vs {len(programColumnNames)}"
            )
            return error

        logger.debug(f'[DEBUG splitHklout] Splitting {inFile}')
        logger.debug(f'[DEBUG splitHklout] Output objects: {miniMtzsOut}')
        logger.debug(f'[DEBUG splitHklout] Input column names: {programColumnNames}')
        if outputColumnNames:
            logger.debug(f'[DEBUG splitHklout] Output column names (relabeling): {outputColumnNames}')

        # Build outfiles list for splitMtz
        outfiles = []
        for i, (obj_name, input_col_string) in enumerate(zip(miniMtzsOut, programColumnNames)):
            # Get output file object from container.outputData
            if not hasattr(self.container.outputData, obj_name):
                error.append(
                    klass=self.__class__.__name__,
                    code=302,
                    details=f"Output object '{obj_name}' not found in container.outputData"
                )
                continue

            file_obj = getattr(self.container.outputData, obj_name)

            # Get output file path
            output_path = file_obj.getFullPath()
            pass  # DEBUG: print(f"[DEBUG splitHklout] {obj_name}.getFullPath() returned: '{output_path}'")

            # Convert relative paths to absolute paths in work directory
            # This handles legacy wrappers that set baseName to just a filename
            from pathlib import Path
            if output_path:
                path_obj = Path(output_path)
                if not path_obj.is_absolute():
                    # Relative path - make it absolute relative to workDirectory
                    output_path = str(self.workDirectory / output_path)
                    pass  # DEBUG: print(f"[DEBUG splitHklout]   Converted to absolute: {output_path}")
            else:
                # Path not set at all - construct default path
                pass  # DEBUG: print(f"[DEBUG splitHklout] Auto-setting path for {obj_name}")

                # Get file extension
                extension = '.mtz'  # Default
                if hasattr(file_obj, 'fileExtensions') and callable(file_obj.fileExtensions):
                    try:
                        extensions = file_obj.fileExtensions()
                        if extensions and len(extensions) > 0:
                            extension = f'.{extensions[0]}'
                    except Exception:
                        pass

                # Construct path: workDirectory/OBJECTNAME.ext
                filename = f"{obj_name}{extension}"
                output_path = str(self.workDirectory / filename)
                pass  # DEBUG: print(f"[DEBUG splitHklout]   Created path: {output_path}")

                # Set the baseName so getFullPath() will work next time
                if hasattr(file_obj.baseName, 'set'):
                    file_obj.baseName.set(output_path)
                elif hasattr(file_obj.baseName, 'value'):
                    file_obj.baseName.value = output_path
                else:
                    file_obj.baseName = output_path

            # Build outfile spec for splitMtz
            # Format: [output_path, input_columns] or [output_path, input_columns, output_columns]

            # Determine output column names (with auto-inference from CONTENT_SIGNATURE_LIST)
            output_col_string = None

            if outputColumnNames is not None:
                # Explicit relabeling provided
                output_col_string = outputColumnNames[i]
                logger.debug(f'[DEBUG splitHklout]   Using explicit output columns: {output_col_string}')
            else:
                # Try to infer standard column names from file object's CONTENT_SIGNATURE_LIST
                # Use contentFlag to select the correct signature (contentFlag is 1-indexed)
                if hasattr(file_obj, 'CONTENT_SIGNATURE_LIST'):
                    sig_list = file_obj.CONTENT_SIGNATURE_LIST
                    if sig_list and len(sig_list) > 0:
                        # Get contentFlag from file object (1-indexed, default to 1)
                        content_flag = 1
                        if hasattr(file_obj, 'contentFlag'):
                            cf = file_obj.contentFlag
                            # Handle CInt wrapper
                            while hasattr(cf, 'value'):
                                cf = cf.value
                            if cf is not None:
                                try:
                                    content_flag = int(cf)
                                except (ValueError, TypeError):
                                    pass

                        # CONTENT_SIGNATURE_LIST is 0-indexed, contentFlag is 1-indexed
                        sig_index = max(0, content_flag - 1)
                        if sig_index < len(sig_list):
                            standard_cols = sig_list[sig_index]
                        else:
                            # Fall back to first signature if index out of range
                            standard_cols = sig_list[0]
                            logger.debug(f'[DEBUG splitHklout]   contentFlag={content_flag} out of range for CONTENT_SIGNATURE_LIST (len={len(sig_list)}), using [0]')

                        if standard_cols:
                            output_col_string = ','.join(standard_cols)
                            logger.debug(f'[DEBUG splitHklout]   Inferred standard column names from {file_obj.__class__.__name__}.CONTENT_SIGNATURE_LIST[{sig_index}] (contentFlag={content_flag}): {output_col_string}')

            # Add to outfiles list
            if output_col_string is not None:
                # Relabeling (explicit or inferred) - use 3-element format
                outfiles.append([output_path, input_col_string, output_col_string])
                logger.debug(f'[DEBUG splitHklout]   {obj_name} -> {output_path}')
                logger.debug(f'[DEBUG splitHklout]      Input columns: {input_col_string} -> Output columns: {output_col_string}')
            else:
                # No relabeling - use 2-element format
                outfiles.append([output_path, input_col_string])
                logger.debug(f'[DEBUG splitHklout]   {obj_name} -> {output_path} (columns: {input_col_string})')

        # Return early if there were errors
        if error.count() > 0:
            return error

        # Call splitMtz to do the actual work
        result = self.splitMtz(inFile, outfiles, logFile)

        if result == self.FAILED:
            error.append(
                klass=self.__class__.__name__,
                code=304,
                details="splitMtz failed"
            )

        return error

    def splitHkloutList(
        self,
        miniMtzsOut: list,
        programColumnNames: list,
        outputBaseName: list = None,
        outputContentFlags: list = None,
        infileList: list = None,
        logFile: str = None,
        **kwargs
    ) -> 'CErrorReport':
        """
        Split multiple HKLOUT files into mini-MTZ files (list processing version).

        This method processes a list of input MTZ files, splitting each one into
        multiple output files. For each input file and each miniMtzsOut entry,
        it creates a corresponding output file with auto-generated paths.

        Args:
            miniMtzsOut: List of output object names in container.outputData
                        (e.g., ['MAPOUT', 'ABCDOUT'])
            programColumnNames: List of comma-separated column name strings
                              (e.g., ['FWT,PHWT', 'HLA,HLB,HLC,HLD'])
            outputBaseName: List of base names for output files (one per miniMtzsOut)
                          If None, uses miniMtzsOut names as base names
            outputContentFlags: List of contentFlag values (one per miniMtzsOut)
                              If provided, sets contentFlag on each output file object
            infileList: List of input MTZ file paths or file objects
                       If None, uses [workDirectory/hklout.mtz]
            logFile: Log file path (default: workDirectory/splitmtz.log)
            **kwargs: Legacy compatibility

        Returns:
            CErrorReport with any errors

        Example:
            >>> # Split 2 input files, each into MAPOUT and ABCDOUT
            >>> error = self.splitHkloutList(
            ...     miniMtzsOut=['MAPOUT', 'ABCDOUT'],
            ...     programColumnNames=['FWT,PHWT', 'HLA,HLB,HLC,HLD'],
            ...     outputBaseName=['MAPOUT', 'ABCDOUT'],
            ...     infileList=self.container.outputData.HKLOUT
            ... )
            >>> # Creates: MAPOUT_1.mtz, ABCDOUT_1.mtz, MAPOUT_2.mtz, ABCDOUT_2.mtz
        """
        from ccp4i2.core.base_object.error_reporting import CErrorReport
        from pathlib import Path
        error = CErrorReport()

        # Default parameters
        if outputBaseName is None:
            outputBaseName = miniMtzsOut[:]  # Copy the list

        if infileList is None:
            infileList = [str(self.workDirectory / 'hklout.mtz')]

        # Validate arguments
        if len(miniMtzsOut) != len(programColumnNames):
            error.append(
                klass=self.__class__.__name__,
                code=305,
                details=f"miniMtzsOut and programColumnNames must have same length. "
                        f"Got {len(miniMtzsOut)} vs {len(programColumnNames)}"
            )
            return error

        if len(miniMtzsOut) != len(outputBaseName):
            error.append(
                klass=self.__class__.__name__,
                code=306,
                details=f"miniMtzsOut and outputBaseName must have same length. "
                        f"Got {len(miniMtzsOut)} vs {len(outputBaseName)}"
            )
            return error

        if outputContentFlags and len(miniMtzsOut) != len(outputContentFlags):
            error.append(
                klass=self.__class__.__name__,
                code=307,
                details=f"miniMtzsOut and outputContentFlags must have same length. "
                        f"Got {len(miniMtzsOut)} vs {len(outputContentFlags)}"
            )
            return error

        logger.debug(f'[DEBUG splitHkloutList] Processing {len(infileList)} input files')
        logger.debug(f'[DEBUG splitHkloutList] Output objects: {miniMtzsOut}')
        logger.debug(f'[DEBUG splitHkloutList] Output base names: {outputBaseName}')

        # Process each input file
        for file_index, infile in enumerate(infileList):
            # Handle file objects vs. string paths
            if hasattr(infile, '__str__'):
                infile_path = str(infile)
            elif hasattr(infile, 'getFullPath'):
                infile_path = infile.getFullPath()
            else:
                infile_path = infile

            logger.debug(f'[DEBUG splitHkloutList] Processing file {file_index + 1}: {infile_path}')

            # Validate input file exists
            if not Path(infile_path).exists():
                error.append(
                    klass=self.__class__.__name__,
                    code=308,
                    details=f"Input file not found: {infile_path}"
                )
                continue

            # Prepare output file objects and paths
            outfiles = []
            for obj_index, (obj_name, col_string, base_name) in enumerate(
                zip(miniMtzsOut, programColumnNames, outputBaseName)
            ):
                # Get or create output list object
                if not hasattr(self.container.outputData, obj_name):
                    error.append(
                        klass=self.__class__.__name__,
                        code=309,
                        details=f"Output object '{obj_name}' not found in container.outputData"
                    )
                    continue

                obj_list = getattr(self.container.outputData, obj_name)

                # Get file type label for output filename suffix
                file_type_label = ''
                if hasattr(obj_list, 'subItemQualifiers') and callable(obj_list.subItemQualifiers):
                    try:
                        label = obj_list.subItemQualifiers('label')
                        if label and label is not NotImplemented and len(label) > 0:
                            file_type_label = '_' + label
                    except Exception:
                        pass

                # Get GUI label for annotation
                gui_label = 'Output '
                if hasattr(obj_list, 'subItemQualifiers') and callable(obj_list.subItemQualifiers):
                    try:
                        label = obj_list.subItemQualifiers('guiLabel')
                        if label and label is not NotImplemented and len(label) > 0:
                            gui_label = label + ' '
                    except Exception:
                        pass

                # Ensure the list has enough items (expand if necessary)
                while file_index >= len(obj_list):
                    obj_list.append(obj_list.makeItem())

                # Construct output file path: baseName_fileIndex_typeLabel.mtz
                output_filename = f"{base_name}_{file_index + 1}{file_type_label}.mtz"
                output_path = str(self.workDirectory / output_filename)

                # Set the path on the output file object
                obj_list[file_index].setFullPath(output_path)

                # Set annotation
                annotation_text = f"{gui_label}{file_index + 1}"
                if hasattr(obj_list[file_index], 'annotation'):
                    if hasattr(obj_list[file_index].annotation, 'set'):
                        obj_list[file_index].annotation.set(annotation_text)
                    else:
                        obj_list[file_index].annotation = annotation_text

                # Set contentFlag if provided
                if outputContentFlags and obj_index < len(outputContentFlags):
                    if hasattr(obj_list[file_index], 'contentFlag'):
                        if hasattr(obj_list[file_index].contentFlag, 'set'):
                            obj_list[file_index].contentFlag.set(outputContentFlags[obj_index])
                        else:
                            obj_list[file_index].contentFlag = outputContentFlags[obj_index]

                # Add to outfiles for splitMtz: [path, input_columns, output_columns]
                # Input and output columns are the same
                outfiles.append([output_path, col_string, col_string])

                logger.debug(f'[DEBUG splitHkloutList]   {obj_name}[{file_index}] -> {output_path} (columns: {col_string})')

            # Call splitMtz to split this input file
            if outfiles:
                result = self.splitMtz(infile_path, outfiles, logFile)

                if result == self.FAILED:
                    error.append(
                        klass=self.__class__.__name__,
                        code=310,
                        details=f"splitMtz failed for input file: {infile_path}"
                    )

        return error

    def appendCommandScript(self, text=None, fileName=None, oneLine=False, clear=False):
        """
        Add text to the command script (list of lines for stdin/script file).

        Args:
            text: String or list of strings to add
            fileName: Path to file whose contents should be added
            oneLine: If text is a list, join into single line
            clear: Clear existing script before adding

        Returns:
            CErrorReport with any errors encountered
        """
        from ccp4i2.core.CCP4ErrorHandling import CErrorReport

        error = CErrorReport()

        if clear:
            self.commandScript = []

        # Load from file if provided
        if fileName is not None:
            fileName = str(fileName)
            if not os.path.exists(fileName):
                error.append(
                    klass=self.__class__.__name__,
                    code=16,
                    details=f"File not found: {fileName}",
                    name=fileName
                )
                return error
            try:
                with open(fileName, 'r') as f:
                    text = f.read()
            except Exception as e:
                error.append(
                    klass=self.__class__.__name__,
                    code=16,
                    details=f"Error reading file: {e}",
                    name=fileName
                )
                return error

        # Process text
        if text is not None:
            # Handle list inputs
            if isinstance(text, list):
                if not oneLine:
                    # Add each item separately
                    for item in text:
                        sub_error = self.appendCommandScript(item)
                        if sub_error.count() > 0:
                            error.extend(sub_error)
                    return error
                else:
                    # Join into single line
                    text = ' '.join(str(item) for item in text)

            # Convert to string
            try:
                text_str = str(text)
            except Exception:
                error.append(
                    klass=self.__class__.__name__,
                    code=5,
                    details="Could not convert text to string"
                )
                return error

            # Ensure newline at end
            if text_str and not text_str.endswith('\n'):
                text_str += '\n'

            self.commandScript.append(text_str)

        return error

    def writeCommandFile(self, qualifier=None):
        """
        Write the command script to a file.

        Args:
            qualifier: Optional qualifier for filename (creates com_{qualifier}.txt)

        Returns:
            Path to written file, or None on error
        """
        if not self.commandScript:
            return None

        # Prepend header comment
        script_lines = [f'# Task {self.TASKNAME} command script\n'] + self.commandScript

        # Generate filename
        fileName = self.makeFileName('COM', qualifier=qualifier)

        try:
            with open(fileName, 'w') as f:
                f.writelines(script_lines)
            return fileName
        except Exception as e:
            print(f'[ERROR] Writing command file {fileName}: {e}')
            self.errorReport.append(
                klass=self.__class__.__name__,
                code=7,
                details=f'Error writing command file: {e}',
                name=fileName
            )
            return None

    def makeFileName(self, format='COM', ext='', qualifier=None):
        """
        Generate consistent names for output files.

        Args:
            format: File type (COM, LOG, STDOUT, STDERR, etc.)
            ext: Custom extension (overrides format default)
            qualifier: Optional qualifier to add to filename

        Returns:
            Full path to file in work directory
        """
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
            'COM': 'com.txt',
            'MGPICDEF': 'report.mgpic.py',
            'PIC': 'report.png',
            'RVAPIXML': 'i2.xml'
        }

        fileName = defNames.get(format, 'unknown.unk')

        # Add qualifier if provided
        if qualifier is not None:
            base, extension = fileName.rsplit('.', 1)
            fileName = f'{base}_{qualifier}.{extension}'

        # Use custom extension if provided
        if ext:
            base = fileName.rsplit('.', 1)[0]
            fileName = f'{base}.{ext}'

        return str(self.workDirectory / fileName)

    def createWarningsXML(self, logfiles=[]):
        """
        Parse log files for WARNING, ERROR, and ATTENTION messages and add to XML tree.

        This is primarily used by prosmart_refmac to extract warnings from refmac logfiles.
        It searches for lines containing ERROR, WARNING, or ATTENTION and adds them to
        self.xmlroot as structured XML.

        Args:
            logfiles: List of log file paths to parse
        """
        from lxml import etree
        import re

        if not hasattr(self, 'xmlroot'):
            print("Warning: createWarningsXML called but self.xmlroot doesn't exist")
            return

        warningsNode = etree.SubElement(self.xmlroot, "Warnings")
        for f in logfiles:
            if not os.path.exists(f):
                print(f"Warning: Log file does not exist: {f}")
                continue

            fileNode = etree.SubElement(warningsNode, "logFile")
            fileNameNode = etree.SubElement(fileNode, "fileName")
            fileNameNode.text = f

            try:
                with open(f) as fh:
                    lines = fh.readlines()
                    il = 1
                    for l in lines:
                        errors = re.findall("ERROR", l)
                        if len(errors) > 0:
                            warningNode = etree.SubElement(fileNode, "warning")
                            typeNode = etree.SubElement(warningNode, "type")
                            warningTextNode = etree.SubElement(warningNode, "text")
                            typeNode.text = "ERROR"
                            lineNumberNode = etree.SubElement(warningNode, "lineNumber")
                            warningTextNode.text = l.rstrip("\n")
                            lineNumberNode.text = str(il)

                        warnings = re.findall("WARNING", l)
                        if len(warnings) > 0:
                            warningNode = etree.SubElement(fileNode, "warning")
                            typeNode = etree.SubElement(warningNode, "type")
                            warningTextNode = etree.SubElement(warningNode, "text")
                            typeNode.text = "WARNING"
                            lineNumberNode = etree.SubElement(warningNode, "lineNumber")
                            warningTextNode.text = l.rstrip("\n")
                            lineNumberNode.text = str(il)

                        attentions = re.findall("ATTENTION", l)
                        if len(attentions) > 0:
                            warningNode = etree.SubElement(fileNode, "warning")
                            typeNode = etree.SubElement(warningNode, "type")
                            warningTextNode = etree.SubElement(warningNode, "text")
                            typeNode.text = "ATTENTION"
                            lineNumberNode = etree.SubElement(warningNode, "lineNumber")
                            warningTextNode.text = l.rstrip("\n")
                            lineNumberNode.text = str(il)

                        il += 1
            except Exception as e:
                print(f"Error reading log file {f}: {e}")

    def getProcessId(self):
        """
        Get the process ID for this plugin instance.

        In synchronous execution, this is just the Python object ID.
        In async execution with Qt, this would be the QProcess pid.

        Returns:
            Process identifier (int)
        """
        return id(self)

    def appendErrorReport(self, code=0, details='', name=None, label=None, cls=None,
                         recordTime=False, stack=True, exc_info=None, severity=None):
        """
        Append an error to the plugin's error report.

        This is a simplified version for legacy compatibility.

        Args:
            code: Error code number
            details: Error message details
            name: Error name
            label: Error label
            cls: Class where error occurred
            recordTime: Whether to record timestamp
            stack: Whether to include stack trace
            exc_info: Exception info tuple
            severity: Error severity (defaults to ERROR if 'exception' in details, otherwise WARNING)
        """
        if cls is None:
            cls = self.__class__

        # Determine severity
        # If severity not explicitly provided, use ERROR for exceptions, WARNING otherwise
        if severity is None:
            # If details mention 'exception', treat as ERROR
            if 'exception' in details.lower():
                severity = SEVERITY_ERROR
            else:
                # Legacy code often uses appendErrorReport for non-fatal issues
                severity = SEVERITY_WARNING

        self.errorReport.append(
            klass=cls.__name__,
            code=code,
            details=details,
            severity=severity
        )
