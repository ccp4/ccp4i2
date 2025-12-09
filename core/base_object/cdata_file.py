"""
CDataFile - Base class for file-related CData classes.

Provides file path management with:
- setFullPath() / getFullPath() for file path operations
- fullPath property for convenient access
- __str__() returns the file path
- Automatic contentFlag introspection
- File I/O placeholder methods (load_from_file, save_to_file)

File path handling:
- In standalone mode: Uses baseName attribute
- In database mode (future): Uses project + relPath + baseName
"""

from typing import Optional, Any
from pathlib import Path
import logging

from .cdata import CData
from .class_metadata import cdata_class, attribute, AttributeType

# Configure logger
logger = logging.getLogger(__name__)


@cdata_class(
    attributes={
        "project": attribute(AttributeType.CUSTOM, custom_class="CUUID"),
        "baseName": attribute(AttributeType.CUSTOM, custom_class="CFilePath"),
        "relPath": attribute(AttributeType.CUSTOM, custom_class="CFilePath"),
        "dbFileId": attribute(AttributeType.CUSTOM, custom_class="CUUID"),
        "annotation": attribute(AttributeType.CUSTOM, custom_class="CString"),
        "subType": attribute(AttributeType.CUSTOM, custom_class="CInt"),
        "contentFlag": attribute(AttributeType.CUSTOM, custom_class="CInt"),
    },
    qualifiers={
        "allowUndefined": True,
        "mustExist": False,
        "fromPreviousJob": False,
        "jobCombo": True,
        "mimeTypeName": "",
        "mimeTypeDescription": "",
        "fileLabel": None,
        "fileExtensions": [],
        "fileContentClassName": None,
        "isDirectory": False,
        "saveToDb": True,
        "requiredSubType": None,
        "requiredContentFlag": None
    }
)
class CDataFile(CData):
    """Base class for file-related CData classes.

    Attributes are automatically created from embedded metadata:
    - baseName: CString - Base filename
    - relPath: CString - Relative path to file
    - annotation: CString - File annotation
    - subType: CInt - File subtype (optional)
    - contentFlag: CInt - Content flag (min=0, optional)

    CDataFile has special behavior for file path handling:
    - Assigning a string to a CDataFile sets the baseName attribute
    - __str__() returns the file path (baseName when not database-connected)
    - fullPath property and getFullPath() method return the full path
    """

    def __init__(self, file_path: str = None, parent=None, name=None, **kwargs):
        # Pass per-instance metadata overrides to base
        super().__init__(parent=parent, name=name, **kwargs)

        # Legacy compatibility
        self.file_path = file_path

        # File content instance (instantiated by loadFile())
        self.content = None

        # Set baseName if file_path provided
        if file_path is not None:
            self.setFullPath(file_path)

    @property
    def fileContent(self):
        """
        Get or auto-create the file content object.

        This property provides fault-tolerant access to file content:
        - If content already exists, return it
        - If content is None, auto-create it based on fileContentClassName qualifier

        This supports legacy code patterns like:
            file.fileContent.loadFile(path)

        Returns:
            CDataFileContent instance or None if content class not configured
        """
        # If content already exists but appears unloaded, try to load it
        if self.content is not None and hasattr(self.content, 'loadFile'):
            # Check if the content appears to be unloaded (e.g., resolutionRange.high is not set for MTZ files)
            if hasattr(self.content, 'resolutionRange') and self.content.resolutionRange is not None:
                high_val = self.content.resolutionRange.high
                is_unloaded = (high_val is None or
                             (hasattr(high_val, 'value') and high_val.value is None) or
                             (hasattr(high_val, 'isSet') and not high_val.isSet()))
                if is_unloaded:
                    try:
                        error = self.content.loadFile()
                        if error and hasattr(error, 'count') and error.count() > 0:
                            logger.debug(
                                "[fileContent property] loadFile() returned errors: %s",
                                error
                            )
                        else:
                            logger.debug(
                                "[fileContent property] Successfully loaded existing content for %s",
                                self.name if hasattr(self, 'name') else 'unnamed'
                            )
                    except Exception as e:
                        logger.warning(
                            "[fileContent property] Failed to load existing content: %s",
                            e
                        )

        if self.content is None:
            # Auto-create content object
            content_class_name = self.get_qualifier('fileContentClassName')
            if content_class_name:
                try:
                    # Import the content class
                    content_class = None
                    for module_name in ['core.CCP4XtalData', 'core.CCP4ModelData', 'core.CCP4CootData']:
                        try:
                            module = __import__(module_name, fromlist=[content_class_name])
                            content_class = getattr(module, content_class_name, None)
                            if content_class is not None:
                                break
                        except (ImportError, AttributeError):
                            continue

                    if content_class is not None:
                        # Create content instance
                        self.content = content_class(parent=self, name='fileContent')
                        logger.debug(
                            "[fileContent property] Auto-created content for %s: %s",
                            self.name if hasattr(self, 'name') else 'unnamed',
                            content_class_name
                        )

                        # Auto-load by calling loadFile() without a path
                        # Let loadFile() get the path from the parent CDataFile
                        if hasattr(self.content, 'loadFile'):
                            try:
                                logger.debug("[fileContent property] Calling loadFile() to get path from parent")
                                error = self.content.loadFile()  # loadFile will call self.get_parent().getFullPath()
                                if error and hasattr(error, 'count') and error.count() > 0:
                                    logger.warning(
                                        "[fileContent property] loadFile() returned errors: %s",
                                        error
                                    )
                                else:
                                    logger.debug(
                                        "[fileContent property] Successfully auto-loaded file content"
                                    )
                            except Exception as load_error:
                                logger.warning(
                                    "[fileContent property] Failed to auto-load: %s",
                                    load_error
                                )
                except Exception as e:
                    logger.warning(
                        "[fileContent property] Failed to auto-create content '%s': %s",
                        content_class_name, e
                    )

        return self.content

    @fileContent.setter
    def fileContent(self, value):
        """Set the file content object."""
        self.content = value

    def setFullPath(self, path: str):
        """Set the full path of the file.

        In non-database-connected mode, this sets the baseName attribute.
        In database-connected mode, this parses the path and sets project,
        relPath, baseName, and dbFileId appropriately.

        Args:
            path: Full file path as a string
        """
        from pathlib import Path

        logger.debug(
            "[setFullPath] Called for %s, input path: %s, hasattr baseName: %s",
            self.name if hasattr(self, 'name') else 'unnamed',
            path,
            hasattr(self, 'baseName')
        )
        if hasattr(self, 'baseName'):
            logger.debug("  baseName type: %s, has .value: %s", type(self.baseName), hasattr(self.baseName, 'value'))

        # Check if we're in database-aware mode FIRST, before setting baseName
        plugin = self._find_plugin_parent()
        is_db_mode = False
        dbHandler = self._get_db_handler()
        if plugin and dbHandler and hasattr(plugin, '_dbProjectId') and plugin._dbProjectId:
            is_db_mode = True
            logger.debug("  DB-aware mode: will parse path into project/relPath/baseName")
            try:
                # Parse path into project/relPath/baseName for output files
                # This will set baseName to just the filename
                self._parse_output_path_database(path, plugin)
            except Exception as e:
                # Don't fail on database errors - fall back to simple mode
                logger.debug("Failed to parse output path for database: %s", e)
                is_db_mode = False

        # If not in database mode, or if parsing failed, use simple path storage
        if not is_db_mode:
            logger.debug("  Non-DB mode: setting baseName to full path")
            if hasattr(self, 'baseName'):
                # If baseName is a CFilePath or similar, set its value
                if hasattr(self.baseName, 'value'):
                    self.baseName.value = path
                    logger.debug("  Set baseName.value = %s", path)
                else:
                    # Otherwise set it directly
                    self.baseName = path
                    logger.debug("  Set baseName = %s", path)
            else:
                # baseName doesn't exist yet - store in file_path for now
                object.__setattr__(self, 'file_path', path)
                logger.debug("  Set file_path = %s (baseName doesn't exist)", path)

        # After path is set, attempt to auto-detect content flag if file exists
        try:
            final_path = self.getFullPath()
            if final_path and Path(final_path).exists():
                logger.debug("[setFullPath] File exists at %s, attempting to detect content flag", final_path)
                self.setContentFlag()  # Auto-detect based on file introspection
        except Exception as e:
            # Don't fail setFullPath if content detection fails
            # This is an optional enhancement, not a critical operation
            logger.debug("[setFullPath] Failed to auto-detect content flag: %s", e)

    def _find_plugin_parent(self):
        """Walk up the parent hierarchy to find the CPluginScript parent."""
        from core.CCP4PluginScript import CPluginScript

        # Check for temporary plugin reference (set during template expansion or pipeline finishUp)
        if hasattr(self, '_temp_plugin_ref'):
            return self._temp_plugin_ref

        # Walk up parent hierarchy
        current = self.parent()
        depth = 0
        while current is not None:
            if isinstance(current, CPluginScript):
                return current
            # parent is now a method, not a property - must call it
            current = current.parent() if hasattr(current, 'parent') and callable(current.parent) else None
            depth += 1
            if depth > 10:  # Safety limit
                break
        return None

    def _get_db_handler(self):
        """Get the database handler from the plugin parent, if available."""
        plugin = self._find_plugin_parent()
        if plugin and hasattr(plugin, '_dbHandler'):
            return plugin._dbHandler
        return None

    def _update_from_database(self, path: str, plugin):
        """Update file attributes from database if the file matches an existing database record.

        Uses the database handler attached to the plugin, if available. This keeps
        CDataFile decoupled from Django and allows for different database backends.

        Args:
            path: Full file path
            plugin: Parent CPluginScript instance
        """
        from pathlib import Path
        import re

        # Get database handler from plugin
        db_handler = self._get_db_handler()
        if not db_handler:
            return  # No database handler available

        abs_path = Path(path).resolve()

        # Try to extract job number from path
        # Path structure: .../project_dir/CCP4_JOBS/job_1/job_2/filename.ext
        # Should extract "1.2" from job_1/job_2
        path_str = str(abs_path)

        # Look for CCP4_JOBS directory pattern
        jobs_match = re.search(r'CCP4_JOBS/(job_\d+(?:/job_\d+)*)', path_str)
        if not jobs_match:
            return  # Not a job directory path

        # Extract job path like "job_1/job_2"
        job_path = jobs_match.group(1)
        # Convert to job number format: "job_1/job_2" → "1.2"
        job_numbers = re.findall(r'job_(\d+)', job_path)
        job_number = '.'.join(job_numbers)

        filename = abs_path.name

        # Query database via handler (returns dict, not Django model)
        try:
            file_info = db_handler.find_file_by_path_sync(
                file_path=str(abs_path),
                job_number=job_number,
                filename=filename
            )

            if file_info:
                # Update CDataFile attributes from database record
                if hasattr(self, 'dbFileId') and hasattr(self.dbFileId, 'value'):
                    self.dbFileId.value = file_info['uuid']

                # Set relPath if available
                if file_info.get('relative_path'):
                    if hasattr(self, 'relPath') and hasattr(self.relPath, 'value'):
                        self.relPath.value = file_info['relative_path']

                # Set project ID
                project_id = getattr(plugin, '_dbProjectId', None)
                if project_id:
                    if hasattr(self, 'project') and hasattr(self.project, 'value'):
                        self.project.value = str(project_id)

        except Exception as e:
            # Silently ignore errors - we're in a file setter
            import logging
            logger = logging.getLogger(__name__)
            logger.debug(f"Failed to update from database: {e}")

    def _parse_output_path_database(self, path: str, plugin):
        """Parse output file path and set project/relPath/baseName for database context.

        This method parses the full path provided to setFullPath() and extracts:
        - project: UUID of the project (from plugin._dbProjectId)
        - relPath: Relative path from project root (e.g., "CCP4_JOBS/job_17" or "CCP4_IMPORTED_FILES")
        - baseName: Just the filename (e.g., "XYZOUT.mtz")

        For external files (paths without CCP4_JOBS or CCP4_IMPORTED_FILES markers):
        - If the input is a relative path, preserve it as-is in baseName
        - If the input is an absolute path, store the full path in baseName

        This distinction is important for ASU content files where source/baseName
        is a user-provided filename like "beta.seq" that should NOT be converted
        to an absolute path using CWD.

        Args:
            path: Full file path (e.g., "/path/to/project/CCP4_JOBS/job_17/XYZOUT.mtz")
            plugin: Parent CPluginScript instance
        """
        from pathlib import Path
        import re

        input_path = Path(path)

        # Only resolve to absolute path if the input is already absolute
        # or if it contains known project markers (CCP4_JOBS, CCP4_IMPORTED_FILES)
        # This preserves relative external paths like "beta.seq" as-is
        is_project_path = 'CCP4_JOBS' in path or 'CCP4_IMPORTED_FILES' in path

        if input_path.is_absolute() or is_project_path:
            abs_path = input_path.resolve()
            path_str = str(abs_path)
        else:
            # Keep relative paths as-is - don't convert to absolute using CWD
            path_str = path

        # Get project ID from plugin
        project_id = getattr(plugin, '_dbProjectId', None)
        if not project_id:
            logger.debug("  No _dbProjectId on plugin, cannot parse for database")
            return

        # Determine relPath by looking for known directory patterns FIRST
        # Pattern 1: CCP4_JOBS/job_17 or CCP4_JOBS/job_17/job_1 (nested jobs)
        jobs_match = re.search(r'(CCP4_JOBS/job_\d+(?:/job_\d+)*)', path_str)
        if jobs_match:
            rel_path = jobs_match.group(1)
            # For project paths, set project UUID and extract baseName
            if hasattr(self, 'project') and hasattr(self.project, 'set'):
                self.project.set(str(project_id))
                logger.debug("  Set project = %s", project_id)
            if hasattr(self, 'baseName') and hasattr(self.baseName, 'set'):
                self.baseName.set(Path(path_str).name)
                logger.debug("  Set baseName = %s", Path(path_str).name)
            if hasattr(self, 'relPath') and hasattr(self.relPath, 'set'):
                self.relPath.set(rel_path)
                logger.debug("  Set relPath = %s (from CCP4_JOBS pattern)", rel_path)
            return

        # Pattern 2: CCP4_IMPORTED_FILES
        if 'CCP4_IMPORTED_FILES' in path_str:
            # For imported files, set project UUID and extract baseName
            if hasattr(self, 'project') and hasattr(self.project, 'set'):
                self.project.set(str(project_id))
                logger.debug("  Set project = %s", project_id)
            if hasattr(self, 'baseName') and hasattr(self.baseName, 'set'):
                self.baseName.set(Path(path_str).name)
                logger.debug("  Set baseName = %s", Path(path_str).name)
            if hasattr(self, 'relPath') and hasattr(self.relPath, 'set'):
                self.relPath.set('CCP4_IMPORTED_FILES')
                logger.debug("  Set relPath = CCP4_IMPORTED_FILES")
            return

        # Pattern 3: External path (outside project structure)
        # Keep the path as provided (relative or absolute) in baseName
        # Do NOT set project for external files - they are not part of the project
        logger.debug("  Path is external (outside project structure), keeping as-is in baseName")
        if hasattr(self, 'baseName') and hasattr(self.baseName, 'set'):
            self.baseName.set(path_str)
            logger.debug("  Set baseName = %s (external path)", path_str)

    def setOutputPath(self, jobName: str = "", projectId: str = None, relPath: str = None):
        """Set output file path using project/relPath/baseName structure (legacy API).

        This is called by set_output_file_names() before job execution to prepare
        output file paths. It sets project, relPath, and baseName fields.

        Args:
            jobName: Optional prefix for the filename (usually empty)
            projectId: Project UUID
            relPath: Relative path from project root (e.g., "CCP4_JOBS/job_17")
        """
        from pathlib import Path

        logger.debug(f"[setOutputPath] Called for {self.name if hasattr(self, 'name') else 'unknown'}")
        logger.debug(f"  jobName={jobName}, projectId={projectId}, relPath={relPath}")

        # Set project ID
        if projectId:
            logger.debug(f"  Has project attr: {hasattr(self, 'project')}")
            if hasattr(self, 'project'):
                logger.debug(f"  Project type: {type(self.project).__name__}")
                logger.debug(f"  Project has set: {hasattr(self.project, 'set')}")
                logger.debug(f"  Project has value: {hasattr(self.project, 'value')}")
                if hasattr(self.project, 'set'):
                    self.project.set(projectId)
                    logger.debug(f"  Called project.set({projectId})")
                elif hasattr(self.project, 'value'):
                    self.project.value = projectId
                    logger.debug(f"  Set project.value = {projectId}")
            logger.debug(f"[setOutputPath] Set project = {projectId}")

        # Set relPath
        if relPath:
            logger.debug(f"  Has relPath attr: {hasattr(self, 'relPath')}")
            if hasattr(self, 'relPath'):
                logger.debug(f"  relPath has set: {hasattr(self.relPath, 'set')}")
                if hasattr(self.relPath, 'set'):
                    self.relPath.set(str(relPath))
                    logger.debug(f"  Called relPath.set({relPath})")
                elif hasattr(self.relPath, 'value'):
                    self.relPath.value = str(relPath)
                    logger.debug(f"  Set relPath.value = {relPath}")
            logger.debug(f"[setOutputPath] Set relPath = {relPath}")

        # Generate baseName from object name if not already set
        if hasattr(self, 'baseName'):
            # Check if baseName is already explicitly set
            basename_is_set = False
            if hasattr(self.baseName, 'isSet'):
                basename_is_set = self.baseName.isSet('value') and bool(str(self.baseName).strip())

            logger.debug(f"  baseName is already set: {basename_is_set}")
            if not basename_is_set:
                # Generate from object name
                obj_name = self.objectName() if hasattr(self, 'objectName') else (self.name if hasattr(self, 'name') else 'output')
                logger.debug(f"  obj_name from objectName/name: {obj_name}")

                # Add job prefix if provided
                if jobName:
                    filename = f"{jobName}{obj_name}"
                else:
                    filename = obj_name

                # Add extension if missing
                if not any(filename.endswith(ext) for ext in ['.mtz', '.pdb', '.cif', '.log', '.xml', '.txt']):
                    # Get appropriate extension from class metadata fileExtensions
                    default_ext = '.mtz'  # Fallback default

                    # Try getting from class metadata
                    from core.base_object.class_metadata import get_class_metadata_by_type
                    try:
                        meta = get_class_metadata_by_type(type(self))
                        if meta and meta.qualifiers and 'fileExtensions' in meta.qualifiers:
                            file_exts = meta.qualifiers['fileExtensions']
                            if file_exts and isinstance(file_exts, list) and len(file_exts) > 0:
                                default_ext = '.' + file_exts[0].lstrip('.')
                                logger.debug(f"  Using extension {default_ext} from class metadata fileExtensions={file_exts}")
                    except Exception as e:
                        logger.debug(f"  Failed to get class metadata: {e}")

                    filename += default_ext

                logger.debug(f"  Generated filename: {filename}")
                if hasattr(self.baseName, 'set'):
                    self.baseName.set(filename)
                    logger.debug(f"  Called baseName.set({filename})")
                elif hasattr(self.baseName, 'value'):
                    self.baseName.value = filename
                    logger.debug(f"  Set baseName.value = {filename}")
                else:
                    self.baseName = filename
                    logger.debug(f"  Set baseName directly = {filename}")

                logger.debug(f"[setOutputPath] Set baseName = {filename}")

        logger.debug(f"[setOutputPath] Done!")

    def getFullPath(self) -> str:
        """Get the full file path as a string.

        Path construction logic:
        1. If baseName contains an absolute path → return it directly
        2. If dbFileId is set (database-aware) → retrieve from database
        3. If relPath is set → construct: workDirectory/relPath/baseName
           (This is the standard for output files: relPath = "CCP4_JOBS/job_n/job_m/job_y")
        4. If only baseName is set → return as-is (may be relative to CWD)
        5. Fallback to legacy file_path attribute

        For output files computed in checkOutputData():
        - relPath: "CCP4_JOBS/job_n/job_m/job_y"
        - baseName: just the filename (e.g., "refined.pdb")
        - fullPath: workDirectory + relPath + baseName

        For user-specified input files:
        - baseName: absolute or relative path (e.g., "/path/to/file.mtz" or "demo_data/test.mtz")
        - relPath: typically empty
        - fullPath: baseName as-is

        Returns:
            Full path to the file, or empty string if not set
        """
        from pathlib import Path

        # First, get baseName value if it exists
        basename_value = None
        if hasattr(self, 'baseName') and self.baseName is not None:
            if hasattr(self.baseName, 'value'):
                basename_value = self.baseName.value
            else:
                basename_value = self.baseName

        # If baseName is an absolute path, check if it's stale or valid
        if basename_value:
            basename_str = str(basename_value)
            if basename_str and Path(basename_str).is_absolute():
                abs_path = Path(basename_str)
                is_temp_file = str(abs_path).startswith(('/tmp/', '/var/tmp/', 'C:\\Temp\\', 'C:\\Windows\\Temp\\'))
                has_relpath = hasattr(self, 'relPath') and self.relPath is not None and \
                             (hasattr(self.relPath, 'value') and self.relPath.value) or self.relPath

                pass  # DEBUG: print(f"[DEBUG getFullPath] baseName is absolute: {basename_str}, exists={abs_path.exists()}, is_temp_file={is_temp_file}, has_relpath={has_relpath}")

                # If file exists at the absolute path, use it
                if abs_path.exists():
                    pass  # DEBUG: print(f"[DEBUG getFullPath] File exists at absolute path, returning: {basename_str}")
                    return basename_str

                # If we have relPath set, the absolute path in baseName is likely stale
                # (from CWD-relative path expansion). Strip to just filename and use relPath logic.
                if has_relpath:
                    pass  # DEBUG: print(f"[DEBUG getFullPath] Have relPath, treating absolute baseName as stale, stripping to filename")
                    basename_value = abs_path.name
                    basename_str = basename_value
                elif not is_temp_file:
                    # No relPath, not a temp file, use the absolute path as-is
                    pass  # DEBUG: print(f"[DEBUG getFullPath] No relPath, using absolute baseName as-is: {basename_str}")
                    return basename_str
                else:
                    # Temp file that doesn't exist - strip to filename
                    pass  # DEBUG: print(f"[DEBUG getFullPath] Stale temp file, stripping to filename")
                    basename_value = abs_path.name
                    basename_str = basename_value

        # Database-aware mode: Check if dbFileId is set
        if hasattr(self, 'dbFileId') and self.dbFileId is not None:
            db_file_id = None
            if hasattr(self.dbFileId, 'value') and self.dbFileId.value is not None:
                db_file_id = self.dbFileId.value
            elif isinstance(self.dbFileId, str):
                db_file_id = self.dbFileId

            if db_file_id:
                db_handler = self._get_db_handler()
                if db_handler:
                    try:
                        import uuid
                        file_uuid = uuid.UUID(str(db_file_id))
                        path = db_handler.get_file_path_sync(file_uuid)
                        if path:
                            logger.debug(f"Retrieved path from database via dbFileId: {path}")
                            return path
                    except Exception as e:
                        logger.debug(f"Failed to retrieve path from database: {e}")

        # Standard path construction: workDirectory + relPath + baseName
        # This is used for output files where relPath = "CCP4_JOBS/job_n/job_m/job_y"
        relpath_value = None
        if hasattr(self, 'relPath') and self.relPath is not None:
            if hasattr(self.relPath, 'value'):
                relpath_value = self.relPath.value
            else:
                relpath_value = self.relPath

        if relpath_value:
            relpath_str = str(relpath_value).strip()
            if relpath_str:
                # Handle environment variable expansion in relPath
                # Patterns like "$CCP4/path/to/file" or "$CCP4I2_ROOT/demo_data/file"
                # should expand the environment variable and construct the full path
                import os
                import re
                env_var_match = re.match(r'^\$([A-Za-z_][A-Za-z0-9_]*)(.*)', relpath_str)
                if env_var_match:
                    env_var_name = env_var_match.group(1)
                    rest_of_path = env_var_match.group(2).lstrip('/')  # Remove leading slash
                    env_value = os.environ.get(env_var_name)
                    if env_value:
                        # Construct full path: env_value / rest_of_path / baseName
                        base_dir = Path(env_value)
                        if rest_of_path:
                            base_dir = base_dir / rest_of_path
                        if basename_value:
                            full_path = base_dir / str(basename_value)
                            logger.debug("Constructed path from env var $%s: %s", env_var_name, full_path)
                            return str(full_path)
                        else:
                            logger.debug("Constructed dir path from env var $%s: %s", env_var_name, base_dir)
                            return str(base_dir)
                    else:
                        logger.warning("Environment variable $%s not set for relPath: %s", env_var_name, relpath_str)

                # Determine base directory: PROJECT root for imported files, workDirectory for job outputs
                plugin = self._find_plugin_parent()
                if plugin:
                    # For CCP4_IMPORTED_FILES, use project directory as base
                    if relpath_str == 'CCP4_IMPORTED_FILES' or relpath_str.startswith('CCP4_IMPORTED_FILES/'):
                        # Need to get project directory
                        # Get project UUID from self.project or plugin._dbProjectId
                        project_id = None
                        if hasattr(self, 'project') and self.project is not None:
                            if hasattr(self.project, 'value'):
                                project_id = self.project.value
                            else:
                                project_id = str(self.project)
                        elif hasattr(plugin, '_dbProjectId'):
                            project_id = plugin._dbProjectId

                        if project_id:
                            # Try to get project directory from database handler
                            db_handler = self._get_db_handler()
                            if db_handler and hasattr(db_handler, 'getProjectDirectory'):
                                try:
                                    project_dir = db_handler.getProjectDirectory(project_id)
                                    if project_dir:
                                        base_dir = Path(project_dir)
                                        if basename_value:
                                            full_path = base_dir / relpath_str / str(basename_value)
                                            logger.debug("Constructed imported file path from db: %s", full_path)
                                            return str(full_path)
                                        else:
                                            full_path = base_dir / relpath_str
                                            logger.debug("Constructed imported dir path from db: %s", full_path)
                                            return str(full_path)
                                except Exception as e:
                                    logger.debug(f"Failed to get project directory from db: {e}")

                            # Fallback: infer project directory from workDirectory
                            # workDirectory is like: /path/to/project/CCP4_JOBS/job_X/job_Y
                            # We need to go up to the project root
                            if hasattr(plugin, 'workDirectory'):
                                work_dir = Path(plugin.workDirectory)
                                # Walk up until we find a directory containing CCP4_IMPORTED_FILES
                                current = work_dir
                                for _ in range(10):  # Max 10 levels up
                                    imported_dir = current / 'CCP4_IMPORTED_FILES'
                                    if imported_dir.exists():
                                        # Found project root
                                        if basename_value:
                                            full_path = current / relpath_str / str(basename_value)
                                            logger.debug("Constructed imported file path (inferred): %s", full_path)
                                            return str(full_path)
                                        else:
                                            full_path = current / relpath_str
                                            logger.debug("Constructed imported dir path (inferred): %s", full_path)
                                            return str(full_path)
                                    current = current.parent
                                    if current == current.parent:  # Reached root
                                        break

                    # For CCP4_JOBS paths: relPath is ALWAYS relative to PROJECT directory, not workDirectory
                    # Get project directory and construct: project_dir / relPath / baseName
                    project_id = None
                    if hasattr(self, 'project') and self.project is not None:
                        if hasattr(self.project, 'value'):
                            project_id = self.project.value
                        else:
                            project_id = str(self.project)
                    elif hasattr(plugin, '_dbProjectId'):
                        project_id = plugin._dbProjectId

                    if project_id:
                        # Get project directory from database
                        db_handler = self._get_db_handler()
                        if db_handler and hasattr(db_handler, 'getProjectDirectory'):
                            try:
                                project_dir = db_handler.getProjectDirectory(project_id)
                                if project_dir:
                                    base_dir = Path(project_dir)
                                    if basename_value:
                                        full_path = base_dir / relpath_str / str(basename_value)
                                        logger.debug("Constructed path from project + relPath + baseName: %s", full_path)
                                        return str(full_path)
                                    else:
                                        full_path = base_dir / relpath_str
                                        logger.debug("Constructed path from project + relPath: %s", full_path)
                                        return str(full_path)
                            except Exception as e:
                                logger.debug(f"Failed to get project directory: {e}")

                    # Fallback: if we can't get project directory, try using workDirectory
                    # and walking up to find project root
                    if hasattr(plugin, 'workDirectory'):
                        work_dir = Path(plugin.workDirectory)
                        # Walk up to find project root (directory containing CCP4_JOBS)
                        current = work_dir
                        for _ in range(10):
                            if (current / 'CCP4_JOBS').exists():
                                # Found project root
                                if basename_value:
                                    full_path = current / relpath_str / str(basename_value)
                                    logger.debug("Constructed path from inferred project + relPath + baseName: %s", full_path)
                                    return str(full_path)
                            current = current.parent
                            if current == current.parent:  # Reached filesystem root
                                break

        # FINAL FALLBACK: Try to construct path without plugin
        # This handles cases where files were imported/autopopulated but plugin hierarchy
        # isn't available (e.g., after loading from XML during job execution)
        # NOTE: Only returns paths to EXISTING files - safe for both input and output files
        # because the exists() check ensures we don't return phantom paths.
        if relpath_value and basename_value:
            relpath_str = str(relpath_value).strip()
            if relpath_str:
                # Try to find project directory from CCP4I2_PROJECTS_DIR + project UUID
                import os
                if 'CCP4I2_PROJECTS_DIR' in os.environ:
                    projects_dir = Path(os.environ['CCP4I2_PROJECTS_DIR'])
                    # Get project UUID from self.project
                    project_uuid = None
                    if hasattr(self, 'project') and self.project is not None:
                        if hasattr(self.project, 'value'):
                            project_uuid = str(self.project.value).strip()
                        else:
                            project_uuid = str(self.project).strip()

                    if project_uuid:
                        # Find project directory by UUID (might have different case in name)
                        if projects_dir.exists():
                            for project_dir in projects_dir.iterdir():
                                if project_dir.is_dir():
                                    # Check if this directory belongs to our project
                                    # Project dirs typically named like "tmp_XyZ123" or custom names
                                    # We need to match by checking if CCP4_IMPORTED_FILES or similar exists
                                    test_path = project_dir / relpath_str / str(basename_value)
                                    if test_path.exists():
                                        logger.debug(f"Found file via CCP4I2_PROJECTS_DIR search: {test_path}")
                                        return str(test_path)

        # Return baseName as-is (may be relative path, or empty)
        if basename_value is not None:
            return str(basename_value) if basename_value else ""

        # Fallback to legacy file_path
        if hasattr(self, 'file_path') and self.file_path is not None:
            return str(self.file_path)

        return ""

    def exists(self) -> bool:
        """Check if the file exists on disk.

        Returns:
            True if the file exists, False otherwise
        """
        from pathlib import Path
        path = self.getFullPath()
        if path:
            return Path(path).exists()
        return False

    def getExt(self) -> str:
        """Get the file extension including the leading dot.

        Returns:
            File extension (e.g., '.pdb', '.mtz', '.cif') or empty string if no extension

        Example:
            >>> file = CPdbDataFile()
            >>> file.setFullPath('/path/to/model.pdb')
            >>> file.getExt()
            '.pdb'
        """
        from pathlib import Path
        path = self.getFullPath()
        if path:
            return Path(path).suffix
        return ''

    def isSet(self, field_name: str = None, allowUndefined: bool = False,
              allowDefault: bool = False, allSet: bool = True) -> bool:
        """Check if the file has been set.

        For CDataFile, a file is considered "set" if its baseName attribute is set.
        This overrides the base CData.isSet() to provide file-specific semantics.

        Args:
            field_name: If None, checks baseName. Otherwise delegates to parent.
            allowUndefined: If True, allow None/undefined values to be considered "set"
            allowDefault: If False, consider values that equal the default as "not set"
            allSet: For container types - if True, all children must be set

        Returns:
            True if baseName is set (has a non-None, non-empty value)
        """
        if field_name is None:
            # For files, check if baseName is set
            if hasattr(self, 'baseName'):
                if hasattr(self.baseName, 'value'):
                    # baseName is a CData wrapper - check its value
                    return self.baseName.value is not None and self.baseName.value != ""
                else:
                    # baseName is a plain value
                    return self.baseName is not None and self.baseName != ""
            # Fallback to legacy file_path
            if hasattr(self, 'file_path') and self.file_path is not None:
                return self.file_path != ""
            return False
        else:
            # For specific field names, delegate to parent
            return super().isSet(field_name, allowUndefined=allowUndefined,
                               allowDefault=allowDefault, allSet=allSet)

    @property
    def fullPath(self):
        """Property to access the full file path as a string.

        Returns a string representing the full path to the file.
        For backward compatibility with legacy Qt code, returns a QtStringCompat
        object that supports both string operations and .get() method.

        Returns:
            full path string with .get() method for Qt compatibility
        """
        from ..base_object.fundamental_types import CString
        from ..base_object.qt_compat import QtStringCompat
        path_string = self.getFullPath()
        # Return a string-like object with .get() method for Qt API compatibility
        # Always return QtStringCompat for consistent API (even for empty strings)
        return QtStringCompat(path_string) if path_string else QtStringCompat("")

    @fullPath.setter
    def fullPath(self, path: str):
        """Set the full file path.

        Args:
            path: Full file path as a string
        """
        self.setFullPath(path)

    def __str__(self) -> str:
        """Return string representation of the file (its path).

        Returns:
            Full path to the file, or empty string if not set
        """
        return self.getFullPath()

    def get(self) -> dict:
        """Get file attributes as a dict. Compatible with old CCP4i2 API.

        For file objects, this returns the fullPath as the primary value.
        Overrides base CData.get() to provide file-specific behavior.

        Returns:
            Dict with fullPath and other file attributes
        """
        # Get base attributes
        result = super().get()
        # Add fullPath as a convenience
        result['fullPath'] = self.getFullPath()
        return result

    def set(self, value={}, **kw):
        """Set file attributes. Compatible with old CCP4i2 API.

        Args:
            value: Can be:
                - None: Clear the file (unset path)
                - str: File path (calls setFullPath)
                - dict: Dict of attributes to set
                - CDataFile: Another file object to copy from
            **kw: Additional keyword arguments passed to parent
        """
        if value is None:
            # None argument: clear the file path
            self.setFullPath('')
            self.unSet()
            return
        if isinstance(value, str):
            # String argument: set as file path
            self.setFullPath(value)
        elif isinstance(value, CDataFile):
            # Another CDataFile: copy its attributes
            super().set(value.get())
            # IMPORTANT: Also copy non-metadata attributes like 'content'
            # that are set directly in __init__ but not tracked by metadata
            if hasattr(value, 'content') and value.content is not None:
                # Deep-copy the content and re-parent it to self
                # A shallow reference copy (self.content = value.content) doesn't work
                # because the content's parent would still be the source file,
                # and it wouldn't be serialized as a child of this file
                source_content = value.content
                content_type = type(source_content)
                try:
                    # Create a new content object with self as parent
                    new_content = content_type(parent=self, name='fileContent')
                    # Deep copy the data from source to new content
                    if hasattr(new_content, '_deep_copy_from'):
                        new_content._deep_copy_from(source_content)
                    self.content = new_content
                except Exception:
                    # Fallback to reference copy if deep copy fails
                    self.content = source_content
            if hasattr(value, 'file_path'):
                self.file_path = value.file_path
        elif isinstance(value, dict):
            # Dict: handle special cases for file paths
            if 'fullPath' in value:
                self.setFullPath(value['fullPath'])
            else:
                # Regular dict - pass to parent
                super().set(value)
        else:
            # Fallback to parent implementation
            super().set(value)

    def load_from_file(self, file_path: str):
        """Load data from file."""
        self.setFullPath(file_path)
        # TODO: Implement file loading logic

    def save_to_file(self, file_path: str = None):
        """Save data to file."""
        path = file_path or self.getFullPath()
        if not path:
            raise ValueError("No file path specified")
        # TODO: Implement file saving logic

    def loadFile(self, initialise: bool = False):
        """
        Load file content by instantiating fileContentClassName and calling its loadFile().

        This base method:
        1. Gets fileContentClassName from qualifiers
        2. Instantiates that class if not already created
        3. Delegates to content class's loadFile() method
        4. Handles errors and emits dataLoaded signal

        Args:
            initialise: If True, recreate content object even if it exists

        Returns:
            CErrorReport with any errors encountered

        Example:
            >>> mtz_file = CMtzDataFile()
            >>> mtz_file.setFullPath('/path/to/data.mtz')
            >>> error = mtz_file.loadFile()
            >>> if error.count() == 0:
            ...     print(f"Cell: {mtz_file.content.cell}")
        """
        from core.base_object.error_reporting import CErrorReport, CException, SEVERITY_ERROR
        from pathlib import Path

        error = CErrorReport()

        # Check if file loading is blocked
        if hasattr(self, 'BLOCK_LOAD_FILES') and self.BLOCK_LOAD_FILES:
            return error

        # Get content class name from qualifiers
        content_class_name = self.get_qualifier('fileContentClassName')
        if content_class_name is None:
            return error  # No content class specified - not an error

        # Instantiate content object if needed
        if self.content is None or initialise:
            try:
                # Import the content class
                # Try local imports first (CCP4XtalData, CCP4ModelData, etc.)
                content_class = None
                for module_name in ['core.CCP4XtalData', 'core.CCP4ModelData', 'core.CCP4CootData']:
                    try:
                        module = __import__(module_name, fromlist=[content_class_name])
                        content_class = getattr(module, content_class_name, None)
                        if content_class is not None:
                            break
                    except (ImportError, AttributeError):
                        continue

                if content_class is None:
                    error.append(
                        klass=self.__class__.__name__,
                        code=105,
                        details=f"Content class '{content_class_name}' not found",
                        name=self.object_name() if hasattr(self, 'object_name') else ''
                    )
                    return error

                # Create content instance
                self.content = content_class(parent=self, name='fileContent')

            except Exception as e:
                error.append(
                    klass=self.__class__.__name__,
                    code=106,
                    details=f"Error instantiating content class '{content_class_name}': {e}",
                    name=self.object_name() if hasattr(self, 'object_name') else ''
                )
                return error

        # Get file path
        file_path = self.getFullPath()
        if not file_path:
            # No path set - unset content
            if self.content is not None:
                self.content.unSet()
            return error

        # Check file exists
        path_obj = Path(file_path)
        if path_obj.exists() and path_obj.is_file():
            try:
                # Delegate to content class's loadFile()
                content_error = self.content.loadFile(file_path)
                if content_error is not None:
                    error.extend(content_error)

            except CException as e:
                error.extend(e)

            except Exception as e:
                error.append(
                    klass=self.__class__.__name__,
                    code=320,
                    details=f"Unrecognised error loading file '{file_path}': {e}",
                    name=self.object_name() if hasattr(self, 'object_name') else ''
                )
        else:
            # File doesn't exist - unset content
            if self.content is not None:
                self.content.unSet()

        # Emit signal (if event system available)
        if hasattr(self, 'dataLoaded'):
            try:
                self.dataLoaded.emit()
            except:
                pass  # Signal system may not be available

        return error

    def _get_conversion_output_path(
        self,
        target_content_type: str,
        target_extension: Optional[str] = None,
        work_directory: Optional[Any] = None
    ) -> str:
        """
        Calculate output path for converted file (generic for all file types).

        This is a base class method used by as_CONTENTTYPE() conversion methods
        in subclasses. It generates output paths following the pattern:
        {inputroot}_as{CONTENT_TYPE}{extension}

        Naming pattern: {inputroot}_as{CONTENT_TYPE}{extension}
        Location: Input file's directory (if writable), else work_directory

        Args:
            target_content_type: Name of target content type (e.g., 'FMEAN', 'MMCIF')
            target_extension: Optional extension override (e.g., '.cif').
                            If None, uses input file's extension.
            work_directory: Fallback directory if input dir not writable

        Returns:
            Full path to output file

        Examples:
            >>> # MTZ conversion (same extension)
            >>> obs_file._get_conversion_output_path('FMEAN')
            '/data/input_asFMEAN.mtz'

            >>> # PDB to mmCIF (different extension)
            >>> pdb_file._get_conversion_output_path('MMCIF', '.cif')
            '/data/model_asMMCIF.cif'
        """
        from pathlib import Path

        input_path = Path(self.getFullPath())
        input_dir = input_path.parent
        input_stem = input_path.stem  # Filename without extension

        # Use target_extension if provided, otherwise preserve input extension
        extension = target_extension if target_extension else input_path.suffix

        # Calculate output filename (one underscore before 'as', none between 'as' and TYPE)
        output_name = f"{input_stem}_as{target_content_type}{extension}"

        # Determine if source file is within the project hierarchy
        # Priority: 1) Check dbFileId + database, 2) Check path patterns
        is_within_project = False

        # Method 1: Check if file has dbFileId (indicates it's been imported into project)
        if hasattr(self, 'dbFileId') and self.dbFileId is not None:
            db_file_id = None
            if hasattr(self.dbFileId, 'value') and self.dbFileId.value is not None:
                db_file_id = self.dbFileId.value
            elif isinstance(self.dbFileId, str) and self.dbFileId:
                db_file_id = self.dbFileId

            if db_file_id:
                # File has a dbFileId, so it's been imported into the project
                # Try to get its path from the database to confirm
                db_handler = self._get_db_handler()
                if db_handler:
                    try:
                        import uuid
                        file_uuid = uuid.UUID(str(db_file_id))
                        db_path = db_handler.get_file_path_sync(file_uuid)
                        if db_path:
                            # Successfully retrieved path from DB
                            # Check if it's within project directories
                            is_within_project = 'CCP4_JOBS' in str(db_path) or 'CCP4_IMPORTED_FILES' in str(db_path)
                    except Exception:
                        pass  # Fall through to path-based check

        # Method 2: Fallback to path-based detection
        if not is_within_project:
            # Simple path pattern matching as fallback
            is_within_project = 'CCP4_JOBS' in str(input_path) or 'CCP4_IMPORTED_FILES' in str(input_path)

        if is_within_project:
            # Source is within project: try to place derived file adjacent to source
            if input_dir.exists() and input_dir.is_dir():
                output_path = input_dir / output_name
                # Check if we can write there
                try:
                    test_file = input_dir / f".write_test_{id(self)}"
                    test_file.touch()
                    test_file.unlink()
                    return str(output_path)
                except (PermissionError, OSError):
                    pass  # Fall through to work_directory

        # Source is external to project OR input directory not writable:
        # Use work_directory (typically points to CCP4_IMPORTED_FILES or can derive it)
        if work_directory:
            work_path = Path(work_directory)

            # If work_directory is in CCP4_JOBS, derive CCP4_IMPORTED_FILES path
            if 'CCP4_JOBS' in str(work_path):
                # Walk up to find CCP4_JOBS
                while work_path.name.startswith('job_') or work_path.name in ['ctruncate', 'refmac', 'prosmart', 'aimless']:
                    work_path = work_path.parent

                if work_path.name == 'CCP4_JOBS':
                    # Go to project root and use CCP4_IMPORTED_FILES
                    project_root = work_path.parent
                    imported_files_dir = project_root / 'CCP4_IMPORTED_FILES'
                    imported_files_dir.mkdir(parents=True, exist_ok=True)
                    return str(imported_files_dir / output_name)

            # Direct use of work_directory
            work_path.mkdir(parents=True, exist_ok=True)
            return str(work_path / output_name)

        # Last resort: same as input (may fail at write time if not writable)
        return str(input_dir / output_name)

    def setContentFlag(self, content_flag: Optional[int] = None, reset: bool = False):
        """
        Set or auto-detect the content flag for this file.

        For backward compatibility, this method can be called without arguments
        to trigger automatic introspection of the file's content.

        Args:
            content_flag: If provided, sets contentFlag to this value directly.
                         If None, attempts to auto-detect by inspecting the file
                         (behavior depends on subclass implementation).
            reset: If True, forces re-detection of content flag (legacy parameter,
                   currently equivalent to calling without content_flag argument).

        Examples:
            >>> # Explicit assignment
            >>> obs_file.setContentFlag(4)  # Set to FMEAN

            >>> # Auto-detection (for MTZ files)
            >>> obs_file.setContentFlag()  # Introspects file, sets based on columns

            >>> # Force re-detection (legacy API)
            >>> obs_file.setContentFlag(reset=True)
        """
        # If reset=True is passed, treat as auto-detection request
        if reset and content_flag is None:
            content_flag = None  # Trigger auto-detection path below
        if content_flag is not None:
            flag_value = content_flag
        else:
            # Auto-detection: delegate to subclass-specific implementation
            flag_value = self._introspect_content_flag()
            if flag_value is None:
                return  # No detection, don't change contentFlag

        # Set the value, handling both CData wrappers and plain values
        if hasattr(self, 'contentFlag') and hasattr(self.contentFlag, 'value'):
            # contentFlag is a CData wrapper - set its .value attribute
            self.contentFlag.value = flag_value
        else:
            # contentFlag doesn't exist or isn't a proper CData wrapper
            # Use simple assignment and hope smart assignment works
            # (If contentFlag is CInt with _is_value_type(), smart assign will work)
            # (If contentFlag is generic CData without _is_value_type(), it will be replaced)
            # TODO: Fix metadata system to create CInt instead of CData
            self.contentFlag = flag_value

    def _introspect_content_flag(self) -> Optional[int]:
        """
        Auto-detect the content flag by inspecting the file.

        This is the base implementation that does nothing. Subclasses
        (like CMiniMtzDataFile) override this to provide file-type-specific
        introspection logic.

        Returns:
            Detected content flag value, or None if cannot be determined
        """
        # Base implementation: no introspection capability
        return None

    def getFileContent(self):
        """
        Load and return the file content object.

        This method:
        1. Calls loadFile() to ensure the file is loaded
        2. Returns the fileContent property (the CDataFileContent instance)

        Returns:
            The loaded CDataFileContent object (e.g., CPdbData for CPdbDataFile)

        Example:
            >>> pdb_file = CPdbDataFile()
            >>> pdb_file.setFullPath("/path/to/structure.pdb")
            >>> content = pdb_file.getFileContent()  # Returns CPdbData instance
            >>> content.sequences  # Access sequences from the loaded content
        """
        self.loadFile()
        return self.fileContent

    def checksum(self) -> str:
        """
        Compute MD5 checksum of the file.

        Returns:
            32-character hexadecimal MD5 checksum string, or empty string if file doesn't exist.
        """
        import hashlib
        from pathlib import Path

        file_path = self.getFullPath()
        if not file_path or not Path(file_path).exists():
            return ""

        md5_hash = hashlib.md5()
        with open(file_path, "rb") as f:
            # Read in chunks to handle large files efficiently
            for chunk in iter(lambda: f.read(4096), b""):
                md5_hash.update(chunk)

        return md5_hash.hexdigest()

    def stripedName(self):
        """
        Get the base filename without path or extension.

        Legacy compatibility method for getting a simple filename for display.

        Returns:
            str: Base filename without directory path or extension

        Example:
            >>> file_obj.setFullPath('/path/to/myfile.mtz')
            >>> file_obj.stripedName()
            'myfile'
        """
        file_path = self.getFullPath()
        if not file_path:
            return ""

        from pathlib import Path
        # Get filename without directory path and extension
        return Path(file_path).stem

    def to_metadata_dict(self) -> dict:
        """
        Extract complete metadata as dict for DB storage/serialization.

        Replaces server/ccp4x/lib/cdata_utils.extract_file_metadata() as a core method.

        This method collects all file metadata including:
        - File paths (fullPath, baseName, relPath)
        - Database identifiers (dbFileId, project)
        - Content information (contentFlag, subType, annotation)
        - File type metadata (mimeTypeName, fileExtensions)

        Returns:
            Dict with all file metadata suitable for database storage or JSON serialization

        Example:
            >>> pdb_file = CPdbDataFile()
            >>> pdb_file.setFullPath('/path/to/structure.pdb')
            >>> metadata = pdb_file.to_metadata_dict()
            >>> metadata['fullPath']
            '/path/to/structure.pdb'
            >>> metadata['mimeTypeName']
            'application/pdb'
        """
        from typing import Any, Dict
        from pathlib import Path

        metadata: Dict[str, Any] = {}

        # Basic file paths
        full_path = self.getFullPath()
        if full_path:
            metadata['fullPath'] = full_path
            metadata['exists'] = Path(full_path).exists()

        # Extract attribute values (baseName, relPath, etc.)
        for attr_name in ['baseName', 'relPath', 'annotation', 'subType', 'contentFlag', 'dbFileId', 'project']:
            if hasattr(self, attr_name):
                attr_obj = getattr(self, attr_name)
                if attr_obj is not None:
                    # Handle CData wrapped values
                    if hasattr(attr_obj, 'value'):
                        value = attr_obj.value
                    else:
                        value = attr_obj

                    # Only include if value is set and not None
                    if value is not None:
                        # Convert UUIDs to strings for JSON serialization
                        if attr_name in ['dbFileId', 'project']:
                            metadata[attr_name] = str(value)
                        else:
                            metadata[attr_name] = value

        # Class metadata from qualifiers
        qualifier_keys = ['mimeTypeName', 'mimeTypeDescription', 'fileLabel',
                         'fileExtensions', 'fileContentClassName', 'isDirectory']
        for key in qualifier_keys:
            value = self.get_qualifier(key)
            if value is not None:
                metadata[key] = value

        # Object metadata
        if hasattr(self, 'object_name') and callable(self.object_name):
            metadata['objectName'] = self.object_name()
        elif hasattr(self, 'name'):
            metadata['objectName'] = self.name

        # Class name
        metadata['className'] = self.__class__.__name__

        return metadata

    def validity(self):
        """Validate this file object and return an error report.

        Checks the following conditions:
        - If allowUndefined=False and no file path is set → ERROR
        - If mustExist=True and file path is set but file doesn't exist → ERROR
        - If file path is not set (even when allowUndefined=True) → WARNING
          This helps draw the user's eye to unset optional files in the GUI

        This validation is part of the CData hierarchy - CDataFile knows its
        own validation rules (allowUndefined, mustExist qualifiers) and reports
        issues via CErrorReport.

        Returns:
            CErrorReport containing any validation errors/warnings
        """
        from .error_reporting import CErrorReport, SEVERITY_ERROR, SEVERITY_WARNING

        report = CErrorReport()

        # Get file path - may be None/empty if not set
        full_path = self.getFullPath()
        has_path = bool(full_path and str(full_path).strip())

        # Get qualifiers for this file object
        allow_undefined = self.get_qualifier('allowUndefined')
        must_exist = self.get_qualifier('mustExist')

        # Check allowUndefined
        # If allowUndefined is False (or not set), the file path is required
        if allow_undefined is False and not has_path:
            # Required file not set
            obj_path = self.object_path()
            report.append(
                klass=self.__class__.__name__,
                code=101,
                details=f'Required input file not set: {obj_path}',
                name=obj_path,
                severity=SEVERITY_ERROR
            )
        elif allow_undefined is True and not has_path:
            # Optional file not set - add WARNING to draw user's attention in GUI
            obj_path = self.object_path()
            report.append(
                klass=self.__class__.__name__,
                code=103,
                details=f'Optional input file not set: {obj_path}',
                name=obj_path,
                severity=SEVERITY_WARNING
            )

        # Check mustExist
        # If mustExist=True and path is set, file must exist
        if must_exist and has_path:
            from pathlib import Path
            if not Path(full_path).exists():
                obj_path = self.object_path()
                report.append(
                    klass=self.__class__.__name__,
                    code=102,
                    details=f'File does not exist: {full_path}',
                    name=obj_path,
                    severity=SEVERITY_ERROR
                )

        return report
