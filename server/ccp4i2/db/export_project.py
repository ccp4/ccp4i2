import logging
import zipfile
import shutil
from pathlib import Path
from xml.etree import ElementTree as ET
from xml.dom import minidom
from datetime import datetime
from typing import List, Dict, Optional, Any, TypedDict, Set

from .models import (
    Project,
    Job,
    File,
    FileUse,
    FileImport,
    JobFloatValue,
    JobCharValue,
    ProjectTag,
    JobValueKey,
    FileType,
)
from .ccp4i2_static_data import FILETYPELIST, KEYTYPELIST

logger = logging.getLogger(f"ccp4i2:{__name__}")


def generate_project_xml_tree(
    project: Project, job_selection: Set[str] = None
) -> ET.Element:
    """
    Generate an XML ElementTree for a CCP4 project.

    Args:
        project (Project): The Django Project model instance to export
        job_selection (Set[str], optional): Set of job numbers to export. If None, exports all jobs.

    Returns:
        ET.Element: The root XML element containing all project data
    """
    root = ET.Element("database")

    # Add header information
    header = ET.SubElement(root, "ccp4i2_header")
    project_id_elem = ET.SubElement(header, "projectId")
    project_id_elem.text = _format_uuid_for_xml(project.uuid)

    export_time = ET.SubElement(header, "exportTime")
    export_time.text = str(int(datetime.now().timestamp()))

    # Add body containing all data
    body = ET.SubElement(root, "ccp4i2_body")

    # Export project table (always export full project info)
    _export_project_table(body, project)

    # Export job table (filtered by selection)
    _export_job_table(body, project, job_selection)

    # Export file table (filtered by job selection)
    _export_file_table(body, project, job_selection)

    # Export file use table (filtered by job selection)
    _export_file_use_table(body, project, job_selection)

    # Export file import table (filtered by job selection)
    _export_file_import_table(body, project, job_selection)

    # Export job key value tables (filtered by job selection)
    _export_job_key_value_tables(body, project, job_selection)

    # Export tag tables (always export full tags)
    _export_tag_tables(body, project)

    return root


def write_xml_tree_to_file(root: ET.Element, output_path: Path) -> Path:
    """
    Write an XML ElementTree to a formatted file.

    Args:
        root (ET.Element): The root XML element to write
        output_path (Path): Path where the XML file will be saved

    Returns:
        Path: The path to the created XML file
    """
    # Write formatted XML to file
    xml_string = ET.tostring(root, encoding="unicode")
    pretty_xml = minidom.parseString(xml_string).toprettyxml(indent="  ")

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(pretty_xml)

    logger.info(f"XML tree written to {output_path}")
    return output_path


def export_project_to_xml(project: Project, output_path: Path) -> Path:
    """
    Export a CCP4 project to XML format.

    Args:
        project (Project): The Django Project model instance to export
        output_path (Path): Path where the XML file will be saved

    Returns:
        Path: The path to the created XML file
    """
    # Generate the XML tree
    root = generate_project_xml_tree(project)

    # Write to file
    write_xml_tree_to_file(root, output_path)

    logger.info(f"Exported project {project.name} to {output_path}")
    return output_path


def export_project_to_zip(
    project: Project, output_path: Path, job_selection: Set[str] = None
) -> Path:
    """
    Export a CCP4 project to a ZIP archive containing XML and project files.

    Args:
        project (Project): The Django Project model instance to export
        output_path (Path): Path where the ZIP file will be saved
        job_selection (Set[str], optional): Set of job numbers to export. If None, exports all jobs.

    Returns:
        Path: The path to the created ZIP file
    """
    # Create temporary directory for XML
    temp_dir = Path(output_path).parent / f"temp_{project.uuid}"
    temp_dir.mkdir(exist_ok=True)

    try:
        # Generate XML tree and write to temporary file
        xml_path = temp_dir / "DATABASE.db.xml"
        root = generate_project_xml_tree(project, job_selection=job_selection)
        write_xml_tree_to_file(root, xml_path)

        # Create ZIP archive
        with zipfile.ZipFile(output_path, "w", zipfile.ZIP_DEFLATED) as zip_archive:
            # Add XML file
            zip_archive.write(xml_path, "DATABASE.db.xml")

            # Add project directories and files
            _add_project_files_to_zip(zip_archive, project, job_selection)

        logger.info(f"Exported project {project.name} to ZIP: {output_path}")
        return output_path

    finally:
        # Clean up temporary directory
        if temp_dir.exists():
            shutil.rmtree(temp_dir)


def _format_uuid_for_xml(uuid_value) -> str:
    """Remove hyphens from UUID for XML serialization."""
    return str(uuid_value).replace("-", "")


def _export_project_table(body: ET.Element, project: Project) -> None:
    """Export project table data."""
    project_table = ET.SubElement(body, "projectTable")
    project_elem = ET.SubElement(project_table, "project")

    project_elem.set("projectid", _format_uuid_for_xml(project.uuid))
    project_elem.set("projectname", project.name)
    project_elem.set("lastjobnumber", str(project.last_job_number))
    project_elem.set("projectdirectory", str(project.directory))
    project_elem.set("projectcreated", str(int(project.creation_time.timestamp())))


def _get_jobs_from_selection(project: Project, job_selection: Set[str]) -> Set[Job]:
    """
    Get Job objects from job number selection, including descendants and jobs
    that have files used by the selected jobs.

    Args:
        project (Project): The project containing the jobs
        job_selection (Set[str]): Set of job numbers as strings

    Returns:
        Set[Job]: Set of Job objects including descendants and file dependency jobs
    """
    if job_selection is None:
        return None

    # Convert string numbers to integers for lookup
    job_numbers = []
    for job_num_str in job_selection:
        try:
            job_numbers.append(int(job_num_str))
        except ValueError:
            logger.warning(f"Invalid job number format: {job_num_str}")
            continue

    # Get the specified top-level jobs
    top_level_jobs = set(Job.objects.filter(project=project, number__in=job_numbers))

    if not top_level_jobs:
        return set()

    # Collect all descendants
    selected_jobs = set(top_level_jobs)
    for job in top_level_jobs:
        selected_jobs.update(_get_job_descendants(job))

    # Add jobs that own files used by any of the selected jobs
    file_dependency_jobs = _get_jobs_with_files_used_by_jobs(selected_jobs)
    selected_jobs.update(file_dependency_jobs)

    return selected_jobs


def _get_jobs_with_files_used_by_jobs(jobs: Set[Job]) -> Set[Job]:
    """
    Get jobs that own files used by the given jobs through FileUse relationships.

    Args:
        jobs (Set[Job]): Set of jobs to check for file usage

    Returns:
        Set[Job]: Set of jobs that own files used by the input jobs
    """
    file_owner_jobs = set()

    # Get all FileUse entries for the selected jobs
    job_uuids = [job.uuid for job in jobs]
    file_uses = FileUse.objects.filter(job__uuid__in=job_uuids).select_related(
        "file__job"
    )

    # Collect the jobs that own the used files
    for file_use in file_uses:
        file_owner_job = file_use.file.job
        file_owner_jobs.add(file_owner_job)

    return file_owner_jobs


def _export_job_table(
    body: ET.Element, project: Project, job_selection: Set[str] = None
) -> None:
    """Export job table data."""
    job_table = ET.SubElement(body, "jobTable")

    if job_selection is not None:
        jobs = _get_jobs_from_selection(project, job_selection)
        if not jobs:
            return  # No valid jobs found
    else:
        jobs = Job.objects.filter(project=project)

    # Sort jobs by number for consistent output
    sorted_jobs = sorted(jobs, key=lambda j: j.number)

    for job in sorted_jobs:
        job_elem = ET.SubElement(job_table, "job")

        job_elem.set("jobid", _format_uuid_for_xml(job.uuid))
        job_elem.set("projectid", _format_uuid_for_xml(project.uuid))
        job_elem.set("status", str(job.status))
        job_elem.set("jobnumber", str(job.number))
        job_elem.set("taskname", job.task_name or "")
        job_elem.set("creationtime", str(int(job.creation_time.timestamp())))
        job_elem.set("title", job.title or job.task_name or "")

        if job.evaluation is not None:
            job_elem.set("evaluation", str(job.evaluation))

        if job.finish_time:
            job_elem.set("finishtime", str(int(job.finish_time.timestamp())))

        if job.parent:
            job_elem.set("parentjobid", _format_uuid_for_xml(job.parent.uuid))


def _export_file_table(
    body: ET.Element, project: Project, job_selection: Set[str] = None
) -> None:
    """Export file table data."""
    file_table = ET.SubElement(body, "fileTable")

    if job_selection is not None:
        jobs = _get_jobs_from_selection(project, job_selection)
        if not jobs:
            return  # No valid jobs found
        # Filter files by selected jobs (now includes file dependency jobs)
        job_uuids = [job.uuid for job in jobs]
        files = File.objects.filter(job__uuid__in=job_uuids)
    else:
        files = File.objects.filter(job__project=project)

    for file_obj in files:
        file_elem = ET.SubElement(file_table, "file")

        file_elem.set("fileid", _format_uuid_for_xml(file_obj.uuid))
        file_elem.set("jobid", _format_uuid_for_xml(file_obj.job.uuid))
        file_elem.set("filename", file_obj.objectName() or "")

        # Use the path property to get directory information
        if file_obj.path:
            file_elem.set("pathflag", str(file_obj.directory or ""))
        else:
            file_elem.set("pathflag", "")

        file_elem.set("jobparamname", file_obj.job_param_name or "")

        # Find file type ID from static data
        file_type_id = _get_file_type_id(file_obj.type.name if file_obj.type else "")
        file_elem.set("filetypeid", str(file_type_id))

        if file_obj.annotation:
            file_elem.set("annotation", file_obj.annotation)
        if file_obj.sub_type:
            file_elem.set("filesubtype", str(file_obj.sub_type))
        if file_obj.content:
            file_elem.set("filecontent", str(file_obj.content))


def _export_file_use_table(
    body: ET.Element, project: Project, job_selection: Set[str] = None
) -> None:
    """Export file use table data, including cross-job file usage."""
    fileuse_table = ET.SubElement(body, "fileuseTable")

    if job_selection is not None:
        # Get the originally selected jobs (top-level + descendants)
        selected_jobs = _get_jobs_from_selection_without_file_dependencies(
            project, job_selection
        )
        if not selected_jobs:
            return  # No valid jobs found

        # Get FileUse entries for the selected jobs
        job_uuids = [job.uuid for job in selected_jobs]
        file_uses = FileUse.objects.filter(job__uuid__in=job_uuids)

        # Also get FileUse entries where selected jobs use files from other jobs
        # This is already covered by the above query since we're filtering by the using job
    else:
        file_uses = FileUse.objects.filter(job__project=project)

    for file_use in file_uses:
        fileuse_elem = ET.SubElement(fileuse_table, "fileuse")

        fileuse_elem.set("fileid", _format_uuid_for_xml(file_use.file.uuid))
        fileuse_elem.set("jobid", _format_uuid_for_xml(file_use.job.uuid))
        fileuse_elem.set("roleid", str(file_use.role))
        fileuse_elem.set("jobparamname", file_use.job_param_name or "")


def _get_jobs_from_selection_without_file_dependencies(
    project: Project, job_selection: Set[str]
) -> Set[Job]:
    """
    Get Job objects from job number selection, including descendants but NOT file dependency jobs.
    This is used when we want the original selection for FileUse filtering.
    """
    if job_selection is None:
        return None

    # Convert string numbers to integers for lookup
    job_numbers = []
    for job_num_str in job_selection:
        try:
            job_numbers.append(int(job_num_str))
        except ValueError:
            logger.warning(f"Invalid job number format: {job_num_str}")
            continue

    # Get the specified top-level jobs
    top_level_jobs = set(Job.objects.filter(project=project, number__in=job_numbers))

    if not top_level_jobs:
        return set()

    # Collect all descendants (but not file dependency jobs)
    selected_jobs = set(top_level_jobs)
    for job in top_level_jobs:
        selected_jobs.update(_get_job_descendants(job))

    return selected_jobs


def _export_file_import_table(
    body: ET.Element, project: Project, job_selection: Set[str] = None
) -> None:
    """Export file import table data."""
    import_table = ET.SubElement(body, "importfileTable")

    if job_selection is not None:
        jobs = _get_jobs_from_selection(project, job_selection)
        if not jobs:
            return  # No valid jobs found
        # Filter file imports by selected jobs
        job_uuids = [job.uuid for job in jobs]
        file_imports = FileImport.objects.filter(file__job__uuid__in=job_uuids)
    else:
        file_imports = FileImport.objects.filter(file__job__project=project)

    for file_import in file_imports:
        import_elem = ET.SubElement(import_table, "importfile")

        import_elem.set("fileid", _format_uuid_for_xml(file_import.file.uuid))
        import_elem.set("creationtime", str(int(file_import.time.timestamp())))
        import_elem.set("sourcefilename", file_import.name or "")
        import_elem.set("checksum", str(file_import.checksum) or "")


def _export_job_key_value_tables(
    body: ET.Element, project: Project, job_selection: Set[str] = None
) -> None:
    """Export job key value tables (both float and char values)."""
    keyvalue_table = ET.SubElement(body, "jobkeyvalueTable")

    if job_selection is not None:
        jobs = _get_jobs_from_selection(project, job_selection)
        if not jobs:
            return  # No valid jobs found
        job_uuids = [job.uuid for job in jobs]
        float_values = JobFloatValue.objects.filter(job__uuid__in=job_uuids)
        char_values = JobCharValue.objects.filter(job__uuid__in=job_uuids)
    else:
        float_values = JobFloatValue.objects.filter(job__project=project)
        char_values = JobCharValue.objects.filter(job__project=project)

    # Export float values
    for float_value in float_values:
        keyvalue_elem = ET.SubElement(keyvalue_table, "jobkeyvalue")

        keyvalue_elem.set("jobid", _format_uuid_for_xml(float_value.job.uuid))
        key_type_id = str(
            _get_key_type_id(float_value.key.name if float_value.key else "")
        )
        keyvalue_elem.set("keytypeid", str(key_type_id))
        keyvalue_elem.set("value", str(float_value.value))

    # Export char values
    for char_value in char_values:
        keyvalue_elem = ET.SubElement(keyvalue_table, "jobkeycharvalue")

        keyvalue_elem.set("jobid", _format_uuid_for_xml(char_value.job.uuid))
        key_type_id = str(
            _get_key_type_id(char_value.key.name if char_value.key else "")
        )
        keyvalue_elem.set("keytypeid", str(key_type_id))
        keyvalue_elem.set("value", str(char_value.value))


def _export_tag_tables(body: ET.Element, project: Project) -> None:
    """Export tag and project tag tables."""
    # Get all tags associated with this project
    project_tags = ProjectTag.objects.filter(projects=project)

    if project_tags.exists():
        # Export tag table
        tag_table = ET.SubElement(body, "tagTable")
        tag_id_map = {}

        for i, project_tag in enumerate(project_tags):
            tag_elem = ET.SubElement(tag_table, "tag")
            tag_id = i + 1  # Simple sequential ID
            tag_id_map[project_tag.id] = tag_id

            tag_elem.set("tagid", str(tag_id))
            tag_elem.set("text", project_tag.text or "")

        # Export project tag table
        projecttag_table = ET.SubElement(body, "projecttagTable")

        for project_tag in project_tags:
            projecttag_elem = ET.SubElement(projecttag_table, "projecttag")

            projecttag_elem.set("projectid", _format_uuid_for_xml(project.uuid))
            projecttag_elem.set("tagid", str(tag_id_map[project_tag.id]))


def _get_file_type_id(file_type_name: str) -> int:
    """Get file type ID from static data."""
    if not file_type_name:
        return 0
    for file_type_id, name, description in FILETYPELIST:
        if name == file_type_name:
            return file_type_id
    return 0  # Default if not found


def _get_key_type_id(key_name: str) -> int:
    """Get key type ID from static data."""
    if not key_name:
        return 0
    for key_type_id, name, description in KEYTYPELIST:
        if name == key_name:
            return key_type_id
    return 0  # Default if not found


def _add_directory_to_zip(
    zip_archive: zipfile.ZipFile, source_dir: Path, archive_dir: str
) -> None:
    """Recursively add directory contents to ZIP archive under the given archive_dir, ensuring all parent directories are present."""
    # Ensure parent directories are added as empty entries
    parts = Path(archive_dir).parts
    for i in range(1, len(parts) + 1):
        parent_dir = Path(*parts[:i])
        zip_info = zipfile.ZipInfo(str(parent_dir) + "/")
        # Only add if not already present
        if zip_info.filename not in zip_archive.namelist():
            zip_archive.writestr(zip_info, "")

    for item in source_dir.rglob("*"):
        if item.is_file():
            # Calculate relative path within the archive under archive_dir
            relative_path = Path(archive_dir) / item.relative_to(source_dir)
            zip_archive.write(item, str(relative_path))
        elif item.is_dir():
            # Add empty directories under archive_dir
            relative_path = Path(archive_dir) / item.relative_to(source_dir)
            zip_info = zipfile.ZipInfo(str(relative_path) + "/")
            if zip_info.filename not in zip_archive.namelist():
                zip_archive.writestr(zip_info, "")


def _add_project_files_to_zip(
    zip_archive: zipfile.ZipFile, project: Project, job_selection: Set[str] = None
) -> None:
    """Add project directories and files to ZIP archive, filtered by job selection."""
    project_dir = Path(project.directory)
    if not project_dir.exists():
        return

    # Get the set of files to export
    if job_selection is not None:
        # Get ALL jobs (including file dependency jobs) for file export
        jobs = _get_jobs_from_selection(project, job_selection)
        if not jobs:
            return  # No valid jobs found
        exported_files = _get_files_for_jobs(jobs)

        # Only get directories for TOP-LEVEL jobs from the original selection
        top_level_jobs = _get_top_level_jobs_from_selection(project, job_selection)
        job_directories = _get_job_directories_from_job_objects(
            top_level_jobs, project_dir
        )

        # Also add directories for file dependency jobs that are outside the descendant tree
        file_dependency_jobs = _get_jobs_with_files_used_by_jobs(
            _get_jobs_from_selection_without_file_dependencies(project, job_selection)
        )
        dependency_job_directories = _get_job_directories_from_job_objects(
            file_dependency_jobs, project_dir
        )
        job_directories.update(dependency_job_directories)

    else:
        # Export all files
        exported_files = File.objects.filter(job__project=project)
        job_directories = set()
        ccp4_jobs_dir = project_dir / "CCP4_JOBS"
        if ccp4_jobs_dir.exists():
            job_directories = {d for d in ccp4_jobs_dir.iterdir() if d.is_dir()}

    # Standard project directories to include
    standard_dirs = [
        "CCP4_COOT",
        "CCP4_DOWNLOADED_FILES",
        "CCP4_PROJECT_FILES",
        "CCP4_TMP",
    ]
    # Special handling for CCP4_IMPORTED_FILES: only add empty placeholder directory
    imported_files_dir = project_dir / "CCP4_IMPORTED_FILES"
    if imported_files_dir.exists():
        zip_info = zipfile.ZipInfo("CCP4_IMPORTED_FILES/")
        if zip_info.filename not in zip_archive.namelist():
            zip_archive.writestr(zip_info, "")
        # Do NOT add files from this directory here; files will be added individually below

    # Add standard directories
    for subdir_name in standard_dirs:
        subdir_path = project_dir / subdir_name
        if subdir_path.exists():
            _add_directory_to_zip(zip_archive, subdir_path, subdir_name)

    # Add selected job directories and file dependency job directories
    for job_dir in job_directories:
        relative_path = job_dir.relative_to(project_dir)
        _add_directory_to_zip(zip_archive, job_dir, str(relative_path))

    # Add individual files referenced by exported File objects using their path property
    file_paths_added = set()  # Track to avoid duplicates

    for file_obj in exported_files:
        if file_obj.path:
            file_path = Path(file_obj.path)
            if (
                file_path.exists()
                and file_path.is_file()
                and file_path not in file_paths_added
            ):
                # Calculate relative path from project directory
                try:
                    relative_path = file_path.relative_to(project_dir)
                    # Only add the file if it hasn't already been added at this relative path
                    if str(relative_path) not in zip_archive.namelist():
                        zip_archive.write(file_path, str(relative_path))
                    file_paths_added.add(file_path)
                except ValueError:
                    # File is outside project directory, skip or handle as needed
                    logger.warning(
                        f"File {file_path} is outside project directory, skipping"
                    )


def _get_top_level_jobs_from_selection(
    project: Project, job_selection: Set[str]
) -> Set[Job]:
    """
    Get only the top-level Job objects from job number selection (no descendants).

    Args:
        project (Project): The project containing the jobs
        job_selection (Set[str]): Set of job numbers as strings

    Returns:
        Set[Job]: Set of top-level Job objects only
    """
    if job_selection is None:
        return None

    # Convert string numbers to integers for lookup
    job_numbers = []
    for job_num_str in job_selection:
        try:
            job_numbers.append(int(job_num_str))
        except ValueError:
            logger.warning(f"Invalid job number format: {job_num_str}")
            continue

    # Get ONLY the specified top-level jobs (no descendants)
    top_level_jobs = set(Job.objects.filter(project=project, number__in=job_numbers))

    return top_level_jobs


def _get_files_for_jobs(jobs: Set[Job]) -> List[File]:
    """Get all File objects associated with the selected jobs."""
    job_uuids = [job.uuid for job in jobs]
    return list(File.objects.filter(job__uuid__in=job_uuids))


def _get_job_directories_from_job_objects(
    jobs: Set[Job], project_dir: Path
) -> Set[Path]:
    """Get job directories directly from Job objects using their directory property."""
    job_directories = set()

    for job in jobs:
        if job.directory:
            job_dir_path = Path(job.directory)
            # Verify the directory exists and is within the project directory
            if job_dir_path.exists() and job_dir_path.is_dir():
                try:
                    # Ensure it's within the project directory
                    job_dir_path.relative_to(project_dir)
                    job_directories.add(job_dir_path)
                except ValueError:
                    # Directory is outside project directory
                    logger.warning(
                        f"Job {job.number} directory {job_dir_path} is outside project directory, skipping"
                    )

    return job_directories


def _get_job_descendants(job: Job) -> Set[Job]:
    """Recursively get all descendant jobs."""
    descendants = set()

    # Get direct children
    children = Job.objects.filter(parent=job)
    for child in children:
        descendants.add(child)
        # Recursively get descendants of children
        descendants.update(_get_job_descendants(child))

    return descendants


# Convenience functions for different export formats
def export_project_xml(project_id: str, output_path: Path) -> Path:
    """Export project by ID to XML format."""
    project = Project.objects.get(uuid=project_id)
    return export_project_to_xml(project, output_path)


def export_project_zip(project_id: str, output_path: Path) -> Path:
    """Export project by ID to ZIP format."""
    project = Project.objects.get(uuid=project_id)
    return export_project_to_zip(project, output_path)


# Example usage
if __name__ == "__main__":
    # Example: Export project to XML
    try:
        project = Project.objects.get(name="my_project")
        xml_output = Path("exported_project.xml")
        zip_output = Path("exported_project.zip")

        export_project_to_xml(project, xml_output)
        export_project_to_zip(project, zip_output)

        print(f"Project exported to {xml_output} and {zip_output}")

    except Project.DoesNotExist:
        print("Project not found")
    except Exception as e:
        print(f"Export failed: {e}")
