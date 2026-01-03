import logging

"""
This module provides functions to import CCP4 project data from XML and ZIP files into a Django-based database.

Functions:
    job_number_hash(dotted_number: str) -> str:
        Generates a hash for a job number by padding each element and concatenating them.

    import_ccp4_project_zip(zip_path: Path, relocate_path: Path = None):
        Imports a CCP4 project from a ZIP file, extracting files and handling job remapping.

    renumber_top_job(job_node: ET.Element, root_node: ET.Element) -> str:
        Renumbers top-level jobs to avoid conflicts with existing jobs in the database.

    import_i2xml_from_file(xml_path: Path, relocate_path: Path = None):
        Imports CCP4 project data from an XML file.

    import_i2xml(root_node: ET.Element, relocate_path: Path) -> dict:
        Imports CCP4 project data from an XML root node and returns a job map.

    import_project(node: ET.Element, relocate_path: Path = None):
        Imports a project from an XML node, creating or updating the project in the database.

    import_job(node: ET.Element):
        Imports a job from an XML node, creating or updating the job in the database.

    import_file(node: ET.Element):
        Imports a file from an XML node, creating or updating the file in the database.

    import_file_use(node: ET.Element):
        Imports file usage information from an XML node, creating or updating the file usage in the database.

    import_file_import(node: ET.Element):
        Imports file import information from an XML node, creating or updating the file import in the database.

    import_job_key_value(node: ET.Element):
        Imports job key-value pairs from an XML node, creating or updating the job float value in the database.

    import_job_key_char_value(node: ET.Element):
        Imports job key-character value pairs from an XML node, creating or updating the job char value in the database.

    import_tag(node: ET.Element):
        Imports a tag from an XML node and stores it in a tag map.

    import_project_tag(node: ET.Element):
        Imports a project tag from an XML node, creating or updating the project tag in the database.
"""
import datetime
import zipfile

from pathlib import Path
from xml.etree import ElementTree as ET
from ..api.serializers import (
    ProjectSerializer,
    JobSerializer,
    FileSerializer,
    FileUseSerializer,
    FileImportSerializer,
    JobCharValueSerializer,
    JobFloatValueSerializer,
    ProjectTagSerializer,
)
from .models import (
    Project,
    Job,
    FileType,
    File,
    JobValueKey,
    ProjectTag,
    FileUse,
    FileImport,
    JobFloatValue,
    JobCharValue,
)
from .ccp4i2_static_data import FILETYPELIST, KEYTYPELIST

logger = logging.getLogger(f"ccp4i2:{__name__}")

tag_map = {}


def job_number_hash(dotted_number: str):
    job_elements = dotted_number.split(".")
    # left_pad each element.  NB: this limits the job number count at every level!
    job_elements = [job_element.ljust(8, "0") for job_element in job_elements]
    # left_pad agregate element.  NB: this limits the number of job levels!
    return "".join(job_elements).ljust(32 * 8, "0")


def import_ccp4_project_zip(zip_path: Path, relocate_path: Path = None):
    """
    Imports a CCP4 project from a zip archive.
    This function extracts the contents of a CCP4 project zip file into the appropriate
    directories and handles any necessary remapping of job numbers.
    Args:
        zip_path (Path): The path to the zip file containing the CCP4 project.
        relocate_path (Path, optional): The path to relocate the project files. Defaults to None.
    Raises:
        Project.DoesNotExist: If the project specified in the XML does not exist in the database.
    """

    with zipfile.ZipFile(zip_path, "r") as zip_archive:
        with zip_archive.open("DATABASE.db.xml", "r") as database_file:
            root_node = ET.parse(database_file).getroot()
            import_i2xml_result = import_i2xml(root_node, relocate_path=relocate_path)
            # print(import_i2xml_result)
            all_archive_files = zip_archive.namelist()
            this_project_node = root_node.findall("ccp4i2_header/projectId")
            this_project = Project.objects.get(uuid=this_project_node[0].text.strip())
            for subdir in [
                "CCP4_COOT",
                "CCP4_DOWNLOADED_FILES",
                "CCP4_PROJECT_FILES",
                "CCP4_IMPORTED_FILES",
                "CCP4_TMP",
            ]:
                subdir_files = [
                    filename
                    for filename in all_archive_files
                    if filename.startswith(f"{subdir}/")
                ]
                zip_archive.extractall(
                    str(Path(this_project.directory)), subdir_files, None
                )

            # Special handling for JOBS, since there may have been job remapping on import
            # Same could apply to imported files...
            top_level_job_numbers = [
                job_number
                for job_number in import_i2xml_result["job_map"]
                if "." not in job_number
            ]
            (Path(this_project.directory) / "CCP4_JOBS").mkdir(exist_ok=True)
            for top_level_job_number in top_level_job_numbers:
                job_files = [
                    zip_entry
                    for zip_entry in all_archive_files
                    if zip_entry.startswith(f"CCP4_JOBS/job_{top_level_job_number}/")
                ]
                for src in job_files:
                    new_job_number = import_i2xml_result["job_map"][
                        top_level_job_number
                    ]
                    destination = Path(this_project.directory) / src.replace(
                        f"CCP4_JOBS/job_{top_level_job_number}/",
                        f"CCP4_JOBS/job_{new_job_number}/",
                        1,
                    )
                    if zip_archive.getinfo(src).is_dir():
                        destination.mkdir(exist_ok=True)
                    else:
                        with zip_archive.open(src, "r") as src_file:
                            with open(destination, "wb") as destination_file:
                                destination_file.write(src_file.read())


def renumber_top_job(job_node: ET.Element, root_node: ET.Element):
    """
    Renumber the top-level job node if necessary to avoid conflicts with existing job numbers.
    Args:
        job_node (ET.Element): The job node to be renumbered.
        root_node (ET.Element): The root node containing all job nodes.
    Returns:
        str: The (possibly new) job number of the job node.
    Notes:
        - If the job node is not a top-level job (i.e., its job number contains a dot),
          the function returns the original job number without changes.
        - If the job number already exists in the database and belongs to a different job,
          the job node is assigned the next available job number.
        - The function also updates the job numbers of all descendant job nodes to reflect
          the new top-level job number.
    """

    # null operation if not top evel job
    if "." in job_node.attrib["jobnumber"]:
        return job_node.attrib["jobnumber"]
    # Look out for case that jobs with matching job numbers are proposed for import
    original_job_number = job_node.attrib["jobnumber"]
    all_import_job_nodes = root_node.findall("ccp4i2_body/jobTable/job")
    top_level_job_nodes = [
        node
        for node in all_import_job_nodes
        if (("parentjobid" not in node.attrib) or node.attrib["parentjobid"] is None)
    ]
    top_level_job_numbers = [node.attrib["jobnumber"] for node in top_level_job_nodes]
    matching_project = Project.objects.get(uuid=job_node.attrib["projectid"])
    all_db_job_numbers = [
        job.number
        for job in Job.objects.filter(parent__isnull=True, project=matching_project)
    ]
    all_existing_job_numbers = list(set(top_level_job_numbers + all_db_job_numbers))
    if original_job_number in all_db_job_numbers:
        old_job = Job.objects.get(
            number=job_node.attrib["jobnumber"],
            project=matching_project.pk,
        )
        if str(old_job.uuid) == job_node.attrib["jobid"]:
            return job_node.attrib["jobnumber"]
        next_free_job_number = (
            max(int(job_number) for job_number in all_existing_job_numbers) + 1
        )
        job_node.attrib["jobnumber"] = str(next_free_job_number)
        descendent_import_job_nodes = [
            node
            for node in all_import_job_nodes
            if node.attrib["jobnumber"].startswith(f"{original_job_number}.")
        ]
        for descendent_import_job_node in descendent_import_job_nodes:
            job_number_elements = descendent_import_job_node.attrib["jobnumber"].split(
                "."
            )
            new_job_number = ".".join(
                [job_node.attrib["jobnumber"]] + job_number_elements[1:]
            )
            descendent_import_job_node.attrib["jobnumber"] = new_job_number
    # print(original_job_number, job_node.attrib["jobnumber"])
    return job_node.attrib["jobnumber"]


def import_i2xml_from_file(xml_path: Path, relocate_path: Path = None):
    root_node = ET.parse(xml_path).getroot()
    import_i2xml(root_node, relocate_path=relocate_path)


def import_i2xml(root_node: ET.Element, relocate_path: Path):
    """
    Imports data from an XML structure into the system, processing various elements
    such as projects, jobs, files, and tags.
    Args:
        root_node (ET.Element): The root element of the XML structure to import.
        relocate_path (Path): The path to which certain elements may need to be relocated.
    Returns:
        dict: A dictionary containing a mapping of original job numbers to their new numbers.
    """

    job_map = {}
    for node in root_node.findall("ccp4i2_body/projectTable/project"):
        import_project(node, relocate_path)
    job_nodes = root_node.findall("ccp4i2_body/jobTable/job")
    job_nodes = sorted(
        job_nodes,
        key=lambda val: job_number_hash(val.attrib["jobnumber"]),
        reverse=False,
    )
    for node in job_nodes:
        original_job_number = node.attrib["jobnumber"]
        job_map[original_job_number] = renumber_top_job(node, root_node)
        import_job(node)
    for node in root_node.findall("ccp4i2_body/fileTable/file"):
        import_file(node)
    for node in root_node.findall("ccp4i2_body/fileuseTable/fileuse"):
        import_file_use(node)
    for node in root_node.findall("ccp4i2_body/importfileTable/importfile"):
        import_file_import(node)
    for node in root_node.findall("ccp4i2_body/jobkeyvalueTable/jobkeyvalue"):
        import_job_key_value(node)
    for node in root_node.findall("ccp4i2_body/jobkeyvalueTable/jobkeycharvalue"):
        import_job_key_char_value(node)
    for node in root_node.findall("ccp4i2_body/tagTable/tag"):
        import_tag(node)
    for node in root_node.findall("ccp4i2_body/projecttagTable/projecttag"):
        import_project_tag(node)
    return {"job_map": job_map}


def import_project(node: ET.Element, relocate_path: Path = None):

    create_dict = {}
    create_dict["uuid"] = node.attrib["projectid"]
    create_dict["name"] = node.attrib["projectname"]
    create_dict["last_job_number"] = node.attrib["lastjobnumber"]
    directory = node.attrib["projectdirectory"]

    if relocate_path is not None:
        directory = relocate_path / Path(directory).name
    node.attrib["projectdirectory"] = directory
    create_dict["directory"] = str(directory)
    create_dict["creation_time"] = datetime.datetime.fromtimestamp(
        float(node.attrib["projectcreated"])
    )
    create_dict["creation_time"] = datetime.datetime.fromtimestamp(
        float(node.attrib["projectcreated"])
    )

    try:
        instance = Project.objects.get(uuid=create_dict["uuid"])
        item_form = ProjectSerializer(data=create_dict, instance=instance)
    except Project.DoesNotExist as err:
        logging.info(f"Attempting to create new project {err}")
        item_form = ProjectSerializer(data=create_dict)

    if item_form.is_valid():
        new_item = item_form.save()
        new_project = new_item.save()
        logging.info(f"Created new project {new_project}")
        return new_project
    else:
        logging.error(f"Issues creating new project {item_form.errors}")
        return item_form.errors


def import_job(node: ET.Element):
    create_dict = {}
    create_dict["uuid"] = node.attrib["jobid"]
    create_dict["status"] = node.attrib["status"]
    if "evaluation" in node.attrib:
        create_dict["evaluation"] = node.attrib["evaluation"]
    create_dict["job_number"] = node.attrib["jobnumber"]
    create_dict["task_name"] = node.attrib["taskname"]
    create_dict["creation_time"] = datetime.datetime.fromtimestamp(
        float(node.attrib["creationtime"])
    )
    if "finishtime" in node.attrib and node.attrib["finishtime"] is not None:
        create_dict["finish_time"] = datetime.datetime.fromtimestamp(
            float(node.attrib["finishtime"])
        )
    create_dict["number"] = node.attrib["jobnumber"]
    if "title" in node.attrib:
        create_dict["title"] = node.attrib["title"]
    else:
        create_dict["title"] = node.attrib["taskname"]
    create_dict["project"] = Project.objects.get(uuid=node.attrib["projectid"]).pk
    if "parentjobid" in node.attrib and node.attrib["parentjobid"] is not None:
        create_dict["parent"] = Job.objects.get(uuid=node.attrib["parentjobid"]).pk

    try:
        instance = Job.objects.get(uuid=create_dict["uuid"])
        item_form = JobSerializer(data=create_dict, instance=instance)
    except Job.DoesNotExist as err:
        logging.info(f"Attempting to create new job {err}")
        item_form = JobSerializer(data=create_dict)

    if item_form.is_valid():
        new_item = item_form.save()
        new_job = new_item.save()
        logging.info(f"Created new job {new_job}")
        return new_job
    else:
        logging.error(f"Issues creating new job {create_dict} {item_form.errors}")
        return item_form.errors


def import_file(node: ET.Element):
    create_dict = {}
    create_dict["uuid"] = node.attrib["fileid"]
    create_dict["name"] = node.attrib["filename"]
    create_dict["directory"] = node.attrib["pathflag"]
    create_dict["job"] = Job.objects.get(uuid=node.attrib["jobid"]).pk
    create_dict["job_param_name"] = node.attrib["jobparamname"]
    filetypeid = int(node.attrib["filetypeid"])
    file_types = [fileType for fileType in FILETYPELIST if filetypeid == fileType[0]]
    file_type = file_types[0]
    create_dict["type"] = FileType.objects.get(name=file_type[1])

    if "annotation" in node.attrib:
        create_dict["annotation"] = node.attrib["annotation"]
    if "filesubtype" in node.attrib:
        create_dict["sub_type"] = node.attrib["filesubtype"]
    if "filecontent" in node.attrib:
        create_dict["content"] = node.attrib["filecontent"]

    try:
        instance = File.objects.get(uuid=create_dict["uuid"])
        item_form = FileSerializer(data=create_dict, instance=instance)
    except File.DoesNotExist as err:
        logging.info(f"Attempting to create new File {err}")
        item_form = FileSerializer(data=create_dict)

    if item_form.is_valid():
        new_item = item_form.save()
        new_file = new_item.save()
        logging.info(f"Created new File {new_file}")
        return new_file
    else:
        logging.error(f"Issues creating new File {item_form.errors}")
        return item_form.errors


def import_file_use(node: ET.Element):
    create_dict = {}
    create_dict["file"] = File.objects.get(uuid=node.attrib["fileid"]).pk
    create_dict["job"] = Job.objects.get(uuid=node.attrib["jobid"]).pk
    create_dict["role"] = node.attrib["roleid"]
    create_dict["job_param_name"] = node.attrib["jobparamname"]

    try:
        instance = FileUse.objects.get(
            file=create_dict["file"],
            job=create_dict["job"],
            role=create_dict["role"],
            job_param_name=create_dict["job_param_name"],
        )
        item_form = FileUseSerializer(data=create_dict, instance=instance)
    except FileUse.DoesNotExist as err:
        logging.info(f"Attempting to create new FileUse {err}")
        item_form = FileUseSerializer(data=create_dict)

    if item_form.is_valid():
        new_item = item_form.save()
        new_file_use = new_item.save()
        logging.info(f"Created new FileUse {new_file_use}")
        return new_file_use
    else:
        logging.error(f"Issues creating new FileUse {item_form.errors}")
        return item_form.errors


def import_file_import(node: ET.Element):
    create_dict = {}
    create_dict["file"] = File.objects.get(uuid=node.attrib["fileid"]).pk
    create_dict["time"] = datetime.datetime.fromtimestamp(
        float(node.attrib["creationtime"])
    )
    create_dict["name"] = node.attrib["sourcefilename"]
    create_dict["checksum"] = node.attrib["checksum"]

    try:
        instance = FileImport.objects.get(file=create_dict["file"])
        item_form = FileImportSerializer(data=create_dict, instance=instance)
    except FileImport.DoesNotExist as err:
        logging.info(f"Attempting to create new FileImport {err}")
        item_form = FileImportSerializer(data=create_dict)

    if item_form.is_valid():
        new_item = item_form.save()
        new_file_import = new_item.save()
        logging.info(f"Created new FileImport {new_file_import}")
        return new_file_import
    else:
        logging.error("Issues creating new FileImport {item_form.errors}")
        return item_form.errors


def import_job_key_value(node: ET.Element):
    create_dict = {}
    create_dict["job"] = Job.objects.get(uuid=node.attrib["jobid"]).pk
    key_type_id = int(node.attrib["keytypeid"])
    key_types = [key_type for key_type in KEYTYPELIST if key_type_id == key_type[0]]
    key_type = key_types[0]
    create_dict["key"] = JobValueKey.objects.get(name=key_type[1]).pk
    create_dict["value"] = float(node.attrib["value"])

    try:
        instance = JobFloatValue.objects.get(
            job=create_dict["job"], key=create_dict["key"]
        )
        item_form = JobFloatValueSerializer(data=create_dict, instance=instance)
    except JobFloatValue.DoesNotExist as err:
        logging.info(f"Attempting to create new project {err}")
        item_form = JobFloatValueSerializer(data=create_dict)

    if item_form.is_valid():
        new_item = item_form.save()
        new_job_float_value = new_item.save()
        logging.info(f"Created new JobFloatValue {new_job_float_value}")
        return new_job_float_value
    else:
        logging.error(f"Issues creating new JobFloatValue {item_form.errors}")
        return item_form.errors


def import_job_key_char_value(node: ET.Element):
    create_dict = {}
    create_dict["job"] = Job.objects.get(uuid=node.attrib["jobid"]).pk
    key_type_id = int(node.attrib["keytypeid"])
    key_types = [key_type for key_type in KEYTYPELIST if key_type_id == key_type[0]]
    key_type = key_types[0]
    create_dict["key"] = JobValueKey.objects.get(name=key_type[1]).pk
    create_dict["value"] = str(node.attrib["value"])

    try:
        instance = JobCharValue.objects.get(
            job=create_dict["job"], key=create_dict["key"]
        )
        item_form = JobCharValueSerializer(data=create_dict, instance=instance)
    except JobCharValue.DoesNotExist as err:
        logging.info(f"Attempting to create new JobCharValue {err}")
        item_form = JobCharValueSerializer(data=create_dict)

    if item_form.is_valid():
        new_item = item_form.save()
        new_job_char_value = new_item.save()
        logging.info(f"Created new JobCharValue {new_job_char_value}")
        return new_job_char_value
    else:
        logging.error(f"Issues creating new JobCharValue {item_form.errors}")
        return item_form.errors


def import_tag(node: ET.Element):
    tag_map[node.attrib["tagid"]] = node.attrib["text"]


def import_project_tag(node: ET.Element):
    # Check if tg with this text already exists if so,  retrieve and add project to is projects otherwise create with project
    project_tag = None
    try:
        tag_text = tag_map[node.attrib["tagid"]]
        project_tag = ProjectTag.objects.get(text=tag_text)
        project_tag.projects.add(Project.objects.get(uuid=node.attrib["projectid"]))
        project_tag.save()
        return project_tag
    except KeyError as err:
        logging.error(f"Cannot determine text of tag  {node} {err}")
    except ProjectTag.DoesNotExist as err:
        logging.info(f"Creating new ProjectTag with text {tag_text} {err}")
        create_dict = {}
        create_dict["text"] = tag_map[node.attrib["tagid"]]
        create_dict["parent"] = None
        create_dict["projects"] = [
            Project.objects.get(uuid=node.attrib["projectid"]).pk
        ]
        item_form = ProjectTagSerializer(data=create_dict)
        if item_form.is_valid():
            new_item = item_form.save()
            new_project_tag = new_item.save()
            logging.info(f"Created new JobCharValue {new_project_tag}")
            return new_project_tag
        else:
            logging.error(f"Issues creating new JobCharValue {item_form.errors}")
            return item_form.errors
