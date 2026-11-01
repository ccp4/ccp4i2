"""Trace a refined model back through its refinement job to the scaling job,
collecting the files and settings the deposition needs.

Django models are imported lazily: this runs on the server, from a
``plugin_method`` the task interface calls when the user picks a model.
"""

REFINEMENT_TASKS = ("prosmart_refmac", "servalcat_pipe", "buster")
FLAG_NAMES = ("USEANOMALOUS", "USE_TWIN")


def refinement_flags(job_dir):
    """``{"USEANOMALOUS": bool, "USE_TWIN": bool}`` (whichever are found)
    from the refinement job's saved parameters, never from its plugin.

    A pipeline (prosmart_refmac, servalcat_pipe) keeps the refmac settings
    in a sub-container or in a sub-job's parameters, so the files are
    searched for the elements wherever they sit: the job's own
    input_params.xml / params.xml first, then its sub-jobs' (job_N_1 ...).
    Never raises: a shape it does not know yields nothing.
    """
    from pathlib import Path

    from lxml import etree

    found = {}
    job_dir = Path(str(job_dir))
    candidates = [job_dir / "params.xml", job_dir / "input_params.xml"]
    try:
        for sub in sorted(job_dir.iterdir()):
            if sub.is_dir() and sub.name.startswith("job_"):
                candidates += [sub / "params.xml", sub / "input_params.xml"]
    except OSError:
        pass
    for path in candidates:
        if len(found) == len(FLAG_NAMES):
            break
        try:
            root = etree.parse(str(path)).getroot()
        except (OSError, etree.XMLSyntaxError):
            continue
        for name in FLAG_NAMES:
            if name in found:
                continue
            node = next(iter(root.iter(name)), None)
            if node is not None and node.text is not None and node.text.strip():
                found[name] = node.text.strip().lower() in ("true", "1", "yes")
    return found


def trace_xyzin_lineage(file_uuid):
    """For a coordinates file: the refinement job that produced it, that job's
    inputs (reflections, free set, TLS, dictionaries), its parameters, and the
    scaling job's statistics and unmerged data.

    Returns a JSON-serialisable dict: ``success``, ``files`` (parameter name
    to file uuid), ``paths``, ``params``, ``warnings``, ``refinement_job_id``,
    ``refinement_task_name``, ``scaling_job_id``, ``software`` (program
    versions found along the way).
    """
    from pathlib import Path

    from lxml import etree

    from ccp4i2.db.models import File, FileUse, Job

    result = {"success": False, "files": {}, "paths": {}, "params": {},
              "warnings": [], "software": []}
    try:
        the_file = File.objects.get(uuid=file_uuid)
    except File.DoesNotExist:
        result["warnings"].append("File not found in database")
        return result
    ref_job = the_file.job
    if ref_job is None:
        result["warnings"].append("File has no producing job")
        return result
    if ref_job.task_name not in REFINEMENT_TASKS:
        result["warnings"].append(
            "Coordinates should come from a refinement job, not '%s'" % ref_job.task_name)
        return result
    if ref_job.status != Job.Status.FINISHED:
        result["warnings"].append("Source refinement job did not finish successfully")
        return result

    result["success"] = True
    result["refinement_job_id"] = ref_job.id
    result["refinement_task_name"] = ref_job.task_name

    for use in FileUse.objects.filter(job=ref_job, role=FileUse.Role.IN).select_related("file"):
        uuid = str(use.file.uuid)
        param = use.job_param_name
        if param in ("F_SIGF", "HKLIN"):
            result["files"]["F_SIGF"] = uuid
        elif param in ("FREERFLAG", "TLSIN"):
            result["files"][param] = uuid
        elif param in ("DICT", "DICT_LIST"):
            result["files"].setdefault("DICT_LIST", []).append(uuid)
    for f in File.objects.filter(job=ref_job):
        if f.job_param_name == "TLSOUT":
            result["files"]["TLSIN"] = str(f.uuid)
        elif f.job_param_name in ("FPHIOUT", "DIFFPHIOUT"):
            result["files"][f.job_param_name] = str(f.uuid)

    if ref_job.task_name != "buster":
        found = refinement_flags(ref_job.directory)
        result["params"].update(found)
        if not found:
            result["warnings"].append(
                "Could not find the refinement's anomalous/twin settings; check them on the "
                "Files page")

    f_sigf_uuid = result["files"].get("F_SIGF")
    if not f_sigf_uuid:
        return result
    try:
        scaling_job = File.objects.select_related("job").get(uuid=f_sigf_uuid).job
    except File.DoesNotExist:
        result["warnings"].append("Could not trace F_SIGF to its producing job")
        return result
    if scaling_job is None:
        return result
    result["scaling_job_id"] = scaling_job.id

    cifstats = File.objects.filter(job=scaling_job, job_param_name="CIFSTATSOUT")
    if cifstats.exists():
        result["files"]["CIFSTATSOUT"] = str(cifstats.first().uuid)
        result["params"]["USEAIMLESSXML"] = True

    xml_paths = []
    xmlout = File.objects.filter(job=scaling_job, job_param_name="XMLOUT").first()
    if xmlout is not None:
        xml_paths.append((Path(scaling_job.directory) / xmlout.name, str(xmlout.uuid)))
    xml_paths.append((Path(scaling_job.directory) / "program.xml", None))
    for xml_path, uuid in xml_paths:
        if not xml_path.exists():
            continue
        try:
            root = etree.parse(str(xml_path)).getroot()
        except etree.XMLSyntaxError:
            continue
        for program in ("AIMLESS", "POINTLESS", "CTRUNCATE"):
            node = root.find(".//%s" % program)
            if node is not None and node.get("version"):
                result["software"].append({"name": program.capitalize(), "version": node.get("version"),
                                           "classification": "data scaling"})
        if "CIFSTATSOUT" not in result["files"] and \
                (root.findall(".//AIMLESS") or root.findall(".//AIMLESS_PIPE")) and \
                root.tag != "IMPORT_MERGED":
            if uuid:
                result["files"]["AIMLESSXML"] = uuid
            else:
                result["paths"]["AIMLESSXML"] = str(xml_path)
            result["params"]["USEAIMLESSXML"] = True
        break

    unmerged = File.objects.filter(job=scaling_job, job_param_name__startswith="UNMERGEDOUT")
    if unmerged.count() >= 1:
        result["files"]["SCALEDUNMERGED"] = str(unmerged[0].uuid)
        result["params"]["INCLUDEUNMERGED"] = True
        if unmerged.count() > 1:
            result["warnings"].append("Multiple unmerged outputs found; using the first of %d"
                                      % unmerged.count())
    else:
        result["warnings"].append("Scaling job did not output unmerged data")
    return result
