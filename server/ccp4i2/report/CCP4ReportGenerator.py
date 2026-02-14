import logging

from ccp4i2.core.task_manager.metadata import TITLES
from ccp4i2.db.ccp4i2_static_data import FILETYPES_CLASS


logger = logging.getLogger(f"ccp4i2:{__name__}")


def getReportJobInfo(jobId=None, projectName=None, jobNumber=None):
    """
    Get job information for report generation.

    This function retrieves comprehensive job information needed for generating
    reports. It uses the database abstraction layer (PROJECTSMANAGER().db())
    which supports both database-attached and database-free operation.

    Args:
        jobId: Job UUID string
        projectName: Optional project name
        jobNumber: Optional job number

    Returns:
        dict: Job information dictionary with keys like 'taskname', 'status',
              'inputfiles', 'outputfiles', 'filenames', etc.
    """
    from ccp4i2.core import CCP4Data, CCP4Modules

    # File role constants (matching legacy CCP4DbApi constants)
    FILE_ROLE_OUT = 0
    FILE_ROLE_IN = 1

    db = CCP4Modules.PROJECTSMANAGER().db()

    if projectName is None:
        projectNameInfo = db.getJobInfo(jobId=jobId)
        projectId = projectNameInfo['projectid']
        projectName = db.getProjectInfo(projectId, mode='projectname')

    if jobId is None:
        jobId = db.getJobId(projectName=projectName, jobNumber=jobNumber)

    jobInfo = db.getJobInfo(
        jobId=jobId,
        mode=['runtime', 'status', 'taskname', 'taskversion', 'jobnumber',
              'descendentjobs', 'projectid', 'projectname', 'jobtitle',
              'creationtime']
    )

    jobInfo['jobid'] = jobId
    jobInfo['tasktitle'] = TITLES.get(jobInfo['taskname'])
    jobInfo['fileroot'] = CCP4Modules.PROJECTSMANAGER().makeFileName(
        jobId=jobId, mode='ROOT'
    )

    # Get imported files (if supported by database handler)
    jobInfo['importedfiles'] = []
    try:
        importedfiles = db.getJobImportFiles(jobId=jobId)
        # Returned list JobId, FileID, ImportId, FileTypeId, Filename,
        # Annotation
        for imp in importedfiles:
            jobInfo['importedfiles'].append({
                'fileId': imp[1],
                'filetypeid': imp[3],
                'filetype': '',
                'filetypeclass': FILETYPES_CLASS[imp[3]] if imp[3] < len(FILETYPES_CLASS) else '',
                'filename': imp[4],
                'relpath': '',
                'projectname': projectName,
                'projectid': projectId,
                'annotation': imp[5]
            })
    except (AttributeError, NotImplementedError):
        # Database handler doesn't support getJobImportFiles
        pass

    # Get input and output files
    for key, role in [['inputfiles', FILE_ROLE_IN],
                      ['outputfiles', FILE_ROLE_OUT]]:
        jobInfo[key] = []
        try:
            fileIdList = db.getJobFiles(jobId=jobId, role=role)
            for fileId in fileIdList:
                fileInfo = db.getFileInfo(
                    fileId=fileId,
                    mode=[
                        'filetypeid',
                        'filetype',
                        'filename',
                        'relpath',
                        'projectname',
                        'projectid',
                        'annotation',
                        'jobparamname'])
                filetypeid = fileInfo.get('filetypeid', 0)
                fileInfo['filetypeclass'] = FILETYPES_CLASS[filetypeid] if filetypeid < len(
                    FILETYPES_CLASS) else ''
                fileInfo['fileId'] = fileId
                jobInfo[key].append(fileInfo)
        except (AttributeError, NotImplementedError):
            # Database handler doesn't support this operation
            pass
    jobInfo['filenames'] = {}

    # Try to get params container if database handler supports it
    try:
        container = db.getParamsContainer(jobId=jobId)
    except (AttributeError, NotImplementedError, Exception):
        # Database handler doesn't support getParamsContainer or job not found
        pass
    else:
        # Get the data keyed by task parameter name
        # Not necessarily all files but should all have __str__ methods
        try:
            for key in container.inputData.dataOrder():
                if isinstance(container.inputData.find(key), CCP4Data.CList):
                    jobInfo['filenames'][key] = []
                    for item in container.inputData.find(key):
                        jobInfo['filenames'][key].append(str(item))
                else:
                    jobInfo['filenames'][key] = str(
                        container.inputData.find(key))

            for key in container.outputData.dataOrder():
                if isinstance(container.outputData.find(key), CCP4Data.CList):
                    jobInfo['filenames'][key] = []
                    for item in container.outputData.find(key):
                        jobInfo['filenames'][key].append(str(item))
                else:
                    try:
                        jobInfo['filenames'][key] = str(
                            container.outputData.find(key))
                    except Exception:
                        jobInfo['filenames'][key] = ""
        except Exception as e:
            logger.debug(f"Error extracting filenames from container: {e}")

    return jobInfo
