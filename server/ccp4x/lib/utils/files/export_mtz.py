import logging

from ccp4i2.core import CCP4TaskManager
from ccp4i2.core import CCP4XtalData
from .get_source_reflection_file import get_source_reflection_file
from ccp4x.db import models
from ccp4x.db.ccp4i2_static_data import FILETYPES_TEXT

logger = logging.getLogger(f"ccp4x:{__name__}")


def export_job_mtz_file(job_uuid):
    """
    Export the merged mtz file for a job.
    If the job has not been run yet, or the mtz file does not exist, return None.
    If the job has been run, but no mtz file was created, return None.
    If the job has been run and an mtz file was created, return the path to the mtz file.
    The mtz file is created by running CAD on the input reflection files.
    The input reflection files are either from aimless or an imported file.
    The column labels in the output mtz file are set by the LABOUT command in CAD.
    The mtz file is created in the job directory with the name 'exportMtz.mtz'.
    The function returns the path to the exported mtz file.
    If the mtz file already exists, it is returned without running CAD again.
    """
    the_job = models.Job.objects.get(uuid=job_uuid)
    job_dir = the_job.directory
    export_file = job_dir / "exportMtz.mtz"
    if export_file.exists():
        return str(export_file)
    export_params = CCP4TaskManager.TASKMANAGER().getTaskAttribute(
        the_job.task_name, "EXPORTMTZPARAMS", default=[]
    )
    if len(export_params) == 0:
        logger.warning("No export parameters found")
        return None
    refln_info = get_source_reflection_file(
        jobId=the_job.uuid, jobParamNameList=export_params[0]
    )
    logger.warning(
        "CProjectViewer.exportJobMtzFile getSourceReflectionFile %s %s",
        the_job,
        refln_info,
    )
    if refln_info.get("fileName", None) is None:
        return None  # No source reflection file found
    print(
        "CProjectViewer.exportJobMtzFile getSourceReflectionFile", the_job, refln_info
    )

    file_info: list[models.File] = []
    param_name_list: list[str] = []
    label_list: list[str] = []
    for ep in export_params[1:]:
        if isinstance(ep, list):
            param, lab = ep
        else:
            param = ep
            lab = None
        try:
            the_file = models.File.objects.get(job__uuid=job_uuid, job_param_name=param)
            file_info.append(the_file)
            param_name_list.append(param)
            label_list.append(
                the_file.job.number
                + "_"
                + CCP4TaskManager.TASKMANAGER().getTaskLabel(the_file.job.task_name)
            )
            try:
                _ = models.FileImport.objects.get(file=the_file)
                label_list[-1] = label_list[-1] + "_import"
            except models.FileImport.DoesNotExist:
                pass
            if lab is not None:
                label_list[-1] = label_list[-1] + "_" + lab
        except models.File.DoesNotExist:
            logger.warning("File not found for job %s and param %s", job_uuid, param)
    print("file_info", file_info)
    col_tag_list = CCP4TaskManager.TASKMANAGER().exportMtzColumnLabels(
        taskName=the_job.task_name,
        jobId=str(job_uuid),
        paramNameList=param_name_list,
        sourceInfoList=file_info,
    )
    if len(col_tag_list) > 0:
        for ii in range(len(col_tag_list)):
            label_list[ii] = col_tag_list[ii]
    # print 'CProjectViewer.exportJobMtzFile fileInfo',file_info
    if len(file_info) == 0:
        return None
    # Create CAD inputFiles and command lines
    file_no = 1
    com_lines = []
    # Dont need reflnInfo['fileName'] in inputFiles as it is used for the CMtzDataFile object that runs CAD
    input_files = []
    for ifile, f_info in enumerate(file_info):
        file_type = FILETYPES_TEXT.index(f_info.type.name) if f_info.type else 0
        print(
            {field.name: getattr(f_info, field.name) for field in f_info._meta.fields},
            f_info.type.pk,
            file_type,
        )
        if file_type in [11, 12]:
            # Beware alternative fileContent
            if f_info.content is not None:
                f_c = f_info.content
            else:
                if file_type == 11:
                    p = CCP4XtalData.CObsDataFile(str(f_info.path))
                else:
                    p = CCP4XtalData.CPhsDataFile(str(f_info.path))
                p.setContentFlag()
                f_c = int(p.contentFlag)
        # Use CAD LABOUT line to set column labels in export file
        file_no += 1
        input_files.append(str(f_info.path))
        if f_info.type == 10:
            com_lines.append(
                "LABOUT FILENUMBER " + str(file_no) + " E1=FREER_" + label_list[ifile]
            )
        elif file_type == 11:
            if f_info.content == 1:
                com_lines.append(
                    "LABOUT FILENUMBER "
                    + str(file_no)
                    + " E1=I(+)_"
                    + label_list[ifile]
                    + " E2=SIGI(+)_"
                    + label_list[ifile]
                    + " E3=I(-)_"
                    + label_list[ifile]
                    + " E4=SIGI(-)_"
                    + label_list[ifile]
                )
            elif f_info.content == 2:
                com_lines.append(
                    "LABOUT FILENUMBER "
                    + str(file_no)
                    + " E1=F(+)_"
                    + label_list[ifile]
                    + " E2=SIGF(+)_"
                    + label_list[ifile]
                    + " E3=F(-)_"
                    + label_list[ifile]
                    + " E4=SIGF(-)_"
                    + label_list[ifile]
                )
            elif f_info.content == 3:
                com_lines.append(
                    "LABOUT FILENUMBER "
                    + str(file_no)
                    + " E1=I_"
                    + label_list[ifile]
                    + " E2=SIGI_"
                    + label_list[ifile]
                )
            else:
                com_lines.append(
                    "LABOUT FILENUMBER "
                    + str(file_no)
                    + " E1=F_"
                    + label_list[ifile]
                    + " E2=SIGF_"
                    + label_list[ifile]
                )
        elif file_type == 12:
            if f_info.content == 1:
                com_lines.append(
                    "LABOUT FILENUMBER "
                    + str(file_no)
                    + " E1=HLA_"
                    + label_list[ifile]
                    + " E2=HLB_"
                    + label_list[ifile]
                    + " E3=HLC_"
                    + label_list[ifile]
                    + " E4=HLD_"
                    + label_list[ifile]
                )
            else:
                com_lines.append(
                    "LABOUT FILENUMBER "
                    + str(file_no)
                    + " E1=PHI_"
                    + label_list[ifile]
                    + "  E2=FOM_"
                    + label_list[ifile]
                )
        elif file_type == 13:
            com_lines.append(
                "LABOUT FILENUMBER "
                + str(file_no)
                + " E1=F_"
                + label_list[ifile]
                + "  E2=PHI_"
                + label_list[ifile]
            )
    # print 'CProjectViewer.exportJobMtzFile',input_files
    # print 'CProjectViewer.exportJobMtzFile',com_lines
    #  Create an CMtzDataFile object and initialise with the refln data file
    m = CCP4XtalData.CMtzDataFile(refln_info["fileName"])
    # print m.runCad.__doc__   #Print out docs for the function
    print("com_lines", com_lines)
    outfile, err = m.runCad(str(export_file), input_files, com_lines)
    # print 'CProjectViewer.exportJobMtzFile',outfile,err.report()
    return outfile


""" def exportJobMtzFile(self, jobId):
        from ccp4i2.core import CCP4XtalData
        from ccp4i2.core import CCP4TaskManager
        # Devise name for the merged file and check if it has already been created
        jobDir = self.jobDirectory(jobId=jobId, create=False)
        jobInfo = self.db().getJobInfo(jobId, ['taskname', 'jobnumber'])
        exportFile = os.path.join(jobDir, 'exportMtz.mtz')
        if os.path.exists(exportFile):
            return exportFile
        exportParams = CCP4TaskManager.TASKMANAGER().getTaskAttribute(jobInfo['taskname'], 'EXPORTMTZPARAMS', default=[])
        if len(exportParams) == 0:
            return None
        # Get the source reflection data either from aimless or an imported file
        # getSourceReflectionFile() returns a dict with elements: fileName, source, jobNumber
        reflnInfo = self.getSourceReflectionFile(jobId = jobId, jobParamNameList=exportParams[0])
        #print 'CProjectViewer.exportJobFile getSourceReflectionFile',jobInfo,reflnInfo
        if reflnInfo.get('fileName', None) is None:
            return None
        # Query database for filenames and job info for the input and ouptput objects
        
        fileInfo = []
        paramNameList = []
        for ep in exportParams[1:]:
            if isinstance(ep, list):
                param, lab = ep
            else:
                param = ep
                lab = None
            fInfo = self.db().getJobFilesInfo(jobId=jobId, jobParamName=param)
            if len(fInfo) == 0:
                fInfo = self.db().getJobFilesInfo(jobId=jobId, jobParamName=param, input=True)
            if len(fInfo) > 0:
                fileInfo.append(fInfo[0])
                paramNameList.append(param)
                fileInfo[-1]['label'] = fileInfo[-1]['jobnumber'] + '_' + CCP4TaskManager.TASKMANAGER().getTaskLabel(fileInfo[-1]['taskname'])
                if fileInfo[-1]['importId'] is not None:
                    fileInfo[-1]['label'] = fileInfo[-1]['label'] + '_import'
                if lab is not None:
                    fileInfo[-1]['label'] = fileInfo[-1]['label'] + '_' + lab
        colTagList = CCP4TaskManager.TASKMANAGER().exportMtzColumnLabels(taskName=jobInfo['taskname'], jobId=jobId,
                                                                         paramNameList=paramNameList, sourceInfoList=fileInfo)
        if len(colTagList) > 0:
            for ii in range(len(colTagList)):
                fileInfo[ii]['label'] = colTagList[ii]
        #print 'CProjectViewer.exportJobMtzFile fileInfo',fileInfo
        if len(fileInfo) == 0:
            return None
        # Create CAD inputFiles and command lines
        fileNo = 1
        comLines = []
        # Dont need reflnInfo['fileName'] in inputFiles as it is used for the CMtzDataFile object that runs CAD
        inputFiles = []
        for fInfo in fileInfo:
            if fInfo['fileTypeId'] in [11, 12]:
                # Beware alternative fileContent
                if fInfo.get('fileContent',None) is not None:
                    fC = fInfo['fileContent']
                else:
                    if fInfo['fileTypeId'] == 11:
                        p = CCP4XtalData.CObsDataFile(fInfo['fullPath'])
                    else:
                        p = CCP4XtalData.CPhsDataFile(fInfo['fullPath'])
                    p.setContentFlag()
                    fC = int(p.contentFlag)
            # Use CAD LABOUT line to set column labels in export file
            fileNo += 1
            inputFiles.append(fInfo['fullPath'])
            if  fInfo['fileTypeId'] == 10:
                comLines.append('LABOUT FILENUMBER ' + str(fileNo) + ' E1=FREER_' + fInfo['label'])
            elif fInfo['fileTypeId'] == 11:
                if fC == 1:
                    comLines.append('LABOUT FILENUMBER ' + str(fileNo) + ' E1=I(+)_' + fInfo['label'] + ' E2=SIGI(+)_' + fInfo['label'] + ' E3=I(-)_' + fInfo['label'] + ' E4=SIGI(-)_' + fInfo['label'])
                elif fC == 2:
                    comLines.append('LABOUT FILENUMBER ' + str(fileNo) + ' E1=F(+)_' + fInfo['label'] + ' E2=SIGF(+)_' + fInfo['label'] + ' E3=F(-)_' + fInfo['label'] + ' E4=SIGF(-)_' + fInfo['label'])
                elif fC == 3:
                    comLines.append('LABOUT FILENUMBER ' + str(fileNo) + ' E1=I_' + fInfo['label'] + ' E2=SIGI_' + fInfo['label'])
                else:
                    comLines.append('LABOUT FILENUMBER ' + str(fileNo) + ' E1=F_' + fInfo['label'] + ' E2=SIGF_' + fInfo['label'])
            elif fInfo['fileTypeId'] == 12:
                if fC == 1:
                    comLines.append('LABOUT FILENUMBER ' + str(fileNo) + ' E1=HLA_' + fInfo['label'] + ' E2=HLB_' + fInfo['label'] + ' E3=HLC_'+fInfo['label']+' E4=HLD_'+fInfo['label'] )
                else:
                    comLines.append('LABOUT FILENUMBER ' + str(fileNo) + ' E1=PHI_' + fInfo['label'] + '  E2=FOM_' + fInfo['label'])
            elif fInfo['fileTypeId'] == 13:
                comLines.append('LABOUT FILENUMBER ' + str(fileNo) + ' E1=F_' + fInfo['label'] + '  E2=PHI_' + fInfo['label'])
        #print 'CProjectViewer.exportJobMtzFile',inputFiles
        #print 'CProjectViewer.exportJobMtzFile',comLines
        #  Create an CMtzDataFile object and initialise with the refln data file
        m = CCP4XtalData.CMtzDataFile(reflnInfo['fileName'])
        #print m.runCad.__doc__   #Print out docs for the function
        outfile, err = m.runCad(exportFile, inputFiles ,comLines)
        #print 'CProjectViewer.exportJobMtzFile',outfile,err.report()
        return outfile

    def getSourceReflectionFile(self, jobId=None, jobParamNameList=None):
        from ccp4i2.core import CCP4XtalData
        from ccp4i2.core import CCP4TaskManager
        exportTaskName = self.db().getJobInfo(jobId=jobId, mode='taskname')
        #print 'into getSourceReflectionFile',exportTaskName
        if not isinstance(jobParamNameList,list):
            jobParamNameList = [jobParamNameList]
        reflnFileList = []
        pList = []
        fList= []
        iList = []
        for jobParam in jobParamNameList:
            reflnFile = None
            fileInfoList = self.db().getJobFilesInfo(jobId=jobId, jobParamName=jobParam, input=True)
            #print 'getSourceReflectionFile fileInfoList',jobParam,fileInfoList
            if len(fileInfoList) > 0:
                taskName = self.db().getJobInfo(jobId=fileInfoList[0]['jobId'], mode='taskname')
                jobNumber =  self.db().getJobInfo(jobId=fileInfoList[0]['jobId'], mode='jobnumber')
                if taskName == 'aimless_pipe':
                    mode = 'Data reduction'
                    from pipelines.aimless_pipe.script import aimless_pipe
                    reflnFile = aimless_pipe.exportJobFile(jobId=fileInfoList[0]['jobId'], mode='complete_mtz')
                elif taskName == 'import_merged':
                    mode = 'Imported by import merged'
                    from pipelines.import_merged.script import import_merged
                    reflnFile = import_merged.exportJobFile(jobId=fileInfoList[0]['jobId'], mode='complete_mtz')
                elif taskName == 'AlternativeImportXIA2':
                    mode = 'Imported from XIA2'
                    import AlternativeImportXIA2
                    reflnFile = AlternativeImportXIA2.exportJobFile(jobId=fileInfoList[0]['jobId'], mode='complete_mtz', fileInfo=fileInfoList[0])
                else:
                    mode = 'Imported file'
                    importInfo = self.db().getImportFileInfo(fileId=fileInfoList[0]['fileId'])
                    #print 'getSourceReflectionFile importInfo',importInfo
                    reflnFile = importInfo['sourcefilename']
                if reflnFile is not None:
                    reflnFileList.append([jobParam, reflnFile, fileInfoList[0]])
                    pList.append(jobParam)
                    fList.append(reflnFile)
                    iList.append(fileInfoList[0])
        colTagList = CCP4TaskManager.TASKMANAGER().exportMtzColumnLabels(taskName=exportTaskName, jobId=jobId, paramNameList=pList, sourceInfoList=iList)
        if len(colTagList) == 0:
            colTagList = pList
        if len(fList) == 0:
            return {}
        if len(fList) == 1:
            return {'fileName': fList[0], 'source': mode, 'jobNumber' : jobNumber}
        else:
            # Multiple reflection input need to be cadded and colmn labels sorted
            # Likely to have redundant freer?
            #print 'CProjectManager.getSourceReflectionFile reflnFileList',reflnFileList
            jobDir = self.jobDirectory(jobId=jobId, create=False)
            exportFile = os.path.join(jobDir, 'exportMtz_reflns.mtz')
            fileNo = 0
            comLines = []
            inputFiles = []
            for ii in range(len(fList)):
                if ii == 0:
                    exObj = CCP4XtalData.CMtzDataFile(fList[ii])
                else:
                    inputFiles.append(fList[ii])
                refnObj = CCP4XtalData.CMtzDataFile(fList[ii])
                comLine = 'LABOUT FILENUMBER ' + str(ii + 1)
                jj = 0
                for column in refnObj.fileContent.listOfColumns:
                    jj += 1
                    comLine += ' E' + str(jj)+'=' + str(column.columnLabel) + '_' + colTagList[ii]
                comLines.append(comLine)
            #print 'CProjectManager.getSourceReflectionFile input files',inputFiles
            outfile, err = exObj.runCad(exportFile, inputFiles ,comLines )
            #print 'CProjectManager.getSourceReflectionFile',outfile,err.report()
            return {'fileName': outfile, 'source': mode, 'jobNumber' : jobNumber}
"""
