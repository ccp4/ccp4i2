from __future__ import print_function
"""
    adding_stats_to_mmcif_i2_gui.py: CCP4 GUI Project
    
    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.
    
    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    """

from qtgui.CCP4TaskWidget import CTaskWidget
from ccp4i2.baselayer import QtGui, QtWidgets
from ccp4i2.baselayer import QtCore
import os
from dbapi import CCP4DbApi
import functools
from ccp4i2.core.CCP4Modules import PROJECTSMANAGER
from ccp4i2.core.CCP4ErrorHandling import CException
import xml.etree.ElementTree as ET

# -------------------------------------------------------------------


class adding_stats_to_mmcif_i2_gui(CTaskWidget):
    # -------------------------------------------------------------------

    # Subclass CTaskWidget to give specific task window
    # this has to match the pluginName given in the corresponding .def.xml
    TASKNAME = 'adding_stats_to_mmcif_i2'
    TASKVERSION = 0.1
    # Section in the task list where this task will be listed e.g. 'refinement','model_building' for full list see MODULE_ORDER in core/CCP4TaskManager.py
    TASKMODULE = ['export']
    SHORTTASKTITLE = 'New deposition task'
    TASKTITLE = 'Prepare and validate files for deposition'
    DESCRIPTION = '''Add data reduction statistics and sequence information into coordinates for deposition'''
    MGDISPLAYFILES = ['XYZIN']
    WHATNEXT = ['coot_rebuild']

    def __init__(self, parent):
        CTaskWidget.__init__(self, parent)

    def drawContents(self):
        self.openFolder(folderFunction='inputData', followFrom=False)

        self.createLine(
            ['subtitle', 'Asymmetric unit content (i.e. sequences)'])
        self.openSubFrame(frame=True)
        self.createLine(['tip', 'ASU Content', 'widget', 'ASUCONTENT'])
        self.closeSubFrame()

        self.openSubFrame(frame=True, title='Coordinates')
        self.createLine(['tip', 'Input model', 'widget',
                        '-browseFiles', False, 'XYZIN'])
        self.container.inputData.XYZIN.dataChanged.connect(self.xyzinSet)
        self.closeSubFrame()
        self.openSubFrame(frame=True, title='Task control')
        self.createLine(['label', 'Use validation server', 'tip',
                        'Send to validation server - requires internet access', 'widget', 'SENDTOVALIDATIONSERVER'])
        self.createLine(['label', 'Add statistics from Aimless', 'tip',
                        "If you didn't run data reduction in CCP4i2 you won't have this", 'widget', 'USEAIMLESSXML'])
        self.container.controlParameters.USEAIMLESSXML.dataChanged.connect(
            self.UseAimlessXmlChanged)
        self.createLine(['label', 'Include unmerged from scaling job', 'tip',
                        "If you didn't run data reduction in CCP4i2 you won't have this", 'widget', 'INCLUDEUNMERGED'])
        self.container.controlParameters.INCLUDEUNMERGED.dataChanged.connect(
            self.IncludeUnmergedChanged)
        self.closeSubFrame()

        self.openSubFrame(frame=True, title='Related files')
        self.createLine(['widget', 'USEANOMALOUS', 'label',
                        'Anomalous used in  refmac job'])
        self.createLine(['widget', 'USE_TWIN', 'label',
                        'Twinning (i.e. intensities) used in refmac job'])
        self.createLine(['tip', 'Input reflections', 'widget',
                        '-browseFiles', False, 'F_SIGF'])
        self.createLine(['widget', '-browseFiles', True, 'SCALEDUNMERGED'])
        self.createLine(['widget', '-guiLabel', 'Aimless XML',
                        '-browseFiles', False, 'AIMLESSXML'])
        self.createLine(['widget', '-guiLabel', 'Refmac task XML',
                        '-browseFiles', False, 'REFMACINPUTPARAMSXML'])
        self.createLine(['tip', 'FreeR flag', 'widget',
                        '-browseFiles', False, 'FREERFLAG'])
        self.createLine(['tip', 'Input TLS', 'widget',
                        '-browseFiles', False, 'TLSIN'])
        self.createLine(['tip', 'Dictionary list', 'widget', 'DICT_LIST'])
        self.createLine(
            ['widget', '-guiLabel', 'Map coefficients', '-browseFiles', False, 'FPHIOUT'])
        self.createLine(['widget', '-guiLabel', 'Map coefficients',
                        '-browseFiles', False, 'DIFFPHIOUT'])
        for widgetName in ['TLSIN', 'DICT_LIST', 'FREERFLAG', 'AIMLESSXML', 'REFMACINPUTPARAMSXML', 'F_SIGF', 'USEANOMALOUS', 'USE_TWIN', 'FPHIOUT', 'DIFFPHIOUT', 'SCALEDUNMERGED']:
            self.getWidget(widgetName).setDisabled(True)
        self.getWidget('DICT_LIST').setListVisible(True)
        self.closeSubFrame()

    def updateFromRefmacJob(self, jobId, includeXYZIN=False):
        # Start off by unsetting
        for paramName in ['TLSIN', 'FREERFLAG', 'AIMLESSXML', 'REFMACINPUTPARAMSXML', 'F_SIGF', 'FPHIOUT', 'DIFFPHIOUT', 'SCALEDUNMERGED']:
            getattr(self.container.inputData, paramName).unSet()

        jobDir = PROJECTSMANAGER().db().jobDirectory(jobId)
        inputParamsPath = os.path.join(jobDir, 'input_params.xml')
        self.container.inputData.REFMACINPUTPARAMSXML.set(inputParamsPath)

        # Empty existing dictionary list
        while len(self.container.inputData.DICT_LIST) > 0:
            self.container.inputData.DICT_LIST.remove(
                self.container.inputData.DICT_LIST[-1])

        from report.CCP4ReportGenerator import getReportJobInfo
        # print getReportJobInfo(jobId=jobId)

        paramsContainer = PROJECTSMANAGER().db().getParamsContainer(jobId=jobId)
        self.container.controlParameters.USEANOMALOUS = paramsContainer.controlParameters.USEANOMALOUS
        self.container.controlParameters.USE_TWIN = paramsContainer.controlParameters.USE_TWIN

        jobFiles = PROJECTSMANAGER().db().getJobFilesInfo(jobId=jobId, input=True)
        for inputFile in jobFiles:
            print(inputFile, "\n")
            # Here identify inputmonomer dictionary entries
            if inputFile['fileTypeId'] in [8]:
                self.container.inputData.DICT_LIST.append(
                    self.container.inputData.DICT_LIST.makeItem())
                self.container.inputData.DICT_LIST[-1].setDbFileId(
                    inputFile['fileId'])
            elif inputFile['fileTypeId'] in [9]:  # Here identify input TLS file
                self.container.inputData.TLSIN.setDbFileId(inputFile['fileId'])
            elif inputFile['fileTypeId'] in [10]:  # Here identify input FREER Flag file
                self.container.inputData.FREERFLAG.setDbFileId(
                    inputFile['fileId'])
            elif inputFile['fileTypeId'] in [11]:  # Here identify input Reflections
                self.container.inputData.F_SIGF.setDbFileId(
                    inputFile['fileId'])
                self.container.inputData.F_SIGF = self.container.inputData.F_SIGF
                jobDir = PROJECTSMANAGER().db().jobDirectory(
                    jobId=inputFile['jobId'])
                xmlPath = os.path.join(jobDir, "program.xml")
                xmlStats = False
                if os.path.isfile(xmlPath):
                    tree = ET.parse(xmlPath)
                    root = tree.getroot()
                    aimlessNodes = root.findall('.//AIMLESS')
                    if len(aimlessNodes) > 0 and root.tag != 'IMPORT_MERGED':
                        self.container.inputData.AIMLESSXML.set(xmlPath)
                        self.getWidget('AIMLESSXML').setDisabled(True)
                        xmlStats = True
                if not xmlStats:
                    self.warn("Data reduction statistics needed",
                              "The input reflections were not apparently generated in a tracked AIMLESS job",
                              "Where data reduction was done in a job of this project, merging statistics are automatically tracked. For the set used here, you will have to locate XML output from the program AIMLESS or uncheck the 'Add statistics from Aimless option'")
                    self.container.inputData.AIMLESSXML.unSet()
                    self.getWidget('AIMLESSXML').setDisabled(False)
                # Look to see if there was a scaled unmerged file in this job
                scalingJobFiles = PROJECTSMANAGER().db(
                ).getJobFilesInfo(jobId=inputFile['jobId'])
                outputUnmergeds = [
                    scalingJobFile for scalingJobFile in scalingJobFiles if scalingJobFile['jobParamName'] == 'UNMERGEDOUT']
                if len(outputUnmergeds) == 1:
                    self.warn('Scaled unmerged warning',
                              'Including scaled unmerged data from scaling job',
                              'To exclude from deposition uncheck "Include unmerged from scaling job"')
                    self.container.inputData.SCALEDUNMERGED.setDbFileId(
                        outputUnmergeds[0]['fileId'])
                    self.container.controlParameters.INCLUDEUNMERGED.set(True)
                    self.updateViewFromModel()
                elif len(outputUnmergeds) == 0:
                    self.warn('No scaled unmerged warning',
                              'Job that generated data used in refinement did not output scaled unmerged',
                              'Either uncheck "Include unmerged from scaling job" or manually set "Scaled unmerged datafile')
                    self.getWidget('SCALEDUNMERGED').setDisabled(False)
                    self.updateViewFromModel()
                else:
                    self.warn('Multiple scaled unmerged warning',
                              'Job that generated data used in refinement made multiple scaled unmerged',
                              'Either uncheck "Include unmerged from scaling job" or manually set "Scaled unmerged datafile"')
                    self.getWidget('SCALEDUNMERGED').setDisabled(False)
                    self.updateViewFromModel()

        # Override if an updated version was output
        jobFiles = PROJECTSMANAGER().db().getJobFilesInfo(jobId=jobId, input=False)
        for file in jobFiles:
            if file['jobParamName'] == 'TLSOUT':
                self.container.inputData.TLSIN.setDbFileId(file['fileId'])
            if file['jobParamName'] == 'XYZOUT' and includeXYZIN:
                self.container.inputData.XYZIN.setDbFileId(file['fileId'])
            # Get output map coefficients from REFMAC
            if file['jobParamName'] == 'FPHIOUT':
                self.container.inputData.FPHIOUT.setDbFileId(file['fileId'])
            if file['jobParamName'] == 'DIFFPHIOUT':
                self.container.inputData.DIFFPHIOUT.setDbFileId(file['fileId'])

    @QtCore.Slot()
    def xyzinSet(self):
        print("dbFIleId ", self.container.inputData.XYZIN.dbFileId)
        # pull inputs from the relevant job
        try:
            creatingJob = PROJECTSMANAGER().db().getFileInfo(
                fileId=self.container.inputData.XYZIN.dbFileId, mode="jobid")
            jobInfo = PROJECTSMANAGER().db().getJobInfo(
                jobId=creatingJob, mode=["status", "taskname"])
            print(jobInfo)
            if jobInfo["status"] == 'Finished' and jobInfo["taskname"] == "prosmart_refmac":
                self.updateFromRefmacJob(creatingJob)
            elif jobInfo["taskname"] != "prosmart_refmac":
                #from PyQt4.QtCore import *
                self.warn("Wrong coordinate source warning",
                          "Coordinates to be deposited should be the output of a prosmart_refmac job",
                          "For deposition, the input coordinates have to be output from a standard CCP4i2 refinement job.  Please run such a job and retry deposition using the output")
                self.container.inputData.XYZIN.unSet()
        except CException as err:
            print("Exception retrieving creatingJob or jobInfo")
            self.updateViewFromModel()

    @QtCore.Slot()
    def UseAimlessXmlChanged(self):
        if not self.container.controlParameters.USEAIMLESSXML:
            self.container.inputData.AIMLESSXML.setQualifier(
                'allowUndefined', True)
        else:
            self.container.inputData.AIMLESSXML.setQualifier(
                'allowUndefined', False)
        self.updateViewFromModel()

    @QtCore.Slot()
    def IncludeUnmergedChanged(self):
        if not self.container.controlParameters.INCLUDEUNMERGED:
            self.container.inputData.SCALEDUNMERGED.setQualifier(
                'allowUndefined', True)
            self.container.inputData.SCALEDUNMERGED.unSet()
        else:
            self.container.inputData.SCALEDUNMERGED.setQualifier(
                'allowUndefined', False)
        self.updateViewFromModel()

    def warn(self, windowTitle, informativeText, detailedText):
        from ccp4i2.baselayer.QtWidgets import QMessageBox
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Warning)
        msg.setText(windowTitle)
        msg.setInformativeText(informativeText)
        msg.setWindowTitle(windowTitle)
        msg.setDetailedText(detailedText)
        msg.setStandardButtons(QMessageBox.Ok)
        # msg.buttonClicked.connect(msgbtn)
        retval = msg.exec_()
