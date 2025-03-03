from __future__ import print_function


"""
     CCP4ErrorReportViewer.py: CCP4 GUI Project
     Copyright (C) 2011 University of York

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

"""
   Liz Potterton Jan 2011 - List program error log
"""

##@package CCP4ProjectWidget View a project
import os
import time
import tempfile
from PySide2 import QtGui, QtWidgets,QtCore
from qtgui import CCP4AbstractViewer
from qtgui import CCP4TextViewer
#from CCP4ErrorHandling import *
from core.CCP4ProgramLog import *


class CErrorReportSelection(QtWidgets.QFrame):

    def __init__(self,parent):
        QtWidgets.QFrame.__init__(self,parent)
        self.setLayout(QtWidgets.QHBoxLayout())
        self.layout().addWidget(QtWidgets.QLabel('Select reports from class:'))
        self.classCombo = QtWidgets.QComboBox(self)
        self.layout().addWidget(self.classCombo)
        self.layout().addWidget(QtWidgets.QLabel('after'))
        self.startTime = QtWidgets.QDateTimeEdit(self)
        self.layout().addWidget(self.startTime)
        self.layout().addWidget(QtWidgets.QLabel('and before'))
        self.endTime = QtWidgets.QDateTimeEdit(self)
        self.layout().addWidget(self.endTime)

    def loadClassCombo(self,model=None):
        if model is None:
            return
        self.classCombo.clear()
        for item in model.classesInReport():
            self.classCombo.addItem(item)

    def initialiseTime(self,model=None):
        if model is None:
            return
        timeRange = model.getTimeRange()
        #print 'CErrorReportSelection.initialiseTime',timeRange
        start = QtCore.QDateTime()
        start.setTime_t(int(timeRange[0]))
        start.addMSecs(int((timeRange[0]-int(timeRange[0]))*1000))
        self.startTime.setDateTime(start)
        self.startTime.setMinimumDateTime(start)
        local = time.localtime()
        midnight = QtCore.QDateTime(QtCore.QDate(local.tm_year, local.tm_mon,local.tm_mday), QtCore.QTime(23, 59, 0))
        self.endTime.setDateTime(midnight)
        self.endTime.setMinimumDateTime(start)

class CErrorReportViewer(CCP4TextViewer.CTextViewer):

    TIME_FORMAT = '%H:%M:%S'

    def __init__(self,parent):
        CCP4TextViewer.CTextViewer.__init__(self, parent)
        self.setFont(style='fixed_width')
        self._model = None
        self.lastErrorIndex = -1

    def setModel(self,model=None):
        #print 'CErrorReportViewer.setModel',model,type(model)
        self._model = model

    def open(self,filename=None):
        self.loadText()
        self.fileName = None
        self.lastModTime= None

    def clear(self):
        self.viewer.clear()

    def loadText(self,start=0):
        text = ''
        indx = -1
        for indx in range(start, len(self._model)):
            report = self._model[indx]
            #print 'CErrorReportViewer.loadText',report
            if 'time' in report:
                timeText = time.strftime(CErrorReportViewer.TIME_FORMAT,time.localtime(report['time']))
            else:
                timeText = ''
            if hasattr(report['class'],'ERROR_CODES'):
                description = report['class'].ERROR_CODES[report['code']].get('description',' ')
                severity = report['class'].ERROR_CODES[report['code']].get('severity',SEVERITY_ERROR)
                if severity == SEVERITY_OK:
                    severityColour = '00FF00'
                elif severity == SEVERITY_ERROR:
                    severityColour = 'FF0000'
                else:
                    severityColour = 'FFC000'
                line = '''{0!s:12} {1:20} <font color=#{5:6}>{2:3}</font> {3:60} {4:}'''.format(timeText,report['class'].__name__, 
                                                                                                report['code'], description, report['details'],
                                                                                                severityColour)
                text = text + line + '<br>'
        self.viewer.append(text)
        self.lastErrorIndex = indx

    def update(self):
        self.loadText(start=self.lastErrorIndex + 1)


class CSendJobError(QtWidgets.QDialog):
    def __init__(self, parent=None, projectId=None, projectName=None):
        from qtgui import CCP4Widgets
        QtWidgets.QDialog.__init__(self, parent=parent)
        self.projectId = projectId
        self.projectName = projectName
        self.setWindowTitle('Send bug report to CCP4')
        self.setModal(True)
        self.setLayout(QtWidgets.QVBoxLayout())
        line = QtWidgets.QHBoxLayout()
        line.addWidget(QtWidgets.QLabel('Send', self))
        self.sendLog = QtWidgets.QCheckBox('log files', self)
        self.sendLog.setChecked(True)
        line.addWidget(self.sendLog)
        self.sendData = QtWidgets.QCheckBox('data files', self)
        self.sendData.setChecked(False)
        line.addWidget(self.sendData)
        self.sendee = QtWidgets.QComboBox(self)
        for addr in  ['ccp4@ccp4.ac.uk']:
            self.sendee.addItem(addr)
        self.sendee.setEditable(True)
        line.addWidget(QtWidgets.QLabel('to',self))
        line.addWidget(self.sendee)
        self.layout().addLayout(line)
        self.layout().addWidget(QtWidgets.QLabel('for the job(s)..',self))
        self.jobListWidget = CCP4Widgets.CFinishedJobsListWidget(self,projectId)
        self.jobListWidget.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.jobListWidget.load()
        self.layout().addWidget(self.jobListWidget)
        self.layout().addWidget(QtWidgets.QLabel('Add comments..'))
        self.comments = QtWidgets.QTextEdit(self)
        self.layout().addWidget(self.comments)
        buttonBox = QtWidgets.QDialogButtonBox(self)
        button = buttonBox.addButton('Send',QtWidgets.QDialogButtonBox.ApplyRole)
        button.clicked.connect(self.applySend)
        button = buttonBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
        button.clicked.connect(self.close)
        line = QtWidgets.QHBoxLayout()
        line.addWidget(buttonBox)
        self.layout().addLayout(line)

    @QtCore.Slot()
    def applySend(self):
        from core import CCP4Utils
        from dbapi import CCP4DbUtils
        from core import CCP4Modules
        print('applySend')
        comments = self.comments.toPlainText().__str__()
        sendee = self.sendee.currentText().__str__()
        sendLog = self.sendLog.isChecked()
        sendData = self.sendData.isChecked()
        jobList = self.jobListWidget.selectedJobs()
        # Put everything to send in a temp directory
        dirName = self.projectName + '_' + str(int((time.time())))
        tmpDir = tempfile.mkdtemp()
        tarDir = os.path.join(tmpDir,dirName)
        print('Saving job(s) to ',tarDir)
        if sendLog or sendData:
            copyJobDir = CCP4DbUtils.CCopyJobDirectories(projectId=self.projectId, jobIdList=jobList, targetDir=tarDir, copyData=sendData)
            copyJobDir.copy()
        projectTmpDir = os.path.join(CCP4Modules.PROJECTSMANAGER().getProjectDirectory(projectId=self.projectId), 'CCP4_TMP')
        if not os.path.exists(projectTmpDir):
            os.mkdir(projectTmpDir)
        tarFile = os.path.join(projectTmpDir, dirName + '.tar.gz')
        rv = CCP4Utils.writeTarGzip(tarDir, tarFile=tarFile)
        message = 'mailto:' + sendee + '?subject=CCP4i2 Job Error Report&body= \n' + comments + '\n'
        if tarFile is not None:
            message = message + '\n\nPLEASE ATTACH FILE CONTAINING JOB DIAGNOSTIC:  ' + str(tarFile)
            if not sendData:
                message = message + '\nThe compressed file does not contain any data files but will contain the names of data files.'
        rv =  QtGui.QDesktopServices.openUrl(QtCore.QUrl(message))
        self.close()

class CPrintLogViewer(CCP4TextViewer.CTextViewer):

    def __init__(self, parent):
        CCP4TextViewer.CTextViewer.__init__(self, parent)
        self.setFont(style='fixed_width')

    def openThread(self, thread = 'main_thread'):
        from core import CCP4Modules
        from core import CCP4Utils
        versionText = CCP4Utils.versionLogHeader()
        ph = CCP4Modules.PRINTHANDLER()
        text = ph.getContent(thread)
        self.viewer.append(versionText + text)

