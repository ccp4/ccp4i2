"""
Copyright (C) 2011 University of York
Liz Potterton Jan 2011 - List program error log
"""

##@package CCP4ProjectWidget View a project

import os
import sys
import tempfile
import time

from PySide2 import QtCore, QtGui, QtWidgets

from . import CCP4TextViewer
from .. import __version__
from ..core import CCP4Utils
from ..core.CCP4Modules import PROJECTSMANAGER
from ..core.CCP4PrintHandler import PRINTHANDLER
from ..core.CCP4Version import CCP4_VERSION


class CSendJobError(QtWidgets.QDialog):
    def __init__(self, parent=None, projectId=None, projectName=None):
        from . import CCP4Widgets
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
        from ..dbapi import CCP4DbUtils
        print('applySend')
        comments = str(self.comments.toPlainText())
        sendee = str(self.sendee.currentText())
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
        projectTmpDir = os.path.join(PROJECTSMANAGER().getProjectDirectory(projectId=self.projectId), 'CCP4_TMP')
        if not os.path.exists(projectTmpDir):
            os.mkdir(projectTmpDir)
        tarFile = os.path.join(projectTmpDir, dirName + '.tar.gz')
        CCP4Utils.writeTarGzip(tarDir, tarFile=tarFile)
        message = 'mailto:' + sendee + '?subject=CCP4i2 Job Error Report&body= \n' + comments + '\n'
        if tarFile is not None:
            message = message + '\n\nPLEASE ATTACH FILE CONTAINING JOB DIAGNOSTIC:  ' + str(tarFile)
            if not sendData:
                message = message + '\nThe compressed file does not contain any data files but will contain the names of data files.'
        QtGui.QDesktopServices.openUrl(QtCore.QUrl(message))
        self.close()


class CPrintLogViewer(CCP4TextViewer.CTextViewer):

    def __init__(self, parent):
        CCP4TextViewer.CTextViewer.__init__(self, parent)
        self.setFont(style='fixed_width')

    def openThread(self, thread = 'main_thread'):
        versionText = (
            f"CCP4i2 version: {__version__}\n"
            f"Running CCP4 version: {CCP4_VERSION}\n"
            f"Using Python version: {sys.version}\n"
            f"Using Qt version: {QtCore.qVersion()}\n"
        )
        ph = PRINTHANDLER()
        text = ph.getContent(thread)
        self.viewer.append(versionText + text)
