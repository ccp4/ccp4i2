from __future__ import print_function


"""
     CCP4MessageBox.py: CCP4 GUI Project
     Copyright (C) 2012 STFC

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
   Liz Potterton Feb 2012 -wrap QMessageBox
"""

##@package CCP4ProjectWidget View a project
                            
from PySide6 import QtGui, QtWidgets,QtCore

class CMessageBox(QtWidgets.QDialog):

  DEVELOPERS = [ [ 'Liz' , 'liz.potterton@york.ac.uk' ],
                 ['Andrey' , 'andrey.lebedev@stfc.ac.uk'] ,
                 [ 'Stuart', 'stuart.mcnicholas@york.ac.uk'] ]


  def __init__(self,parent=None,title=None,message='',exception=None,details=None,jobId=None,openJob=None):
    QtWidgets.QDialog.__init__(self,parent)
    if title is None: title = parent.window().windowTitle()
    self.message = message
    self.setWindowTitle(title)
    if openJob is not None and openJob.jobId is not None:
      self.openJob = openJob
    elif jobId is not None:
      from dbapi import CCP4DbUtils
      self.openJob = CCP4DbUtils.COpenJob(jobId=jobId)
    else:
      self.openJob = None

    self.setLayout(QtWidgets.QVBoxLayout())
    self.textWidget = QtWidgets.QTextEdit(self)
    self.textWidget.setObjectName('messageBox')
    self.textWidget.setPlainText(self.message)
    self.textWidget.setReadOnly(True)
    self.layout().addWidget(self.textWidget)

    buttons = QtWidgets.QDialogButtonBox(self)
    b = buttons.addButton('Report',QtWidgets.QDialogButtonBox.ActionRole)
    b.clicked.connect(self.handleReport)
    self.showButton = buttons.addButton('Show details',QtWidgets.QDialogButtonBox.ActionRole)
    self.showButton.clicked.connect(self.toggleDetails)
    b = buttons.addButton(QtWidgets.QDialogButtonBox.Close)
    b.clicked.connect(self.close)
    self.layout().addWidget(buttons)

    if exception is not None:
      import sys,traceback
      try:
        stack =  traceback.format_exc()
      except:
        stack = 'No traceback info recovered'
      self.details = str(exception) + '\n' + stack
    elif details is not None:
      self.details = details
    else:
      self.details = 'No details'
    

    self.setModal(True)
    self.show()
    self.raise_()

  @QtCore.Slot()
  def handleReport(self):
    self.reportDialog = QtWidgets.QDialog(self)
    self.reportDialog.setLayout(QtWidgets.QVBoxLayout())
    if self.openJob is not None:
      self.sendJobWidget = QtWidgets.QCheckBox('Send files for job number '+self.openJob.jobnumber,self)
      self.reportDialog.layout().addWidget(self.sendJobWidget)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Use archive format',self))
    self.archiveWidget = QtWidgets.QComboBox(self)
    for item in ['zip','gztar','tar']: self.archiveWidget.addItem(item)
    line.addWidget(self.archiveWidget)
    self.reportDialog.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Send to',self))
    self.devWidget = QtWidgets.QComboBox(self)
    for dev,adr in CMessageBox.DEVELOPERS:
      self.devWidget.addItem(dev)
    line.addWidget(self.devWidget)
    self.reportDialog.layout().addLayout(line)
    buttons = QtWidgets.QDialogButtonBox(self)
    self.reportDialog.layout().addWidget(buttons)
    b = buttons.addButton('Send',QtWidgets.QDialogButtonBox.ActionRole)
    b.clicked.connect(self.send)
    b = buttons.addButton(QtWidgets.QDialogButtonBox.Cancel)
    b.clicked.connect(self.reportDialog.close)
    self.reportDialog.show()
    self.reportDialog.raise_()
                 
    
  @QtCore.Slot()
  def send(self):
    import os,sys,shutil
    from core import CCP4Modules
    from dbapi import CCP4DbApi
    selectedAdr = None
    selectedDev = self.devWidget.currentText().__str__()
    for dev,adr in  CMessageBox.DEVELOPERS:
      if selectedDev == dev : selectedAdr = adr
    if self.openJob is not None and self.sendJobWidget.isChecked():
      #print 'CMessageBox.send jobId',self.openJob.jobId
      jobDir =  os.path.split(CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=self.openJob.jobId,mode='ROOT'))[0]
      #Copy input files into the same directory
      fileList = CCP4Modules.PROJECTSMANAGER().db().getJobFiles(jobId=self.openJob.jobId,role=CCP4DbApi.FILE_ROLE_IN,mode='fullPath')
      inputFilesDir = os.path.join(jobDir,'INPUT_FILES')
      if not os.path.exists(inputFilesDir): os.mkdir(inputFilesDir)
      for inpFile in fileList:
        dstFile = os.path.join(inputFilesDir,os.path.split(inpFile)[1])
        if not os.path.exists(dstFile):
          try:
            shutil.copyfile(inpFile,dstFile)
          except:
            print('ERROR copying file '+inpFile+' to '+dstFile)
      # Make a zip archive
      archFormat = self.archiveWidget.currentText().__str__()
      if archFormat == 'gztar':
        zipFile = jobDir + '.tar.gz'
      else:
        zipFile = jobDir + '.' + archFormat
      try:
        shutil.make_archive(jobDir,archFormat,jobDir)
      except:
        print('ERROR making archive of job directory '+zipFile)
      # Should we delete the inputFilesDir?
    else:
      zipFile = None
    # Get version / os info
    from core import CCP4File,CCP4Utils
    version = CCP4File.CI2XmlHeader()
    version.loadFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(),'core','version.params.xml'))
    print('version',version)
    versionText = 'CCP4i2 version: ' + str(version.ccp4iVersion) + '\n' + \
                  'SVN version:  ' +  str(version.pluginVersion) + '\n' + \
                  'Creation date: ' + str(version.creationTime) + '\n' + \
                  'Platform: ' + sys.platform +  '\n'
    # Stick it all in a mailto
    message = 'mailto:'+selectedAdr+'?subject=CCP4i2 Problem&body= \n' + versionText + '\n' + \
              self.message + '\n' + self.details + '\n'
    if zipFile is not None: message = message + '\n\nPLEASE ATTACH FILE CONTAINING JOB DATA:'+zipFile + \
        '\nThe compressed file does not contain any data files but will contain the names of data files.'
    rv =  QtGui.QDesktopServices.openUrl( QtCore.QUrl(message) )
    self.reportDialog.close()
    self.close()
    

    
  @QtCore.Slot()
  def toggleDetails(self):
    if str(self.showButton.text())[0:4] == 'Show':
      self.showButton.setText('Hide details')
      self.textWidget.setReadOnly(False)
      self.textWidget.setPlainText(self.details)
      self.textWidget.setReadOnly(True)
    else:
      self.showButton.setText('Show details')
      self.textWidget.setReadOnly(False)
      self.textWidget.setPlainText(self.message)
      self.textWidget.setReadOnly(True)

  
