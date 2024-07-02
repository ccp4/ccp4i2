from __future__ import print_function

"""
     qtgui/CProcessManagerWidget.py: CCP4 GUI Project
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2010 University of York

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
    May 2010 Liz Potterton - copied from graphical part of ccp4mg/qtgui/JobControl.py
"""
##@package CCP4ProcessManagerWidget (QtGui) Placeholder for manager of running processes

from PySide6 import QtCore,QtGui, QtWidgets
from core.CCP4Modules import *
from core import CCP4ProcessManager

class CProcessManagerWidget(QtWidgets.QWidget):

  def __init__(self,parent=None,processManager=None):    
    QtWidgets.QWidget.__init__(self,parent)
    self.processManager = processManager
    '''
    MAINWINDOW().addMenuDefinition('job_control', { 'text' : 'Review jobs',
                                                    'tip' : 'Show details of external program runs',
                                                    'slot' : self.openReview },
                                   'Tools')
    '''


#------------------------------------------------------------------------------
  def handle_progress_gui(self,**keywords):
#------------------------------------------------------------------------------
    #print "in handle_progress_gui"
    pass
  
  

#--------------------------------------------------------------------
  def openReview(self):
#--------------------------------------------------------------------
    import MGSimpleDialog
    if not self.reviewWindow:
      self.reviewWindow = jobReviewGui(MAINWINDOW())
      self.reviewWindow.showSignal.connect(self.showStatus)
    self.reviewWindow.show()
    self.reviewWindow.raise_()

#--------------------------------------------------------------------
  @QtCore.Slot()
  def showStatus(self):
#--------------------------------------------------------------------
    root = self.reviewWindow.selectedJob()
    if root not in self.pending_jobs: return
    cmd = self.pending_jobs[root]['command'] + ' '
    for item in self.pending_jobs[root]['arglist']: cmd = cmd + ' ' + item
    self.reviewWindow.loadText( root, self.pending_jobs[root]['name'],
                                'Run command: ' + cmd +'\n\n' +
                                'Finish status: '+str(self.pending_jobs[root]['status'])+'\n\n' +
                                'Standard error:\n'+self.pending_jobs[root]['error']+'\n\n' +
                                 'Standard out:\n'+self.pending_jobs[root]['output'],
                                self.pending_jobs[root]['logfile'])
    
    
#--------------------------------------------------------------------
class CProgress( QtWidgets.QDialog ):
#--------------------------------------------------------------------
  
#--------------------------------------------------------------------
  def __init__(self,text,title='Progress',help=''):
#--------------------------------------------------------------------

    QtWidgets.QDialog.__init__(self,MAINWINDOW())

    self.cancelCall = ''
    
    self.setWindowTitle(title)
    self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
    
    layout = QtWidgets.QVBoxLayout()

    for line in text.split('\n'):
      widget = QtWidgets.QLabel(line,self)
      layout.addWidget(widget)
    
    widget = QtWidgets.QProgressBar(self)
    widget.setRange(0,0)
    widget.setObjectName('progress')
    layout.addWidget(widget)

    widget = QtWidgets.QPushButton('Cancel job',self)
    widget.clicked.connect(self.cancelProgress)
    layout.addWidget(widget)
    
    
    self.setLayout(layout)
    self.show()


#--------------------------------------------------------------------
  def setCancel(self,cancelCall):
#--------------------------------------------------------------------
    self.cancelCall = cancelCall

#--------------------------------------------------------------------
  def setProgress(self,value):
#--------------------------------------------------------------------
    widget = self.findChild(QtWidgets.QProgressBar,'progress')
    if widget: widget.setValue(value)

#--------------------------------------------------------------------
  @QtCore.Slot()
  def cancelProgress(self):
#--------------------------------------------------------------------
    if  self.cancelCall: self.cancelCall()
    self.close()
    
#--------------------------------------------------------------------    
class mgProgressDialog(QtWidgets.QProgressDialog):
#--------------------------------------------------------------------
  def __init__(self,labelText='',cancelButtonText='',minimum=1,maximum=100,parent=None,title='',help=''):
    if not parent: parent = MAINWINDOW()
    QtWidgets.QProgressDialog.__init__(self,labelText,cancelButtonText,minimum,maximum,parent)
    self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
    self.setWindowTitle(title)
    

#--------------------------------------------------------------------
class jobReviewGui(QtWidgets.QDialog):
#--------------------------------------------------------------------

  showSignal = QtCore.Signal()

  def __init__(self,parent):
    QtWidgets.QDialog.__init__(self,parent)
    self.setWindowTitle('CCP4mg: Review jobs')
    layout = QtWidgets.QVBoxLayout()
    frame_margin = 1

    self.select_frame =  QtWidgets.QGroupBox('Select job to review',self)
    frame_layout =  QtWidgets.QVBoxLayout()
    frame_layout.setContentsMargins(frame_margin,frame_margin,frame_margin,frame_margin)
    frame_layout.setSpacing(frame_margin)
    self.selectJob= QtWidgets.QListWidget(self)
    self.selectJob.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
    self.selectJob.setMinimumHeight(30)
    #self.selectJob.setDefaultFont(QtGui.QFont('Courier',12))
    font = self.selectJob.viewOptions().font
    font.setFamily('Courier')
    import guiUtils
    self.selectJob.currentRowChanged.connect(self.showSignal.emit)
    frame_layout.addWidget(self.selectJob)
    self.select_frame.setLayout(frame_layout)
    
    self.error_frame =  QtWidgets.QGroupBox('Program output',self)
    frame_layout =  QtWidgets.QVBoxLayout()
    frame_layout.setContentsMargins(frame_margin,frame_margin,frame_margin,frame_margin)
    frame_layout.setSpacing(frame_margin)
    self.error = QtWidgets.QTextEdit(self)
    self.error.setReadOnly(1)
    self.error.setMinimumWidth(600)
    self.error.setMinimumHeight(100)
    self.error.setFontFamily('Courier')
    frame_layout.addWidget(self.error)
    self.error_frame.setLayout(frame_layout)

    self.log_frame =  QtWidgets.QGroupBox('Log file',self)
    frame_layout =  QtWidgets.QVBoxLayout()
    frame_layout.setContentsMargins(frame_margin,frame_margin,frame_margin,frame_margin)
    frame_layout.setSpacing(frame_margin)
    self.logfile = QtWidgets.QTextEdit(self)    
    self.logfile.setReadOnly(1)
    self.logfile.setFontFamily('Courier')
    frame_layout.addWidget(self.logfile)
    self.log_frame.setLayout(frame_layout)

    self.splitter = QtWidgets.QSplitter(QtCore.Qt.Vertical,self)
    layout.addWidget(self.splitter)
    self.splitter.addWidget(self.select_frame)
    self.splitter.addWidget(self.error_frame)
    self.splitter.addWidget(self.log_frame)
    self.splitter.setStretchFactor(2,3)
    self.splitter.setSizes([2,2,3])

    dialog_buttons = QtWidgets.QDialogButtonBox(self)
    cancel_button = dialog_buttons.addButton('Close',QtWidgets.QDialogButtonBox.RejectRole)
    #self.showButton = dialog_buttons.addButton('Show job details',QtWidgets.QDialogButtonBox.ApplyRole)
    cancel_button.clicked.connect(self.close)
    #import guiUtils
    layout.addWidget(dialog_buttons)

    
    self.loadSelectJob()
    self.selectJob.show()
    self.setLayout(layout)
    

  def loadSelectJob(self):
    self.selectJob.blockSignals(True)
    self.selectJob.clear()
    for item in JOBCONTROL().get_jobs_list():
      widget = self.selectJob.addItem(item)
    self.selectJob.blockSignals(False)
      
  def addJob(self,item):
     self.selectJob.blockSignals(True)
     self.selectJob.addItem(item)
     self.selectJob.blockSignals(False)
     
  def selectedJob(self):
    irow = self.selectJob.currentRow()
    if irow<0: return ''
    text = str(self.selectJob.item(irow).text())
    return int(text.split()[0])
    
  def loadText(self,root=-1,name='',error='',logfile=''):
    self.error_frame.setTitle('Output from job number '+str(root)+' '+name)
    self.error.setText(error)
    log_text = ''
    import os
    if logfile:
      self.log_frame.setTitle('Log file: '+logfile)
      if not os.path.exists(logfile):
        log_text = log_text + 'Error - log file does not exists'
      else:
        try:
          f = open(logfile,'r')
          text = f.read()
          f.close()
          log_text = log_text + '\n' + text
        except:
          log_text = log_text + 'Error reading log file'
    else:
      self.log_frame.setTitle('No log file')

    self.logfile.setText(log_text)
    sizes = self.splitter.sizes()
    if sizes[0]>0.9:
      l_sel = self.selectJob.count()
      l_log = len(log_text.split('\n'))
      l_err = len(error.split('\n'))

  def updateSelectJob(self,root,item_text):
    row = root-1
    #print 'updateSelectJob',row,item_text,'count', self.selectJob.count()
    self.selectJob.blockSignals(True)
    if row< self.selectJob.count():
      item = self.selectJob.item(row)
      if item: item.setText(item_text)
    else:
      self.selectJob.addItem(item_text)
      
    self.selectJob.blockSignals(False)
