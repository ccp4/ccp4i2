from __future__ import print_function


"""
     CCP4ProjectWidget.py: CCP4 GUI Project
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
   Liz Potterton Feb 2010 - Copied from earlier review.py used for Developers meeting demo
"""

##@package CCP4ProjectWidget (QtGui) The demo project view widget
                            
from PySide6 import QtGui, QtWidgets,QtCore
from core.CCP4Modules import WEBBROWSER

FONT_SIZE = 14

STATUS_TEXT = { -1: 'unknown', 0 :'queued',1:'running',2:'cancelled',3:'finished',4:'forwarded' ,5:'terminated'}
STATUS_FINISHED = [3,4,5]
STAGES_TITLES = { 'DP' : 'Data processing','DR' : 'Data reduction','PP' : 'Phasing preparation','PH' : 'Phasing','PR' : 'Phase refinement', 'BB' : 'Building','RF' : 'Refinement'} 
STAGES_CODE_LIST = ['DP' , 'DR' , 'PP', 'PH', 'PR', 'BB', 'RF' ]
COLUMN_ORDER = ['parent_id','label','priority','estimated_time','confidence','status','action']
COLUMN_TITLES = { 'id' : 'ID',
                  'parent_id': 'Parent',
                  'label' : 'Description',
                  'estimated_time' : 'Time',
                  'estimated_confidence' : 'Confidence',
                  'priority' : 'Priority',
                  'confidence' : 'Quality',
                  'status' : 'Status',
                  'action' : 'Action' }

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
class CProjectWidgetDemo(QtWidgets.QFrame):
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

  MARGIN = 1
  
  def __init__(self,parent=None):
      QtWidgets.QFrame.__init__(self)

      layout = QtWidgets.QVBoxLayout()
      layout.setMargin(CProjectWidgetDemo.MARGIN)
      layout.setContentsMargins(CProjectWidgetDemo.MARGIN,CProjectWidgetDemo.MARGIN,
                                CProjectWidgetDemo.MARGIN,CProjectWidgetDemo.MARGIN)
      layout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)

      
      top_line_layout = QtWidgets.QHBoxLayout()
      l = QtWidgets.QLabel('Show jobs for',self)
      l.font().setPointSize(FONT_SIZE)
      top_line_layout.addWidget(l)
      self.stage_selector = QtWidgets.QComboBox(self)
      for item in STAGES_CODE_LIST: self.stage_selector.addItem(STAGES_TITLES[item])
      top_line_layout.addWidget(self.stage_selector)
      
      '''
      l = QtWidgets.QLabel('following from job')
      l.font().setPointSize(FONT_SIZE)
      top_line_layout.addWidget(l)
      self.parent_selector = QtWidgets.QLineEdit(self)
      self.parent_selector.setMaximumWidth(40)
      top_line_layout.addWidget(self.parent_selector)
      self.show_low_scoring = QtWidgets.QCheckBox('Show low scoring',self)
      top_line_layout.addWidget(self.show_low_scoring)
      self.show_low_scoring.setChecked(1)
      '''
  
      top_line_layout.addStretch(5)
      layout.addLayout(top_line_layout)
      
      
      self.tableFrame = QtWidgets.QFrame(self)
      tlayout =  QtWidgets.QVBoxLayout()
      tlayout.setMargin(CProjectWidgetDemo.MARGIN)
      tlayout.setContentsMargins(CProjectWidgetDemo.MARGIN,CProjectWidgetDemo.MARGIN,
                                 CProjectWidgetDemo.MARGIN,CProjectWidgetDemo.MARGIN)
      self.tableFrame.setLayout(tlayout)
      layout.addWidget(self.tableFrame)
      
      self.setLayout(layout)
     
      self.jobsModel = JobTableModel()
      self.table = None


#------------------------------------------------------------------------------------------------------
  def setProjectData(self,filename):
#------------------------------------------------------------------------------------------------------

      global_data = {}
      local_data = {}
      exec(compile(open(filename).read(), filename, 'exec'),global_data,local_data)
      self.DUMMYDATA = {}
      for key,value in list(local_data.get('DUMMYDATA',{}).items()):
        self.DUMMYDATA[key] = value
      # 'CProjectWidgetDemo.setProjectData',self.DUMMYDATA

      self.loadDummy(0)
      
      self.setTable()

      self.stage_selector.currentIndexChanged[int].connect(self.handleStageSelector)
      self.jobsModel.priorityChanged.connect(self.handlePriorityChange)
      


#------------------------------------------------------------------------------------------------------
  def setTable(self):
#------------------------------------------------------------------------------------------------------
      if self.table:
        self.tableFrame.layout().removeWidget(self.table)
        #self.layout().removeWidget(self.table)
        self.table.close()
        self.table = None
      
      self.table = QtWidgets.QTableView(self)
      #self.table.setAttribute(QtCore.Qt.WA_DeleteOnClose)
      self.table.setModel(self.jobsModel)
      delegate = JobDelegate(self)
      self.table.setItemDelegate(delegate)
      delegate.buttonClicked.connect(self.handleShow)

      self.setPersistentEditors()
      self.table.resizeColumnsToContents()
      self.tableFrame.layout().addWidget(self.table)
      #self.layout().addWidget(self.table)
      self.table.show()
      #print 'CProjectReview.setTable',self.table.size().width(),self.table.size().height()
      #self.table.setEditTriggers(self.table.AllEditTriggers)
      #self.tableFrame.show()
      #self.show()


#------------------------------------------------------------------------------------------------------
  def tableSize(self):
#------------------------------------------------------------------------------------------------------
    width = self.table.verticalHeader().width()
    height = self.table.horizontalHeader().height()
    for ic in range(0,self.jobsModel.columnCount()):
      width = width + self.table.columnWidth(ic)

    for ir in range(0,self.jobsModel.rowCount()):
      height = height + self.table.rowHeight(ir)
   
    print('CProjectWidgetDemo.tableSize',width,height)    
    return width,height

#------------------------------------------------------------------------------------------------------
  def Size(self):
#------------------------------------------------------------------------------------------------------
    width,height = self.tableSize()
    width = width + 4*CProjectWidgetDemo.MARGIN
    height = height + 8*CProjectWidgetDemo.MARGIN + self.stage_selector.height()
    return  width,height
    
      
#------------------------------------------------------------------------------------------------------
  def setPersistentEditors(self):
#------------------------------------------------------------------------------------------------------
    for mode in ['priority','action']:
      column = COLUMN_ORDER.index(mode)
      for row in range(self.jobsModel.rowCount()):
        self.table.openPersistentEditor(self.jobsModel.index(row,column))

    
#------------------------------------------------------------------------------------------------------
  def runDemo(self,id):
#------------------------------------------------------------------------------------------------------
    if id != 104: return
    self.loadDummy(1)
    import functools
    QtCore.QTimer.singleShot(1000, functools.partial(self.loadDummy,2))
    QtCore.QTimer.singleShot(3000, functools.partial(self.loadDummy,3))
    QtCore.QTimer.singleShot(8000, functools.partial(self.loadDummy,4))
    
#------------------------------------------------------------------------------------------------------
  def loadDummy(self,step=0):
#------------------------------------------------------------------------------------------------------
    #print 'loadDummy',step
    self.jobsModel.loadJobs(self.DUMMYDATA.get(step))
    self.jobsModel.setVisibleJobs()
    if self.table: self.table.reset()
      
#------------------------------------------------------------------------------------------------------
  @QtCore.Slot(int)
  def handleStageSelector(self,istage):
#------------------------------------------------------------------------------------------------------
    self.jobsModel.setVisibleJobs(stage=STAGES_CODE_LIST[istage])
    self.setTable()

#------------------------------------------------------------------------------------------------------
  def handleLowScoring(self,low):
#------------------------------------------------------------------------------------------------------
    #print 'handleLowScoring',low
    self.jobsModel.setVisibleJobs(filter_low_scoring=1-int(self.show_low_scoring.isChecked()))
    self.setTable()

#------------------------------------------------------------------------------------------------------
  @QtCore.Slot(tuple)
  def handlePriorityChange(self,args):
#------------------------------------------------------------------------------------------------------
    #print 'handlePriorityChange',args
    id,priority = args
    if priority>=0.99:
      self.runDemo(id)

#------------------------------------------------------------------------------------------------------
  @QtCore.Slot('QModelIndex')
  def handleShow(self,index):
#------------------------------------------------------------------------------------------------------
    idx = self.jobsModel.getIdFromIndex(index)
    #print 'CProjectReview.handleShow',idx
    if idx == 104:
      #self.openCCP4i('XVNTGE')
      from core.CCP4Utils import getCCP4I2Dir
      import os
      filename=os.path.join(getCCP4I2Dir(),'test','data','test_job_104.data.xml')
      from core.CCP4TaskManager import TASKMANAGER
      TASKMANAGER().openDataFile(fileName=filename,editableData=False)
    elif idx == 401:
      self.openMG('/Users/lizp/I2/demo/dUTPase.pdb')
      
#------------------------------------------------------------------------------------------------------
  def openCCP4i(self,task):
#------------------------------------------------------------------------------------------------------
    process = QtCore.QProcess(self)
    rv = process.start('ccp4i -t '+task)
    #rv = process.start('ccp4i -t fft')

#------------------------------------------------------------------------------------------------------
  def openMG(self,pdbfile):
#------------------------------------------------------------------------------------------------------
    process = QtCore.QProcess(self)
    process.start('/Users/lizp/Desktop/dev/ccp4mg/bin/ccp4mg  -nore '+pdbfile)



#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
class Job:
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

  def __init__(self,defn):
    #print 'Job defn',defn
    self.id = defn[0]
    self.parent_id = defn[3]
    self.stage = defn[1]
    self.template = defn[2]
    self.label =  defn[4]
    self.estimated_time = defn[5]
    self.estimated_confidence = defn[6]
    self.priority = defn[7]
    self.confidence = defn[8]
    self.status = defn[9]

#------------------------------------------------------------------------------------------------------
  def setValue(self,name,value):
#------------------------------------------------------------------------------------------------------
    setattr(self,name,value)
    
#------------------------------------------------------------------------------------------------------
  def __cmp__(self,other):
#------------------------------------------------------------------------------------------------------
    if self.id > other.id:
      return 1
    elif self.id > other.id:
      return -1
    else:
      return 0



#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
class JobTableModel(QtCore.QAbstractTableModel):
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

  priorityChanged = QtCore.Signal(tuple)

  def __init__(self):
    QtCore.QAbstractTableModel.__init__(self)
    self.jobs = {}
    self.stage = 'PH'
    self.parent_id = 0
    self.filter_low_scoring =0
    self.visibleJobs = []

#------------------------------------------------------------------------------------------------------
  def setVisibleJobs(self,stage='',parent_id=-1, filter_low_scoring=-1):
#------------------------------------------------------------------------------------------------------
     if stage: self.stage = stage
     if parent_id>=0: self.parent_id = parent_id
     if filter_low_scoring>=0: self.filter_low_scoring = filter_low_scoring
     self.visibleJobs = []
     for id,job in list(self.jobs.items()):
       if job.stage == self.stage:
         if  (not self.filter_low_scoring) or (STATUS_FINISHED.count(job.status) and job.confidence>0.3) or  (STATUS_FINISHED.count(job.status)<=0 and job.estimated_confidence>0.3):
           self.visibleJobs.append(id)
     #print 'setVisibleJobs',self.visibleJobs

#------------------------------------------------------------------------------------------------------
  def rowCount(self,index=QtCore.QModelIndex()):
#------------------------------------------------------------------------------------------------------
    #print 'rowCount',len(self.visibleJobs)
    return len(self.visibleJobs)
  
#------------------------------------------------------------------------------------------------------
  def columnCount(self,index=QtCore.QModelIndex()):
#------------------------------------------------------------------------------------------------------
    return len(COLUMN_ORDER)
  
#------------------------------------------------------------------------------------------------------
  def data(self, index, role = QtCore.Qt.DisplayRole):
#------------------------------------------------------------------------------------------------------
    if not index.isValid() or not ( 0 <= index.row() < self.rowCount() ):
#FIXME PYQT - or maybe None? This used to return QVariant.
      return None
    job = self.jobs [ self.visibleJobs[index.row()]]
    name = COLUMN_ORDER[index.column()]
    
    if role == QtCore.Qt.DisplayRole:
      if ['id','parent_id','estimated_time'].count(name):
        if name == 'estimated_time' and STATUS_FINISHED.count(job.status):
#FIXME PYQT - or maybe None? This used to return QVariant.
          return None
        else:
          return str(getattr(job,name,0))
      elif ['estimated_confidence','priority','confidence'].count(name):
        if  name == 'confidence' and not STATUS_FINISHED.count(job.status):
          value = getattr(job,'estimated_confidence',-1.0)
        else:
          value = getattr(job,name,-1.0)
        if value<0.0:
          return ' '
        elif ['estimated_confidence','priority'].count(name) and  STATUS_FINISHED.count(job.status):
#FIXME PYQT - or maybe None? This used to return QVariant.
          return None
        else:
          return "%4.1f"%value
      elif name == 'status':
        return STATUS_TEXT[getattr(job,name,-1)]
      else:
        return getattr(job,name,' ')

    elif role == QtCore.Qt.BackgroundColorRole:
      if job.status == 1:
        return QtGui.QColor('navajowhite')
      elif STATUS_FINISHED.count(job.status):
        return QtGui.QColor('paleturquoise')
      elif  ['confidence','estimated_time'].count(name) and not STATUS_FINISHED.count(job.status):
        return QtGui.QColor('lightgrey')
      else:
        return QtGui.QColor('white')

#FIXME PYQT - or maybe None? This used to return QVariant.
    return None

#------------------------------------------------------------------------------------------------------
  def setData(self, index, value, role=QtCore.Qt.EditRole):
#------------------------------------------------------------------------------------------------------
    id = self.getIdFromIndex(index)
    if id<0: return
    job = self.jobs [ id ]
    #print 'JobTableModel.setData',job.id,value, value.toDouble()
    if COLUMN_ORDER[index.column()] == 'priority':
      job.priority = value
      self.priorityChanged.emit(job.id,job.priority)

#------------------------------------------------------------------------------------------------------
  def getData(self,index,name):
#------------------------------------------------------------------------------------------------------
    job = self.jobs [ self.visibleJobs[index.row()]]
    return getattr(job,name,0.0)


#------------------------------------------------------------------------------------------------------
  def getIdFromIndex(self,index):
#------------------------------------------------------------------------------------------------------
    if index.isValid() and 0 <= index.row() < len(self.visibleJobs):
      id = self.visibleJobs[index.row()]
      return id
    else:
      return -1
  
#------------------------------------------------------------------------------------------------------
  def getIdData(self,id,name):
#------------------------------------------------------------------------------------------------------
    if id in self.jobs:
      return getattr(self.jobs[id],name,0.0)
    else:
      return 0.0

#------------------------------------------------------------------------------------------------------
  def headerData(self,section, orientation, role = QtCore.Qt.DisplayRole):
#------------------------------------------------------------------------------------------------------
    if role ==  QtCore.Qt.DisplayRole:
      if orientation == QtCore.Qt.Horizontal:
        return COLUMN_TITLES[COLUMN_ORDER[section]]
      else:
        if 0 <= section < len(self.visibleJobs):
          job = self.jobs [ self.visibleJobs[section]]
          return str(getattr(job,'id',0))
        else:
          pass
          #print 'headerData section',section

#FIXME PYQT - or maybe None? This used to return QVariant.
    return None
    
#------------------------------------------------------------------------------------------------------
  def loadJobs(self,datalist=[]):
#------------------------------------------------------------------------------------------------------
    for item in datalist:
      self.jobs[item[0]] = Job(item)
    print('Loaded',len(self.jobs),'jobs')

  

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
class JobDelegate(QtWidgets.QItemDelegate):
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

  buttonClicked = QtCore.Signal('QModelIndex')

  def __init__(self,parent=None):
    super(JobDelegate, self).__init__(parent)

#------------------------------------------------------------------------------------------------------
  def createEditor(self,parent, option, index):
#------------------------------------------------------------------------------------------------------
    #print 'JobDelegate.createEditor',parent, option, index,index.column()
    import functools
    if index.column() == COLUMN_ORDER.index('priority'):
      spinbox = QtWidgets.QDoubleSpinBox(parent)
      spinbox.setRange(0.0,1.0)
      spinbox.setSingleStep(0.1)
      spinbox.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignVCenter)
      return spinbox
    if index.column() == COLUMN_ORDER.index('action'):
      button = QtWidgets.QPushButton('Show',parent)
      button.clicked.connect(functools.partial(self.handleButtonClick,index))
      id = index.row()
      return button
    else:
      return QtWidgets.QItemDelegate.createEditor(self,parent, option, index)

#------------------------------------------------------------------------------------------------------
  def setEditorData ( self, editor, index) :
#------------------------------------------------------------------------------------------------------
    if index.column() == COLUMN_ORDER.index('priority'):
      editor.setValue(index.model().getData(index,'priority'))
    else:
      QtWidgets.QItemDelegate.setEditorData( self, editor, index)

#------------------------------------------------------------------------------------------------------
  def setModelData(self,editor,model,index):
#------------------------------------------------------------------------------------------------------
    #print 'setModelData',editor,model,index
    if index.column() == COLUMN_ORDER.index('priority'):
      model.setData(index,editor.value())


#------------------------------------------------------------------------------------------------------
  @QtCore.Slot('QModelIndex')
  def handleButtonClick(self,index):
#------------------------------------------------------------------------------------------------------
    #print 'JobDelegate.handleButtonClick',index
    self.buttonClicked.emit(index)
    
   
#===========================================================================


if __name__ == "__main__":
    from core.CCP4Modules import QTAPPLICATION
    app = QTAPPLICATION()
    window = QtWidgets.QDialog()
    window.setWindowTitle('CCP4 Project: Coseners')
    #print 'font' , window.font().setPointSize(FONT_SIZE)
    layout = QtWidgets.QVBoxLayout()
    widget = CProjectReview()
    layout.addWidget(widget)
    window.setLayout(layout)
    window.show()

    widget.stage_selector.setCurrentIndex(3)
    sys.exit(app.exec_())
