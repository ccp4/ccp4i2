"""
Copyright (C) 2016 STFC
Liz Potterton april 2016 - Gui for remote running
"""

from collections.abc import Callable
import functools
import os
import time

from ccp4mg.qtgui import UtilityThread
from PySide2 import QtCore, QtWidgets

from . import CCP4Widgets
from ..core import CCP4Annotation
from ..core import CCP4JobServer
from ..core import CCP4Utils
from ..core.CCP4ErrorHandling import CException


def JOBCONTROLLERGUI():
    if CServerParamsDialog.insts is None:
        CServerParamsDialog.insts = CServerParamsDialog()
    return  CServerParamsDialog.insts


class CServerParamsDialog(QtWidgets.QDialog):
    insts = None
    guiDef = { 'machine' : [ 'Server' , 'machine', 'Choose a server'],
               'username':  [ 'Username on server' ,  'Your username on server'],
               'password' :  [ 'Password on server','Your password on server' ],
               'keyFilename' : [ 'SSH key filename', 'Local file containing SSH public keys' ],
               'ccp4Dir' :  [ 'CCP4 directory on server', 'Full path of CCP4 on the server machine' ],
               'tempDir' :  [ 'Work directory on server', 'Full path of temporary work directory server machine' ],
               'sge_root' :  [ 'SGE_ROOT (Sun gridengine dir)', 'Full path of Sun grid engine root directory server machine' ],
               'queueOptionsFile' : [ 'Queue option file', 'File containing options for the job queue' ],
               'timeout' :  [ 'Timeout in seconds', 'Wait for connection' ]
             }
      
    def userGuiOptions(self,mechanism,validate=None,customCodeFile=None):
      if mechanism in ['ssh','ssh_shared','qsub_remote','qsub_shared','slurm_remote']:
        if validate is None or validate == 'password':
          ret = ['username','password','machine']
        elif validate == 'key_filename':
          ret = ['username','keyFilename','machine']
        elif validate == 'pass_key_filename':
          ret = ['username','keyFilename','password','machine']
        if mechanism in ['ssh','qsub_remote','slurm_remote']: ret.append('tempDir')
        #if mechanism in ['qsub_shared','qsub_remote']: ret.append('sge_root')
      elif mechanism == 'qsub_local':
        ret = []
      elif mechanism == 'custom':
        if customCodeFile is not None:
          from ..qtcore.CCP4JobController import JOBCONTROLLER
          return JOBCONTROLLER().customHandler(customCodeFile=customCodeFile).SHOW_USER_PARAMS
        else:
          return ['username','password']
      if mechanism in ['qsub_local','qsub_remote','qsub_shared','slurm_remote']:
          ret.append('queueOptionsFile')
      return ret

  
    def getParams(self):
      p = CCP4JobServer.CServerParams(machine = self.get('machine'),
                                             username = self.get('username'),
                                             password = self.get('password'),
                                             keyFilename = self.get('keyFilename'),
                                             mechanism = self.get('mechanism'),
                                             validate= self.get('validate'),
                                             ccp4Dir = self.get('ccp4Dir'),
                                             tempDir = self.get('tempDir'),
                                             sge_root = self.get('sge_root'),
                                             customCodeFile = self.get('customCodeFile'),
                                             queueOptionsFile = self.get('queueOptionsFile'),
                                             serverGroup=self.get('serverGroup'))
      print('CServerParamsDialog.getParams',p)
      return p
                                                 

    def get(self,param):
      if param not in self.widgets:
        serverGroup = self.container.get('SERVERGROUP'+str(self.modeGroup.checkedId()))
        if param == 'serverGroup':
          return str(serverGroup.name)
        elif param == 'machine':
          return str(getattr(serverGroup,'serverList')[0])
        elif param == 'timeout':
          if serverGroup.timeout.isSet():
            return float(serverGroup.timeout)
          else:
            return None
        elif param not in ['username','password']:
          return str(getattr(serverGroup,param))
        else:
          return ''
      
      if param in ['queueOptionsFile']:
        self.widgets['queueOptionsFile'].updateModelFromView()
        return str(self.queueOptionsFile)
      else:
        if isinstance(self.widgets[param],QtWidgets.QLineEdit):
          if param == 'timeout':
            if len(self.widgets[param].text()) == 0:
              return None
            else:
              return float(self.widgets[param].text())
          return str(self.widgets[param].text())
        else:
          return str(self.widgets[param].currentText())

    def setInfo(self,info=None):
      if info is None:
        self.widgets['info'].setText('')
      else:
        self.widgets['info'].setText(info)

    def  valid(self):
      return True
   
    def __init__(self,parent=None,params={}):
        self.widgets = { }
        self.labels = { }
        QtWidgets.QDialog.__init__(self,parent)
        try:
          from ..core.CCP4ServerSetup import SERVERSETUP
          self.container = SERVERSETUP()
        except:
          print('Failed loading server config file from CCP4I2/local_setup/servers_config.params.xml')
        
        self.setWindowTitle('Choose server to run job')
        layout = QtWidgets.QVBoxLayout()
        self.setLayout(layout)
        self.setMinimumWidth(500)


        frame = QtWidgets.QFrame(self)
        frame.setLayout(QtWidgets.QVBoxLayout())
        frame.setFrameStyle(QtWidgets.QFrame.Box | QtWidgets.QFrame.Plain)
        MARGIN = 6
        frame.layout().setSpacing(MARGIN)
        frame.layout().setContentsMargins(MARGIN,MARGIN,MARGIN,MARGIN)
        frame.layout().addWidget(CCP4Widgets.CItalicLabel('Choose group of servers to use..',self))
        self.layout().addWidget(frame)
        
        self.modeGroup = QtWidgets.QButtonGroup(self)
        default = 0
        row = 0
        for idx in range(1,10):
          group = self.container.get('SERVERGROUP'+str(idx))
          #if group is not None: print 'SERVERGROUP',idx,group, group.isSet(), len(group.serverList)
          if group is not None and len(group.serverList)>0:
            but = QtWidgets.QRadioButton(str(group.name),self)
            self.modeGroup.addButton(but,idx)
            frame.layout().addWidget(but)
            row = row + 1
            if default==0: default = idx
        #print 'default modeGroup',default
        try:            
          self.modeGroup.button(default).setChecked(True)
        except:
          print('Failed to set modeGroup',default)

        self.customFrame = QtWidgets.QFrame(self)
        self.customFrame.setLayout(QtWidgets.QVBoxLayout())
        self.layout().addWidget(self.customFrame)

        butBox = QtWidgets.QDialogButtonBox(self)
        but = butBox.addButton("Submit job",QtWidgets.QDialogButtonBox.AcceptRole)
        but.clicked.connect(self.accept)
        but = butBox.addButton("Help",QtWidgets.QDialogButtonBox.HelpRole)
        but.clicked.connect(self.help)
        but = butBox.addButton("Cancel",QtWidgets.QDialogButtonBox.RejectRole)
        but.clicked.connect(self.reject)       
        self.layout().addWidget(butBox)

        self.modeGroup.buttonClicked.connect(self.loadServers)
        self.loadServers()
        
    def drawDetails(self,serverGroup,params={}):
        self.widgets = { }
        self.labels = { }
        #params['ccp4Dir'] = params.get('ccp4Dir',CCP4Utils.getCCP4Dir())
        #params['username'] = params.get('username',CCP4Utils.getUserId())
        
        self.detailsFrame = QtWidgets.QFrame(self)
        layout = QtWidgets.QGridLayout()
        self.detailsFrame.setLayout(layout)
        self.detailsFrame.setFrameStyle(QtWidgets.QFrame.Box | QtWidgets.QFrame.Plain)
        layout.addWidget(CCP4Widgets.CItalicLabel('Details',self),0,0,1,2)
        row = 1
          
        #for param in ['machine', 'ccp4Dir', 'username', 'password','keyFilename','queueOptionsFile']:
        #for param in CServerParamsDialog.guiDef.keys():
        for param in self.userGuiOptions(mechanism=str(serverGroup.mechanism),validate=str(serverGroup.validate),
                                     customCodeFile=str(serverGroup.customCodeFile)):

          self.labels[param] = QtWidgets.QLabel(CServerParamsDialog.guiDef[param][0],self)
          layout.addWidget(self.labels[param] ,row,0)
          if param == 'machine':
            self.widgets[param] = QtWidgets.QComboBox(self)
          elif param == 'queueOptionsFile':
            from ..core import CCP4File
            self.queueOptionsFile = CCP4File.CDataFile(parent=self)
            self.widgets[param] = CCP4Widgets.CDataFileView(self,model=self.queueOptionsFile)
          else:
            self.widgets[param] = QtWidgets.QLineEdit(self)
            self.widgets[param].setText(params.get(param,''))
          self.widgets[param].setToolTip(CServerParamsDialog.guiDef[param][1])
          if param == 'password' : self.widgets[param].setEchoMode(QtWidgets.QLineEdit.Password)
          layout.addWidget(self.widgets[param],row,1)
          row = row + 1

        if serverGroup.mechanism == 'custom':
          for i in range(self.customFrame.layout().count()):
            try:
              w = self.customFrame.layout().takeAt(0).widget()
              w.hide()
              w.deleteLater()
            except:
              pass
          from ..qtcore.CCP4JobController import JOBCONTROLLER
          guiFunction = getattr(JOBCONTROLLER().customHandler(customCodeFile=serverGroup.customCodeFile),'guiFrame',None)
          if guiFunction is not None and isinstance(guiFunction, Callable):
            try:
              guiFunction(parentFrame=self.customFrame)
            except Exception as e:
              print('Failed drawing custom widgets',str(e))
          self.customFrame.show()
        else:
          self.customFrame.hide()


    @QtCore.Slot()
    def loadServers(self):
      serverGroup = self.container.get('SERVERGROUP'+str(self.modeGroup.checkedId()))
      if serverGroup is None:
        warningMess = QtWidgets.QMessageBox.warning(self,'Run on server','''To run on server you must first configure servers.
See under Utilities -> System administrator tools''')
        self.hide()
        return
        
      if hasattr(self,'detailsFrame'):
        idx = self.layout().indexOf(self.detailsFrame)
        self.layout().takeAt(idx)
        self.detailsFrame.hide()
        self.detailsFrame.deleteLater()
      else:
        idx = 1
      self.drawDetails(serverGroup)
      self.layout().insertWidget(idx,self.detailsFrame)
      
      #print 'loadServers',self.modeGroup.checkedId(),serverGroup,len(serverGroup.serverList)
      if 'machine' in self.widgets:
        self.widgets['machine'].setEditable(bool(serverGroup.userExtensible))
        self.widgets['machine'].clear()
        for item in serverGroup.serverList:
          self.widgets['machine'].addItem(str(item))

      for key in ['ccp4Dir','keyFilename','tempDir']:
        if key in self.widgets:
          if getattr(serverGroup,key).isSet():
            self.widgets[key].setText(str(getattr(serverGroup,key)))
          else:
            self.widgets[key].setText('')
            
      if 'password' in self.widgets:
        if serverGroup.validate is None or str(serverGroup.validate) == 'password':
          self.labels['password'].setText('Your password on server')
          self.widgets['password'].setToolTip('Your password on server')
        elif str(serverGroup.validate) == 'key_filename':
          pass
        elif str(serverGroup.validate) == 'pass_key_filename':
          self.labels['password'].setText('Your key file password')
          self.widgets['password'].setToolTip('Your key file password')

    @QtCore.Slot()
    def accept(self):
        QtWidgets.QDialog.accept(self)

    @QtCore.Slot()
    def reject(self):
        QtWidgets.QDialog.reject(self)

    def setInfo(self,text):
      if 'info' not in self.widgets:
         self.widgets['info'] = QtWidgets.QLabel(self)
         self.widgets['info'].setObjectName('warning')
         self.layout().insertWidget(0,self.widgets['info'])
      if text is None:
        self.widgets['info'].setText('')
      else:
        self.widgets['info'].setText(text)

            
    @QtCore.Slot()
    def help(self):
      from .CCP4WebBrowser import WEBBROWSER
      WEBBROWSER().loadWebPage(helpFileName='servers.html',newTab=True)

class CServerGroupView(CCP4Widgets.CComplexLineWidget):
    
  MODEL_CLASS = CCP4Annotation.CServerGroup
  
  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = { 'vboxLayout' : True }
    qualis.update(qualifiers)
    CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualis)
    self.layout().takeAt(0).widget().deleteLater()
    self.frames = {}
    self.widgets['name'] = CCP4Widgets.CStringView(self)
    qualis 
    self.widgets['mechanism'] = CCP4Widgets.CStringView( self,qualifiers=model.mechanism.qualifiers() )
    self.widgets['validate'] = CCP4Widgets.CRadioButtonGroup(self)
    self.widgets['validate'].addRadioButton('password','by password')
    self.widgets['validate'].addRadioButton('key_filename','by key file')
    self.widgets['validate'].addRadioButton('pass_key_filename','by passworded key file')
    self.widgets['mechanism'].widget.currentIndexChanged[int].connect(self.changeMechanism)
    self.widgets['validate'].buttonReleased.connect(self.changeValidate)
    '''
    if model is not None:
      self.widgets['mechanism'].populate(menuItems=model.mechanism.qualifiers('enumerators'),menuText=model.mechanism.qualifiers('menuText'))
    '''
    self.widgets['serverList'] = CCP4Widgets.CListView(self, qualifiers= { 'editorClassName' : 'CStringView' })
    self.widgets['userExtensible'] = CCP4Widgets.CBooleanView(self)
    self.widgets['ccp4Dir'] = CCP4Widgets.CStringView(self,qualifiers={  'allowUndefined' : True } )
    
    self.widgets['tempDir'] = CCP4Widgets.CStringView(self,qualifiers={  'allowUndefined' : True } )
    self.widgets['sge_root'] = CCP4Widgets.CStringView(self,qualifiers={  'allowUndefined' : True } )
    self.widgets['keyFilename'] = CCP4Widgets.CStringView(self,qualifiers={  'allowUndefined' : True } )
    self.widgets['timeout'] = CCP4Widgets.CFloatView(self,qualifiers={  'allowUndefined' : True } )
    self.widgets['maxTries'] = CCP4Widgets.CIntView(self,qualifiers={  'allowUndefined' : False } )
    self.widgets['customCodeFile'] = CCP4Widgets.CDataFileView(self,qualifiers={  'allowUndefined' : True , 'jobCombo' : False, 'browseDb': False } )
    
    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Name for server group',self))
    line.addWidget(self.widgets['name'])
    line.addWidget(QtWidgets.QLabel('using mechanism',self))
    line.addWidget(self.widgets['mechanism'])
    self.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Validate by',self))
    for item in ['password','key_filename','pass_key_filename']:
        line.addWidget(self.widgets['validate'].getButton(item))
    line.addStretch(2)
    self.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    self.widgets['serverList'].setLineWidth(1)
    self.widgets['serverList'].setFrameStyle(QtWidgets.QFrame.Box|QtWidgets.QFrame.Plain)
    line.addWidget(self.widgets['serverList'])
    self.layout().addLayout(line)
    self.widgets['serverList'].setListVisible(True)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(self.widgets['userExtensible'])
    line.addWidget(QtWidgets.QLabel('Allow users to add alternative servers at run time',self))
    line.addStretch(2)
    self.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(CCP4Widgets.CItalicLabel('CCP4 distro directory on these servers',self))
    self.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(self.widgets['ccp4Dir'])
    self.layout().addLayout(line)
    
    frame = self.makeFrame('tempDir')
    line = QtWidgets.QHBoxLayout()
    line.addWidget(CCP4Widgets.CItalicLabel('Temporary directory (eg /tmp/$USER)',self))
    frame.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(self.widgets['tempDir'])
    frame.layout().addLayout(line)

    frame = self.makeFrame('sge_root')
    line = QtWidgets.QHBoxLayout()
    line.addWidget(CCP4Widgets.CItalicLabel('SGE_ROOT (Sun Grid Engine dir)',self))
    frame.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(self.widgets['sge_root'])
    frame.layout().addLayout(line)
   
    frame = self.makeFrame('key_filename')
    frame.layout().addWidget(CCP4Widgets.CItalicLabel('Local SSH key filename (eg $HOME/.ssh/id_rsa)',self))
    frame.layout().addWidget(self.widgets['keyFilename'])

    line = self.makeFrame('timeout',vertical=False)
    line.layout().addWidget(CCP4Widgets.CItalicLabel('Timeout (seconds)',self))
    line.layout().addWidget(self.widgets['timeout'])
    line.layout().addWidget(CCP4Widgets.CItalicLabel('Maximum tries',self))
    line.layout().addWidget(self.widgets['maxTries'])
    self.layout().addWidget(line)

    frame = self.makeFrame('customCodeFile')
    line = QtWidgets.QHBoxLayout()
    line.addWidget(CCP4Widgets.CItalicLabel('Custom server code file',self))
    frame.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(self.widgets['customCodeFile'])
    frame.layout().addLayout(line)
    
    line = QtWidgets.QHBoxLayout()
    self.testButton = QtWidgets.QToolButton(self)
    self.testButton.setText('Test server group')
    self.deleteButton = QtWidgets.QToolButton(self)
    self.deleteButton.setText('Delete server group')
    line.addStretch(5)
    #I am disabling this button for now as I cannot see how this can work - it does not work in 7.1 either, so must have been broken for some time.
    #line.addWidget(self.testButton )
    line.addWidget(self.deleteButton )
    self.layout().addLayout(line)

    self.setModel(model=model)
    self.updateViewFromModel()
    self.widgets['serverList'].updateViewFromModel()
    self.widgets['serverList'].handleRowChange(0,force=True)
    self.widgets['mechanism'].updateViewFromModel()
    self.changeMechanism(0)

  def makeFrame(self,name,vertical=True):
    self.frames[name] = QtWidgets.QFrame()
    if vertical:
      self.frames[name].setLayout(QtWidgets.QVBoxLayout())
    else:
      self.frames[name].setLayout(QtWidgets.QHBoxLayout())
    self.frames[name].layout().setSpacing(1)
    self.frames[name].layout().setContentsMargins(1,1,1,1)
    self.layout().addWidget(self.frames[name])
    return self.frames[name]
  

  def updateViewFromModel(self):
    # Exclude the serverList which is behaving oddly because
    # dataChanged signal from serverList object is relayed by
    # CServerGroup and causes this method to be invoked
    for item in ['name','userExtensible','ccp4Dir','tempDir','keyFilename','customCodeFile','timeout','maxTries']:
        self.widgets[item].updateViewFromModel()
    self.widgets['validate'].setValue(str(self.model.validate))
    self.changeValidate(self.model.validate.qualifiers('enumerators').index(str(self.model.validate)))

  def updateModelFromView(self):
    for item in ['name','userExtensible','ccp4Dir','tempDir','keyFilename','serverList','customCodeFile','timeout','maxTries']:
        self.widgets[item].updateModelFromView()
    self.model.validate.set(self.widgets['validate'].getValue())
    #print 'CServerGroupView.updateModelFromView',self.model.validate

  @QtCore.Slot(int)
  def changeMechanism(self,indx):
    mechanism = self.widgets['mechanism'].getWidgetText()
    #print 'CServerGroupView.changeMechanism',indx,mechanism
    if mechanism in ['custom']:
      self.frames['customCodeFile'].show()
    else:
      self.frames['customCodeFile'].hide()
    if mechanism in ['ssh','ssh_shared','qsub_remote','qsub_shared','slurm_remote']:
      self.frames['timeout'].show()
    else:
      self.frames['timeout'].hide()
    if mechanism in ['ssh','qsub_remote','slurm_remote','custom']:
      self.frames['tempDir'].show()
    else:
      self.frames['tempDir'].hide()
    if mechanism in ['qsub_local','qsub_remote','qsub_shared']:
      self.frames['sge_root'].show()
    else:
      self.frames['sge_root'].hide()
    
  @QtCore.Slot(int)
  def changeValidate(self,id):
    #print 'changeValidate',id
    self.model.validate.set(self.widgets['validate'].getValue())
    if id in [0]:
      self.frames['key_filename'].hide()
    else:
      self.frames['key_filename'].show()
    

class CServerSetupWindow(QtWidgets.QDialog):
  ''' This is the glue to get display the CTaskServerSetup as an independent window '''

  testMessageSignal = QtCore.Signal(str)

  def __init__(self,parent):
    QtWidgets.QDialog.__init__(self,parent)
    self.setWindowTitle('Configure servers for remote running')
    #self.setMaximumWidth(650)
    #self.setMinimumHeight(500)
    self.setLayout(QtWidgets.QVBoxLayout())
    MARGIN = 2
    self.layout().setSpacing(MARGIN)
    self.layout().setContentsMargins(MARGIN,MARGIN,MARGIN,MARGIN)
    self.widgets= {}

    from ..core.CCP4ServerSetup import SERVERSETUP
    self.container = SERVERSETUP()
    self.container.load()
    line = QtWidgets.QHBoxLayout()
    self.sourceLabel = QtWidgets.QLabel('Server setup loaded from ',self)
    line.addWidget(self.sourceLabel)
    self.layout().addLayout(line)
    self.altSourceBut = QtWidgets.QPushButton('Load alternate setup',self)
    line.addWidget(self.altSourceBut)
    self.altSourceBut.clicked.connect(self.handleLoadAlt)                                                  

    self.scrollArea= QtWidgets.QScrollArea(self)
    self.scrollArea.setMaximumWidth(570)
    self.scrollArea.setMinimumWidth(570)
    self.scrollArea.setMinimumHeight(500)
    self.layout().addWidget(self.scrollArea)
    self.contentFrame = QtWidgets.QFrame(self)
    self.scrollArea.setWidget(self.contentFrame)
    self.scrollArea.setWidgetResizable(1)
    self.scrollArea.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
    self.contentLayout = QtWidgets.QVBoxLayout()
    self.contentLayout.setSpacing(MARGIN)
    self.contentLayout.setContentsMargins(MARGIN,MARGIN,MARGIN,MARGIN)
   
    self.contentFrame.setLayout(self.contentLayout)

    self.makeAllWidgets()

    line = QtWidgets.QHBoxLayout()
    self.layout().addLayout(line)
    addButton = QtWidgets.QToolButton(self)
    addButton.setText('Add another server group')
    line.addStretch(4)
    line.addWidget(addButton)
    addButton.clicked.connect(self.handleAdd)

    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel("<b>N.B. You should restart CCP4i2 to make the 'Run on server' button appear<br/>in project toolbars if you have enabled<br/>remote running.</b>"))
    self.layout().addLayout(line)

    self.buttons = QtWidgets.QDialogButtonBox(self)
    but = self.buttons.addButton(QtWidgets.QDialogButtonBox.Help)
    but.clicked.connect(self.help)
    if self.container.writeAccess('installation'):
      but = self.buttons.addButton('Save for CCP4 installation',QtWidgets.QDialogButtonBox.ApplyRole)
      but.clicked.connect(functools.partial(self.doApply,'installation'))
    but = self.buttons.addButton('Save for me only',QtWidgets.QDialogButtonBox.ApplyRole)
    but.clicked.connect(functools.partial(self.doApply,'user'))
    but = self.buttons.addButton(QtWidgets.QDialogButtonBox.Close)
    but.clicked.connect(self.close)
    line = QtWidgets.QHBoxLayout()
    line.addStretch(0.2)
    line.addWidget(self.buttons)
    line.addStretch(0.2)
    self.layout().addLayout(line)

  def makeAllWidgets(self):
    n = 1
    while self.container.get('SERVERGROUP'+str(n)) is not None:
      self.makeWidget(n)
      n = n + 1

  def clear(self):
    for key,w in list(self.widgets.items()):
      w.hide()
      w.deleteLater()
    self.widgets = {}

  def makeWidget(self,n):
    self.widgets['SERVERGROUP'+str(n)] = CServerGroupView(self,self.container.get('SERVERGROUP'+str(n)))
    self.widgets['SERVERGROUP'+str(n)].setMaximumWidth(550)
    self.widgets['SERVERGROUP'+str(n)].setMinimumWidth(550)
    self.widgets['SERVERGROUP'+str(n)].deleteButton.clicked.connect(functools.partial(self.handleDelete,n))
    self.widgets['SERVERGROUP'+str(n)].testButton.clicked.connect(functools.partial(self.handleTest,n))

    self.contentLayout.addWidget(self.widgets['SERVERGROUP'+str(n)])

  def show(self):
    self.setAltSourceButton()
    QtWidgets.QDialog.show(self)
    
  def setAltSourceButton(self):
    if self.container.source is not None:
      self.sourceLabel.setText('Server setup loaded from '+str(self.container.source))
      altSource =  ['user','installation'][1 -  ['user','installation'].index( self.container.source )]
      if self.container.preferencesFile(altSource)[0].exists():
        self.altSourceBut.setText('Load '+altSource+' setup')
        self.altSourceBut.setEnabled(True)
        return
    self.altSourceBut.setText('Load alternate setup')
    self.altSourceBut.setEnabled(False)
    
  @QtCore.Slot()
  def handleLoadAlt(self):
    altSource =  ['user','installation'][1 -  ['user','installation'].index( self.container.source )]
    #print 'handleAltSource',altSource
    self.clear()
    self.container.clear()
    self.container.load(altSource)
    self.makeAllWidgets()
    self.setAltSourceButton()

  def getNServerGroups(self):
    n = 0
    while self.container.get('SERVERGROUP'+str(n+1)) is not None:
      n = n + 1
    return n

  @QtCore.Slot()
  def handleAdd(self):
    n = self.getNServerGroups() + 1
    self.container.setContents( { 'SERVERGROUP'+str(n) : { 'class' :CCP4Annotation.CServerGroup }} )  
    self.container.get('SERVERGROUP'+str(n)).setObjectName('SERVERGROUP'+str(n))
    self.makeWidget(n)
                                
    
  @QtCore.Slot(int)
  def handleDelete(self,n):
    #print 'CServerSetupWindow.handleDelete',n,self.getNServerGroups()
    nTot = self.getNServerGroups()
    if nTot==1: return
    self.widgets['SERVERGROUP'+str(n)].hide()
    self.widgets['SERVERGROUP'+str(n)].deleteLater()
    del self.widgets['SERVERGROUP'+str(n)]
    self.container.get('SERVERGROUP'+str(n)).blockSignals(True)
    self.container.deleteObject('SERVERGROUP'+str(n))

    if n < nTot:
      for ii in range(n,nTot):
        self.widgets['SERVERGROUP'+str(ii+1)].setObjectName('SERVERGROUP'+str(ii))
        self.widgets['SERVERGROUP'+str(ii)] = self.widgets.pop('SERVERGROUP'+str(ii+1))
        self.container.renameObject('SERVERGROUP'+str(ii+1),'SERVERGROUP'+str(ii))
        
    

  @QtCore.Slot(int)
  def handleTest(self,indx):
    if hasattr(self,'testDialog'):
      self.testDialog.hide()
      self.testDialog.deleteLater()
    self.testDialog = QtWidgets.QDialog(self)
    self.testDialog.setWindowTitle('Test remote server')
    self.testDialog.setLayout(QtWidgets.QVBoxLayout())
    self.testDialog.layout().addWidget(QtWidgets.QLabel('Test existance of CCP4 directory and temporary directory on each machine',self))
    line = QtWidgets.QHBoxLayout()
    self.testDialog.layout().addLayout(line)
    self.testDialogUsername = QtWidgets.QLineEdit(self.testDialog)
    line.addWidget(QtWidgets.QLabel('Test with username',self.testDialog))
    line.addWidget(self.testDialogUsername)
    line = QtWidgets.QHBoxLayout()
    self.testDialog.layout().addLayout(line)
    self.testDialogPassword = QtWidgets.QLineEdit(self.testDialog)
    self.testDialogPassword.setEchoMode(QtWidgets.QLineEdit.Password)
    line.addWidget(QtWidgets.QLabel('Test with password',self.testDialog))
    line.addWidget(self.testDialogPassword)
    buttons = QtWidgets.QDialogButtonBox(self.testDialog)
    but = buttons.addButton(QtWidgets.QDialogButtonBox.Cancel)
    but.clicked.connect(self.testDialog.close)
    but = buttons.addButton(QtWidgets.QDialogButtonBox.Apply)
    but.clicked.connect(functools.partial(self.showTests,indx))
    line = QtWidgets.QHBoxLayout()
    line.addStretch(0.2)
    line.addWidget(buttons)
    line.addStretch(0.2)
    self.testDialog.layout().addLayout(line)
    self.testDialog.show()
    self.testDialog.raise_()
                                       
    
  @QtCore.Slot(int)
  def showTests(self,indx):
    self.testDialog.close()
    if hasattr(self,'testMessage'):
      self.testMessage.hide()
      self.testMessage.deleteLater()
    else:
      self.testMessageSignal.connect(self.updateTestReport)
    message = ''
    self.testMessage = QtWidgets.QDialog(self)
    self.testMessage.setWindowTitle('Testing remote server')
    self.testMessage.setLayout(QtWidgets.QVBoxLayout())
    self.testMessage.layout().addWidget(QtWidgets.QLabel('Test existance of CCP4 directory and temp directory',self.testMessage))
    self.testReport = QtWidgets.QTextEdit(self.testDialog)
    self.testMessage.layout().addWidget(self.testReport)
    buttons = QtWidgets.QDialogButtonBox(self.testDialog)
    but = buttons.addButton(QtWidgets.QDialogButtonBox.Cancel)
    self.testMessage.layout().addWidget(buttons)
    but.clicked.connect(functools.partial(self.cancelTests,indx,True))
    self.testMessage.show()
    self.testMessage.raise_()

    if not hasattr(self,'testThreads'): self.testThreads = {}
    self.testThreads[indx] = UtilityThread.UtilityThread(functools.partial(self.runTests,indx))
    self.testThreads[indx].finished.connect(functools.partial(self.cancelTests,indx))
    self.testThreads[indx].start()
    #self.runTests(indx)

  @QtCore.Slot(str)
  def updateTestReport(self,message):
      self.testReport.setReadOnly(False)
      mess=self.testReport.toPlainText() + '\n' + message
      self.testReport.setPlainText(mess)
      self.testReport.setReadOnly(True)
      self.testReport.repaint()
      
  @QtCore.Slot(int,bool)
  def cancelTests(self,indx,force=False):
    if getattr(self,'testThreads',None) is not None and len(self.testThreads)>indx: 
      if force: self.testThreads[indx].exit()
      del self.testThreads[indx]
    if force: self.testMessage.close()

  def runTests(self,indx):
    sG = self.widgets['SERVERGROUP'+str(indx)]
    remoteFileList=[sG.widgets['ccp4Dir'].model.__str__()]
    if sG.widgets['mechanism'].model.__str__() in ['ssh','qsub_remote','slurm_remote']:
      remoteFileList.append(sG.widgets['tempDir'].model.__str__())
    if sG.widgets['timeout'].model.isSet():
      timeout = float(sG.widgets['timeout'].model)
    else:
      timeout = None
    maxTries = int(sG.widgets['maxTries'].model)
    for server in sG.widgets['serverList'].model:
      #print 'testing server',server
      self.testMessageSignal.emit('Testing '+str(server)+'\n')
      try:
        from ..qtcore.CCP4JobController import JOBCONTROLLER
        ret = JOBCONTROLLER().testRemoteFiles(machine=str(server),username=str(self.testDialogUsername.text()),
           password = str(self.testDialogPassword.text()), remoteFileList=remoteFileList,timeout=timeout,maxTries=maxTries)
        message = ''
        labelList = ['CCP4 distro','Tmp dir']
        for ii in range(len(ret)):
          message =  message + labelList[ii] + ': ' + str(ret[ii]) +'\n'
      except CException as e:
        message =  '  Failed connection\n   ' + e.description()[0] + '\n   ' + str(e._reports[0]['details']) +'\n'
      except Exception as e:
        message =  '  Failed connection \n'+str(e)
      
      self.testMessageSignal.emit(message)


  @QtCore.Slot(str)
  def doApply(self,source):
    invalidList = []
    for key in list(self.widgets.keys()):
      self.widgets[key].updateModelFromView()
      self.widgets[key].validate()
      if self.widgets[key].isValid is not None and not self.widgets[key].isValid:
        invalidList.append(key)

    if len(invalidList)>0:
      QtWidgets.QMessageBox.warning(self,'Server setup','Invalid data - not saved')
    else:
      self.container.save(source)
      self.setAltSourceButton()
      from ..qtcore.CCP4JobController import JOBCONTROLLER
      JOBCONTROLLER().resetServersEnabled()

  @QtCore.Slot()
  def close(self):
    QtWidgets.QDialog.close(self)

  @QtCore.Slot()
  def help(self):
    from .CCP4WebBrowser import WEBBROWSER
    WEBBROWSER().loadWebPage(helpFileName='servers.html',newTab=True)


class CPasswordEntry(QtWidgets.QDialog):

  passwordEntered = QtCore.Signal(str)

  def __init__(self,parent,label='Please reenter password',closeMessage='Ignore remote jobs'):
    QtWidgets.QDialog.__init__(self,parent)
    self.setModal(False)
    self.setLayout(QtWidgets.QVBoxLayout())
    self.layout().addWidget(QtWidgets.QLabel(label))
    self.passwordEntry = QtWidgets.QLineEdit(self)
    self.passwordEntry.setEchoMode(QtWidgets.QLineEdit.Password)
    self.layout().addWidget(self.passwordEntry)
    butLayout = QtWidgets.QHBoxLayout()
    for label,connect in [['Apply',self.sendPass],[closeMessage,self.close]]:
      but = QtWidgets.QPushButton(label,self)
      butLayout.addWidget(but)
      but.clicked.connect(connect)
    self.layout().addLayout(butLayout)

  def sendPass(self):
    passw = str(self.passwordEntry.text())
    #print 'sendPass',passw
    self.passwordEntered.emit(passw)
    self.close()


class CListProcesses(QtWidgets.QDialog):
  COLUMNHEADERS = [ 'Machine/Project','Job','Taskname','Running process','Pid','Running time','Executable' ]
  def __init__(self,parent):
    QtWidgets.QDialog.__init__(self,parent)
    self.setWindowTitle('CCP4i2 Running jobs and processes')
    self.setLayout(QtWidgets.QVBoxLayout())
    self.treeWidget = QtWidgets.QTreeWidget(self)
    self.treeWidget.setColumnCount(len(self.COLUMNHEADERS))
    self.treeWidget.setMinimumWidth(400)    
    self.treeWidget.setColumnWidth(1,50)
    self.treeWidget.setColumnWidth(4,50)
    self.treeWidget.setColumnWidth(6,250)
    self.layout().addWidget(self.treeWidget)
    self.treeWidget.itemSelectionChanged.connect(self.handleSelectionChanged)
    
    butLayout = QtWidgets.QHBoxLayout()
    butLayout.addStretch(3)
    for label,connect in [['Update',self.load],
                          ['Help',self.help],
                          ['Close',self.close],
                          ['Kill job',self.killJob],
                          ['Mark as finished',functools.partial(self.markFinished,False)],
                          ['Mark as failed',self.markFinished]]:
      but = QtWidgets.QPushButton(label,self)
      but.clicked.connect(connect)
      butLayout.addWidget(but)
    butLayout.addStretch(3)
    self.killButton = butLayout.itemAt(4).widget()
    self.markFinishedButton = butLayout.itemAt(5).widget()
    self.markFailedButton = butLayout.itemAt(6).widget()
    self.killButton.setEnabled(False)
    self.markFinishedButton.setEnabled(False)
    self.markFailedButton.setEnabled(False)
    self.layout().addLayout(butLayout)
    self.projectNameCache = {}
    from ..qtcore.CCP4JobController import JOBCONTROLLER
    JOBCONTROLLER().remoteProcessesList.connect(self.handleRemoteProcessesList)

  def projectName(self,projectId):
    if projectId not in self.projectNameCache:
      from ..core.CCP4ProjectsManager import PROJECTSMANAGER
      self.projectNameCache[projectId] = PROJECTSMANAGER().db().getProjectInfo(projectId,'projectname')
    return self.projectNameCache[projectId]

  @QtCore.Slot()
  def load(self):
    self.treeWidget.clear()
    #self.tableWidget.setRowCount(0)
    self.treeWidget.setHeaderLabels(self.COLUMNHEADERS)
    # load local info
    now = time.time()
    from ..qtcore.CCP4JobController import JOBCONTROLLER
    procDict = JOBCONTROLLER().listLocalProcesses(containsList=['ccp4'])
    #print 'CListProcesses.load',procDict
    from ..core.CCP4ProjectsManager import PROJECTSMANAGER
    runningJobs = PROJECTSMANAGER().db().getRunningJobs()
    self.drawTree(CCP4Utils.getHostName(),now,procDict,runningJobs)
    
    JOBCONTROLLER().listRemoteProcesses()


  @QtCore.Slot(tuple)
  def handleRemoteProcessesList(self,args):
    machine,jobDict,procDict,atTime = args
    #print 'handleRemoteProcessesList',machine,procDict,procDict,atTime
    #print 'handleRemoteProcessesList processes keys',procDict.keys()
    from ..core.CCP4ProjectsManager import PROJECTSMANAGER
    runningJobs = PROJECTSMANAGER().db().getRunningJobs(remote=True)
    #print 'handleRemoteProcessesList runningJobs',runningJobs
    self.drawTree(machine,atTime,procDict,runningJobs)

  def drawTree(self,machine,now,procDict,runningJobs):
    machineItem = QtWidgets.QTreeWidgetItem([machine,'','','','','',''])
    machineItem.setFlags(QtCore.Qt.ItemIsEnabled)
    font = machineItem.font(0)
    font.setBold(True)
    machineItem.setFont(0,font)
    self.treeWidget.addTopLevelItem(machineItem)
    self.treeWidget.expandItem(machineItem)
    for jobId,jobNumber,taskName,projectId,processId,parentJobId in runningJobs:
      if parentJobId is None and processId is not None:
        from ..core.CCP4TaskManager import TASKMANAGER
        taskTitle = TASKMANAGER().getTitle(taskName)
        #print 'CListProcesses.load',jobId,jobNumber,processId
        projectName = self.projectName(projectId)
        if processId in procDict:
          procDict[processId]['used'] = True
          taskItem = self.makeItem(now,projectName,jobNumber,taskTitle,procDict[processId],italic=True,jobId=jobId)
          taskItem.setFlags(QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEnabled)
          machineItem.addChild(taskItem)
          for child in procDict[processId]['children']:
            if child in procDict:
              procDict[child]['used'] = True
              item = self.makeItem(now,'','','',procDict[child])
              item.setFlags(QtCore.Qt.ItemIsEnabled)
              taskItem.addChild(item)
          self.treeWidget.expandItem(taskItem)
        else:
          machineItem.addChild(self.makeItem(now,projectName,jobNumber,taskTitle,None,bold=True,jobId=jobId))

    # List any unaccounted ccp4 processes
    mypid = os.getpid()
    for key,pDict in list(procDict.items()):
      if 'used' not in pDict and not pDict['pid'] == mypid:
        item = self.makeItem(now,'Not i2 process','','',pDict)
        item.setFlags(QtCore.Qt.ItemIsEnabled)
        machineItem.addChild(item)

  def makeItem(self,now,projectName,jobNumber,taskName,procDict,bold=False,italic=False,jobId=None):
     if procDict is None:
       treeItem = QtWidgets.QTreeWidgetItem([projectName,jobNumber,taskName,'No process found','','',''])
     else:
       treeItem = QtWidgets.QTreeWidgetItem([projectName,jobNumber,taskName,str(procDict['name']),
                  str(procDict['pid']),'%.1f' % (now-procDict['create_time']),str(procDict['exe'])])
     if bold or italic:
       font = treeItem.font(0)
       if bold: font.setBold(True)
       if italic: font.setItalic(True)
       for ii in range(7): treeItem.setFont(ii,font)
     if jobId is not None:
       treeItem.setData(0,QtCore.Qt.UserRole,jobId)
         
     return treeItem

  @QtCore.Slot()
  def help(self):
    from .CCP4WebBrowser import WEBBROWSER
    WEBBROWSER().loadWebPage(helpFileName='tutorial',target='running')

  @QtCore.Slot()
  def handleSelectionChanged(self):
    selItems = self.treeWidget.selectedItems()
    if len(selItems) == 0:
      self.markFinishedButton.setEnabled(False)
      self.markFailedButton.setEnabled(False)
      self.killButton.setEnabled(False)
    else:
      ifRunning = selItems[0].data(3,QtCore.Qt.DisplayRole).__str__() != 'No process found'
      self.markFinishedButton.setEnabled((not ifRunning))
      self.markFailedButton.setEnabled((not ifRunning))
      self.killButton.setEnabled(ifRunning)
    
  @QtCore.Slot()
  def killJob(self):
    from ..dbapi import CCP4DbApi
    jobId = self.treeWidget.selectedItems()[0].data(0,QtCore.Qt.UserRole).__str__()
    print('CListProcesses.killJob',jobId)
    from ..qtcore.CCP4JobController import JOBCONTROLLER
    err = JOBCONTROLLER().killJobProcess(jobId=jobId)
    from ..core.CCP4ProjectsManager import PROJECTSMANAGER
    PROJECTSMANAGER().db().updateJobStatus(jobId,CCP4DbApi.JOB_STATUS_FAILED)

    self.load()

  @QtCore.Slot(bool)
  def markFinished(self,failed=True):
    from ..dbapi import CCP4DbApi
    jobId = self.treeWidget.selectedItems()[0].data(0,QtCore.Qt.UserRole).__str__()
    from ..core.CCP4ProjectsManager import PROJECTSMANAGER
    if failed:
      PROJECTSMANAGER().db().updateJobStatus(jobId,CCP4DbApi.JOB_STATUS_FAILED)
    else:
      PROJECTSMANAGER().db().updateJobStatus(jobId,CCP4DbApi.JOB_STATUS_FINISHED)
    self.load()
