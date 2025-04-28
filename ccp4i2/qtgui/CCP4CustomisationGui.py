"""
Copyright (C) 2013 STFC
Liz Potterton July 2013 - create and manage customisations
"""

import functools

from PySide2 import QtCore, QtWidgets

from ..core.CCP4ErrorHandling import CException, Severity
from ..core.CCP4Modules import DUMMYMAINWINDOW, WEBBROWSER
from ..core.CCP4WarningMessage import warningMessage


class CCustomisationGui(QtWidgets.QDialog):

  insts = []
  
  def __init__(self,parent=None,mode=None,title=None,ifEdit=True,ifClone=True):
    if parent is None: parent = DUMMYMAINWINDOW()
    QtWidgets.QDialog.__init__(self,parent)
    CCustomisationGui.insts.append(self)
    self.mode = mode
    self.setWindowTitle(title)
    self.createWidget = None

    layout = QtWidgets.QGridLayout()
    self.setLayout(layout)
  
    self.customListView = CCustomListWidget(self)
    self.customListView.itemDoubleClicked.connect(self.handleDoubleClick)
    self.customListView.populate()
    self.layout().addWidget(self.customListView,1,0)

    butDef = [          ['New',self.handleNew,True,False],
                        ['Edit',self.handleEdit,True,False],
                        ['Clone',self.handleClone,True,False],
                        ['Export',self.handleExport,True,False],
                        ['Import',self.handleImport,True,False],
                        ['Delete',self.handleDelete,True,False],
                        ['Help',self.handleHelp,True,False]]
    if not ifClone: del butDef[2]
    if not ifEdit: del butDef[1]

    buttonLayout = QtWidgets.QVBoxLayout()
    for label,slot,enabled,default in butDef:
                       
      button =QtWidgets.QPushButton(self,text=label)
      button.setEnabled(enabled)
      if default: button.setDefault(True)
      buttonLayout.addWidget(button)
      button.clicked.connect(slot)
      
    self.layout().addLayout(buttonLayout,1,1)

    
  def done(self,r):
    try:
      CCustomisationGui.insts.remove(self)
    except:
      print('Error in CCustomisationGui.done removing from insts list')
    print('CCustomisationGui.done',r,CCustomisationGui.insts)
    QtWidgets.QDialog.done(self,r)

  def manager(self):
    # Reimplement in sub-class
    return None

  @QtCore.Slot()
  def handleNew(self):
    pass

  @QtCore.Slot(str)
  def handleEdit(self,selected):
    print('CCustomisationGui.handleEdit',selected)
    pass

  @QtCore.Slot()
  def handleClone(self):
    selected = self.customListView.selectedItem()
    if selected is None: return
    if not hasattr(self,'cloneDialog'):
      self.cloneDialog = CCustomCloneDialog(self)
      self.cloneDialog.clone.connect(self.handleClone2)
    self.cloneDialog.setSelected(selected)
    self.cloneDialog.show()
    self.cloneDialog.raise_()

  @QtCore.Slot(tuple)
  def handleClone2(self,original_new):
    original,new = original_new
    try:
      self.manager().clone(original,new)
    except CException as e:
      warningMessage(e, 'Clone '+self.mode,'Error cloning '+original,parent=self)
      return
    self.handleEdit(new)
  
  @QtCore.Slot('QListWidgetItem')
  def handleDoubleClick(self,item):
    nameVar = item.data(QtCore.Qt.UserRole)
    name = nameVar.__str__()
    self.handleEdit(selected=name)

  @QtCore.Slot()
  def handleExport(self):
    selected = self.customListView.selectedItem()
    if selected is None: return
    from . import CCP4FileBrowser
    self.browser = CCP4FileBrowser.CFileDialog(self,
           title='Save '+self.mode+' to compressed file',
           filters= ['CCP4 '+self.mode+' (*.ccp4_'+self.mode+'.tar.gz)'],
           defaultSuffix='ccp4_'+self.mode+'.tar.gz',
           fileMode=QtWidgets.QFileDialog.AnyFile  )
    self.browser.selectFile.connect(functools.partial(self.handleExport2,selected))
    self.browser.show()

  @QtCore.Slot(str,str)
  def handleExport2(self,selected,fileName):
    self.browser.hide()
    self.browser.deleteLater()
    del self.browser
    err = self.manager().export(selected,fileName)
    if err.maxSeverity()>Severity.WARNING:
      warningMessage(err, 'Export '+self.mode,'Error creating compressed file',parent=self)
      
  @QtCore.Slot()
  def handleImport(self):
    from . import CCP4FileBrowser
    self.browser = CCP4FileBrowser.CFileDialog(self,
           title='Import '+self.mode+' compressed file',
           filters= ['CCP4 compressed '+self.mode+' (*.ccp4_'+self.mode+'.tar.gz)'],
           defaultSuffix='ccp4_'+self.mode+'.tar.gz',
           fileMode=QtWidgets.QFileDialog.ExistingFile  )
    self.browser.selectFile.connect(self.handleImport1)
    self.browser.show()

  @QtCore.Slot(str)
  def handleImport1(self,fileName):
    self.browser.hide()
    self.browser.deleteLater()
    del self.browser
    #print 'handleImport1',fileName
    try:
      self.manager().testImport(fileName)
    except CException as e:
      if e.count(code=110) == 0:
        warningMessage(e, 'Error opening compressed '+self.mode+' file',parent=self)
      else:
        name = e[0]['details']
        self.importQuery = CCustomImportDialog(self,name,fileName)
        self.importQuery.importSignal.connect(self.handleImport2)
        self.importQuery.show()
    else:
      try:
        self.manager().uncompress(fileName,self.manager().getDirectory())
      except CException as e:
        warningMessage(e, 'Error importing '+self.mode,'Failed to write to '+self.mode+' directory',parent=self)

  @QtCore.Slot(str,str,bool,bool)
  def handleImport2(self,name,fileName,overwrite,rename):
    #print 'handleImport2',name,overwrite,rename
    if overwrite:
      self.manager().delete(name)
    self.manager().uncompress(fileName, self.manager().getDirectory(),rename=rename)
    
  @QtCore.Slot()
  def handleDelete(self):
    selected = self.customListView.selectedItem()
    if selected is None: return
    ret = QtWidgets.QMessageBox.question(self,self.windowTitle(),'Delete '+self.mode+' '+selected,QtWidgets.QMessageBox.Ok|QtWidgets.QMessageBox.Cancel)
    if ret == QtWidgets.QMessageBox.Cancel: return
    ret = self.manager().delete(selected)
    if ret:
      QtWidgets.QMessageBox.warning(self,self.windowTitle(),'Error deleting '+self.mode+' '+selected)

  @QtCore.Slot()
  def handleHelp(self):
    WEBBROWSER().loadWebPage(helpFileName='customisation')


class CCustomListWidget(QtWidgets.QListWidget):

  def __init__(self,parent):
    QtWidgets.QListWidget.__init__(self,parent)
    self.parent().manager().listChanged.connect(self.populate)

  @QtCore.Slot()
  def populate(self):
    self.clear()
    titleList =self.parent().manager().getTitleList()
    for name,title in titleList:
      obj = QtWidgets.QListWidgetItem(title)
      obj.setData(QtCore.Qt.UserRole,name)
      self.addItem(obj)


  def selectedItem(self):
    seleList = self.selectedItems()
    #print 'contentsTree.selectedItem',seleList
    if len(seleList) == 0:
      return None
    else:
      qvar = seleList[0].data(QtCore.Qt.UserRole)
      if qvar.isNull() or (not qvar.isValid()):
        return seleList[0].text().__str__()
      else:
        return qvar.__str__()
    

class CCustomImportDialog(QtWidgets.QDialog):

  importSignal = QtCore.Signal(str,str,int,str)

  def __init__(self,parent,name=None,fileName=None):
    QtWidgets.QDialog.__init__(self,parent)
    self.setWindowTitle("Import a "+self.parent().mode)
    self.setLayout(QtWidgets.QVBoxLayout())
    self.name = name
    self.fileName = fileName
    self.layout().addWidget(QtWidgets.QLabel('There is already a '+self.parent().mode+' called '+str(name)))
    
    self.bg = QtWidgets.QButtonGroup(self)
    self.bg.setExclusive(True)
    but =  QtWidgets.QRadioButton('Overwrite existing '+self.parent().mode,self)
    self.bg.addButton(but,1)
    self.layout().addWidget(but)
    line = QtWidgets.QHBoxLayout()
    but =  QtWidgets.QRadioButton('Rename to ',self)
    line.addWidget(but)
    self.rename = QtWidgets.QLineEdit(self)
    line.addWidget(self.rename)
    self.layout().addLayout(line)
    self.bg.addButton(but,2)
    but.setChecked(True)
    
    butBox = QtWidgets.QDialogButtonBox(self)
    but = butBox.addButton('Import '+self.parent().mode,QtWidgets.QDialogButtonBox.ActionRole)
    but.setDefault(False)
    but.clicked.connect(self.apply)
    but = butBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
    but.clicked.connect(self.close)
    self.layout().addWidget(butBox)
    
                                                          
  @QtCore.Slot()
  def close(self):
    self.hide()
    self.deleteLater()
    
  @QtCore.Slot()
  def apply(self):
    if self.bg.checkedId() == 2 and len(self.rename.text().__str__())<=0 or self.rename.text().__str__() in self.parent().manager().getList(): return
    self.importSignal.emit(self.name,self.fileName,2-self.bg.checkedId(),self.rename.text().__str__())
    self.close()

class CCustomCloneDialog(QtWidgets.QDialog):
  clone = QtCore.Signal(tuple)

  def __init__(self,parent,selected=None):
    QtWidgets.QDialog.__init__(self,parent)
    self.setWindowTitle("Clone a "+self.parent().mode)
    self.setLayout(QtWidgets.QVBoxLayout())
    self.label = QtWidgets.QLabel(self)
    self.layout().addWidget(self.label)
    
    line = QtWidgets.QHBoxLayout()
    line.addWidget( QtWidgets.QLabel('Name',self))
    self.nameWidget = QtWidgets.QLineEdit(self)
    line.addWidget(self.nameWidget)
    self.layout().addLayout(line)

    butBox = QtWidgets.QDialogButtonBox(self)
    but = butBox.addButton('Clone '+self.parent().mode,QtWidgets.QDialogButtonBox.ActionRole)
    but.setDefault(False)
    but.clicked.connect(self.apply)
    but = butBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
    but.clicked.connect(self.close)
    self.layout().addWidget(butBox)

    self.setSelected(selected)

  def setSelected(self,selected):
    self.selected = selected
    self.label.setText('Clone '+str(selected)+' to new '+self.parent().mode+':')

  @QtCore.Slot()
  def apply(self):
    name = self.nameWidget.text().__str__()
    #title = self.titleWidget.text().__str__()

    if len(name)<=0 or name in self.parent().manager().getList():
      QtWidgets.QMessageBox.warning(self,'Clone '+self.parent().mode,'Please enter unique name for new '+self.parent().mode)
      return

    self.clone.emit((self.selected,name))
    self.close()
