"""
Martin Noble
"""

import tempfile

from PySide2 import QtCore, QtWidgets

from ....qtgui.CCP4TaskWidget import CTaskWidget
from .ProvideAlignment import importAlignment


class CAlignmentSelectorModel(QtCore.QAbstractItemModel):
    
  def __init__(self,parent):
    self.parent_=parent
    QtCore.QAbstractItemModel.__init__(self)
  
  def parent(self,index):
    if not index.isValid():
      return QtCore.QModelIndex()
    else:
      return QtCore.QModelIndex()

  def child(self,row):
    if self.parent_.container.controlParameters.PASTEORREAD == 'HHPREDIN':
      return self.parent_.container.inputData.HHPREDIN.fileContent.alignmentList[row]
    elif self.parent_.container.controlParameters.PASTEORREAD == 'BLASTIN':
      return self.parent_.container.inputData.BLASTIN.fileContent.alignmentList[row]
    else:
      return None
  def rowCount(self,mI):
    if self.parent_.container.controlParameters.PASTEORREAD == 'HHPREDIN':
      return len(self.parent_.container.inputData.HHPREDIN.fileContent.alignmentList)
    elif self.parent_.container.controlParameters.PASTEORREAD == 'BLASTIN':
      return len(self.parent_.container.inputData.BLASTIN.fileContent.alignmentList)
    else:
      return 0
  def columnCount(self,mI):
    return 1

  def data(self,index,role):
    if not index.isValid():
#FIXME PYQT - or maybe None? This used to return QVariant.
      return None
    item = index.internalPointer()
    return item.data(role)
    
  def index(self,row,column,parent):
    try:
      if self.parent_.container.controlParameters.PASTEORREAD == 'HHPREDIN':
        return self.createIndex(row,column,self.parent_.container.inputData.HHPREDIN.fileContent.alignmentList[row])
      elif self.parent_.container.controlParameters.PASTEORREAD == 'BLASTIN':
        return self.createIndex(row,column,self.parent_.container.inputData.BLASTIN.fileContent.alignmentList[row])

    except:
      return self.createIndex(0,0,parent)

class CAlignmentSelectorView(QtWidgets.QListView):

    currentChangedSignal = QtCore.Signal(tuple)

    def __init__(self,parent,model):
      QtWidgets.QListView.__init__(self,parent)
      self.setModel(model)
      self.setSelectionMode(self.SingleSelection)
      
    def currentChanged(self,current,previous):
      print('CAlignmentSelectorWidget.currentChanged',current.row(),previous.row())
      self.currentChangedSignal.emit((current.row(),previous.row()))

class CTaskProvideAlignment(CTaskWidget):
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'ProvideAlignment'
    TASKVERSION = 0.0
    TASKMODULE='data_entry'
    TASKTITLE="Import an alignment"
    WHATNEXT = []
    DESCRIPTION = '''Enter an alignment from a file or by cut and paste'''
    
    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)
    
    def drawContents(self):
        self.selectorModel = CAlignmentSelectorModel(self)
        self.selector = CAlignmentSelectorView(self,self.selectorModel)
        self.selector.currentChangedSignal.connect(self.handleCurrentChanged)

        self.setProgramHelpFile('ProvideAlignment')
        
        folder = self.openFolder(folderFunction='inputData',title='Optional objects from which to start definition',followFrom=False)
        
        self.createLine( ['label','Paste or read alignment, or extract from HHPred or Blast search','widget','PASTEORREAD'])
        self.container.controlParameters.PASTEORREAD.dataChanged.connect(self.handleModeChanged)

        self.createLine( [ 'widget', '-browseDb', True, 'ALIGNIN', 'tip', 'Alignment object or file' ] ,toggle=['PASTEORREAD','open',['ALIGNIN']])
        self.createLine( [ 'widget', '-browseDb', True, 'HHPREDIN', 'tip', 'HHPred results' ] ,toggle=['PASTEORREAD','open',['HHPREDIN']])
        self.createLine( [ 'widget', '-browseDb', True, 'BLASTIN', 'tip', 'Blast results' ] ,toggle=['PASTEORREAD','open',['BLASTIN']])
        self.createLine( [ 'widget', '-guiMode','multiLine','SEQUENCETEXT' ] ,toggle=['PASTEORREAD','open',['PASTE']])

        line = self.createLine(['advice','Choose one alignment from the search file'],toggle=['PASTEORREAD','open',['HHPREDIN','BLASTIN']])
        line = self.createLine(toggle=['PASTEORREAD','open',['HHPREDIN','BLASTIN']])
        line.layout().addWidget(self.selector)
        self.container.inputData.HHPREDIN.dataChanged.connect(self.handleFileSelected)
        self.container.inputData.BLASTIN.dataChanged.connect(self.handleFileSelected)

        self.createLine( ['advice','Annotation for the alignment'])
        line = self.createLine( [ 'widget', 'ANNOTATION' ] )
        self.getWidget('ANNOTATION').setMinimumWidth(600)

        if self.container.controlParameters.PASTEORREAD == 'HHPREDIN':
          if self.container.inputData.HHPREDIN.isSet(): self.loadSelector()
        elif self.container.controlParameters.PASTEORREAD == 'BLASTIN':
          if self.container.inputData.BLASTIN.isSet(): self.loadSelector()

    def handleCurrentChanged(self,row_previous):
      row,previous = row_previous
      self.container.inputData.ALI_INDEX.set(row)
      if self.container.controlParameters.PASTEORREAD == 'HHPREDIN':
        anno = str(self.container.inputData.HHPREDIN.fileContent.alignmentList[row].annotation)
        if anno.count('}'): anno = anno.split('}')[0]+'}'
      else:
        anno = str(self.container.inputData.BLASTIN.fileContent.alignmentList[row].hitId)
      self.container.controlParameters.ANNOTATION.set(anno)
      

    @QtCore.Slot()
    def handleFileSelected(self):
      print('handleFileSelected',self.container.controlParameters.PASTEORREAD == 'BLASTIN',self.container.inputData.BLASTIN)
      self.container.inputData.ALI_INDEX.set(0)
      if self.container.controlParameters.PASTEORREAD == 'HHPREDIN':
        self.container.inputData.HHPREDIN.loadFile()
      elif self.container.controlParameters.PASTEORREAD == 'BLASTIN':
        self.container.inputData.BLASTIN.loadFile()
      self.selector.reset()

      
    @QtCore.Slot()
    def handleModeChanged(self):
      mode = self.container.controlParameters.PASTEORREAD.__str__()
      self.container.controlParameters.ANNOTATION.unSet()
      for item in ['ALIGNIN','HHPREDIN','BLASTIN']:
        self.container.inputData.get(item).unSet()
        self.container.inputData.get(item).setQualifiers({'allowUndefined':True})
      if mode in ['ALIGNIN','HHPREDIN','BLASTIN']:
        self.container.inputData.get(mode).setQualifiers({'allowUndefined':False})
      self.container.inputData.ALI_INDEX.set(0)
      self.selector.update()
      self.validate()

    def loadSelector(self):
      #print 'loadSelector',self.container.controlParameters.PASTEORREAD
      if self.container.controlParameters.PASTEORREAD.__str__() == 'HHPREDIN':
        if self.container.inputData.HHPREDIN.isSet():
          self.selector.update(QtCore.QModelIndex())
      elif self.container.controlParameters.PASTEORREAD.__str__() == 'BLASTIN':
        print('loadSelector',self.container.inputData.BLASTIN)
        if self.container.inputData.BLASTIN.isSet():
          self.selector.update(QtCore.QModelIndex())
      else:
        return
    
      if not self.container.inputData.ALI_INDEX.isSet(): self.container.inputData.ALI_INDEX.set(0)
      self.selector.selectionModel().select(self.selectorModel.index(int(self.container.inputData.ALI_INDEX),0,QtCore.QModelIndex()),
                                                QtCore.QItemSelectionModel.ClearAndSelect)

    def isValid(self):
        #Here override logic of whether there is validinput to task
        invalidElements = super(CTaskProvideAlignment,self).isValid()
        if self.container.controlParameters.PASTEORREAD.__str__() ==  'PASTE':
            if self.container.inputData.ALIGNIN in invalidElements:
                invalidElements.remove(self.container.inputData.ALIGNIN)
            # Create a temporary data object to test for validity (==recognisability) of pasted text
            tempFile = tempfile.NamedTemporaryFile(suffix='.txt',delete=False)
            tempFile.write(bytes(self.container.controlParameters.SEQUENCETEXT.__str__(),"utf-8"))
            tempFile.close()
            alignment, format, commentary = importAlignment(tempFile.name)
            print('isValid',alignment, format, commentary)
            if alignment is not None:
                self.container.inputData.ALIGNIN.setFullPath(tempFile.name)
                if self.container.controlParameters.SEQUENCETEXT in invalidElements:
                    invalidElements.remove(self.container.controlParameters.SEQUENCETEXT)
            else:
                if not self.container.controlParameters.SEQUENCETEXT in invalidElements:
                    invalidElements.append(self.container.controlParameters.SEQUENCETEXT)
        else:
            if self.container.controlParameters.SEQUENCETEXT in invalidElements:
                invalidElements.remove(self.container.controlParameters.SEQUENCETEXT)
        return invalidElements
