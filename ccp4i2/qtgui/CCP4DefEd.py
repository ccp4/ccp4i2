"""
Copyright (C) 2010 University of York

Liz Potterton Oct 2010 Create widget to edit def.xml files
"""

import functools
import os
import re
import shutil
import sys

from PySide2 import QtCore, QtGui, QtWidgets

from .. import __version__
from ..core import CCP4Data
from ..core import CCP4Utils
from ..core.CCP4DataManager import DATAMANAGER
from ..core.CCP4ErrorHandling import CErrorReport, CException
from ..core.CCP4WarningMessage import warningMessage
from . import CCP4Widgets


HELPTEXT = { 'Introduction' : """<html>
<h3>Introduction</h3>
<p>This GUI is intended to help developers create a wrapper. Each CCP4i2 wrapper must have a DEF XML file which defines the input and output data and any control parameters. The main window will help create a  DEF XML file. There are also options under the <i>Tools</i. menu to create the appropriate directory structures with example template scripts.</p>
<p>The first step is to create a wrapper directory using the <i>Tools</i>-><i>Make wrapper plugin</i> option. In the new window you must enter a unique name for the wrapper and choose whether it goes in the general wrappers or a particular pipeline directory. You have a choice of basic templates and of license to prepend to the template. You can add alternative licenses to the directory <i>$CCP4I2/data/plugin_templates/license_templates</i>.</p>
<p>Once the wrapper directory is created you should use the main defEd editor to create the DEF XML file. But, beware, if you are lucky there may already be COM template file for generating the input command file for the program that you which to run. If you want to use an existing COM template file then you must ensure that you parameter names are consistent with those used in the COM file. The DEF file should be saved with the path <i>$CCP4I2/wrappers/myWrapper/script/myWrapper.def.xml</i>.</p>

<h3>Editing a DEF file</h3>
<p>The main functionality of this <b>defEd</b> editor is to create and edit DEF XML files for CCP4i2 GUI and pipelines. The editor is split vertically into three panes.  The top pane shows the current <i>contents tree</i>, the middle pane is the <i>object editor</i> with tools to edit the currently selected item in the data definition and the bottom pane has help text. A right mouse click in any widget in the gui gives a menu with the option of help.</p>
<a name="Contents tree"/><h5>Contents tree</h5>
<p>The current data contents tree is shown in the top pane along with buttons to..
<br/><b>Append container</b> - a container is added at end of the 'top level' of the tree
<br/><b>Insert container</b> - a container is inserted after current item at same level in the tree
<br/><b>Add object</b> - add a data object in either the currently selected container or after the currently selected data object
<br/><b>Delete object</b> - delete the currently selected container or data object
<br/>Click on any item in the tree to make it editable in the <i>object editor</i>.
</p>
<p>
You should start by appending a container and using the <i>Name</i> widget in the <i>object editor</i> to set a suitable
name.  The recommended def file structure (open for discussion) is to have three containers for <i>inputData</i>, <i>outputData</i> and
<i>controlParameters</i>.  The <i>controlParameters</i> could (should?) have sub-containers
for the different components of a pipeline.
Any additional parameters needed for gui function should be in a separate <i>guiControls</i> container.
</p>
<a name="Data objects"/><a name="Data names"/><h5>The object editor</h5>
<p>When you add an object it is created with type <i>CString</i>. Use the table on the left side of the <i>object editor</i> to select the appropriate data class. A right mouse click on any item in this table will give documentation in the help pane. Also edit the <i>Name</i> for the data object. On the right side of the pane the appropriate qualifiers for each data class are shown.  These enable simple customisations of the basic class so avoiding the need to write a Python subclass for trivial data.</p>
<p>The qualifiers are described in the class documentation but some general points:
<ul>
<li>Any data class inherits the qualifiers of its base class. It can have a different default qualifier value than the base class and can have additional qualifiers.</li>
<li>For now please assume that there is no validity checking of the qualifier values that you enter.</li>
<li>The <i>enumerators</i> qualifier can be either obligatory (<i>onlyEnumerators=True</i>) or advisory ( (<i>onlyEnumerators=False</i>). In either case it is good to provide a <i>menuText</i> list of strings, equivalent to the enumerators, to appear in the gui. These should inform the user of the significance of the data value, so something like 'few','normal','many' is good.</li>
<li>Please provide as many qualifiers as possible - these are used to check input data validity and will help to make the system robust. Also provide the 'advisory' enumerators and menuText where possible.</li>
<li>The qualifiers that are lists should be specified with valid Python syntax. For example a menuText might be <i>['few','normal','many']</i>.</li>
<li>Suggestions for other general classes or qualifiers are welcome to lizp@ysbl.york.ac.uk</li>
</ul>
</p>
<a name="Collections"/><h5>Collections: Lists and Dictionaries</h5>
<p>The top line of the <i>object editor</i> enables defining a list or dictionary of the chosen data type. Note that, unlike Python list or dictionary a CList or CDict is of only one data type. A CListist has three qualifiers: <i>listMinLength</i>, <i>listMaxLength</i> and <i>compare</i>. Compare should be set to 1 or -1 to switch on validity checking that successive values in the list are ascending or descending in value.  Note that the list item class must have a __cmp__() method.</p>

<a name="Column group"/><h5>Program column groups</h5>
<p>The class to handle program input of MTZ columns, CProgramColumnGroup, is a special case for which the plugin developer can specify the contents of the class in the XML file.  (It would be conceptually easier to insist on sub-classing!). The 'qualifiers' pane for this class includes a table in which you should enter the 'program' name of the column and the required column type. You can enter a comma-separated list of preferred column names, to be used by default if they appear in the MTZ file. The additional columns for <i>Partner</i> and <i>Offset</i> can probably be ignored.

<p>The <i>Save</i> option on the <i>File</i> menu can be used to save the contents tree to an XML file. You will get an additional dialog box requiring the header information for the XML file.</p>

<p>Beware that <b>defEd</b> currently has no tools for Undo/Redo or moving data items. These are not an immediate priority (c.f. more data classes, command templates etc.)  - please let me know if you find this a serious problem.</p>
""",
             'header' : """<html>
<h4>The DEF File Header</h4>
<p>The header info for the DEF xml file. You only need to
set the name (a single word) and version (in form n.m.i) for the plugin.
Optionally add short title that will be used in the GUI.</p>
</html>
"""
           }

CONTAINER_NAMES = ['inputData','controlParameters','outputData','guiControls']


class CDefEd(QtWidgets.QMainWindow):

  insts = None

  MARGIN = 2

  ERROR_CODES = { 101 : { 'description' : 'Unknown error reading def.xml file' },
                  102 : { 'description' : 'Unknown error displaying contents of def.xml file' },
                  103 : { 'description' :  'No data to save' },
                  104 : { 'description' :  'Unknown error trying to set qualifier'},
                  105 : { 'description' :  'Error attempting to create new object of class' },
                  106 : { 'description' :  'Attempting to change name of data object to name of existing object:' },
                  107 : { 'description' :  'Error finding container for new object' },
                  108 : { 'description' :  'Error creating new CString object' },
                  109 : { 'description' :  'Unknown error adding new object to data container:' },
                  110 : { 'description' :  'Need to select or create object' },
                  111 : { 'description' :  'Error making backup copy of current file' },
                  112 : { 'description' :  'Unknown error saving to file' }
                  }


  def __init__(self,parent=None):
    QtWidgets.QMainWindow.__init__(self,parent)
    self.setWindowTitle('CCP4i2 defEd - DEF file editor')
    if not CDefEd.insts: CDefEd.insts = self

    self.openFileDialog = None
    self.saveFileDialog = None
    self.headerDialog = None
    self.container = None
    self.edited = False
    self.helpDocument = {}
    self.initialiseSave()
    self.newNameCount = 0
    self.nameClashWarned = []
    self.infoHistory = []
    self.infoHistoryPosition = -1
    self.loadedFileName = None

    self.menuDefinitions = {}
    self.menuDefinitions['File'] = ['open','save','saveAs','print','help','close']
    self.menuDefinitions['Tools'] = ['makeWrapper','makePipeline','makeDataXml']

    self.contextMenu = QtWidgets.QMenu(self)
    
    for menuName,menuTitle in [['File','&File'],['Tools','&Tools']]:
      m = self.menuBar().addMenu(menuTitle)
      m.setObjectName(menuName)
      self.connect(m,QtCore.SIGNAL("aboutToShow()"),functools.partial(self.updateMenu,menuName))

    self.initialiseActionDefinitions()

    mainWidget = QtWidgets.QSplitter(QtCore.Qt.Vertical,self)
    
    topFrame = QtWidgets.QFrame()
    topFrame.setLayout(QtWidgets.QGridLayout())
    self.contentsTree = CContentsTree(self)
    self.contentsTree.setToolTip('The current data definition.')
    topFrame.layout().addWidget(self.contentsTree,0,0)
    self.buttonLayout = QtWidgets.QVBoxLayout()
    topFrame.layout().addLayout(self.buttonLayout,0,1)
    mainWidget.addWidget(topFrame)
    #self.connect(self.contentsTree,QtCore.SIGNAL('currentItemChanged(QTreeWidgetItem*,QTreeWidgetItem*)'),self.handleCurrentItemChanged)
    self.connect(self.contentsTree,QtCore.SIGNAL('itemSelectionChanged()'),self.handleSelectionChanged)


    # Create buttons
    self.appendContainerButton = CCP4Widgets.CPushButton(self,text='Append container')
    self.buttonLayout.addWidget(self.appendContainerButton)
    self.connect(self.appendContainerButton,QtCore.SIGNAL('clicked()'),functools.partial(self.handleNewContainer,True))
    self.insertContainerButton = CCP4Widgets.CPushButton(self,text='Insert container')
    self.buttonLayout.addWidget(self.insertContainerButton)
    self.connect(self.insertContainerButton,QtCore.SIGNAL('clicked()'),self.handleNewContainer)
    self.newButton = CCP4Widgets.CPushButton(self,text='Add object')
    self.buttonLayout.addWidget(self.newButton)
    self.connect(self.newButton,QtCore.SIGNAL('clicked()'),self.handleNew)
    self.deleteButton = CCP4Widgets.CPushButton(self,text='Delete object')
    self.buttonLayout.addWidget(self.deleteButton)
    self.connect(self.deleteButton,QtCore.SIGNAL('clicked()'),self.handleDelete)


    self.defEditor = CDefEditor(self)
    self.defEditor.setEditMode()
    mainWidget.addWidget(self.defEditor)
    self.connect(self.defEditor.classBrowser,QtCore.SIGNAL('leftMousePress'),self.handleClassBrowserClick)
    self.connect(self.defEditor.classBrowser,QtCore.SIGNAL('rightMousePress'),self.handleClassInfoRequest)
    self.connect(self.defEditor,QtCore.SIGNAL('qualifierEdited'),self.handleQualifierEdited)
    self.connect(self.defEditor,QtCore.SIGNAL('nameChanged'),self.handleNameEdit)
    self.connect(self.defEditor.collectionFrame.classCombo,QtCore.SIGNAL('currentIndexChanged(int)'),self.handleListStateChanged)
    self.connect(self.defEditor.collectionFrame,QtCore.SIGNAL('qualifierEdited'),self.handleCollectionQualifierEdited)
    self.connect(self.defEditor.columnGroupWidget,QtCore.SIGNAL('edited'),self.handleColumnGroupEdited)
    self.setCentralWidget(mainWidget)

    frame = QtWidgets.QFrame(self)
    frame.setFrameShape(frame.Box)
    frame.setLayout(QtWidgets.QVBoxLayout())
    frame.layout().setContentsMargins(CDefEd.MARGIN,CDefEd.MARGIN,CDefEd.MARGIN,CDefEd.MARGIN)
    frame.layout().setSpacing(CDefEd.MARGIN)
    self.classInfoToolBar = QtWidgets.QToolBar('Class Information',self)
    self.classInfoToolBar.setMaximumSize(9999,25)
    from . import CCP4GuiUtils
    CCP4GuiUtils.populateToolBar(self,toolBarWidget=self.classInfoToolBar,definition=['back','forward'])
    frame.layout().addWidget(self.classInfoToolBar)

    self.classInfoWidget = QtWidgets.QTextBrowser(self)
    self.classInfoWidget.setOpenLinks(False)
    self.connect(self.classInfoWidget,QtCore.SIGNAL('anchorClicked(const QUrl &)'),self.handleAnchorClicked)
    frame.layout().addWidget(self.classInfoWidget)
    mainWidget.addWidget(frame)

    self.connect(self.contentsTree,QtCore.SIGNAL('rightMousePress'),functools.partial(self.updateContextMenu,'Contents tree'))
    self.connect(self.defEditor.columnGroupWidget,QtCore.SIGNAL('rightMousePress'),functools.partial(self.updateContextMenu,'Column group'))
    self.connect(self.defEditor.collectionFrame,QtCore.SIGNAL('rightMousePress'),functools.partial(self.updateContextMenu,'Lists'))
    for w in [self.appendContainerButton,self.insertContainerButton,self.newButton,self.deleteButton]:
      self.connect(w,QtCore.SIGNAL('rightMousePress'),functools.partial(self.updateContextMenu,'Data objects'))
    self.connect(self.defEditor.nameEditor,QtCore.SIGNAL('rightMousePress'),functools.partial(self.updateContextMenu,'Data names'))
    self.connect(self.defEditor.nameCombo,QtCore.SIGNAL('rightMousePress'),functools.partial(self.updateContextMenu,'Data names'))
    self.connect(self.defEditor,QtCore.SIGNAL('rightMousePress'),self.updateContextMenu)

    self.showHelp()
    self.show()
    
    if self.container is None:
      from ..core import CCP4Container
      self.container = CCP4Container.CContainer(name='NEWCONTAINER')

    QtCore.QTimer.singleShot(0,self.handleCommandLine)

  def handleCommandLine(self):
    if len(sys.argv)>1: self.loadDefFile(sys.argv[1])


  def updateMenu(self,menuName):
    widget = self.menuBar().findChild(QtWidgets.QMenu,menuName)
    if widget is None:
      pass
    else:
      widget.clear()
      from . import CCP4GuiUtils
      CCP4GuiUtils.populateMenu(self,widget,self.menuDefinitions[menuName],default_icon='')

  def initialiseActionDefinitions(self):
    
    self.actionDefinitions = {}
    self.actionDefinitions['close'] = dict ( 
        text = 'Close editor',
        tip = 'Close this def file viewer',
        slot = self.close  )

    self.actionDefinitions['open'] = dict (
          text = "Open file",
          tip = "Open def.xml file",
          slot = self.handleOpenFile,
          icon = 'fileopen',
          shortcut = self.tr('cmd+O')
          )
    self.actionDefinitions['print'] = dict (
          text = "Print",
          tip = "Print list of data objects",
          slot = self.handlePrint,
          shortcut = self.tr('cmd+P')
          )
    
    def e(): return  self.edited and (self.loadedFileName is not None)
    self.actionDefinitions['save'] = dict (
          text = "Save",
          tip = "Backup current def.xml and save to same file name",
          slot = self.handleSaveFile,
          enabled = e
          )

    def e(): return self.edited
    self.actionDefinitions['saveAs'] = dict (
          text = "Save as",
          tip = "Save to def.xml file",
          slot = self.handleSaveAsFile,
          enabled = e
          )
    self.actionDefinitions['help'] = dict (
          text = "Help",
          tip = "Show main help",
          slot = self.showHelp
          )
    self.actionDefinitions['makeWrapper'] = dict (
          text = "Make wrapper plugin",
          tip = "Create directories and template files",
          slot = self.makeWrapper
          )
    self.actionDefinitions['makePipeline'] = dict (
          text = "Make pipeline project",
          tip = "Create directories and template files for a pipeline",
          slot = self.makePipeline,
          enabled = 1
          )
    self.actionDefinitions['makeDataXml'] = dict (
          text = "Make data XML file",
          tip = "Specify data for testing pipeline/wrapper",
          slot = self.makeDataXml,
          enabled = 0
          )
    self.actionDefinitions['back'] = dict (
          text = "Back",
          tip = "Go back to previous page",
          slot = self.handleBack,
          enabled = 0
          )
    self.actionDefinitions['forward'] = dict (
          text = "Forward",
          tip = "Go to next page",
          slot = self.handleForward,
          enabled = 0
          )

  def getActionDef(self,name,**info):
    return self.actionDefinitions.get(name,dict(text=name))

  def updateContextMenu(self,mode=None,event=None):
    if event is None or mode is None: return
    self.contextMenu.clear()
    a = self.contextMenu.addAction('Help: '+mode)
    self.connect(a,QtCore.SIGNAL('triggered(bool)'),functools.partial(self.updateHelp,mode))
    self.contextMenu.popup(QtCore.QPoint(event.globalX(),event.globalY()))
                  
   
  def handleOpenFile(self):
    rv = self.inviteSaveFile()
    if rv: return
    fileName = QtWidgets.QFileDialog.getOpenFileName(self,
          "Open def file", '', "Def files (*.def.xml)")                 
    if fileName is not None: self.loadDefFile(str(fileName))

  def handlePrint(self):
    fileName = QtWidgets.QFileDialog.getSaveFileName(self,
          "Save to text file", '', "Text files (*.txt)")
    if fileName is not None:
      text = ''
      root = self.contentsTree.topLevelItem(0)
      for indx in range(self.contentsTree.topLevelItemCount()):
        t = self.treeToText( self.contentsTree.topLevelItem(indx) ,0)
        text = text + t
    CCP4Utils.saveFile(fileName=fileName,text=text)

  def treeToText(self,parent,level):
    d = []
    for i in range(3):
      try:
        v = parent.data(i,0)
        d.append(v.__str__())
      except:
        d.append('')
    text = ''
    for i in range(level): text = text + '       '
    text = text  + "%-25s %-25s %s\n" % (d[0],d[1],d[2])
    level += 1
    for indx in range(parent.childCount()):
      t = self.treeToText( parent.child(indx) ,level)
      text = text + t
    return text
        
  
  def loadDefFile(self,fileName):
    from ..core.CCP4Container import CContainer
    self.container = CContainer(parent=self)
    #Expect loadContentsFromXml to return the error but put in try/except
    # just in case
    #try
    e = self.container.loadContentsFromXml(fileName)
    '''
    except CException as e:
      pass
    except:
      e = CException(self.__class__,101,fileName)
    '''

    if len(e)>0: 
      warningMessage(e)

    e = self.contentsTree.populate(container=self.container)
    if len(e)>0:
      warningMessage(e, 'Error drawing contents tree for def file')
    self.initialiseSave()
    if self.contentsTree.firstDataItem is None:
      self.defEditor.unSet()
    else:
      self.contentsTree.setCurrentItem(self.contentsTree.firstDataItem)
      self.contentsTree.scrollToItem(self.contentsTree.firstDataItem)

    
      '''
      firstDataObject = self.container.firstDataObject()
      if firstDataObject is not None:
        treeItemList = self.contentsTree.findItems(str(firstDataObject.objectName()),QtCore.Qt.MatchExactly)
        if len(treeItemList)>0:
          self.contentsTree.setCurrentItem(treeItemList[0])
          self.contentsTree.scrollToItem(treeItemList[0])
      '''
    self.loadedFileName = fileName
    
    self.setWindowTitle('CCP4i2 defEd: '+os.path.basename(fileName))
    self.edited = False

  def handleSaveFile(self):
    if self.loadedFileName is None: return
    if os.path.exists(self.loadedFileName):
      bak = 1
      bakFile = self.loadedFileName + '.'+str(bak)+'.bak'
      while os.path.exists(bakFile):
        bak = bak+1
        bakFile = self.loadedFileName + '.'+str(bak)+'.bak'
      try:
        shutil.copyfile(self.loadedFileName,bakFile)
      except:
        e = CException(self.__class__,111,bakFile)
        warningMessage(e)

    #print 'CDefEd.handleSaveFile',self.loadedFileName
    try:
      self.container.saveContentsToXml(self.loadedFileName)
    except CException as e:
      warningMessage(e)
    except Exception as e:
      e = CException(self.__class__,112,self.loadedFileName)
      warningMessage(e)
    self.edited = False
 
  def handleSaveAsFile(self):
    if not self.edited:
      QtWidgets.QMessageBox.warning(self,'No data to save','No changes to data to save')
      return
    fileName = QtWidgets.QFileDialog.getSaveFileName (self, 'Save def file',CCP4Utils.getCCP4I2Dir() , "Def files (*.xml)" )
    if len(fileName) ==0: return
    base,ext0 = os.path.splitext(str(fileName))
    ext1 = os.path.splitext(base)[1]
    #print 'handleSaveAsFile',base,ext0,'*',ext1
    if len(ext1) == 0 and ext0 == '.xml':
      fileName = base + '.def' + '.xml'
    self.saveDef(fileName)
    
    '''
    if self.saveFileDialog is None:
      self.saveFileDialog = QtWidgets.QFileDialog(self,"Save def file", '', "Def files (*.def.xml)")
      self.saveFileDialog.setDefaultSuffix('def.xml')
      self.saveFileDialog.setDirectory(CCP4Utils.getCCP4I2Dir())
      self.saveFileDialog.setConfirmOverwrite(False)
      self.saveFileDialog.setAcceptMode(QtWidgets.QFileDialog.AcceptSave)
      self.saveFileDialog.setModal(True)
      self.connect(self.saveFileDialog, QtCore.SIGNAL('filesSelected(const QStringList&)'),self.saveDef)
    self.saveFileDialog.exec_()
    '''

  def initialiseHeader(self,pluginName=''):
    if self.container is None: return
    h = self.container.addHeader()
    h.function='DEF'
    h.creationTime.setCurrentTime()
    h.userId.setCurrentUser()
    h.pluginName = pluginName
    h.ccp4iVersion = __version__

  def saveDef(self,fileName):
    #print 'CDefEd.saveData',fileName
    basename = os.path.splitext(fileName)[0]
    if isinstance(fileName,list):
      if len(fileName) == 0: return
      fileName = str(fileName[0])
    if self.container is None:
      e = CException(self.__class__,103,fileName)
      warningMessage(e)
      return
    
    if str(self.container.objectName()) == 'NEWCONTAINER':
      self.container.setObjectName(basename)
    self.setWindowTitle('CCP4i2 defEd: '+basename)

    self.showHelp('header')

    self.initialiseHeader(pluginName=basename)
    # Recreate headerDialog every time so the call to saveDef2 uses the right fileName
    self.headerDialog = CHeaderDialog(self,model=self.container.header)
    self.connect(self.headerDialog,QtCore.SIGNAL('closed'),functools.partial(self.saveDef2,fileName))
    done = self.headerDialog.show()
    
  def saveDef2(self,fileName):
    self.container.saveContentsToXml(fileName)
    self.edited = False
    self.loadedFileName = fileName

  def handleSelectionChanged(self):
    from ..core.CCP4Container import CContainer
    current = self.contentsTree.selectedItem()
    #print 'handleSelectionChanged',current
    dataObj = self.dataObjectFromTreeWidgetItem(current)
    if dataObj is None:
      pass
    ifContainer = isinstance(dataObj,CContainer)
      #print 'handleCurrentItemChanged Selected CContainer'
    self.defEditor.set(dataObj)
    self.defEditor.setEditMode(ifContainer)
    self.defEditor.setColumnGroupWidget(dataObj)
    self.nameClashWarned = []
      

  def currentDataObject(self):
    treeItem = self.contentsTree.selectedItem()
    if treeItem is None: return None
    dataObj = self.dataObjectFromTreeWidgetItem(treeItem)
    #print 'currentDataObject',treeItem,repr(dataObj)
    return dataObj

  def dataObjectFromTreeWidgetItem(self,item):
    nameList = []
    while isinstance(item,QtWidgets.QTreeWidgetItem):
      nameList.insert(0,str(item.text(0)))
      item = item.parent()
    #print 'dataObjectFromTreeWidgetItem',nameList
    obj = self.container
    try:
      for name in nameList:
        obj = obj.get(name)
    except:
      return None
    #print 'dataObjectFromTreeWidgetItem',repr(obj),obj.__class__
    return obj

  def redrawCurrentItem(self):
    # Update the text of the current item in the content browser
    dataObject = self.currentDataObject()
    if dataObject is None: return
    currentItem = self.contentsTree.currentItem()
    if currentItem is None: return
    self.drawContentTreeItem(currentItem,dataObject)

  def drawContentTreeItem(self,item,dataObject):
    item.setText(0,str(dataObject.objectName()))
    item.setText(1,dataObject.__class__.__name__)

    qText = ''
    qualis = dataObject.qualifiers(default=False)
    for key,value in list(qualis.items()):
      if value is None or value is NotImplemented:
        pass
      elif isinstance(value,list) and len(value)>0 and isinstance(value[0],CCP4Data.CData):
        # This is to deal withthe CProgramColumnGroup columnGroups qualifiers
        l = []
        for obj in value: l.append(str(obj.get()))
        qText = str(l)
      elif isinstance(value,CCP4Data.CData):
        qText = qText+key+':'+str(value.get())
      else:
        try:
          qText = qText+key+':'+str(value)+ ' '
        except:
          pass    
    if isinstance(dataObject,CCP4Data.CCollection):
      qText = qText + 'subItem:'+dataObject.subItemClassName()+' '
      qualis = dataObject.subItemQualifiers(default=False)
      for key,value in list(qualis.items()):
        if value is not None and value is not NotImplemented:
          try:
            qText = qText+key+':'+str(value)+ ' '
          except:
            pass
    #print 'drawContentTreeItem',repr(dataObject),dataObject.objectName(),qualis,qText
    item.setText(2,qText)


  def handleDelete(self):
    dataObject = self.currentDataObject()
    if dataObject is None: return

    container = dataObject.parent()
    if container is not None:
      container.deleteObject(dataObject.objectName())
      self.contentsTree.deleteCurrentItem()
 
    #print 'handleDelete',repr(container),container.contents().keys()
    self.edited= True
  

  def insertPositionInContainer(self):
    from ..core.CCP4Container import CContainer
    # Figure out which container and which afterObject
    dataObject = self.currentDataObject()
    if dataObject is None:
      #Put new item after last item at top level
      container = self.container
      afterObject = None
    elif isinstance(dataObject,CContainer):
      #put new item at end of container
      container = dataObject
      if len(container.dataOrder())>0:
        afterObject = container.dataOrder()[-1]
      else:
        afterObject = None
    else:
      #Assume is a CData in a container
      container = dataObject.parent()
      if container is None or not isinstance(container,CContainer):
        # Very odd!
        e = CException(self.__class__,107)
        warningMessage(e)
        return
      afterObject = str(dataObject.objectName())
    #print 'handleNew',repr(container),repr(afterObject)
    return container,afterObject

  def newName(self,rootName):
    if self.newNameCount == 0:
      newName = rootName
    else:
      newName = rootName+'_'+str(self.newNameCount)
    self.newNameCount = self.newNameCount + 1
    return newName
  
  def handleNew(self):
    from ..core.CCP4Data import CString
    container,afterObject = self.insertPositionInContainer()   
    # unique name
    newName = self.newName('NONAME')

    # Make a default CString object
    try:
      newObject = CString(name=newName,parent=container)
    except:
      e = CException(self.__class__,108)
      warningMessage(e)
      return
    #print 'newObject',repr(newObject)
    # Add object to container
    try:
      container.addObject(newObject,name=newName,afterObject=afterObject)
    except CException as e:
      warningMessage(e)
      return
    except:
      e = CException(self.__class__,109,str(container.objectName()))
      warningMessage(e)
      return

    # Add to contentsTree
    newItem = self.contentsTree.insertObject(name=newName,defn=container.contents(newName))
    self.contentsTree.setCurrentItem(newItem)
    self.defEditor.classBrowser.setSelectedClass('CString')
    self.edited= True
    self.nameClashWarned = []
    
    
  def handleNewContainer(self,append=False):
    from ..core.CCP4Container import CContainer
    if append:
      container = self.container
      afterObject = None
    else:
      container,afterObject = self.insertPositionInContainer()
    
    # unique name
    newName = self.newName('CONTAINER')

    # Make default object
    try:
      newObject = CContainer(name=newName,parent=container)
    except:
      e = CException(self.__class__,108)
      warningMessage(e)
      return
    #print 'newObject',repr(newObject)
    # Add object to container
    try:
      container.addObject(newObject,name=newName,afterObject=afterObject)
    except CException as e:
      warningMessage(e)
      return
    except:
      e = CException(self.__class__,109,str(container.objectName()))
      warningMessage(e)
      return

    # Add to contentsTree
    if append:
      newItem = self.contentsTree.appendObject(name=newName,defn=container.contents(newName))
    else:
      newItem = self.contentsTree.insertObject(name=newName,defn=container.contents(newName))

    self.defEditor.setNameComboDefault(newName)
      
    self.contentsTree.setCurrentItem(newItem)
    self.edited= True
    self.nameClashWarned = []

  def handleClassInfoRequest(self,modelItem=None,className=None):
    #print 'handleClassInfoRequest',modelItem,modelItem.text()
    if className is None: className = str(modelItem.text())
    self.setClassInfoDocument(className)
    self.infoHistory = [className]
    self.infoHistoryPosition = len(self.infoHistory)-1
    self.setInfoHistoryActions()

  def setClassInfoDocument(self,className):
    classInfo = CClassInfo(className=className,parent=self)
    self.classInfoWidget.setDocument(classInfo)
    
  def handleAnchorClicked(self,url):
    className = (str(url)).split('#')[-1]
    self.setClassInfoDocument(className=className)
    if self.infoHistoryPosition < len(self.infoHistory)-1:
      del self.infoHistory[self.infoHistoryPosition+1:]
    self.infoHistory.append(className)
    self.infoHistoryPosition = len(self.infoHistory)-1
    self.setInfoHistoryActions()

  def handleBack(self):
    self.infoHistoryRestore(-1)
  def handleForward(self):
    self.infoHistoryRestore(1)

  def infoHistoryRestore(self,increment=0):
    if increment<=0:
      new=max(0,self.infoHistoryPosition+increment)
    else:
      new = min(len(self.infoHistory)-1,self.infoHistoryPosition+increment)
    if new != self.infoHistoryPosition:
      self.infoHistoryPosition = new
      self.setClassInfoDocument(self.infoHistory[self.infoHistoryPosition])
    self.setInfoHistoryActions()
    
  def updateHelp(self,mode=None):
    #print 'updateHelp',mode
    if mode in HELPTEXT:
      self.showHelp(mode)
    elif self.helpDocument.get('Introduction',None) is not None and self.helpDocument['Introduction'].find('name="'+mode+'"'):
      self.showHelp('Introduction')
      self.classInfoWidget.scrollToAnchor(mode)      
    else:
      cls = DATAMANAGER().getClass(mode)
      if cls is not None: self.handleClassInfoRequest(className=mode)

  def showHelp(self,name=None):
    if name is None: name = 'Introduction'
    if self.helpDocument.get(name,None) is None:
      self.helpDocument[name] = QtGui.QTextDocument(self)
      self.helpDocument[name].setHtml(HELPTEXT[name])
    self.infoHistory = []
    self.infoHistoryPosition = -1
    self.classInfoWidget.setDocument(self.helpDocument[name])
    self.setInfoHistoryActions()
    

  def setInfoHistoryActions(self):
    nS = len(self.infoHistory)
    nP = self.infoHistoryPosition
    self.findChild(QtWidgets.QAction,'forward').setEnabled(nS>0 and nP<(nS-1))
    self.findChild(QtWidgets.QAction,'back').setEnabled(nS>0 and nP>0)


  def close(self):
    if self.edited:
      rv = self.inviteSaveFile()
      if rv: return
    QtWidgets.QMainWindow.close(self)

  def initialiseSave(self):
    self.nAutoSaved = 0
    self.autoSavedStatus = []
    self.autoSavedPosition = -1

  def autoSave(self):
    self.nAutoSaved = self.nAutoSaved + 1
    self.autoSavedStatus.append(self.getSaveStatus())
    self.autoSavedPosition = len(self.autoSavedStatus)-1

  def getSaveStatus(self):
    current = self.currentDataObject()
    if current is not None:
      currentName = str(current.objectName())
    else:
      currentName = None
    return { 'container' : self.container.saveContentsToEtree(), 'current': currentName }
  
  def updateStatus(self,status={}):
    self.container.clear()
    self.container.loadContentsFromEtree(status['container'])
    if status['current'] is None:
      dataObj = None
    else:
      dataObj = self.container.get(status['current'])
    if dataObj is None:
      self.defEditor.unSet()
    else:
      self.defEditor.set(dataObj)

  def handleQualifierEdited(self,qualifierName):
    #print 'handleQualifierEdited',qualifierName
    textValue = self.defEditor.getQualifierValue(qualifierName)
    if textValue is None:
      print('Error getting qualifier value from widget for',qualifierName)
      return
    try:
      exec('value='+textValue)
    except:
      value = textValue
    #print 'handleQualifierEdited',textValue,value
    dataObject = self.currentDataObject()
    if dataObject is not None:
      if isinstance(dataObject,CCP4Data.CCollection):
        subObj = dataObject.subItemObject()
        if subObj is not None:
          try:
            subObj.setQualifier(qualifierName,value)
          
          except CException as e:
            warningMessage(e)
          except:
            e = CException(self.__class__,104,qualifierName)
            warningMessage(e)
          
      else:
        try:
          dataObject.setQualifier(qualifierName,value)        
        except CException as e:
          warningMessage(e)
        except:
          e = CException(self.__class__,104,qualifierName)
          warningMessage(e)
      
      self.redrawCurrentItem()
      self.edited= True

  def handleCollectionQualifierEdited(self,qualifierName):
    value = self.defEditor.collectionFrame.getQualifiers().get(qualifierName)
    #print 'CDefEd.handleCollectionQualifierEdited',qualifierName,value,type(value)
    dataObject = self.currentDataObject()
    #print 'handleCollectionQualifierEdited dataObject',repr(dataObject)
    if dataObject is not None:
      dataObject.setQualifier(qualifierName,value)
      self.redrawCurrentItem()
      self.edited= True

  def handleClassBrowserClick(self,item):
    
    className = str(item.text())
    if className == self.defEditor.currentClassName: return
    if self.defEditor.currentClassName == 'CContainer': return
        
    dataObject = self.currentDataObject()
    if dataObject is None:
      # We dont have any item selected for editing but let the user see
      # The qualifiers editor for the chosen class
      self.defEditor.changeDataClass(className)
      self.defEditor.setDefaultQualifiers()
      return
    
    name = str(dataObject.objectName())
    newObject = self.makeNewObject(className,self.defEditor.collectionFrame.getContainerClassName(),dataObject)
    #print 'handleClassBrowserClick newObject',repr(newObject)
    container = dataObject.parent()
    container.replaceObject(newObject,name=name)
    
    self.defEditor.changeDataClass(className)
    self.defEditor.setQualifiers(newObject.qualifiers())
    self.defEditor.setColumnGroupWidget(newObject)
    self.redrawCurrentItem()
    self.edited= True

  def handleListStateChanged(self,indx):
    # User has toggled the list button so must convert to/from a list
    dataObject = self.currentDataObject()
    if dataObject is None:
      e = CException(self.__class__,110)
      warningMessage(e)
      return
    #print 'handleListStateChanged', state,isinstance(dataObject,CList)
    if isinstance(dataObject,CCP4Data.CCollection):
      className = dataObject.subItemClassName()
    else:
      className = dataObject.__class__.__name__
    newObject = self.makeNewObject(className,self.defEditor.collectionFrame.getContainerClassName(),dataObject)
    if newObject is None: return
    container = dataObject.parent()
    name = str(dataObject.objectName())
    container.replaceObject(newObject,name=name)
    self.redrawCurrentItem()
    self.edited= True
    
  def makeNewObject(self,className=None,containerClassName=None,oldObject=None):
    newObject = None
    name = str(oldObject.objectName())
    if containerClassName is None:
      try:
        newObject = DATAMANAGER().getClass(className)(name=name)
      except:
        e = CException(self.__class__,105,className)
        warningMessage(e)
        return
      if isinstance(oldObject,CCP4Data.CCollection) and oldObject.subItemClassName() == className:
        try:
          newObject.setQualifiers(oldObject.subItemQualifiers())
        except:
          pass
    else:
      if isinstance(oldObject, DATAMANAGER().getClass(containerClassName)):
        #Change the item class in a list/dict
        if className == oldObject.subItemClassName(): return None
        oldObject.setSubItem( { 'class' : DATAMANAGER().getClass(className) , 'qualifiers' : {} } )
        newObject = oldObject
      else:
        cls =  DATAMANAGER().getClass(className)
        if isinstance(oldObject,cls):
          subItemQualifiers=oldObject.qualifiers(default=False)
        elif isinstance(oldObject,CCP4Data.CCollection) and oldObject.subItemClassName() == className:
          subItemQualifiers=oldObject.subItemQualifiers(default=False)
        else:
          subItemQualifiers = {}
        if ['CList','CDict'].count(containerClassName):
          qualifiers=self.defEditor.collectionFrame.getQualifiers()
        else:
          qualifiers = {}
        newObject =  DATAMANAGER().getClass(containerClassName)(name=name,qualifiers=qualifiers,subItem={ 'class' : cls, 'qualifiers' : subItemQualifiers } )
    #print 'makeNewObject',repr(newObject)
    return newObject

      

  def handleNameEdit(self):
    newName = self.defEditor.getName()
    
    dataObject = self.currentDataObject()
    if dataObject is None: return
    oldName = dataObject.objectName()
    if newName == oldName: return
    
    container = dataObject.parent()
    #print 'handleNameEdit',oldName,newName,container
    try:
      container.renameObject(oldName,newName)
    except CException as e:
      if not self.nameClashWarned.count(newName):
        self.nameClashWarned.append(newName)
        warningMessage(e)

    # Beware until the content tree is updated to reflect new name
    # there is descrepancy between gui and container and dataObjectFromTreeWidgetItem and
    # currentDataObject don't work
    self.drawContentTreeItem(self.contentsTree.selectedItem(),dataObject)
    self.defEditor.currentContentName = newName
    self.edited= True


  def handleColumnGroupEdited(self,iRow,iCol):
    data = self.defEditor.columnGroupWidget.get()
    dataObject = self.currentDataObject()
    #print 'handleColumnGroupEdited',data,dataObject.objectName()
    if dataObject is None: return
    dataObject.setColumnGroupQualifier(columnGroup=data,initialise=True)
    self.drawContentTreeItem(self.contentsTree.selectedItem(),dataObject)
    self.edited = True
    
  def inviteSaveFile(self):
    #print 'inviteSaveFile',self.edited
    if not self.edited: return
    '''
    self.saveMessage = QtWidgets.QMessageBox.question(self)
    self.saveMessage.setWindowTitle('CCP4DefEd Save changes')
    self.saveMessage.setText('Save changes to file?')
    self.saveMessage.setStandardButtons(QtWidgets.QMessageBox.Save|QtWidgets.QMessageBox.Cancel)
    rv = self.saveMessage._exec()
    '''
    rv = QtWidgets.QMessageBox.question(self,'CCP4DefEd Save changes','Save changes to file?',
                            QtWidgets.QMessageBox.Save|QtWidgets.QMessageBox.Cancel|QtWidgets.QMessageBox.Discard,
                                                  QtWidgets.QMessageBox.Cancel)
    if rv == QtWidgets.QMessageBox.Cancel:
      return 1
    elif rv == QtWidgets.QMessageBox.Discard:
      return 0
    else:
      self.handleSaveFile()
      return 0

  def makeWrapper(self):
    if not hasattr(self,'wrapperDialog'):
      self.wrapperDialog = CMakeWrapperPlugin(self)
    self.wrapperDialog.show()


  def makePipeline(self):
    if not hasattr(self,'pipelineDialog'):
      self.pipelineDialog = CMakePipeline(self)
    self.pipelineDialog.show()

  def makeDataXml(self):
    pass

    
    
class CContentsTree(QtWidgets.QTreeWidget):

  ERROR_CODES = { 101 : { 'description' : 'Failed to draw container item in tree' },
                  102 : { 'description' : 'Failed to draw data item in tree' },
                  103 : { 'description' : 'Failed to draw qualifiers for data item in tree' },
                  104 : { 'description' : 'No currently selected item to delete' }
                  }

  def __init__(self,parent):
    QtWidgets.QTreeWidget.__init__(self,parent)
    self.setColumnCount(3)
    self.setHeaderLabels(['Name','Class','Qualifiers'])
    self.setItemsExpandable(True)
    self.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
    self.firstDataItem = None

  def sizeHint(self):
    return QtCore.QSize(400,200)

  
  def mousePressEvent(self,event):
    if event.button() == QtCore.Qt.RightButton:
      self.emit(QtCore.SIGNAL('rightMousePress'),event)
      event.accept()
    else:
      QtWidgets.QTreeWidget.mousePressEvent(self,event)


  def populate(self,container=None):
    e = CErrorReport()
    self.clear()
    if container is None: return
    #top = QtWidgets.QTreeWidgetItem([str(container.objectName())])
    #self.addTopLevelItem(top)
    e.extend(self.loadContainer(container=container,parent=None))
    return e

  def clear(self):
    self.firstDataItem = None
    QtWidgets.QTreeWidget.clear(self)

  def loadContainer(self,container=None,parent=None):
    from ..core import CCP4Container
    e = CErrorReport()
    for name in container.dataOrder():
      defn = container.contents(name)
      if defn is not None:
        try:
          item = self.makeTreeWidgetItem(name,defn)
        except CException as err:
          e.extend(err)
        else:
          if parent is None:
            self.addTopLevelItem(item)
          else:
            parent.addChild(item)
          if defn['class'] == CCP4Container.CContainer:
            e.extend(self.loadContainer(container.__getattr__(name),item))
          else:
            if self.firstDataItem is None: self.firstDataItem = item
    return e

  def makeTreeWidgetItem(self,name,defn):
    from ..core import CCP4Container, CCP4Data
    item = None
    if defn['class'] == CCP4Container.CContainer:
      item = QtWidgets.QTreeWidgetItem([name],1001)
      #item.setFlags(QtCore.Qt.ItemIsEnabled)
      item.setFlags(QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEnabled)
    else:
      qText = ''
      if isinstance(defn['class'],CCP4Data.CCollection):
        qText = 'subItem:'+defn['subItem']['class'].__name__
        qualifiers = defn['subItemQualifiers']
      else:
        qualifiers = defn['qualifiers']
      try:
        for key,value in list(qualifiers.items()):
          if isinstance(value,list) and len(value)>0 and isinstance(value[0],CCP4Data.CData):
            # This is to deal withthe CProgramColumnGroup columnGroups qualifiers
            l = []
            for item in value: l.append(str(item.get()))
            qText = str(l)
          elif isinstance(value,CCP4Data.CData):
            qText = qText+key+':'+str(value.get())
          else:
            qText = qText+key+':'+str(value)+ ' '
      except:
        raise CException(self.__class__,103,name)
      #print 'qualifiers',defn['qualifiers'],qText
      try:
        item = QtWidgets.QTreeWidgetItem([name,defn['class'].__name__,qText],1002)
        item.setFlags(QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEnabled)
      except:
        raise CException(self.__class__,102,name)
    return item


  def insertObject(self,name,defn):
    item = self.makeTreeWidgetItem(name,defn)
    currentItem = self.selectedItem()
    #print 'insertObject',name,currentItem,
    #if currentItem is not None: print int(currentItem.type())
    if currentItem is None:
      self.addTopLevelItem(item)
    elif int(currentItem.type()) == 1001:
      currentItem.addChild(item)
    else:
      containerItem = currentItem.parent()
      #print 'insertObject parent',containerItem
      if containerItem is None:
        indx = self.indexOfTopLevelItem(currentItem)
        self.insertTopLevelItem(indx+1,item)
      else:
        indx = containerItem.indexOfChild(currentItem)
        #print 'insertObject indx',indx
        containerItem.insertChild(indx+1,item)
    return item

  def appendObject(self,name,defn):
    item = self.makeTreeWidgetItem(name,defn)
    self.addTopLevelItem(item)
    return item

  def deleteCurrentItem(self):
    currentItem = self.currentItem()
    if currentItem is None:
      raise CException(self.__class__,104)
    parent = currentItem.parent()
    if parent is None:
      # its a top level item
      indx = self.indexOfTopLevelItem(currentItem)
      nextIndx = min (indx, self.topLevelItemCount()-2)
      removed = self.takeTopLevelItem(indx)
      if nextIndx>=0:
        self.setCurrentItem(self.topLevelItem(nextIndx))
      removed = None
    else:
      indx = parent.indexOfChild(currentItem)
      try:
        nextItem = parent.child(indx+1)
      except:
        nextItem = parent
      parent.removeChild(currentItem)
      self.setCurrentItem(nextItem)
      currentItem = None
    
  def selectedItem(self):
    seleList = self.selectedItems()
    #print 'contentsTree.selectedItem',seleList
    if len(seleList) == 0:
      return None
    elif len(seleList) == 1:
      return seleList[0]
    else:
      #print 'Error in CContentsTree - more than one selected item'
      return seleList[0]
  
    

class CDataClassBrowser(QtWidgets.QTreeView):

  def __init__(self,parent):
    QtWidgets.QTreeView.__init__(self,parent)
    #self.setHeaderLabels(['Class','Description'])
    self.setItemsExpandable(True)
    self.setHeaderHidden(True)
    self.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
    #self.connect(self,QtCore.SIGNAL('doubleCLicked(const QModelIndex &)'),self.doubleClicked)
    self.populate()

  def populate(self):
    model = DATAMANAGER().buildQStandardItemModel(parent=self)
    self.setModel(model)

  def selectedClass(self):
    selList = self.selectedItems()
    #print 'selectedClass',selList
    if len(selList) == 1:
      clsName = str(selList[0].text())
      #print 'changeDataClass clsName',clsName
      return clsName
    else:
      return None

  def setSelectedClass(self,clsName=None):
    if clsName is None: return
    itemList = self.searchModelTree(clsName)
    #print 'CDataClassBrowser.setSelectedClass',clsName,type(clsName),itemList
    self.blockSignals(True)
    if len(itemList)==1:
      modelIndex = itemList[0].index()
      if modelIndex is not None:
        #print 'CDataClassBrowser.setSelectedClass modelIndex', modelIndex
        self.setCurrentIndex(modelIndex)
    self.blockSignals(False)

 
  def searchModelTree(self,text):
    # Surely there is  method to search a model for specific text label???
    # But I cant figure what is is so try this ..
    mrow = 0
    hits = []
    #print 'searchModelTree',self.model(),self.model().rowCount()
    while mrow < self.model().rowCount():
      moduleItem = self.model().item(mrow,0)
      irow = 0
      #print 'moduleItem',mrow,moduleItem,moduleItem.rowCount()
      while irow < moduleItem.rowCount():
        clsItem = moduleItem.child(irow,0)
        #print 'clsItem',str(clsItem.text())
        if str(clsItem.text()) == text:
          hits.append(clsItem)
        irow = irow + 1
      mrow = mrow + 1

    return hits

  
  def mousePressEvent(self,event):
    if [QtCore.Qt.RightButton,QtCore.Qt.LeftButton].count(event.button()):
      modelIndex = self.indexAt(event.pos())
      if modelIndex is not None:
        model = modelIndex.model()
        if not model:
          QtWidgets.QTreeView.mousePressEvent(self,event)
          return
        item = model.itemFromIndex(modelIndex)
        #print 'mousePressEvent',modelIndex,model,item
        if item.rowCount()>0:
          #event.ignore()
          pass
        else:
          if event.button() == QtCore.Qt.RightButton:
            self.emit(QtCore.SIGNAL('rightMousePress'),item)
          else:
            self.emit(QtCore.SIGNAL('leftMousePress'),item)
      
    #print 'mousePressEvent calling QTreeView'
    QtWidgets.QTreeView.mousePressEvent(self,event)
  
    

  '''
  def doubleClicked(self,modelIndex):
    print 'doubleClicked',modelIndex
  '''


class CDefEditor(QtWidgets.QFrame):

  ERROR_CODES = {
    101 : { 'description' : 'Failed to retieve qualifier definition for class' },
    102 : { 'description' : 'Failed to draw qualifier definition widget for class' },
    103 : { 'description' : 'Failed to set default value qualifier' },
    104 : { 'description' : 'Error deleting qualifier widget' }
    }

  def __init__(self,parent=None):
    QtWidgets.QFrame.__init__(self,parent)
    self.setLayout(QtWidgets.QVBoxLayout())
    self.layout().setSpacing(CDefEd.MARGIN)
    self.layout().setContentsMargins(CDefEd.MARGIN,CDefEd.MARGIN,CDefEd.MARGIN,CDefEd.MARGIN)
    self.setFrameShape(self.Box)

    self.currentClassName = None

    # The name line
    lineLayout = QtWidgets.QHBoxLayout()
    lineLayout.addWidget(QtWidgets.QLabel('Name:',self))
    lineLayout.setSpacing(CDefEd.MARGIN)
    lineLayout.setContentsMargins(CDefEd.MARGIN,CDefEd.MARGIN,CDefEd.MARGIN,CDefEd.MARGIN)
    self.nameStack =  QtWidgets.QStackedWidget(self)   
    self.nameEditor = CCP4Widgets.CLineEdit(self)
    self.nameStack.addWidget(self.nameEditor)
    # Beware - this will not work with CCp4Combo widget
    self.nameCombo= QtWidgets.QComboBox(self)
    self.nameCombo.addItem('CONTAINER')
    for item in CONTAINER_NAMES:self.nameCombo.addItem(item)
    self.nameCombo.setEditable(True)      
    self.nameStack.addWidget(self.nameCombo)
    lineLayout.addWidget(self.nameStack)
    self.layout().addLayout(lineLayout)
    

    self.containerFrameStack = QtWidgets.QStackedWidget(self)
    self.collectionFrame = CCollectionFrame(self)
    self.containerFrameStack.addWidget(self.collectionFrame)
    dummyLine = QtWidgets.QFrame(self)
    dummyLine.setLayout(QtWidgets.QHBoxLayout())
    dummyLine.layout().addWidget(QtWidgets.QLabel('The selected item is a container - you can edit the name and add other objects'))
    self.containerFrameStack.addWidget(dummyLine)
    self.layout().addWidget(self.containerFrameStack)

    self.splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal,self)
    self.layout().addWidget(self.splitter)
    self.layout().setStretch(1,10)
    
    leftFrame = QtWidgets.QFrame(self)
    leftFrame.setLayout(QtWidgets.QGridLayout())
    leftFrame.layout().setSpacing(CDefEd.MARGIN)
    leftFrame.layout().setContentsMargins(CDefEd.MARGIN,CDefEd.MARGIN,CDefEd.MARGIN,CDefEd.MARGIN)

    lineLayout = QtWidgets.QHBoxLayout()
    lineLayout.addWidget(QtWidgets.QLabel('Class:',self))
    self.classLabel = QtWidgets.QLabel(self)
    lineLayout.addWidget(self.classLabel)
    lineLayout.addStretch(5)
    leftFrame.layout().addLayout(lineLayout,0,0)   
    self.classBrowser = CDataClassBrowser(self)
    leftFrame.layout().addWidget(self.classBrowser,1,0)
    self.splitter.addWidget(leftFrame)

    self.rightFrame =  QtWidgets.QFrame(self)
    self.rightFrame.setLayout(QtWidgets.QGridLayout())
    self.rightFrame.layout().setSpacing(CDefEd.MARGIN)
    self.rightFrame.layout().setContentsMargins(CDefEd.MARGIN,CDefEd.MARGIN,CDefEd.MARGIN,CDefEd.MARGIN)
    
    self.connect(self.nameEditor,QtCore.SIGNAL('editingFinished()'),functools.partial(self.emit,QtCore.SIGNAL('nameChanged')))
    self.connect(self.nameCombo,QtCore.SIGNAL('editTextChanged()'),functools.partial(self.emit,QtCore.SIGNAL('nameChanged')))
    self.connect(self.nameCombo,QtCore.SIGNAL('currentIndexChanged(int)'),functools.partial(self.emit,QtCore.SIGNAL('nameChanged')))

    vLayout = QtWidgets.QVBoxLayout()
    self.rightFrame.layout().addLayout(vLayout,0,0)
    
    self.qFrame = QtWidgets.QFrame(self)
    self.qFrame.setLayout(QtWidgets.QGridLayout())
    self.rightFrame.layout().addWidget(self.qFrame,1,0)

    self.columnGroupWidget = CProgramColumnGroupQualifierView(self)
    self.columnGroupWidget.setVisible(0)
    self.columnGroupWidget.blockSignals(True)
    self.rightFrame.layout().addWidget(self.columnGroupWidget,2,0)
    
    self.splitter.addWidget(self.rightFrame)
    self.splitter.setSizes([150,300])

  def deleteQualiferWidget(self):
    #print 'deleteQualiferWidget',self.qFrame
    if self.qFrame is None: return
    self.splitter.widget(1).layout().removeWidget(self.qFrame)
    self.qFrame.deleteLater()
    self.qFrame = None


  def drawQualifiers(self,cls):
    # Not using CCP4Widgets 'cos want to keep the standard context menu with paste function
    try:
      obj = cls()
      qDefs = obj.qualifiersDefinition()
      qOrder = obj.qualifiersOrder()
    except:
      e = CException(self.__class__,101,cls.__name__)
      warningMessage(e)
      return
    
    self.deleteQualiferWidget()
    self.qFrame = QtWidgets.QFrame(self)
    self.qFrame.setLayout(QtWidgets.QGridLayout())
    self.rightFrame.layout().addWidget(self.qFrame,1,0)
    
    e = CErrorReport()
    iRow = -1
    for quali in qOrder:
      defn = qDefs.get(quali,{})
      #print 'drawQualifiers',quali,defn
      try:
        widget = None
        qType = defn.get('type',None)
        editable =  defn.get('editable',True)
        if editable:
          if qType is None: qType = obj.pythonType()
          if qType == str:
            #widget = CCP4Widgets.CLineEdit(self.qFrame)
            widget = QtWidgets.QLineEdit(self.qFrame)
          elif qType == int:
            #widget = CCP4Widgets.CLineEdit(self.qFrame)
            widget = QtWidgets.QLineEdit(self.qFrame)
            validator = QtGui.QIntValidator(self.qFrame)
            widget.setValidator(validator)
          elif qType == float:
            #widget = CCP4Widgets.CLineEdit(self.qFrame)
            widget = QtWidgets.QLineEdit(self.qFrame)
            validator = QtGui.QDoubleValidator(self.qFrame)
            widget.setValidator(validator)
          elif qType == bool:
            widget = CCP4Widgets.CCheckBox(self.qFrame)
            if quali == 'default' and qDefs.get('allowUndefined',True):
              widget.setTristate()
          elif qType == 'enumerator':
            menu = defn.get('menu',[])
            widget = CCP4Widgets.CComboBox(self.qFrame)
            for item in menu:
              widget.addItem(item)
          else:
            #widget = CCP4Widgets.CLineEdit(self.qFrame)
            widget = QtWidgets.QLineEdit(self.qFrame)

          if widget is not None:
            widget.setObjectName(str(quali))
            iRow = iRow + 1
            #print 'drawQualifiers',quali,qType,widget,iRow
            self.qFrame.layout().addWidget(widget,iRow,1)
            label = QtWidgets.QLabel(quali,self.qFrame)
            self.qFrame.layout().addWidget(label,iRow,0)
            if qType == bool:
              self.connect(widget,QtCore.SIGNAL('clicked(bool)'),functools.partial(self.emit,QtCore.SIGNAL('qualifierEdited'),quali))
            else:
              self.connect(widget,QtCore.SIGNAL('editingFinished()'),functools.partial(self.emit,QtCore.SIGNAL('qualifierEdited'),quali))
            #self.connect(widget,QtCore.SIGNAL('rightMousePress'),functools.partial(self.emit,QtCore.SIGNAL('rightMousePress'),cls.__name__))
      except:
        e.append(self.__class__,102,quali)

    if len(e)>0: warningMessage(e)

  def setDefaultQualifiers(self):
    cls = DATAMANAGER().getClass(self.currentClassName)
    if cls is None: return
    obj = cls()
    if isinstance(cls,CCP4Data.CCollection):
      self.setQualifiers(obj.subItemObject().qualifiers())
    else:
      self.setQualifiers(obj.qualifiers())

  def setQualifiers(self,qualifiers):
    e = CErrorReport()
    if self.qFrame is None: return
    #print 'setObjectQualifiers',qualifiers.items()
    for key,value in list(qualifiers.items()):
      try:
        widget = self.qFrame.findChild(QtWidgets.QWidget,key)
        #print 'setQualifiers',key,type(widget)
        if widget is not None and value is not NotImplemented:
          if isinstance(widget,QtWidgets.QLineEdit):
            if value is None:
              widget.setText('')
            else:
              widget.setText(str(value))
          elif isinstance(widget,QtWidgets.QCheckBox):
            if value is None or value is NotImplemented:
              if widget.isTristate(): widget.setCheckState(QtCore.Qt.PartiallyChecked)
            else:
              widget.setChecked(value)
          else:
            e.append(self.__class__,103,key)
      except:
         e.append(self.__class__,103,key)
         
    if len(e)>0: warningMessage(e)

        
  def changeDataClass(self,clsName):
    if clsName is None: return
    #print 'changeDataClass',self.currentClassName
    self.classLabel.setText(clsName)
    cls = DATAMANAGER().getClass(clsName)
    if cls is None: return
    self.drawQualifiers(cls)
    self.currentClassName = clsName
    self.setDefaultQualifiers()

  def set(self,dataObj=None):
    if dataObj is None or not isinstance(dataObj,CCP4Data.CData): return
    self.nameEditor.setText(str(dataObj.objectName()))
    if isinstance(dataObj,CCP4Data.CCollection):
      self.collectionFrame.set(dataObj)
      self.changeDataClass(dataObj.subItemClass().__name__)
      self.classBrowser.setSelectedClass(dataObj.subItemClassName())
    else:
      self.collectionFrame.set()
      self.changeDataClass(str(dataObj.__class__.__name__))
      self.classBrowser.setSelectedClass(str(dataObj.__class__.__name__))

    self.setColumnGroupWidget(dataObj)
        
    if isinstance(dataObj,CCP4Data.CCollection):
      self.setQualifiers(dataObj.subItemQualifiers())
    else:
      self.setQualifiers(dataObj.qualifiers())


  def unSet(self):
    self.currentClassName = None
    self.nameEditor.setText('')
    self.nameCombo.clearEditText()
    self.classLabel.setText('')
    self.deleteQualiferWidget()
    
  def getQualifierValue(self,name=None):
    widget = self.qFrame.findChild(QtWidgets.QWidget,name)
    if widget is None: return None
    if isinstance(widget,QtWidgets.QCheckBox):
      if widget.isTristate() and (widget.checkState() == QtCore.Qt.PartiallyChecked):
        return None
      elif widget.isChecked():
        return 'True'
      else:
        return 'False'
    else:
      return str(widget.text())
      
  def setEditMode(self,container=False):
    # Expect that if current object is container then disallow editing
    # List and class type
    if container:
      idx = 1
    else:
      idx = 0
    self.containerFrameStack.setCurrentIndex(idx)
    self.nameStack.setCurrentIndex(idx)

  def setColumnGroupWidget(self,currentObject):
    #print 'deEditor.setColumnGroupWidget',currentObject.objectName(),type(currentObject)
    from ..core.CCP4XtalData import CProgramColumnGroup
    if currentObject is None or not isinstance(currentObject,CProgramColumnGroup):
      self.columnGroupWidget.setVisible(False)
      self.columnGroupWidget.blockSignals(True)
    else:
      self.columnGroupWidget.blockSignals(True)
      self.columnGroupWidget.clear()
      data = currentObject.getColumnGroupQualifier()
      #print 'deEditor.setColumnGroupWidget data',data
      self.columnGroupWidget.set(data)
      self.columnGroupWidget.setVisible(True)
      self.columnGroupWidget.blockSignals(False)

  def setNameComboDefault(self,name):
    self.nameCombo.removeItem(0)
    self.nameCombo.insertItem(0,name)
    self.nameCombo.setCurrentIndex(0)
    

  def getName(self):
    if self.nameStack.currentIndex() == 1:
      return str(self.nameCombo.currentText())
    else:
      return str(self.nameEditor.text())

class CCollectionFrame(QtWidgets.QFrame):

  CONTAINER_CLASS_MENU = ['single object','a list (CList)','a dictionary (CDict)']
  CONTAINER_CLASS_NAME = [ None,'CList','CDict']
  
  def __init__(self,parent):
    QtWidgets.QFrame.__init__(self,parent)
    self.setLayout(QtWidgets.QHBoxLayout())
    self.layout().setSpacing(CDefEd.MARGIN)
    self.layout().setContentsMargins(CDefEd.MARGIN,CDefEd.MARGIN,CDefEd.MARGIN,CDefEd.MARGIN)
    self.layout().addWidget(QtWidgets.QLabel('Make the object a',self))
    self.classCombo =CCP4Widgets.CComboBox(self)
    self.classCombo.setMaximumWidth(200)
    for item in CCollectionFrame.CONTAINER_CLASS_MENU:
      self.classCombo.addItem(item,item)
    self.layout().addWidget(self.classCombo )

    self.stack = QtWidgets.QStackedWidget(self)

    dummy = QtWidgets.QLabel('',self)
    self.stack.addWidget(dummy)

    self.listFrame = QtWidgets.QFrame(self)
    self.listFrame.setLayout(QtWidgets.QVBoxLayout())
    listDefaultLayout = QtWidgets.QHBoxLayout()
    listDefaultLayout.addWidget(QtWidgets.QLabel('List default value',self))
    self.listDefault = CCP4Widgets.CLineEdit(self)
    listDefaultLayout.addWidget(self.listDefault)
    self.listFrame.layout().addLayout(listDefaultLayout)
    listQualifiersLayout = QtWidgets.QHBoxLayout()                           
    listQualifiersLayout.addWidget(QtWidgets.QLabel('List of.. min length:'))
    self.minLength = CCP4Widgets.CLineEdit(self)
    validator = QtGui.QIntValidator(self)
    validator.setBottom(0)
    self.minLength.setValidator(validator)
    listQualifiersLayout.addWidget(self.minLength)
    
    listQualifiersLayout.addWidget(QtWidgets.QLabel('max length:',self))
    self.maxLength = CCP4Widgets.CLineEdit(self)
    validator = QtGui.QIntValidator(self)
    validator.setBottom(0)
    self.maxLength.setValidator(validator)
    listQualifiersLayout.addWidget(self.maxLength)

    self.compare = CCP4Widgets.CComboBox(self)
    self.compare.setEditable(False)
    self.compare.addItems(['No comparison','Increasing values','Decreasing values'])
    listQualifiersLayout.addWidget(self.compare)
    self.listFrame.layout().addLayout(listQualifiersLayout)
    self.stack.addWidget(self.listFrame)

    dictDefaultFrame = QtWidgets.QFrame(self)
    dictDefaultFrame.setLayout(QtWidgets.QHBoxLayout())
    dictDefaultFrame.layout().addWidget(QtWidgets.QLabel('Dictionary default',self))
    self.dictDefault =  CCP4Widgets.CLineEdit(self)
    dictDefaultFrame.layout().addWidget(self.dictDefault)
    self.stack.addWidget(dictDefaultFrame)

    tableDefaultFrame = QtWidgets.QFrame(self)
    tableDefaultFrame.setLayout(QtWidgets.QHBoxLayout())
    tableDefaultFrame.layout().addWidget(QtWidgets.QLabel('Table default',self))
    self.tableDefault =  CCP4Widgets.CLineEdit(self)
    tableDefaultFrame.layout().addWidget(self.tableDefault)
    self.stack.addWidget(tableDefaultFrame)
    
    self.layout().addWidget(self.stack)
   
    self.connect(self.classCombo,QtCore.SIGNAL('currentIndexChanged(int)'),self.updateStack)
    self.connect(self.dictDefault,QtCore.SIGNAL('editingFinished()'),functools.partial(self.emit,QtCore.SIGNAL('qualifierEdited'),'default'))
    self.connect(self.listDefault,QtCore.SIGNAL('editingFinished()'),functools.partial(self.emit,QtCore.SIGNAL('qualifierEdited'),'default'))
    self.connect(self.minLength,QtCore.SIGNAL('editingFinished()'),functools.partial(self.emit,QtCore.SIGNAL('qualifierEdited'),'listMinLength'))
    self.connect(self.maxLength,QtCore.SIGNAL('editingFinished()'),functools.partial(self.emit,QtCore.SIGNAL('qualifierEdited'),'listMaxLength'))
    self.connect(self.compare,QtCore.SIGNAL('currentIndexChanged(int)'),functools.partial(self.emit,QtCore.SIGNAL('qualifierEdited'),'listCompare'))

    self.connect(self.minLength,QtCore.SIGNAL('rightMousePress'),functools.partial(self.emit,QtCore.SIGNAL('rightMousePress')))
    self.connect(self.maxLength,QtCore.SIGNAL('rightMousePress'),functools.partial(self.emit,QtCore.SIGNAL('rightMousePress')))
    self.connect(self.compare,QtCore.SIGNAL('rightMousePress'),functools.partial(self.emit,QtCore.SIGNAL('rightMousePress')))


  def updateStack(self,indx):
    #print 'CCollectionFrame.updateStack',indx
    self.stack.setCurrentIndex(indx)
    
  def mousePressEvent(self,event):
    if event.button() == QtCore.Qt.RightButton:
      self.emit(QtCore.SIGNAL('rightMousePress'),event)
      event.accept()
    else:
      event.ignore()

  def set(self,dataObj=None):
    self.blockSignals(True)
    self.classCombo.blockSignals(True)
    if dataObj is not None:
      if isinstance(dataObj,CCP4Data.CList):
        self.setList(dataObj)
      elif isinstance(dataObj,CCP4Data.CDict):
        self.setDict(dataObj)
      else:
        self.setSingle()
    else:
      self.setSingle()
    self.updateStack(self.classCombo.currentIndex())
    self.blockSignals(False)
    self.classCombo.blockSignals(False)
      
  def setSingle(self):
    self.classCombo.setCurrentIndex(0)
    self.minLength.setText('')
    self.maxLength.setText('')
    self.compare.setCurrentIndex(0)
    self.listDefault.setText('')
    self.dictDefault.setText('')

  def setList(self,dataObj):
    self.classCombo.setCurrentIndex(1)
    qualis = dataObj.qualifiers()
    self.minLength.setValue(qualis.get('listMinLength',None))
    self.maxLength.setValue(qualis.get('listMaxLength',None))
    compare = qualis.get('listCompare',None)
    if compare is None:
      self.compare.setCurrentIndex(0)
    elif compare == 1:
      self.compare.setCurrentIndex(1)
    elif compare == -1:
      self.compare.setCurrentIndex(2)
    self.listDefault.setText(str(qualis.get('default')))
    self.dictDefault.setText('')

  def setDict(self,dataObj):
    self.classCombo.setCurrentIndex(2)
    qualis = dataObj.qualifiers()
    self.minLength.setText('')
    self.maxLength.setText('')
    self.compare.setCurrentIndex(0)
    self.listDefault.setText('')
    self.dictDefault.setText(str(qualis.get('default')))

  def getQualifiers(self):
    qualis = {}
    currentIndex = self.classCombo.currentIndex()
    default = None
    if currentIndex == 1:
      t = str(self.minLength.text())
      if len(t)>0:
        qualis['listMinLength'] = int(t)
      else:
        qualis['listMinLength'] = None
      t = str(self.maxLength.text())
      if len(t)>0:
        qualis['listMaxLength'] = int(t)
      else:
        qualis['listMaxLength'] = None
      if qualis['listMinLength'] is not None and  qualis['listMaxLength'] is not None and qualis['listMinLength'] > qualis['listMaxLength']:
        t = qualis['listMinLength']
        qualis['listMinLength'] = qualis['listMaxLength']
        qualis['listMaxLength'] = t
        print('swapped listMinLength/listMaxLength',qualis['listMinLength'],qualis['listMaxLength'])

      ic = self.compare.currentIndex()
      if ic == 0:
        qualis['compare'] = None
      elif ic == 1:
        qualis['compare'] = 1
      elif ic == 2:
        qualis['compare'] = -1
      #print 'CCollectionFrame getQualifiers',qualis
      default =  str(self.listDefault.text())
    elif currentIndex == 2:
      default = str(self.dictDefault.text())

    if default is not None:
      try:
        exec("qualis['default']="+default)
      except:
        qualis['default'] = NotImplemented

    #print 'CCollectionFrame.getQualifiers',qualis
    return qualis

  def getContainerClassName(self):
    return CCollectionFrame.CONTAINER_CLASS_NAME[self.classCombo.currentIndex()]

      
class CClassInfo(QtGui.QTextDocument):
  def __init__(self,className=None,parent=None):
    QtGui.QTextDocument.__init__(self,parent)
    self._className = None
    self.setClassName(className)

    
  def setClassName(self,clsName):
    from ..core.CCP4DataManager import DATAMANAGER
    cls = DATAMANAGER().getClass(clsName)
    if cls is not None:
      self._className=clsName
      self.setPlainText(self.getInfo())
      self.setHtml(self.getInfo())

  def getClass(self):
    from ..core.CCP4DataManager import DATAMANAGER
    return DATAMANAGER().getClass(self._className)
   

  def getClassPathText(self):
    cls = self.getClass()
    if cls is None: return ''
    from ..core.CCP4Data import baseClassList
    baseClassList = baseClassList(cls)
    text = ''
    for clsItem in baseClassList:
      module = str(clsItem.__module__)
      name = clsItem.__name__
      text = text + '<a href="./'+module+'.html#'+name+'">'+module+'.'+name+'</a> -> '
    return text[0:-4]

  def getContentsText(self):
    from ..core import CCP4Data
    cls = self.getClass()
    if cls is None: return ''
    text = '<h4>Contents of class:</h4>\n'
    if issubclass(cls,CCP4Data.CList):
      try:
        className = cls.SUBITEM.get('class').__name__
        moduleName = cls.SUBITEM.get('class').__module__
      except:
        className = ' NOT FOUND'
        moduleName = ' NOT FOUND'
      text = text + 'List of <a href="./'+ moduleName +'.html#'+className+'">' + className + '</a>\n'
    elif issubclass(cls,CCP4Data.CDict):
      try:
        className = cls.SUBITEM.get('class').__name__
        moduleName = cls.SUBITEM.get('class').__module__
      except:
        className = ' NOT FOUND'
        moduleName = ' NOT FOUND'
      text = text + 'CDict of <a href="./'+ moduleName +'.html#' +className+'">' + className + '</a>\n'
    elif len( cls.CONTENTS) == 0:
      text = text + 'This is a simple class without contents\n'
    else:
      text = text + '<table>\n'
      for key,value in list(cls.CONTENTS.items()):
        clsName = value.get('class').__name__
        #text = text + "{0:20}   {1}\n".format(key,value.get('class').__name__)
        className = value.get('class').__name__
        moduleName = value.get('class').__module__
        text = text + '<tr><td>' + key + '</td><td><a href="./'+ moduleName +'.html#' +className+'">' + className + '</a></td></tr>\n'
      text = text + '</table>\n'
    return text

  def getQualifiersInfo(self):
    cls = self.getClass()
    if cls is None: return ''
    try:
      obj = cls()
      qDefs = obj.qualifiersDefinition()
      qOrder = obj.qualifiersOrder()
      qValues = obj.qualifiers()
    except:
      print('Can not get qualifiers info for class')
      return ''
    text = '<h4>Qualifiers for class:</h4>\n<table>\n'
    for quali in qOrder:
      defn = qDefs.get(quali,{})
      desc = defn.get('description','')
      dType = defn.get('type',None)
      if dType is None: dType = obj.PYTHONTYPE
      s = re.search(r"<type '(.*)'>",str(dType))
      if s is not None:
        dType = s.groups()[0]
      else:
        dType = ''
      val = str(qValues.get(quali))
      #text = text + "{0:20}   {1:20}   {2}\n".format(quali,dType,desc)
      text = text + '<tr><td>'+quali+'</td><td><i>'+dType+'</i></td><td><b>'+val+'</b></td><td>'+desc+'</td></tr>\n'
    text = text + '</table>'
    return text
        
      
  def getInfo(self,html=True):
    if self._className is None: return ''
    doc = self.getClass().__doc__
    if html:
      text = '<html>\n'
    else:
      text = ''
    text = text = '<a name="'+self._className+'">'
    if doc is not None:
      text = text + '<h3>' +  self._className+':  ' + doc +'</h3>\n\n'
    else:
      text = text + '<h3>' + self._className  +'</h3>\n'
    text = text + self.getClassPathText() + '\n\n'
    text = text + self.getContentsText() + '\n\n'
    text = text + self.getQualifiersInfo() + '\n'
    if html: text = text + '</html>\n'
    return text


class CHeaderDialog(QtWidgets.QDialog):


  def __init__(self,parent,model=None):
    QtWidgets.QDialog.__init__(self,parent)
    self.setModal(True)
    self.setWindowTitle('DEF file header details')
    self.setLayout(QtWidgets.QVBoxLayout())

    self.editor = CCP4Widgets.CI2XmlHeaderView(self,model=model)
    self.layout().addWidget(self.editor)

    self.done = QtWidgets.QPushButton('OK',self)
    self.done.setAutoDefault(0)
    self.connect(self.done,QtCore.SIGNAL('released()'),self.close)
    doneLayout = QtWidgets.QHBoxLayout()
    doneLayout.addStretch(5)
    doneLayout.addWidget(self.done)
    doneLayout.addStretch(5)
    self.layout().addLayout(doneLayout)

  def close(self):
    self.editor.updateModelFromView()
    self.emit(QtCore.SIGNAL('closed'))
    QtWidgets.QDialog.close(self)

     
class CProgramColumnGroupQualifierView(QtWidgets.QFrame):

  def __init__(self,parent=None):
    QtWidgets.QFrame.__init__(self,parent=None)
    from ..core.CCP4XtalData import CColumnType
    self.setLayout(QtWidgets.QGridLayout())
    maxColumns=5
    iCol = -1
    for label in ['Name','Type','DefaultList','Partner','Offset']:
      iCol = iCol+1
      self.layout().addWidget(QtWidgets.QLabel(label,self),0,iCol)
      
    for iRow in range(1,maxColumns+1):
      for iCol in range(maxColumns):
        if iCol == 1:
          w =  QtWidgets.QComboBox(self)
          for item in CColumnType.QUALIFIERS['enumerators']:
            w.addItem(item)
          w.setCurrentIndex(2)
          self.connect(w,QtCore.SIGNAL('currentIndexChanged(int)'),functools.partial(self.emit,QtCore.SIGNAL('edited'),iRow,iCol))
        else:
          w = CCP4Widgets.CLineEdit(self)
          self.connect(w,QtCore.SIGNAL('editingFinished()'),functools.partial(self.emit,QtCore.SIGNAL('edited'),iRow,iCol))
        w.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.connect(w,QtCore.SIGNAL('rightMousePress'),functools.partial(self.emit,QtCore.SIGNAL('rightMousePress')))
        self.layout().addWidget(w,iRow,iCol)

      
  def mousePressEvent(self,event):
    if event.button() == QtCore.Qt.RightButton:
      self.emit(QtCore.SIGNAL('rightMousePress'),event)
      event.accept()
    else:
      event.ignore()

  def clear(self):
    for iRow in range(1,self.layout().rowCount()):
      for iCol in [0,2,3,4]:
        self.layout().itemAtPosition(iRow,iCol).widget().clear()
      self.layout().itemAtPosition(iRow,1).widget().setCurrentIndex(2)
        
  def get(self):
    colGroupItemList = []
    for iRow in range(1,self.layout().rowCount()):
      if len(str(self.layout().itemAtPosition(iRow,0).widget().text()).strip())>0:
        colGroupItem = {}
        colGroupItem['columnName'] = (str(self.layout().itemAtPosition(iRow,0).widget().text())).strip()
        colGroupItem['columnType'] = str(self.layout().itemAtPosition(iRow,1).widget().currentText())
        colGroupItem['defaultList'] = str(self.layout().itemAtPosition(iRow,2).widget().text())
        colGroupItem['partnerTo'] = str(self.layout().itemAtPosition(iRow,3).widget().text())
        if len(colGroupItem['partnerTo'].strip())==0: colGroupItem['partnerTo'] = None
        offset = str(self.layout().itemAtPosition(iRow,4).widget().text())
        if len(offset) == 0:    
          colGroupItem['partnerOffset'] = None
        else:
          colGroupItem['partnerOffset'] = int(offset)
        colGroupItemList.append(colGroupItem)
    return colGroupItemList

  def set(self,colGroupItemList=[]):
    from ..core.CCP4XtalData import CColumnType
    self.clear()
    iRow = 0
    columnTypes = CColumnType.QUALIFIERS['enumerators']
    for colGroupItem in colGroupItemList:
      #print 'CProgramColumnGroupQualifierView.set',colGroupItem
      iRow = iRow + 1
      if colGroupItem.get('columnName',None) is not None:
        self.layout().itemAtPosition(iRow,0).widget().setText(colGroupItem['columnName'])
      if colGroupItem.get('columnType',None) is not None:
        if isinstance(colGroupItem['columnType'],(list,CCP4Data.CList)):
          colType = str(colGroupItem['columnType'][0])
        else:
          colType = str(colGroupItem['columnType'])
        #print 'CProgramColumnGroupQualifierView.set',colType,type(colType)
        if columnTypes.count(colType):
          iColType =  columnTypes.index(colType)
          self.layout().itemAtPosition(iRow,1).widget().setCurrentIndex(iColType)
      if colGroupItem.get('defaultList',None) is not None:
        self.layout().itemAtPosition(iRow,2).widget().setText(colGroupItem['defaultList'])
      if colGroupItem.get('partnerTo',None) is not None:
        self.layout().itemAtPosition(iRow,3).widget().setText(colGroupItem['partnerTo'])
      if colGroupItem.get('partnerOffset',None) is not None:
        self.layout().itemAtPosition(iRow,4).widget().setText(str(colGroupItem['partnerOffset']))
      
    
class CMakeWrapperPlugin(QtWidgets.QDialog):

  def __init__(self,parent=None):
    QtWidgets.QDialog.__init__(self,parent)
    from ..core.CCP4ScriptManager import SCRIPTMANAGER

    self.setModal(True)
    self.setWindowTitle('Create Plugin')
    self.setLayout(QtWidgets.QVBoxLayout())

    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('See help text for more information',self))
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Name of plugin',self))
    self.name = QtWidgets.QLineEdit(self)
    line.addWidget(self.name)
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Make in wrappers directory or a pipeline',self))
    self.pipeline = QtWidgets.QComboBox(self)
    self.pipeline.addItem('Wrappers')
    for item in SCRIPTMANAGER().listOfPipelines(): self.pipeline.addItem(item)
    #self.pipeline.addItem('Other..')
    line.addWidget(self.pipeline)
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Use plugin template',self))
    self.template =  QtWidgets.QComboBox(self)
    for item in SCRIPTMANAGER().listOfPluginTemplates(): self.template.addItem(item)
    line.addWidget(self.template)
    self.layout().addLayout(line)
    
    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Use license',self))
    self.license =  QtWidgets.QComboBox(self)
    self.license.addItem('No license')
    for item in SCRIPTMANAGER().listOfPluginLicenses(): self.license.addItem(item)
    line.addWidget(self.license)
    self.layout().addLayout(line)
    
    buttonBox = QtWidgets.QDialogButtonBox(self)
    but = buttonBox.addButton(QtWidgets.QDialogButtonBox.Apply)
    self.connect(but,QtCore.SIGNAL('clicked()'),self.makeWrapper)
    but = buttonBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
    self.connect(but,QtCore.SIGNAL('clicked()'),self.close)
    self.layout().addWidget(buttonBox)

  def makeWrapper(self):
    from ..core.CCP4ScriptManager import SCRIPTMANAGER
    name = str(self.name.text()).strip()
    pipeline = str(self.pipeline.currentText())
    if pipeline == 'Wrappers' : pipeline = None
    template =  str(self.template.currentText())
    license =  str(self.license.currentText())
    if license == 'No license': license = None
    
    if len(name)==0:    
      QtWidgets.QMessageBox.warning(self,'No wrapper name','You must provide a wrapper name')
      return

    wrapperDir = None
    try:
      wrapperDir,report = SCRIPTMANAGER().makeWrapperPlugin(name=name,pipeline=pipeline,template=template,license=license)
    except CException as e:
      warningMessage(e)
    except:
       QtWidgets.QMessageBox.warning(self,'Undefined error','Undefined error trying to create wrapper directory')
    else:
      if len(report)>0:
        warningMessage(report)
      elif wrapperDir is not None:
        QtWidgets.QMessageBox.information(self,'New wrapper plugin','New wrapper plugin created in directory '+wrapperDir)


class CMakePipeline(QtWidgets.QDialog):

  def __init__(self,parent=None):
    QtWidgets.QDialog.__init__(self,parent)

    self.setModal(True)
    self.setWindowTitle('Create Pipeline')
    self.setLayout(QtWidgets.QVBoxLayout())

    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('See help text for more information',self))
    self.layout().addLayout(line)

    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Name of pipeline',self))
    self.name = QtWidgets.QLineEdit(self)
    line.addWidget(self.name)
    self.layout().addLayout(line)
    
    buttonBox = QtWidgets.QDialogButtonBox(self)
    but = buttonBox.addButton(QtWidgets.QDialogButtonBox.Apply)
    self.connect(but,QtCore.SIGNAL('clicked()'),self.makePipeline)
    but = buttonBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
    self.connect(but,QtCore.SIGNAL('clicked()'),self.close)
    self.layout().addWidget(buttonBox)


  def makePipeline(self):
    from ..core.CCP4ScriptManager import SCRIPTMANAGER
    name = str(self.name.text()).strip()
    
    if len(name)==0:    
      QtWidgets.QMessageBox.warning(self,'No pipeline name','You must provide a pipeline name')
      return

    pipelineDir = None
    try:
      pipelineDir = SCRIPTMANAGER().makePipeline(name=name)
    except CException as e:
      warningMessage(e)
    #except:
    #   QtWidgets.QMessageBox.warning(self,'Undefined error','Undefined error trying to create pipeline directory')

    if pipelineDir is not None:
      QtWidgets.QMessageBox.information(self,'New pipeline','New pipeline created in directory '+pipelineDir)
