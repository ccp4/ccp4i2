from __future__ import print_function

"""
     qtgui/CCP4ReportWidgets.py: CCP4 Gui Project
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
'''
     Liz Potterton Mar 2010 - Prototype viewer for data expected to be in reports
'''
from PySide6 import QtGui, QtWidgets,QtCore
from core.CCP4ErrorHandling import *

##@package CCP4ReportWidgets (QtGui) QtWebKit plugins for CCP4 Reports
ICONBUTTONSIZE=24
DRAGICONSIZE=16


def loadTypeDefinitions(typeDefinitions=None):
  if not typeDefinitions: return
  for dataType,defn in list({
    'data_table' : [CReportTable,{}] }.items()):
    try:
      typeDefinitions.addType(dataType=dataType,widgetClass=defn[0],widgetQualifiers=defn[1])
    except CException as e:
      e.report()
    except:
      pass


class CReportTableView(QtWidgets.QTableView):

    def __init__(self,parent=None,model=None,qualifiers={}):
        QtWidgets.QTableView.__init__(self,parent)
        self.horizontalHeader().setVisible(1)
        if model is not None: self.setModel(model)

    def setModel(self,table=None):
        #print 'CReportTable.setModel',table
        QtWidgets.QTableView.setModel(self,table)
        self.resizeColumnsToContents()


class CLauncherButton(QtWidgets.QPushButton):
  def __init__(self,parent=None,dataEtreeList=[],**kw):
    QtWidgets.QPushButton.__init__(self,parent)
    #print 'CLauncherButton',dataEtreeList,kw
    self.dataEtreeList = dataEtreeList
    self.lastClickTime = None
    self.argList = {}
    self.argList.update(kw)
    self.setText(kw.get('label','launch'))
    iconName = kw.get('icon')
    if iconName is not None:
      import os
      from core import CCP4Utils
      icon_path = os.path.join(CCP4Utils.getCCP4I2Dir(),'qticons',iconName)
      self.setIcon(QtGui.QIcon(icon_path))
    toolTip = kw.get('tooltip')
    if toolTip is not None:
      self.setToolTip(toolTip)
    self.clicked.connect(self.handleClick)

  @QtCore.Slot()
  def handleClick(self,checked):
    #print 'CLauncherButton.handleClick',self.argList
    # Beware the unintended double click
    import time,os
    if self.lastClickTime is None:
      self.lastClickTime = time.time()
    else:
      dif = time.time() - self.lastClickTime
      self.lastClickTime = time.time()
      if dif < 1.0: return
    from core import CCP4Modules
    if self.argList.get('exe') == 'CCP4mg':
      if self.argList.get('sceneFile',None) is not None:
        sceneFile = os.path.join(CCP4Modules.PROJECTSMANAGER().jobDirectory(self.argList.get('jobId')),os.path.split(self.argList['sceneFile'])[-1])
        CCP4Modules.LAUNCHER().launch(viewer='CCP4mg',argList=['-nore',sceneFile],projectId=self.argList.get('projectId',None),guiParent=self.parent())
      elif self.argList.get('pictureDefinitionFile',None) is not None:
        CCP4Modules.LAUNCHER().launch(viewer='CCP4mg',argList=['-nore','-pict',self.argList['pictureDefinitionFile']],projectId=self.argList.get('projectId',None),guiParent=self.parent())
      elif self.argList.get('jobId',None) is not None:
        CCP4Modules.LAUNCHER().openInViewer(viewer='CCP4mg',jobId=self.argList.get('jobId'),projectId=self.argList.get('projectId',None),guiParent=self.parent())
    
    elif self.argList.get('exe') == 'Coot':
      #print '\n\nArglist',self.argList
      CCP4Modules.LAUNCHER().openInViewer(viewer='coot_job',jobId=self.argList.get('jobId',None),projectId=self.argList.get('projectId',None),guiParent=self.parent())
      
    elif self.argList.get('exe',None) == 'loggraph':
      plotCurrentIndex = -1
      if hasattr(self.parent(),"page"):

          idString = ""
          if len(self.dataEtreeList)>0 and self.argList.get('fileName',None) is not None:
              for item in self.dataEtreeList:
                if item.get('id',None) is not None:
                    idString += item.get('id')+","

          plotCurrentIndexAsQVariant = self.parent().page().mainFrame().evaluateJavaScript('''
              function getGraphPlots(){
                  var graphElements = $('div[data-renderer="CCP4i2Widgets.CCP4FlotRenderer"]');
                  var idCurrentIdxDict = {};
                  for(var ig=0;ig<graphElements.length;ig++) {
                      try {
                          var buttonDiv = graphElements[ig].querySelector('div[id="theButtonDiv"]');
                          var launchObj = buttonDiv.querySelector('object[type="x-ccp4-widget/CLauncherButton"]')
                          var idParam = launchObj.querySelector('param[name="ccp4_data_id"]')
                          var currentIndexParam = launchObj.querySelector('param[name="ccp4_data_current_index"]')
                          var theID = idParam.getAttribute("value");
                          var currentIndex = currentIndexParam.getAttribute("value");
                          idCurrentIdxDict[theID] = currentIndex;
                      } catch(err) {
                      }
                  }

                  return idCurrentIdxDict;
              }
              function getPlotsCurrentIndex(){
                  var idIdxs = getGraphPlots();
                  var theId = "'''+idString+'''";
                  if(theId in idIdxs){
                      if(idIdxs.hasOwnProperty(theId)){
                          return idIdxs[theId];
                      }
                  }
                  return -1;
              }
              getPlotsCurrentIndex();
          ''')
#FIXME PYQT - toInt removed.
          plotCurrentIndexVal = plotCurrentIndexAsQVariant
          plotCurrentIndex = int(plotCurrentIndexVal)
      if len(self.dataEtreeList)>0 and self.argList.get('fileName',None) is not None:
        launchArgs = []
        for item in self.dataEtreeList:
          if item.get('id',None) is not None:
            launchArgs.extend(['--select',item.get('id')])
        if plotCurrentIndex > -1:
            launchArgs.extend(["-G",str(plotCurrentIndex)])
        launchArgs.append(self.argList['fileName'])
        CCP4Modules.LAUNCHER().launch(viewer='loggraph',argList=launchArgs)
      else:
        CCP4Modules.LAUNCHER().openInViewer(viewer='loggraph',fileName=self.argList.get('fileName',None),
                          jobId=self.argList.get('jobId',None),projectId=self.argList.get('projectId',None))

        
class CDownloadButton(QtWidgets.QPushButton):
  def __init__(self,parent=None,dataEtreeList=[],**kw):
    QtWidgets.QPushButton.__init__(self,parent)
    #print 'CDownloadButton',parent,kw
    self.dataName = kw.get('dataName',None)
    self.jobId = kw.get('jobId',None)
    self.lastClickTime = None
    self.setText(kw.get('label','Download'))
    self.setStyleSheet("QPushButton { font-size:8px; }")
    self.setMaximumWidth(70)
    self.setToolTip('Save table as a csv (comma separated values) file for import to Word, Excel etc')
    self.clicked.connect(self.handleClick)

  def handleClick(self,checked):
    # Beware the unintended double click
    import time,os
    from core import CCP4Modules
    if self.lastClickTime is None:
      self.lastClickTime = time.time()
    else:
      dif = time.time() - self.lastClickTime
      self.lastClickTime = time.time()
      if dif < 1.0: return
    self.parent().emitDownloadRequest(self.jobId,self.dataName)

      
class CLaunchTaskButton(QtWidgets.QPushButton):
  def __init__(self,parent=None,dataEtree=None,**kw):
    QtWidgets.QPushButton.__init__(self,parent)
    #print 'CTaskButton',kw
    self.dataEtree = dataEtree
    self.lastClickTime = None
    self.argList = {}
    self.argList.update(kw)
    self.setText(kw.get('label'))
    self.taskName = kw.get('taskName')
    self.clicked.connect(self.handleClick)

  @QtCore.Slot()
  def handleClick(self):
    #print 'CTaskButton.handleClick',self.parent(),self.argList,self.lastClickTime
    # Beware the unintended double click
    import time
    if self.lastClickTime is None:
      self.lastClickTime = time.time()
    else:
      dif = time.time() - self.lastClickTime
      self.lastClickTime = time.time()
      if dif < 1.0: return
    jobInfo =  self.parent().report.getJobInfo()
    jobId = jobInfo.get('jobId',None)
    if jobId is None: return
    #print 'CTaskButton.handleClick jobId',jobId
    from qtgui import CCP4ProjectViewer
    viewer = CCP4ProjectViewer.PROJECTVIEWER(projectId=jobInfo['projectId'],open=True)
    if viewer is not None:
      viewer.openTask(taskName=self.taskName,followJobId=jobId)
    
class CSubJobButton(QtWidgets.QPushButton):
  def __init__(self,parent=None,dataEtree=None,**kw):
    QtWidgets.QPushButton.__init__(self,parent)
    #print 'CLauncherButton',kw
    self.dataEtree = dataEtree
    self.argList = {}
    self.argList.update(kw)
    self.setText(kw.get('label','Sub job'))
    self.setMinimumWidth(600)
    self.setMinimumHeight(30)
    self.setMaximumHeight(30)
    self.clicked.connect(self.handleClick)

  @QtCore.Slot()
  def handleClick(self):
    print('CSubJobButton.handleClick')
    
class CDisplayObjectIcon(QtWidgets.QToolButton):

  leftMousePress = QtCore.Signal('QMouseEvent')
  leftMouseRelease = QtCore.Signal('QMouseEvent')
  rightMouseRelease = QtCore.Signal('QMouseEvent')

  def __init__(self,parent=None,dataEtree=None,**kw):
    QtWidgets.QToolButton.__init__(self,parent)
    #print 'CDisplayObjectIcon',kw
    self.dataEtree = dataEtree
    self.jobId = kw.get('jobId',None)
    self.iconMenu = QtWidgets.QMenu(self)
    self.iconMenu.triggered.connect(self.handleMenu)

  def sizeHint(self):
    return QtCore.QSize(ICONBUTTONSIZE,ICONBUTTONSIZE)

  def mousePressEvent(self,event):
        #print 'CIconButton.mousePressEvent',event.button(),self.dragType()
        if event.button() == QtCore.Qt.LeftButton:
          self.leftMousePress.emit(event)
          #if self.dragType() is not None:
          #  self.dragStartPosition = event.pos()
          #  event.accept()
        elif event.button() == QtCore.Qt.RightButton:
            self.iconMenu.popup(event.globalPos())
            QtWidgets.QToolButton.mousePressEvent(self,event)
        return
      
  def mouseReleaseEvent(self,event):
        #print 'CIconButton.mouseReleaseEvent'
        if event.button() == QtCore.Qt.RightButton:
            self.rightMouseRelease.emit(event)
        if event.button() == QtCore.Qt.LeftButton:
            self.leftMouseRelease.emit(event)
        QtWidgets.QToolButton.mouseReleaseEvent(self,event)

  '''
  def mouseMoveEvent(self,event):
        #print 'CIconButton.mouseMoveEvent'
        #if event.button() != QtCore.Qt.LeftButton:
        #    print 'wrong button',event.button(),QtCore.Qt.LeftButton
        #   return
        if self.dragType() is None: return
        if self.dragStartPosition is not None and \
          (event.pos() - self.dragStartPosition).manhattanLength() < QTAPPLICATION().startDragDistance():
            return
        self.parent().startDrag()
  '''

  def handleMenu(self,action):
    pass


class CReferenceIcon(CDisplayObjectIcon):
  def __init__(self,parent=None,dataEtree=None,**kw):
    CDisplayObjectIcon.__init__(self,parent,dataEtree=dataEtree,**kw)
    #print 'CReferenceIcon',kw
    self.taskName = kw.get('taskName',None)
    from qtgui import CCP4GuiUtils
    icon =CCP4GuiUtils.createIcon(name='book')
    self.setIcon(icon)
    #action = self.iconMenu.addAction('View')
    import functools
    for sele in ['all','preferred']:
      subMenu = QtWidgets.QMenu('Download '+sele+' references',self)
      subMenu.triggered.connect(functools.partial(self.handleDownloadMenu,sele))
      self.iconMenu.addMenu(subMenu)
      action = subMenu.addAction('PubMed file for EndNote')
      action = subMenu.addAction('BibTeX file for LaTex')
    self.iconMenu.addAction('Help')
    self.fileBrowser = None

  @QtCore.Slot('QAction')
  def handleMenu(self,action):
    label = str(action.text())
    if label == 'Help':
      import os
      from core import CCP4Utils,CCP4Modules
      CCP4Modules.WEBBROWSER().openFile(os.path.join(CCP4Utils.getCCP4I2Dir(),'docs','general','references.html'))

  @QtCore.Slot(str,'QAction')
  def handleDownloadMenu(self,selection,action):
    if self.fileBrowser is not None: return
    import functools
    from core import CCP4TaskManager
    label = str(action.text())
    #print 'CReferenceIcon.handleMenu',label
    format = None
    if label.count('PubMed'):
      format = 'medline'
    elif label.count('BibTeX'):
      format = 'bibtex'
    if format is not None:
      fileNameList = CCP4TaskManager.TASKMANAGER().searchReferenceFile(self.taskName,cformat=format,drillDown=True)
      if len(fileNameList)==0: return
      title = CCP4TaskManager.TASKMANAGER().getTitle(self.taskName)
      from qtgui import CCP4FileBrowser
      fileBrowser = CCP4FileBrowser.CFileDialog(parent=self,
           title='Save '+selection+' bibliographic references for '+title,
           filters= ['Reference file (*.'+format+'.txt)'],
           defaultSuffix=format+'.txt',
           fileMode=QtWidgets.QFileDialog.AnyFile  )
      fileBrowser.selectFile.connect(functools.partial(self.downloadReferences,selection,format,fileNameList[0]))
      fileBrowser.show()

    
  @QtCore.Slot(str,str,str,str)
  def downloadReferences(self,selection,format,sourceFile,targetFile):
    if selection == 'all':
      import shutil  
      try:
        shutil.copyfile(sourceFile,targetFile)
      except:
        QtWidgets.QMessageBox.warning(self,'Downloading references file','Failed to copy reference file to '+str(targetFile))
    elif selection == 'preferred':
      from core import CCP4Utils
      text = CCP4Utils.readFile(sourceFile)
      if format == 'medline':
        firstRef = '\nPMID- ' + text.split('\nPMID- ')[1]
      elif format == 'bibtex':
        firstRef = '\n%' + text.split('\n%')[1]
      CCP4Utils.saveFile(targetFile,firstRef)
