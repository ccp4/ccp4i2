from __future__ import print_function

"""
     CCP4AbstractViewer.py: CCP4 GUI Project
     Copyright (C) 2009-2010 University of York

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
     Liz Potterton Jan 2010 - Create CCP4AbstractViewer
"""

## @package CCP4AbstractViewer (QtGui) Base for file viewers in CCP4WebBrowser
import os
from PySide2 import QtGui, QtWidgets,QtCore
from core import CCP4Modules
from core.CCP4ErrorHandling import *

def handleFileChanged(fileName):
    fileName = str(fileName)
    #print 'CCP4AbstractViewer.handleFileChanged',fileName
    indx = 0
    browser = CCP4Modules.WEBBROWSER(indx)
    while browser is not None:
        tab = browser.fileOpenInTab(fileName)
        #print 'CCP4AbstractViewer.handleFileChanged tab',tab
        if tab >= 0:
            #print 'CCP4AbstractViewer.handleFileChanged reloading',fileName
            browser.tab().widget(tab).reload()
            return
        indx = indx + 1
        browser = CCP4Modules.WEBBROWSER(indx)


#-------------------------------------------------------------------
class CAbstractViewer(QtWidgets.QScrollArea):
#-------------------------------------------------------------------
    MENUTEXT = 'Viewer'
    FILEWATCHER = None
    ERROR_CODES = {1 : {'severity' : SEVERITY_ERROR, 'description' : 'Failed opening file'},
                   2 : {'severity' : SEVERITY_ERROR, 'description' : 'Failed reading file'},
                   3 : {'severity' : SEVERITY_ERROR, 'description' : 'No file name or file does not exist'},}
# Subclassed to display various file types in the CCP4WebBrowser
#-------------------------------------------------------------------
    def __init__(self,parent=None,fileName=None):
#-------------------------------------------------------------------
        #QtWidgets.QWidget.__init__(self,parent)
        QtWidgets.QScrollArea.__init__(self,parent)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        #self.setWidgetResizable(1)
        self.fileName = None
        self.lastModTime = None
        if fileName is not None: self.open(fileName)

#-------------------------------------------------------------------
    def open(self,fileName):
#-------------------------------------------------------------------
        # Expect this to be reimplemented in sub-class
        #print 'CAbstractViewer.open',fileName
        self.fileName=os.path.abspath(fileName)
        self.lastModTime = os.path.getmtime(fileName)
        if self.fileName is not None: 
            self.setObjectName(self.title())

    def isFileModified(self):
        #print 'CAbstractViewer.isFileModified', self.fileName, self.lastModTime
        if self.fileName is None or  self.lastModTime is None: return False
        modTime =  os.path.getmtime(self.fileName)
        if modTime > self.lastModTime:
            return True
        else:
            return False

    def reload(self):
        if self.fileName is not None and os.path.exists(self.fileName):
            self.open(self.fileName)

    def title(self):
        if self.fileName is not None:
            root,base = os.path.split(self.fileName)
            if len(root) > 0:
                root, jobDir = os.path.split(root)
                if jobDir[0:4]=='job_': base = jobDir+'/'+base
            return base
        else:
            return ''

    def watchFile(self,modTime=180.0):
        import time
        if self.fileName is None: return
        # Only watch files last modified in last three minutes
        if time.time() - os.path.getmtime(self.fileName) > modTime: return
        CCP4Modules.FILEWATCHER().addPath(self.fileName)

    def close(self):
        if self.fileName is None: return
        #print 'CAbstractViewer.close removing watch',self.fileName
        try:
            CCP4Modules.FILEWATCHER().removePath(self.fileName)
        except:
            pass
  
    def Size(self):
        return QtCore.QSize()

#-------------------------------------------------------------------
    def mgSize(self):
#-------------------------------------------------------------------
        return QtCore.QSize()

#-------------------------------------------------------------------
    def Print(self,painter):
#-------------------------------------------------------------------
        pass

#-------------------------------------------------------------------
    def Save(self,fileName):
#-------------------------------------------------------------------
        pass

#-------------------------------------------------------------------
    def isPrintable(self):
#-------------------------------------------------------------------
        return 0

#-------------------------------------------------------------------
    def isSaveable(self):
#-------------------------------------------------------------------
        return 0
  
#-------------------------------------------------------------------
    def isSearchable(self):
#-------------------------------------------------------------------
        return 0
  
#-------------------------------------------------------------------
    def isScaleable(self):
#-------------------------------------------------------------------
        return 0

    def isRunable(self):
        return 0

#-------------------------------------------------------------------
    def getFileExt(self):
#-------------------------------------------------------------------
#FIXME
        return None
        #print 'getFileExt',self.__class__
        from qtcore import CCP4CustomMimeTypes
        mimeHandler = CCP4CustomMimeTypes.MimeTypesHandler()
        
        mimetype = mimeHandler.getMimeTypeForViewer(self.__class__)
        if  mimetype:
            exts =  mimeHandler.getMimeTypeInfo(mimetype,'fileExtensions')
        else:
            exts = []
        #print 'getFileExt',exts
        return exts

#-------------------------------------------------------------------
    def getLabel(self):
#-------------------------------------------------------------------
        from qtcore import CCP4CustomMimeTypes
        mimeHandler = CCP4CustomMimeTypes.MimeTypesHandler()
            
        mimetype = mimeHandler.getMimeTypeForViewer(self.__class__)
        if  mimetype:
            label =  mimeHandler.getMimeTypeInfo(mimetype,'description')
        else:
            label = ''
        #print 'getLabel',label
        return label

    def browserWindow(self):
        return self.parent()


    def handleTabbedOpen(self):
        pass

    def handleTabbedClosed(self):
        pass

class CFileWatchTimer(QtCore.QObject):

    doHandleTimeout = QtCore.Signal()
    fileChanged = QtCore.Signal(str)

    insts = None

    def __init__(self,interval=2000):
        parent = CCP4Modules.QTAPPLICATION()
        QtCore.QObject.__init__(self,parent)
        parent.aboutToQuit.connect(self.Exit)
        self.files = {}
        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(interval)
        self.timer.timeout.connect(self._handleTimeout)
        self.doHandleTimeout.connect(self.handleTimeout)
        self.timer.start()
        self.blockExit = False

    @QtCore.Slot()
    def Exit(self):
        import sys
        self.timer.stop()
        sys.__stdout__.write('CFileWatcher.Exit blockExit'+str(self.blockExit)+'\n');sys.__stdout__.flush()

    def addPath(self,fileName):
        from core import CCP4Utils
        for item in list(self.files.keys()):
            if CCP4Utils.samefile(fileName,item,default=False): return False
        self.files[fileName] = os.path.getmtime(fileName)
        return True

    def removePath(self,fileName):
        from core import CCP4Utils
        if fileName in self.files:
            del self.files[fileName]
            return True
        for item in list(self.files.keys()):
            if CCP4Utils.samefile(item,fileName,default=False):
                del self.files[item]
                return True
        return False

    @QtCore.Slot()
    def _handleTimeout(self):
        self.doHandleTimeout.emit()

    @QtCore.Slot()
    def handleTimeout(self):
        self.blockExit = True
        delList = []
        for fileName in list(self.files.keys()):
            try:
                mtime = os.path.getmtime(fileName)
                if mtime > self.files[fileName]:
                    #print 'CFileWatchTimer.handleTimeout emiting',fileName
                    self.fileChanged.emit(fileName)
                    self.files[fileName] = mtime
            except:
                delList.append(fileName)
        # Remove any that screwed up so don't do it again 
        for fileName in delList:
            del self.files[fileName]
        self.blockExit = False
   
