"""
     CCP4Modules.py: CCP4 GUI Project
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

   Liz Potterton Jan 2010 - 
"""

import sys
import time

from . import CCP4ComFilePatchManager
from . import CCP4CustomTaskManager
from . import CCP4ImportedJobManager
from . import CCP4Preferences
from . import CCP4PrintHandler
from . import CCP4ProcessManager
from . import CCP4ProjectsManager
from . import CCP4TaskManager
from . import CCP4WorkflowManager
from ..qtcore import CCP4CustomMimeTypes
from ..qtcore import CCP4DbThread
from ..qtcore import CCP4HTTPServerThread
from ..qtcore import CCP4JobController
from ..qtcore import CCP4Launcher
from ..qtgui import CCP4AbstractViewer
from ..qtgui import CCP4DemoData
from ..qtgui import CCP4JobControlGui
from ..qtgui import CCP4WebBrowser
from ..qtgui import CCP4Widgets
from ..utils import QApp
from .CCP4Config import GRAPHICAL


TIMING=False

def QTAPPLICATION(graphical=None):
    # NB can not use QApplication.instance() as it returns the QApplication object
    # rather than CApplication. Suppose could reimplement CApplication.instance() but
    # would just do the same as this function.
    if graphical is None:
        graphical = GRAPHICAL()
    if QApp.MYAPPLICATION is None:
        #print 'Starting Qt, graphical mode:',graphical
        if graphical:
            QApp.MYAPPLICATION = QApp.CGuiApplication(sys.argv)
        else:
            QApp.MYAPPLICATION = QApp.CApplication(sys.argv)
    return QApp.MYAPPLICATION


def PROJECTSMANAGER():
    if CCP4ProjectsManager.CProjectsManager.insts is None:
        #print 'CCP4Modules.PROJECTSMANAGER setting new'
        #traceback.print_stack()
        if TIMING:t1=time.time()
        CCP4ProjectsManager.CProjectsManager.insts= CCP4ProjectsManager.CProjectsManager()
        if TIMING:
            t2=time.time()
            print("Starting 1",CCP4ProjectsManager.CProjectsManager.insts, t2-t1)
    return CCP4ProjectsManager.CProjectsManager.insts

def PROCESSMANAGER():
    if CCP4ProcessManager.CProcessManager.insts is None:
        if TIMING:t1=time.time()
        t = CCP4ProcessManager.CProcessManager()
        if TIMING:
            t2=time.time()
            print("Starting 2",t.__class__, t2-t1)
    return CCP4ProcessManager.CProcessManager.insts

def LAUNCHER():
    if CCP4Launcher.CLauncher.insts is None:
        parent = QTAPPLICATION()
        if TIMING:t1=time.time()
        t = CCP4Launcher.CLauncher(parent)
        if TIMING:
            t2=time.time()
            print("Starting 3",t.__class__, t2-t1)
    return CCP4Launcher.CLauncher.insts

# Modules below here are graphical

def PIXMAPMANAGER():
    if CCP4Widgets.CPixmapManager.insts is None:
        if TIMING:t1=time.time()
        t = CCP4Widgets.CPixmapManager()
        if TIMING:
            t2=time.time()
            print("Starting 4",t.__class__, t2-t1)
    return CCP4Widgets.CPixmapManager.insts

def WEBBROWSER(index = -1,new=False,mini=False):
    if new or (len(CCP4WebBrowser.CBrowserWindow.Instances)==0 and index<0):
        if TIMING:t1=time.time()
        t = CCP4WebBrowser.CBrowserWindow(welcome=False)
        if TIMING:
            t2=time.time()
            print("Starting 5",t.__class__, t2-t1)
        if mini: t.setMini(True)
        t.show()
        t.raise_()
        return t
    index = min(index,0)
    if index<len(CCP4WebBrowser.CBrowserWindow.Instances):
        win = list(CCP4WebBrowser.CBrowserWindow.Instances)[index]
        print('WEBBROWSER raise_')
        try:
            win.show()
            win.raise_()
        except:
            pass
        return win
    else:
        return None

def DUMMYMAINWINDOW():
    if CCP4WebBrowser.CBrowserWindow.Dummy is None:
        CCP4WebBrowser.CBrowserWindow.Dummy = CCP4WebBrowser.CBrowserWindow()
        CCP4WebBrowser.CBrowserWindow.Dummy.hide()  # KJS: Looks like there is no hide present. Check to see if this is used anywhere.
        try:
            CCP4WebBrowser.CBrowserWindow.Instances.remove(CCP4WebBrowser.CBrowserWindow.Dummy)
        except:
            print('Error in CCP4Modules.DUMMYMAINWINDOW')
    return CCP4WebBrowser.CBrowserWindow.Dummy

def JOBCONTROLLER():
    if not CCP4JobController.CJobController.insts:
        if TIMING: t1=time.time()
        t=CCP4JobController.CJobController()
        if TIMING: 
            t2=time.time()
            print("Starting 6",t.__class__, t2-t1)
    return CCP4JobController.CJobController.insts

def MIMETYPESHANDLER():
#FIXME
    if not CCP4CustomMimeTypes.CCustomMimeTypes.insts:
        if TIMING: t1=time.time()
        t = CCP4CustomMimeTypes.CCustomMimeTypes()
        if TIMING: 
            t2=time.time()
            print("Starting 7",t.__class__, t2-t1)
    return CCP4CustomMimeTypes.CCustomMimeTypes.insts

def TASKMANAGER():
    if CCP4TaskManager.CTaskManager.insts is None:
        if TIMING: t1=time.time()
        CCP4TaskManager.CTaskManager.insts = CCP4TaskManager.CTaskManager()
        if TIMING: 
            t2=time.time()
            print("Starting 8",CCP4TaskManager.CTaskManager.insts.__class__, t2-t1)
    return CCP4TaskManager.CTaskManager.insts

def WORKFLOWMANAGER():
    if CCP4WorkflowManager.CWorkflowManager.insts is None:
        if TIMING: t1=time.time()
        CCP4WorkflowManager.CWorkflowManager.insts = CCP4WorkflowManager.CWorkflowManager()
        if TIMING: 
            t2=time.time()
            print("Starting 9",CCP4WorkflowManager.CWorkflowManager.insts.__class__, t2-t1)
    return CCP4WorkflowManager.CWorkflowManager.insts

def COMFILEPATCHMANAGER():
    if CCP4ComFilePatchManager.CComFilePatchManager.insts is None:
        if TIMING: t1=time.time()
        CCP4ComFilePatchManager.CComFilePatchManager.insts = CCP4ComFilePatchManager.CComFilePatchManager()
        if TIMING: 
            t2=time.time()
            print("Starting 10",CCP4ComFilePatchManager.CComFilePatchManager.insts.__class__, t2-t1)
    return CCP4ComFilePatchManager.CComFilePatchManager.insts

def CUSTOMTASKMANAGER():
    if CCP4CustomTaskManager.CCustomTaskManager.insts is None:
        if TIMING: t1=time.time()
        CCP4CustomTaskManager.CCustomTaskManager.insts = CCP4CustomTaskManager.CCustomTaskManager()
        if TIMING: 
            t2=time.time()
            print("Starting 11",CCP4CustomTaskManager.CCustomTaskManager.insts.__class__, t2-t1)
    return CCP4CustomTaskManager.CCustomTaskManager.insts

def IMPORTEDJOBMANAGER():
    if CCP4ImportedJobManager.CImportedJobManager.insts is None:
        if TIMING: t1=time.time()
        CCP4ImportedJobManager.CImportedJobManager.insts = CCP4ImportedJobManager.CImportedJobManager()
        if TIMING: 
            t2=time.time()
            print("Starting 12",CCP4ImportedJobManager.CImportedJobManager.insts.__class__, t2-t1)
    return CCP4ImportedJobManager.CImportedJobManager.insts

def PREFERENCES():
    if CCP4Preferences.CPreferences.insts is None:
        if TIMING: t1=time.time()
        p = CCP4Preferences.CPreferences()
        if TIMING: 
            t2=time.time()
            print("Starting 13",p.__class__, t2-t1)
    return CCP4Preferences.CPreferences.insts

  
def PRINTHANDLER():
    if CCP4PrintHandler.CPrintHandler.insts is None:
        if TIMING: t1=time.time()
        obj = CCP4PrintHandler.CPrintHandler()
        if TIMING: 
            t2=time.time()
            print("Starting 14",obj.__class__, t2-t1)
    return CCP4PrintHandler.CPrintHandler.insts

def HTTPSERVER(fileName=None):
    if CCP4HTTPServerThread.CHTTPServerThread.insts is None:
        if TIMING: t1=time.time()
        obj = CCP4HTTPServerThread.CHTTPServerThread(fileName=fileName)
        if TIMING: 
            t2=time.time()
            print("Starting 15",obj.__class__, t2-t1)
    return CCP4HTTPServerThread.CHTTPServerThread.insts

def DBSERVER(fileName=None):
    if CCP4DbThread.CDbThread.insts is None:
        if TIMING: t1=time.time()
        obj = CCP4DbThread.CDbThread(fileName)
        obj.start()
        if TIMING: 
            t2=time.time()
            print("Starting 16",obj.__class__, t2-t1)
    return CCP4DbThread.CDbThread.insts

def FILEWATCHER():
    if CCP4AbstractViewer.CFileWatchTimer.insts is None:
        if TIMING: t1=time.time()
        CCP4AbstractViewer.CFileWatchTimer.insts =  CCP4AbstractViewer.CFileWatchTimer()
        if TIMING: 
            t2=time.time()
            print("Starting 17",CCP4AbstractViewer.CFileWatchTimer.insts.__class__, t2-t1)
    return CCP4AbstractViewer.CFileWatchTimer.insts

def DEMODATAMANAGER():
    if CCP4DemoData.CDemoData.insts is None:
        if TIMING: t1=time.time()
        CCP4DemoData.CDemoData.insts =  CCP4DemoData.CDemoData()
        if TIMING: 
            t2=time.time()
            print("Starting 18",CCP4DemoData.CDemoData.insts.__class__, t2-t1)
    return CCP4DemoData.CDemoData.insts

def JOBCONTROLLERGUI():
    if CCP4JobControlGui.CServerParamsDialog.insts is None:
        if TIMING: t1=time.time()
        CCP4JobControlGui.CServerParamsDialog.insts = CCP4JobControlGui.CServerParamsDialog()
        if TIMING: 
            t2=time.time()
            print("Starting 19",CCP4JobControlGui.CServerParamsDialog.insts.__class__, t2-t1)
    return  CCP4JobControlGui.CServerParamsDialog.insts

def SERVERSETUP():
    if CCP4JobController.CServerSetup.insts is None:
        if TIMING: t1=time.time()
        p = CCP4JobController.CServerSetup()
        if TIMING: 
            t2=time.time()
            print("Starting 20",CCP4JobController.CServerSetup.insts.__class__, t2-t1)
    return CCP4JobController.CServerSetup.insts
