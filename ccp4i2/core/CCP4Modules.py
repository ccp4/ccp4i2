# Used because of circular imports


def COMFILEPATCHMANAGER():
    from . import CCP4ComFilePatchManager
    return CCP4ComFilePatchManager.COMFILEPATCHMANAGER()

def CUSTOMTASKMANAGER():
    from . import CCP4CustomTaskManager
    return CCP4CustomTaskManager.CUSTOMTASKMANAGER()

def DBSERVER(fileName=None):
    from ..qtcore import CCP4DbThread
    return CCP4DbThread.DBSERVER(fileName)

def DEMODATAMANAGER():
    from ..qtgui import CCP4DemoData
    return CCP4DemoData.DEMODATAMANAGER()

def DUMMYMAINWINDOW():
    from ..qtgui import CCP4WebBrowser
    return CCP4WebBrowser.DUMMYMAINWINDOW()

def FILEWATCHER():
    from ..qtgui import CCP4AbstractViewer
    return CCP4AbstractViewer.FILEWATCHER()

def HTTPSERVER(fileName=None):
    from ..qtcore import CCP4HTTPServerThread
    return CCP4HTTPServerThread.HTTPSERVER(fileName)

def IMPORTEDJOBMANAGER():
    from . import CCP4ImportedJobManager
    return CCP4ImportedJobManager.IMPORTEDJOBMANAGER()

def JOBCONTROLLER():
    from ..qtcore import CCP4JobController
    return CCP4JobController.JOBCONTROLLER()

def JOBCONTROLLERGUI():
    from ..qtgui import CCP4JobControlGui
    return CCP4JobControlGui.JOBCONTROLLERGUI()

def LAUNCHER():
    from ..qtcore import CCP4Launcher
    return CCP4Launcher.LAUNCHER()

def MIMETYPESHANDLER():
    from ..qtcore import CCP4CustomMimeTypes
    return CCP4CustomMimeTypes.MIMETYPESHANDLER()

def PIXMAPMANAGER():
    from qtgui import CCP4Widgets
    return CCP4Widgets.PIXMAPMANAGER()

def PREFERENCES():
    from . import CCP4Preferences
    return CCP4Preferences.PREFERENCES()

def PROCESSMANAGER():
    from . import CCP4ProcessManager
    return CCP4ProcessManager.PROCESSMANAGER()

def PROJECTSMANAGER():
    from . import CCP4ProjectsManager
    return CCP4ProjectsManager.PROJECTSMANAGER()

def SERVERSETUP():
    from . import CCP4ServerSetup
    return CCP4ServerSetup.SERVERSETUP()

def TASKMANAGER():
    from . import CCP4TaskManager
    return CCP4TaskManager.TASKMANAGER()

def WEBBROWSER(index=-1, new=False, mini=False):
    from ..qtgui import CCP4WebBrowser
    return CCP4WebBrowser.WEBBROWSER(index, new, mini)

def WORKFLOWMANAGER():
    from . import CCP4WorkflowManager
    return CCP4WorkflowManager.WORKFLOWMANAGER()
