"""
     CCP4WebBrowser.py: CCP4 GUI Project
   Liz Potterton Jan 2010 - Copied from MGWebBrowser with changes to make stand-alone and to enable
                            Qt plugins and custom mime types
                 Feb 2010 - Make menus/toolbars customisable, add CCP4ProjectManger as placeholder
"""

##@package CCP4WebBrowser (QtWebKit) The web browser framework

from collections.abc import Callable
import functools
import glob
import os
import re
import sys
import time

from lxml import etree
from PySide2 import QtCore, QtGui, QtWebEngineWidgets, QtWidgets

from qtgui import CCP4WebView
from core import CCP4Modules
from core.CCP4ErrorHandling import *
from qtgui import CCP4WebToolBarButtons

def WEBBROWSER(index = -1,new=False,mini=False):
    if new or (len(CBrowserWindow.Instances)==0 and index<0):
        t = CBrowserWindow(welcome=False)
        if mini: t.setMini(True)
        t.show()
        t.raise_()
        return t
    index = min(index,0)
    if index<len(CBrowserWindow.Instances):
        win = list(CBrowserWindow.Instances)[index]
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
    if CBrowserWindow.Dummy is None:
        CBrowserWindow.Dummy = CBrowserWindow()
        CBrowserWindow.Dummy.hide()  # KJS: Looks like there is no hide present. Check to see if this is used anywhere.
        try:
            CBrowserWindow.Instances.remove(CBrowserWindow.Dummy)
        except:
            print('Error in DUMMYMAINWINDOW')
    return CBrowserWindow.Dummy


def setupWebkit():
    try:   # nice.... not.... also our old friend the try-except-pass (check this is indent'ed ok )
        if sys.platform == 'win32':
            homedir = os.environ.get('USERPROFILE')
        else:
            homedir = os.environ.get('HOME')
    except:
        pass
    if homedir and os.path.exists(homedir):
        webkit_db = os.path.join(homedir,'.qwebkit_idb')
        if os.path.exists(webkit_db):
            QtWebKit.QWebSettings.setIconDatabasePath(webkit_db)
        else:
            try:
                os.mkdir(webkit_db)
                QtWebKit.QWebSettings.setIconDatabasePath(webkit_db)
            except:
                # Favicons disabled
                pass

def checkRunningJobs():
    runningJobs = CCP4Modules.PROJECTSMANAGER().db().getRunningJobs()
    topRunningJobs = {}
    for job in runningJobs:
        if job[5] is None and len(job[1].split('.')) == 1:
            projectName=CCP4Modules.PROJECTSMANAGER().db().getProjectInfo(projectId=job[3],mode='projectname')
            if projectName in topRunningJobs:
                topRunningJobs[projectName].append(job)
            else:
                topRunningJobs[projectName] = [job]
    return topRunningJobs

def exitBrowser():
    if "-dbLocal" in sys.argv:
        topRunningJobs = checkRunningJobs()
        if len(topRunningJobs) > 0:
            message = "The following jobs are still running and you are in 'Local database mode'."
            if bool(CCP4Modules.PREFERENCES().DBLOCAL_QUIT_RUNNING):
                message += "<br/><br/>If you quit, results from these running jobs may be lost <em> but the job will appear to still be running when you restart CCP4i2</em>.<br/><br/>Do you really want to quit?<br/><br/>"
            else:
                message += "<br/><br/>Please wait for these jobs to finish before quitting. Or kill them if you really want to quit now.<br/><br/>"
            for proj,jobs in list(topRunningJobs.items()):
                message += "<b>Project: "+proj+"</b><br>"
                for job in jobs:
                    message += str(job[1]) + " " + job[2] + "<br/>"
                message += "<br/><br/>"
    
            mb = QtWidgets.QMessageBox(QtWidgets.QMessageBox.Question,"Jobs still running",message)

            if bool(CCP4Modules.PREFERENCES().DBLOCAL_QUIT_RUNNING):
                mb.addButton("Really quit",QtWidgets.QMessageBox.AcceptRole)
                mb.addButton("Cancel",QtWidgets.QMessageBox.RejectRole)
                rv = mb.exec_()
                if rv == 1:
                    return
            else:
                mb.addButton("Cancel",QtWidgets.QMessageBox.RejectRole)
                rv = mb.exec_()
                return
    purgeStatusFiles(0)
    from qtgui import CCP4ProjectViewer
    #print '*CCP4WebBrowser exitBrowser'
    CMainWindow.queryClose=False
    rv = saveStatus()
    purgeStatusFiles(2)
    for win in CCP4ProjectViewer.CProjectViewer.Instances:
        win.Exit()
    CCP4Modules.PREFERENCES().save()
    #print 'purgeStatusFiles done'
    #CCP4Modules.QTAPPLICATION().setActiveWindow(CCP4Modules.WEBBROWSER())
    CCP4Modules.QTAPPLICATION().closeAllWindows()

def saveStatus():
    from core import CCP4Utils
    from qtgui import CCP4I1Projects
    from core import CCP4File
    if CMainWindow.STATUS_SAVED:
        return ''
    from qtgui import CCP4ProjectViewer
    body = etree.Element('body')
    root = etree.Element('windows')
    body.append(root)
    #print 'saveStatus',CBrowserWindow.Instances,CCP4ProjectViewer.CProjectViewer.Instances
    #import traceback
    #traceback.print_stack(limit=5)
    for win in CCP4I1Projects.CI1ProjectViewer.Instances:
        win.Exit()
    for win in CBrowserWindow.Instances:
        winEle = etree.Element('browser')
        try:
            size = (win.size().width(),win.size().height())
            sizeEle = etree.Element('windowSize')
            sizeEle.text = str(size[0]) + ',' + str(size[1])
            winEle.append(sizeEle)
        except:
            pass
        #print 'CMainWindow.saveStatus browser',win.tab().count()
        for ii in range(win.tab().count()):
            tabEle = etree.Element('tab')
            winEle.append(tabEle)
            fileEle =  etree.Element('filename')
            try:
                url = win.tab().widget(ii).url()
            except:
                pass
            else:
                #filename = os.path.normpath(str(win.tab().widget(ii).fileName))
                #print 'CMainWindow.saveStatus',url.isLocalFile()
                if url.isLocalFile():
                    filename = url.toLocalFile().__str__().strip('file:')
                    if not os.path.exists(filename) and os.path.splitext(filename)[1] == '.htm':
                        filename = filename + 'l'
                    if os.path.relpath(filename, CCP4Utils.getCCP4I2Dir()).startswith('docs'):
                        fileEle.text = '$CCP4I2/' + os.path.relpath(filename, CCP4Utils.getCCP4I2Dir())
                    else:
                        fileEle.text = filename
                else:
                    fileEle.text = url.toString()
            #print 'CMainWindow.saveStatus browser tab file',win.tab().widget(ii).fileName
            tabEle.append(fileEle)
            title = win.tab().widget(ii).title()
            if title is not None:
                titleEle = etree.Element('title')
                titleEle.text = title
            tabEle.append(titleEle)
        root.append(winEle)
    for win in CCP4ProjectViewer.CProjectViewer.Instances:
        size = None
        try:
            projectName = CCP4Modules.PROJECTSMANAGER().db().getProjectInfo(projectId=win.getProject(), mode='projectname')
            jobNumber = win.getOpenJobNumber()
            size = (win.size().width(),win.size().height())
        except:
            pass
        else:
            winEle = etree.Element('projectViewer')
            projectEle = etree.Element('project')
            #print 'CMainWindow.saveStatus projectName',projectName
            projectEle.text = projectName
            winEle.append(projectEle)
            if jobNumber is not None:
                jobEle = etree.Element('openJobNumber')
                jobEle.text = str(jobNumber)
                winEle.append(jobEle)
                if size is not None:
                    sizeEle = etree.Element('windowSize')
                    sizeEle.text = str(size[0])+','+str(size[1])
                    winEle.append(sizeEle)
            # Issues with handing ByteArray - skip this for now
            #ele = etree.SubElement(projectEle,'headerState')
            #ele.text = win.projectWidget().getHeaderState()
            root.append(winEle)
    #print 'CCP4WebBrowser.saveStatus',etree.tostring(body,pretty_print=True)
    statusFile = os.path.join(CCP4Utils.getDotDirectory(), 'status', 'status_' + str(int(time.time())) + '.ccp4i2_status.xml')
    f= CCP4File.CI2XmlDataFile(statusFile)
    f.header.function.set('STATUS')
    f.header.setCurrent()
    f.saveFile(bodyEtree=body)
    CMainWindow.STATUS_SAVED = True
    print('CCP4i2 status saved to', statusFile)
    return statusFile

def restoreStatus():
    from core import CCP4Utils
    from qtgui import CCP4ProjectViewer
    # For now just try to restore projects
    statusEtree = retrieveStatus()
    if statusEtree is None:
        return
    nBrowserWindows = 0
    for winEle in statusEtree.iterchildren():
        #print 'CCP4WebBrowser.restoreStatus', winEle.tag
        if str(winEle.tag) == 'projectViewer':
            for  projectEle in winEle.iter('project'):
                projectName=projectEle.text
                try:
                    projectId = CCP4Modules.PROJECTSMANAGER().db().getProjectId(projectName=projectName)
                except:
                    projectId = None
                    print('Unable to open project - ' + projectEle.text + '- not found in database')
                # 'restoreStatus ',projectEle.text,projectId
                if projectId is not None:
                    # Do we know the open job?
                    #print 'Opening project:',projectEle.text
                    jobId= None
                    jobEle = winEle.find('openJobNumber')
                    if jobEle is not None:
                        try:
                            jobId = CCP4Modules.PROJECTSMANAGER().db().getJobId(jobNumber=str(jobEle.text), projectName=projectName)
                        except CException as e:
                            #print e.report()
                            jobId = None
                    try:
                        if jobEle is not None:
                            print('Opening Job:', jobEle.text, jobId)
                        proj = CCP4ProjectViewer.CProjectViewer(projectId=projectId, jobId=jobId)
                        try:
                            sizeStr = winEle.find('windowSize').text.split(',')
                            size = (int(sizeStr[0]), int(sizeStr[1]))
                        except:
                            size = CCP4ProjectViewer.DEFAULT_WINDOW_SIZE
                        proj.resize(size[0], size[1])
                        # Issues with handing ByteArray - skip this for now
                        #headerEle = winEle.find('headerState')
                        #if headerEle is not None:
                        #  proj.projectWidget().setHeaderState(str(headerEle.text()))
                        proj.menuBar().show()
                        proj.show()
                        proj.raise_()
                    except:
                        raise
                        print("Failed to open project", projectId, jobId)
        elif str(winEle.tag) == 'browser':
            if nBrowserWindows == 0:
                browser = CCP4Modules.WEBBROWSER()
                for w in CBrowserWindow.Instances:
                    browser.loadFinished.connect(w.handleFindNext)
            else:
                browser = CBrowserWindow()
                for w in CBrowserWindow.Instances:
                    browser.loadFinished.connect(w.handleFindNext)
                browser.show()
                browser.raise_()
            nBrowserWindows = nBrowserWindows + 1
            for tab in winEle.iterchildren():
                if tab.tag == 'windowSize':
                    try:
                        sizeStr = tab.text.split(',')
                        browser.resize(int(sizeStr[0]), int(sizeStr[1]))
                    except:
                        pass
                elif tab.tag == 'tab':
                    fileName = None
                    title = None
                    for ele in tab.iterchildren():
                        if ele.tag == 'filename':
                            fileName = str(ele.text)
                            #print 'restoreStatus fileName',fileName
                            if fileName.startswith('$CCP4I2'):
                                fileName = os.path.join(CCP4Utils.getCCP4I2Dir(), fileName[8:])
                            elif ele.tag == 'title':
                                title = str(ele.text)
                            if fileName is not None:
                                #print 'Opening file:',fileName
                                if fileName.startswith('http'):
                                    url = QtCore.QUrl(fileName)
                                    browser.loadPage(url)
                                else:
                                    browser.openFile(fileName, title=title, internal=True)

def applyCommandLine(args):
    # If there is a file/project on the com line.
    from core import CCP4Utils
    from qtgui import CCP4ProjectViewer
    #print 'CCP4Browser.applyCommandLine',args
    iArg = 0
    while iArg < len(args):
        filepath = None
        if args[iArg][0:2] == '-t':
            top_path = CCP4Utils.getCCP4I2Dir()
            filepath = os.path.join(top_path, 'test', 'data', 'test_plugin.html')
        else:
            projectName = args[iArg]
            try:
                projectId = CCP4Modules.PROJECTSMANAGER().db().getProjectId(projectName=projectName)
            except:
                projectId = None
            #print 'CCP4WebBrowser.applyCommandLine projectId', projectId
            if projectId is None:
                filepath  = os.path.abspath(args[iArg])
                #print 'CCP4Browser.applyCommandLine filepath', filepath
                if os.path.exists(filepath):
                    #print 'Opening file:', filepath
                    CCP4Modules.WEBBROWSER().openFile(filepath, internal=True)
                else:
                    print(args[iArg],'not recognised as project name or file')
            else:
                #print 'Opening project:',projectName
                proj = CCP4ProjectViewer.CProjectViewer(projectId=projectId)
                proj.resize(CCP4ProjectViewer.DEFAULT_WINDOW_SIZE[0], CCP4ProjectViewer.DEFAULT_WINDOW_SIZE[1])
                proj.show()
                proj.raise_()
        iArg = iArg + 1
    if len(CBrowserWindow.Instances) + len(CCP4ProjectViewer.CProjectViewer.Instances) <= 0:
        #print 'Opening default browser window'
        browser = CBrowserWindow()
        browser.show()


def retrieveStatus():
    from core import CCP4File
    from core import CCP4Utils
    statusFileList = glob.glob(os.path.join(CCP4Utils.getDotDirectory(), 'status', 'status_*.ccp4i2_status.xml'))
    if len(statusFileList) == 0:
        return None
    else:
        statusFileList.sort()
        print('Retrieving status file:' + statusFileList[-1])
        f = CCP4File.CI2XmlDataFile(statusFileList[-1])
        body = f.getBodyEtree()
        root = body.find('windows')
        return root

def purgeStatusFiles(leave=1):
    from core import CCP4Utils
    #print 'purgeStatusFiles leave',leave
    statusFileList = sorted(glob.glob(os.path.join(CCP4Utils.getDotDirectory(), 'status', 'status_*.ccp4i2_status.xml')))
    #print 'purgeStatusFiles',statusFileList
    if leave == 0:
        for sFile in statusFileList:
            os.remove(sFile)
    elif len(statusFileList) > leave:
        for sFile in statusFileList[0:-leave]:
            os.remove(sFile)

#-------------------------------------------------------------------
def OPENFILE(fileName=None, cformat=None, title=None, toFront=False):
#-------------------------------------------------------------------
    mimeTypeHandler = CCP4Modules.MIMETYPESHANDLER()
    if os.path.isdir(fileName):
        cformat = 'dir'
    if mimeTypeHandler is None:
        # No handler so just throw it at web browser and hope
        if cformat is None:
            cformat = 'text/html'
    else:
        if cformat is None:
            cformat = mimeTypeHandler.formatFromFileExt(fileName)
        if not mimeTypeHandler.isSupportedFormat(format):
            cformat = 'text/html'
    #print 'CCP4WebBrowser.OPENFILE format',fileName,format
    if mimeTypeHandler.useDesktopServices(format):
        abs_fileName = os.path.abspath(fileName)
        #print 'calling QDesktopServices',abs_fileName
        url = QtCore.QUrl.fromLocalFile(abs_fileName)
        rv = QtGui.QDesktopServices.openUrl(url)
        if not rv:
            QtWidgets.QMessageBox.warning(None, 'CCP4i2 display ' + fileName,
                                      'Attempting to display file ' + os.path.split(abs_fileName)[-1] + '\nusing desktop services failed')
            return None
    widgetClassList = mimeTypeHandler.getViewers(format)
    if len(widgetClassList) > 0:
        widgetClass = widgetClassList[0]
        if isinstance(widgetClass,str):
            # This is a keyword for the launcher
            CCP4Modules.LAUNCHER().openInViewer(viewer=widgetClass, fileName=fileName)
            return None
    CCP4Modules.WEBBROWSER().openFile(fileName=fileName, format=format, title=title, toFront=toFront)


class CTabWidget(QtWidgets.QTabWidget):

    def __init__(self, parent=None, name=None, mini=False):
        QtWidgets.QTabWidget.__init__(self, parent)
        self.setTabsClosable(True)
        self.tabCloseRequested.connect(self.closeTab);
        self.setObjectName(name)
        self.currentChanged.connect(self.handleCurrentChanged)
        self.currentOpen = -1
        #try:
        #  self.setTabsCloseable(True)
        #except:

        """
        if mini:
            self.setCornerWidget(QtWidgets.QPushButton(self.tr("Expand window")))
        else:
            self.setCornerWidget(QtWidgets.QPushButton(self.tr("Close tab")))
        """

    def closeTab(self,i):
        if i == -1:
            return
        tabItem = self.widget(i)
        self.removeTab(i)
        #Is this actually necessary?
        tabItem.close()
        tabItem.deleteLater()

    def deleteTabWidget(self, idx=None, widget=None):
        if idx is None:
            for ii in range(self.count()):
                if self.widget(ii) == widget: 
                    idx = ii
                    break
        if idx is None: return
        self.widget(idx).close()
        self.widget(idx).deleteLater()
        # We are keeping track of currentOpen to be sure to
        # call handleTabbedClosed() when necessary
        if self.currentOpen > idx:
            self.currentOpen = self.currentOpen - 1
        elif self.currentOpen == idx:
            self.currentOpen = -1
        self.removeTab(idx)

    @QtCore.Slot(int)
    def handleCurrentChanged(self,indx):
        #print 'CTabWidget.handleCurrentChanged', indx
        self.currentOpen = indx
        if self.widget(indx) is None:
            return
        self.parent().parent().editSplitter.addressEdit.setText(self.widget(indx).title())
        if hasattr(self.widget(indx),"url") and isinstance(self.widget(indx).url, Callable):
            self.parent().parent().editSplitter.addressEdit.setText(self.widget(indx).url().toString())
        else:
            if hasattr(self.widget(indx),"fileName"):
                self.parent().parent().editSplitter.addressEdit.setText(QtCore.QUrl.fromLocalFile(self.widget(indx).fileName).toString())
        self.widget(indx).handleTabbedOpen()
        self.widget(indx).setFocus(QtCore.Qt.OtherFocusReason)
        self.parent().parent().updateActionEnabled()

class CEditSplitter(QtWidgets.QSplitter):
    def __init__(self,parent):
        QtWidgets.QSplitter.__init__(self,parent)
        self.addressEdit = QtWidgets.QLineEdit()
        self.searchEdit = QtWidgets.QLineEdit()
        searchStyle = """QLineEdit {border: 1px solid gray;border-radius: 6;border-style: inset;padding: 0 3px;}"""
        self.searchEdit.setStyleSheet(searchStyle)
        self.addWidget(self.addressEdit) 
        self.addWidget(self.searchEdit)
        self.setStretchFactor(0, 2)


class CMenuBar(QtWidgets.QMenuBar):

    def __init__(self, parent):
        from core.CCP4Bazaar import bzrlib_exists
        QtWidgets.QMenuBar. __init__(self, parent)
        self.menuDefinitions = {}
        for menuName, menuTitle in [ ['File', '&File/Projects'], ['Edit', '&Edit'], ['View', '&View'],
                                    ['Utilities', '&Utilities'],['Help', 'Help/Tutorials']]:
            print("Adding menu",menuName, menuTitle)
            self.addMenu(menuName, menuTitle)
        #self.menuDefinitions['File'] = ['open_tab','open_browser', 'open', 'sep','close_window', 'quit']
        #self.menuDefinitions['Edit'] = ['find','preferences','ccp4i2_config']

#It is a bit cheeky inserting help_about here, but we know that OS X will intercept this and stick in Application Menu.
#And the 'Edit' menu does not seem to be redrawn (at present!). I guess a special '_dummy' menu as MG does it would allow
#'Edit' to be dynamic in future.
        if sys.platform == "darwin":
            self.menuDefinitions['Edit'] = ['find', 'preferences','help_about']
        else:
            self.menuDefinitions['Edit'] = ['find', 'preferences']
        self.menuDefinitions['View'] = ['zoomIn','zoomOut','resetZoom']
        self.menuDefinitions['Utilities'] = ['listProcesses', 'manageImportFiles', 'sendReport',['System administrator tools','update_core',
                                                                         'serverSetup','pdb_redo_setup','import_task'],
                                             ['Developer tools','redo_report','view_report_source','remake_cached_lookups','view_test_report',
                                              'export_task','list_tasks','grab_task','compress_demo_data','auto_task_docs',['Program log', 'programPrintLog','programHTTPLog'], 'clearStatusFiles', 'editScript', 'backup_db_xml']]
        self.menuDefinitions['Customisation'] = ['manageWorkflows', 'patchComFile', 'customTask', 'importJob']
        self.menuDefinitions['file'] = ['manager_projects',
                                            'new_project',['View old CCP4i projects', 'open_i1projects_default',
                                                           'open_i1projects_select']]
        if sys.platform == "darwin":
            self.menuDefinitions['Help'] = ['demo_data_info', 'help_quickstart', 'help_quickexpert', 'help_youtube',
                                        'task_docs', 'cloud_docs', 'help_ccp4i2',
                                        'tips_of_the_day']
        else:
            self.menuDefinitions['Help'] = ['help_about', 'demo_data_info', 'help_quickstart', 'help_quickexpert', 'help_youtube',
                                        'task_docs', 'cloud_docs', 'help_ccp4i2',
                                        'tips_of_the_day']
        if bzrlib_exists:
            self.menuDefinitions['Utilities'][-1].insert(1, 'update_gui')
        # Beware need to define the 'quit' in File menu at startup for the slot to exitBrowser() to work
        # This is likely a mac-specific thing since the quit gets move to application menu
        from qtgui import CCP4GuiUtils
        CCP4GuiUtils.populateMenu(self.parent(), self.menuWidget('Edit'), self.menuDefinition('Edit'), default_icon='')
        CCP4GuiUtils.populateMenu(self.parent(), self.menuWidget('View'), self.menuDefinition('View'), default_icon='')
        CCP4GuiUtils.populateMenu(self.parent(), self.menuWidget('File'), self.menuDefinition('File'), default_icon='')

    def menuWidget(self,menuName):
        return self.findChild(QtWidgets.QMenu, menuName)

    def updateMenu(self, menuName):
        from qtgui import CCP4GuiUtils
        menuWidget = self.menuWidget(menuName)
        if menuWidget is None:
            pass
        elif menuName == 'History':
            menuDefn = []
            menuWidget.clear()
            menuWidget.setEnabled(False)
            if hasattr(self.parent(),"tab"):
                widget = self.parent().tab().currentWidget()
                items = list(widget.history().items())
                def goToHistory(i):
                    item = widget.history().itemAt(i)
                    widget.history().goToItem(item)
                i = 0
                for item in items:
                    self.parent().setActionDefinition('history_' + str(item.url().toString()), 
                        dict(text = str(item.url().toString()), tip = "Go to history item",
                        slot = functools.partial(goToHistory,i),
                        enabled = 1))
                    menuDefn.append('history_' + str(item.url().toString()))
                    i += 1
            CCP4GuiUtils.populateMenu(self.parent(), menuWidget, menuDefn, default_icon='')
            if len(menuDefn) > 0:
                menuWidget.setEnabled(True)
        elif menuName == 'File':
            recentProjects = CCP4Modules.PROJECTSMANAGER().db().getRecentProjects(order='access', limit=11)
            menuDefn = ['manage_projects', 'new_project', 'export_project',
                        'import_project', 'sep']
            for projectId, projectName, accessTime in recentProjects[:10]:
                if 'open_project_' + str(projectId) not in self.parent().actionDefinitions:
                    self.parent().setActionDefinition('open_project_' + str(projectId), 
                                                      dict(text = projectName, tip = "Open this project",
                                                           slot = functools.partial(self.parent().openProject, projectId, projectName),
                                                           enabled = 1))
                menuDefn.append('open_project_' + str(projectId))

            #SJM 16/03/2018 - This is arguably unneccessary. Manage/open is perhaps better.
            #                 I have however removed 'more_projects' as it is duplicate of Manage/open.
            proj_dir_list0=CCP4Modules.PROJECTSMANAGER().db().getProjectDirectoryList()
            if len(proj_dir_list0)>15:
                submenuMore = ["More projects"]
                chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                binsize = 3
                allTheseChars = []
                submenDict = {}

                for i in range(len(chars)//int(binsize)):
                    if i < len(chars)//int(binsize) -1:
                        theseChars = chars[i*binsize:(i+1)*binsize]
                    else:
                        theseChars = chars[i*binsize:]
                    idx = theseChars[0].upper()+"-"+theseChars[-1].upper()
                    submenuMore.append([idx])
                    allTheseChars.append((theseChars,idx))
                    submenDict[idx] = submenuMore[-1]

                submenuMore.append(["0-9"])
                allTheseChars.append(("0123456789","0-9"))
                submenDict["0-9"] = submenuMore[-1]

                submenuOther = ["other"]
                submenuMore.append(submenuOther)
                
                for pid,pName,pDir,parent,created,lastAccess in proj_dir_list0:
                    if 'open_project_' + str(pid) not in self.parent().actionDefinitions:
                        #If this really is a performance hit, surely we can generate this list at start up.
                        self.parent().setActionDefinition('open_project_' + str(pid), 
                                                      dict(text = pName, tip = "Open this project",
                                                           slot = functools.partial(self.parent().openProject, pid, pName),
                                                           enabled = 1, checkable = 1,
                                                           checked = functools.partial( self.isProjectOpen, pid)))
                    skipped = True
                    for submen in allTheseChars:
                        if pName.upper()[0] in submen[0]:
                            submenDict[submen[1]].append('open_project_' + str(pid))
                            skipped = False
                    if skipped:
                        print("Skipped",pName)
                        submenuOther.append('open_project_' + str(pid))

                menuDefn.append(submenuMore)
            menuDefn.append("sep")
            menuDefn.append(['View old CCP4i projects', 'open_i1projects_default','open_i1projects_select'])
            menuDefn.append(['Browser','open_browser', 'open', 'sep','close_window'])
            menuDefn.append('quit')

            #menuDefn.append('more_projects')

            # Enable more_projects if >10 projects in db 
            if 'more_projects' in self.parent().actionDefinitions:
                self.parent().actionDefinitions['more_projects']['enabled'] = len(recentProjects) > 10
            menuWidget.clear()
            CCP4GuiUtils.populateMenu(self.parent(), menuWidget, menuDefn, default_icon='')
            #self._redrawProjectMenu = False
        elif menuName == 'Edit':
            menuWidget.clear()
            if sys.platform == "darwin":
                CCP4GuiUtils.populateMenu(self.parent(), menuWidget, ['find'], default_icon='')
            else:
                CCP4GuiUtils.populateMenu(self.parent(), menuWidget, ['find', 'preferences'], default_icon='')
        elif menuName == 'Utilities':
            menuWidget.clear()
            menuDefn = [['Copy demo data to project']]
            if self.parent().isProjectViewer():
                testDatasets = CCP4Modules.DEMODATAMANAGER().getTestDatasets()
                for dataset,label in testDatasets:
                    self.parent().setActionDefinition('download_test_' + dataset,
                                                      dict(text=label, tip="Copy this data to project directory",
                                                           slot = functools.partial(CCP4Modules.DEMODATAMANAGER().copyDemoDataToProject,
                                                                                    self, self.parent().taskFrame.openJob.projectId, dataset),
                                                           enabled = self.parent().isProjectViewer))
                    menuDefn[0].append('download_test_' + dataset)
                menuDefn[0].append('sep')
            #menuDefn[0].append('demo_data_info')
            #menuDefn[0].append('download_demo_data')
            menuDefn.extend( self.menuDefinitions['Utilities'] )
            CCP4GuiUtils.populateMenu(self.parent(), menuWidget, menuDefn, default_icon='')
        else:
            menuWidget.clear()
            CCP4GuiUtils.populateMenu(self.parent(), menuWidget, self.menuDefinition(menuName), default_icon='')

    def isProjectOpen(self, projectId):
        from qtgui import CCP4ProjectViewer
        for proj in CCP4ProjectViewer.CProjectViewer.Instances:
            if hasattr(proj,'taskFrame') and hasattr(proj.taskFrame,'openJob') and proj.taskFrame.openJob.projectId == projectId:
                return True
        return False

    def addMenu(self, menuName=None, menuTitle=None, updateFunction=None):
    #def addMenu(self, menuName=None, menuTitle=None):
        if menuName is None or menuTitle is None:
            print('ERROR in CMenuBar.addMenu; no menu name or no menu title given')
            return None
        widget = self.menuWidget(menuName)
        if widget is None:
            widget = QtWidgets.QMenuBar.addMenu(self, menuTitle)
            widget.setObjectName(menuName)
            if updateFunction is None:
                widget.aboutToShow.connect(functools.partial(self.updateMenu, menuName))
                widget.addAction("dum")
            else:
                widget.aboutToShow.connect(functools.partial(updateFunction, menuName))
                widget.addAction("dum")
        return widget

    def removeMenu(self, menuName=''):
        # Save current menus, clear and add all menus again
        currentMenuList = []
        found = False
        currentMenuWidgets = self.findChildren(QtWidgets.QMenu)
        for w in currentMenuWidgets:
            name = str(w.objectName())
            if name == menuName:
                found = True
            else:
                currentMenuList.append([str(w.objectName()), str(w.title())])
        if found:
            self.clear()
            for menuName, menuTitle in currentMenuList:
                self.addMenu(menuName, menuTitle)

    def setMenuTitle(self, menuName='', menuTitle=''):
        widget = self.menuWidget(menuName)
        if widget is not None:
            widget.setTitle(menuTitle)

    def menuDefinition(self, menuName=''):
        return self.menuDefinitions.get(menuName,[])

    def setMenuDefinition(self, menuName='', definition=[]):
        self.menuDefinitions[menuName] = definition

    def appendToMenu(self, menuName='', menuItem=''):
        # menuItem should be the name of an action or 'sep' (separator)
        if menuName not in self.menuDefinitions:
            print('ERROR in CMenuBar.appendToMenu: ' + menuName + ' menu does not exist')
        self.menuDefinitions[menuName].append(menuItem)

class CToolBar(QtWidgets.QToolBar):
    #FIXME And all the use of str should go in PyQt5.

    toolBarPreferencesMapping = {
       "task_menu" : "SHOW_TASK_MENU_TOOLBUTTON",
       "job_search" : "SHOW_JOB_SEARCH_TOOLBUTTON",
       "export_project" : "SHOW_EXPORT_PROJECT_TOOLBUTTON",
       "run" : "SHOW_RUN_TOOLBUTTON",
       "run_remote" : "SHOW_RUN_REMOTE_TOOLBUTTON",
       "clone" : "SHOW_CLONE_TOOLBUTTON",
       "task_help" : "SHOW_TASK_HELP_TOOLBUTTON",
       "references" : "SHOW_REFERENCES_TOOLBUTTON",
       "export_mtz" : "SHOW_EXPORT_MTZ_TOOLBUTTON",
       "view_coot" : "SHOW_VIEW_COOT_TOOLBUTTON",
       "view_ccp4mg" : "SHOW_VIEW_CCP4MG_TOOLBUTTON",
       "show_log" : "SHOW_SHOW_LOG_TOOLBUTTON",
       "show_i2run" : "SHOW_SHOW_I2RUN_TOOLBUTTON",
       "new_project2" : "NEW_PROJECT_TOOLBUTTON"
    }

    def __init__(self,parent,name='',title=''):
        QtWidgets.QToolBar.__init__(self, title, parent)
        self.setObjectName(name)
        if hasattr(CCP4Modules.PREFERENCES(), "TOOLBARBUTTONSSTYLE"):
            if int(CCP4Modules.PREFERENCES().TOOLBARBUTTONSSTYLE) == 0:
                self.setToolButtonStyle(QtCore.Qt.ToolButtonIconOnly)
            elif int(CCP4Modules.PREFERENCES().TOOLBARBUTTONSSTYLE) == 1:
                self.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
            else:
                self.setToolButtonStyle(QtCore.Qt.ToolButtonTextUnderIcon)
        # KJS : This is the place it is set, not sure where this is redrawn.
        if CCP4Modules.PREFERENCES().HD_ICONS:
            self.setIconSize(QtCore.QSize(48, 48))
        else:
            self.setIconSize(QtCore.QSize(24, 24))
        self.definition = []

    def append(self,toolName=''):
        if not self.definition.count(toolName):
            self.definition.append(toolName)
            self.redraw()
    
    def extend(self,toolNameList=[]):
        for toolName in toolNameList:
            if not self.definition.count(toolName) or toolName=="sep":
                self.definition.append(toolName)
        self.redraw()
    
    def remove(self,toolName=''):
        if self.definition.count(toolName):
            self.definition.remove(toolName)
            self.redraw()

    def redraw(self):
        self.clear()
        from qtgui import CCP4GuiUtils
        CCP4GuiUtils.populateToolBar(parent=self.parent(), toolBarWidget=self, definition=self.definition)

        children = self.findChildren(QtWidgets.QToolButton)
        for child in children:
            if child.defaultAction() is not None:
                theName = str(child.defaultAction().objectName())
                if theName in self.toolBarPreferencesMapping:
                    if not CCP4Modules.PREFERENCES().get(self.toolBarPreferencesMapping[theName]):
                        child.defaultAction().setVisible(False)

    def editPreferences(self):
        from . import CCP4ProjectViewer

        prefWidget = QtWidgets.QDialog()
        prefWidget.setWindowTitle("Edit visible tool buttons")
        children = self.findChildren(QtWidgets.QToolButton)
        listWidget = QtWidgets.QListWidget()
        layout = QtWidgets.QVBoxLayout()
        prefWidget.setLayout(layout)
        layout.addWidget(listWidget)

        for child in children:
            if child.defaultAction() is not None:
                item = QtWidgets.QListWidgetItem()
                item.setText(child.defaultAction().text())
                item.setIcon(child.defaultAction().icon())
                listWidget.addItem(item)
                item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
                item.setCheckState(QtCore.Qt.Unchecked)
                if child.defaultAction().isVisible():
                    item.setCheckState(QtCore.Qt.Checked)
                item.setData(QtCore.Qt.UserRole,child)

        def setItemVisibilities(item):
            checkState, tb = item.checkState(), item.data(QtCore.Qt.UserRole)
            theName = tb.defaultAction().objectName()
            for window in CCP4ProjectViewer.CProjectViewer.Instances:
                acts = window.findChildren(QtWidgets.QAction,theName)
                for act in acts: #Should be length 0 or 1
                    if checkState:
                        act.setVisible(True)
                    else:
                        act.setVisible(False)
            mapping = self.toolBarPreferencesMapping[str(theName)]
            if checkState:
                val = True
            else:
                val = False
            #FIXME - Aargh. There must be a nicer way.
            #MN Trying to increase code compactness/readability
            getattr(CCP4Modules.PREFERENCES(), mapping).set(val)
            #Replaces
            '''
            if mapping == "SHOW_TASK_MENU_TOOLBUTTON":
                CCP4Modules.PREFERENCES().SHOW_TASK_MENU_TOOLBUTTON.set(val)
            elif mapping == "SHOW_JOB_SEARCH_TOOLBUTTON":
                CCP4Modules.PREFERENCES().SHOW_JOB_SEARCH_TOOLBUTTON.set(val)
            elif mapping == "SHOW_EXPORT_PROJECT_TOOLBUTTON":
                CCP4Modules.PREFERENCES().SHOW_EXPORT_PROJECT_TOOLBUTTON.set(val)
            elif mapping == "SHOW_RUN_TOOLBUTTON":
                CCP4Modules.PREFERENCES().SHOW_RUN_TOOLBUTTON.set(val)
            elif mapping == "SHOW_RUN_REMOTE_TOOLBUTTON":
                CCP4Modules.PREFERENCES().SHOW_RUN_REMOTE_TOOLBUTTON.set(val)
            elif mapping == "SHOW_CLONE_TOOLBUTTON":
                CCP4Modules.PREFERENCES().SHOW_CLONE_TOOLBUTTON.set(val)
            elif mapping == "SHOW_TASK_HELP_TOOLBUTTON":
                CCP4Modules.PREFERENCES().SHOW_TASK_HELP_TOOLBUTTON.set(val)
            elif mapping == "SHOW_REFERENCES_TOOLBUTTON":
                CCP4Modules.PREFERENCES().SHOW_REFERENCES_TOOLBUTTON.set(val)
            elif mapping == "SHOW_EXPORT_MTZ_TOOLBUTTON":
                CCP4Modules.PREFERENCES().SHOW_EXPORT_MTZ_TOOLBUTTON.set(val)
            elif mapping == "SHOW_VIEW_COOT_TOOLBUTTON":
                CCP4Modules.PREFERENCES().SHOW_VIEW_COOT_TOOLBUTTON.set(val)
            elif mapping == "SHOW_VIEW_CCP4MG_TOOLBUTTON":
                CCP4Modules.PREFERENCES().SHOW_VIEW_CCP4MG_TOOLBUTTON.set(val)
            elif mapping == "SHOW_SHOW_LOG_TOOLBUTTON":
                CCP4Modules.PREFERENCES().SHOW_SHOW_LOG_TOOLBUTTON.set(val)
            elif mapping == "NEW_PROJECT_TOOLBUTTON":
                CCP4Modules.PREFERENCES().NEW_PROJECT_TOOLBUTTON.set(val)
            '''
        listWidget.itemChanged.connect(setItemVisibilities)
        prefWidget.exec_()

    def contextMenuEvent(self,e):
        if self.parent().objectName() != 'projectViewer':
            return
        menu = QtWidgets.QMenu(self)
        custAct = QtWidgets.QAction("Customize",self)
        custAct.triggered.connect(self.editPreferences)
        menu.addAction(custAct);
        menu.exec_(e.globalPos());

def isAlive(qobj):
    import shiboken2
    return shiboken2.isValid(qobj)

def mainWindowIcon():
    if CMainWindow._MAINWINDOWICON is None:
        from core import CCP4Utils
        fileName = os.path.join(CCP4Utils.getCCP4I2Dir(), 'qticons', 'ccp4.png')
        CMainWindow._MAINWINDOWICON = QtGui.QIcon(QtGui.QPixmap(fileName))
    return CMainWindow._MAINWINDOWICON


class CMainWindow(QtWidgets.QMainWindow):
    projectManagerDialog = None
    queryClose = False
    STATUS_SAVED = False
    _MAINWINDOWICON = None
    ERROR_CODES = {201 : {'description' : 'Error opening zip file for write'},
                   202 : {'description' : 'Error saving directory to zip file'},
                   203 : {'description' : 'Error closing zip file from read'},
                   204 : {'description' : 'Error opening zip file to read'},
                   205 : {'description' : 'Error extracting zip file from directory'},
                   206 : {'description' : 'Compressed task file does not contain appropriately named task'},
                   207 : {'description' : 'Can not import task as it requires overwriting exisiting task'},
                   208 : {'description' : 'Error deleting existing task directory'},  # KJS : Changed below to 209 .....
                   209 : {'description' : 'Selected inappropriate directory to save as compressed task file'}}

    def showTipsOfTheDay(self):
        from qtgui.CCP4TipsOfTheDay import CTipsOfTheDay
        tipsOfTheDay = CTipsOfTheDay()
        tipsOfTheDay.exec_()

    def backupDB(self):
        CCP4Modules.PROJECTSMANAGER().backupDB()

    def backupDBXML(self):
        CCP4Modules.PROJECTSMANAGER().backupDBXML()

    def __init__(self, parent=None):
        from qtcore import CCP4UpdateManager
        QtWidgets.QMainWindow.__init__(self,parent)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowIcon(mainWindowIcon())
        # Remove access to bzr as Andrey claims it causes access to repository that is acting badly possibly due to loading
        try:
            from core import CCP4Update
            self.version = CCP4Update.get_ccp4_str() + ' '
        except:
            self.version = 'DEVELOPMENT CCP4i2'
        self.preferenceWindow=None
        self.configWindow=None
        self.fileDialog = None
        self.findFrame= None
        self.um = CCP4UpdateManager.um
        self.initialiseActionDefinitions()
        self.setMenuBar(CMenuBar(self))
        if sys.platform == "darwin":
            self.setUnifiedTitleAndToolBarOnMac(True)
            shortcut = QtWidgets.QShortcut(QtGui.QKeySequence.Quit, self);
            shortcut.setContext(QtCore.Qt.ApplicationShortcut)
            shortcut.activated.connect(exitBrowser)

    def closeEvent(self, event):
        from qtgui import CCP4ProjectViewer
        from qtgui import CCP4I1Projects
        if CBrowserWindow.Dummy is None:
            nOpen = 0
        else:
            nOpen = len(CBrowserWindow.Dummy.findChildren(QtWidgets.QDialog))
        #print 'closeEvent Dummy.findChildren',nOpen
        nOpen = nOpen + len(CBrowserWindow.Instances) + len(CCP4ProjectViewer.CProjectViewer.Instances) + len(CCP4I1Projects.CI1ProjectViewer.Instances)
        #print 'closeEvent Instances',nOpen
        #print 'closingDown',CCP4Modules.QTAPPLICATION().closingDown()
        if nOpen > 1:
            #self.close()
            event.accept()
        else:
            if CMainWindow.queryClose:
                rv = QtWidgets.QMessageBox.question(self,'Close CCP4Browser','Close CCP4 Browser?',
                      QtWidgets.QMessageBox.Close|QtWidgets.QMessageBox.Cancel, QtWidgets.QMessageBox.Cancel)
            else:
                rv = QtWidgets.QMessageBox.Close
            if rv ==  QtWidgets.QMessageBox.Close:
                self.Exit()
                saveStatus()
                CCP4Modules.PREFERENCES().save()
                purgeStatusFiles(2)
                event.accept()
            else:
                event.ignore()

    def Exit(self):
        pass

    def windowTitle(self):
        return str(QtWidgets.QMainWindow.windowTitle(self))

    def initialiseActionDefinitions(self):
        self.actionDefinitions = {}
        self.actionDefinitions['quit'] = dict(text='Quit CCP4i2', tip='Close all CCP4i2 windows', slot=exitBrowser)
        self.actionDefinitions['help_about'] = dict(text="About", tip="CCP4i2 background",
                                                    slot=self.showAbout)
        self.actionDefinitions['welcome_ccp4i2'] = dict(text="Welcome - options to get started",
                                                        tip="Options to set up CCP4i2 Project",
                                                        slot=functools.partial(self.showHelp, 'ccp4i2'))
        self.actionDefinitions['help_quickstart'] = dict(text="Quickstart - 10 minute intro",
                                                         tip="Quick introduction to using i2",
                                                         slot = functools.partial(self.showHelp, 'quickstart'))
        self.actionDefinitions['help_quickexpert'] = dict(text = "Quick expert - more quick hints",
                                                          tip = "Part 2 of quick introduction to i2",
                                                          slot = functools.partial(self.showHelp, 'quickexpert'))
        self.actionDefinitions['help_youtube'] = dict(text = "View YouTube video",
                                                      tip = "Quick introduction in youtube (6 mins)",
                                                      slot = functools.partial(self.showHelp, 'youtube'), icon = 'youtube')
        self.actionDefinitions['help_ccp4i2'] = dict(text = "CCP4i2 Help", tip = "Documentation and 'where to start'",
                                                     slot = functools.partial(self.showHelp, 'tutorial'),
                                                     shortcut = QtGui.QKeySequence.HelpContents, icon = 'help')
        self.actionDefinitions['task_docs'] = dict(text = "Task documentation", tip = "Documentation for the individual tasks",
                                                   slot = functools.partial(self.showHelp, 'task_docs'))
        self.actionDefinitions['cloud_docs'] = dict(text = "Task Info from CCP4Cloud", tip = "Additional Documentation (from cloud project)",
                                                   slot = functools.partial(self.showHelp, 'cloud_docs'))
        self.actionDefinitions['help_ccp4_home'] = dict(text = "CCP4 Home Page", tip = "For further information on CCP4",
                                                        slot = functools.partial(self.showHelp, 'ccp4_home'), icon = 'ccp4')
        self.actionDefinitions['help_updates'] = dict(text = "CCP4 Updates", tip = "",
                                                      slot = functools.partial(self.showHelp, 'updates'))
        self.actionDefinitions['help_license'] = dict(text = "CCP4 Licence", tip = "",
                                                      slot = functools.partial(self.showHelp, 'license'),)
        self.actionDefinitions['tips_of_the_day'] = dict(text = "Tips of the day", tip = "",
                                                      slot = self.showTipsOfTheDay)
        self.actionDefinitions['backup_db_xml'] = dict(text = "Backup to XML files", tip = "Backup the database to XML files",
                                                      slot = self.backupDBXML)
        self.actionDefinitions['open'] = dict(text = "Open file", tip = "Open html, text, image or log file", 
                                              slot = self.openFileDialog, icon = 'fileopen', shortcut = QtGui.QKeySequence.Open, enabled = self.isNotProjectViewer)
        self.actionDefinitions['close_window'] = dict(text = "Close window", tip = "Close this CCP4i2 window",
                                                      slot = self.close, shortcut = QtGui.QKeySequence.Close)
        self.actionDefinitions['save'] = dict(text = "Save", tip = "Save to file", slot = self.handleSave,
                                              enabled = self.widgetIsSaveable)
        self.actionDefinitions['print'] = dict(text = "Print", tip = "Print file", slot = self.handlePrint,
                                               enabled = self.widgetIsPrintable)
        #FIXME - What does this do?
        self.actionDefinitions['run'] = dict(text = "Run", tip = "Run script", slot = self.handleRun,
                                             enabled = self.widgetIsRunable)
        self.actionDefinitions['find'] = dict(text = "Find", tip = "Find in file", slot = self.openFind, icon = 'search',
                                              checkable = 1, checked = self.isFindFrameOpen, enabled = self.widgetIsSearchable, shortcut = QtGui.QKeySequence.Find)
        self.actionDefinitions['preferences'] = dict(text = "Preferences", slot = self.openPreferences,
                                                     checkable = 1, checked = self.isPreferencesOpen,)
        self.actionDefinitions['ccp4i2_config'] = dict(text = "Configure CCP4i2", slot = self.openConfig,
                                                       checkable = 1, checked = self.isConfigOpen,)
        self.actionDefinitions['open_tab'] = dict(text = "New tab", tip = "Open a new web browser tab",
                                                      slot = self.createWindow, enabled = self.isNotProjectViewer, shortcut = "Ctrl+T",icon="newtab")
        self.actionDefinitions['open_browser'] = dict(text = "New browser window", tip = "Open a new web browser window",
                                                      slot = self.openBrowserWindow, enabled = 1, shortcut = QtGui.QKeySequence.New)
        self.actionDefinitions['open_i1projects_default'] = dict(text = "View default CCP4i projects", tip = "Import and display projects from previous CCP4 interface",
                                                                 slot =  functools.partial(self.openI1Projects, 'default'), enabled = 1)
        self.actionDefinitions['open_i1projects_select'] = dict(text = "View selected CCP4i projects", tip = "Import and display projects from previous CCP4 interface",
                                                                slot = functools.partial(self.openI1Projects, 'select'), enabled = 1)
        self.actionDefinitions['back'] = dict(text = "Back", tip = "Go back one page",
                                              slot = self.historyBack, shortcut = QtGui.QKeySequence.Back, enabled = 1)
        self.actionDefinitions['forward'] = dict(text = "Forward",tip = "Go forward one page",
                                                 slot = self.historyForward, shortcut = QtGui.QKeySequence.Forward, enabled = 1)
        self.actionDefinitions['reload'] = dict(text = "Reload",tip = "Reload current page",
                                                slot = self.reloadPage, shortcut = "Ctrl+R", enabled = 1)
        self.actionDefinitions['listProcesses'] = dict(text = "Running jobs and processes", tip = "Check jobs are still running",
                                                       slot = self.listProcesses)
        self.actionDefinitions['manageImportFiles'] = dict(text = "Manage imported files", tip = "Record provenance or remove imported files",
                                                           slot = self.openManageImportFiles, enabled = self.isProjectViewer)
        self.actionDefinitions['manageWorkflows'] = dict(text = "Workflows", tip = "Create & manage user defined work flows",
                                                         slot = self.openWorkflow)
        self.actionDefinitions['patchComFile'] = dict(text = "Task parameters and patches", tip = "Save/restore task parameters or customise a program command file",
                                                      slot = self.openPatchComFile)
        self.actionDefinitions['customTask'] = dict(text = "Custom tasks", tip = "Define a task",
                                                    slot = self.openCustomTasks)
        self.actionDefinitions['importJob'] = dict(text = "Import job", tip = "Report job run outside CCP4i2",
                                                   slot = self.openImportJobs)
        self.actionDefinitions['programPrintLog'] = dict(text = "Program print log", tip = "Show printed diagnostic for this run of the program",
                                                         slot = functools.partial(self.openApplication,'programPrintLog'), enabled = 1)
        self.actionDefinitions['programHTTPLog'] = dict(text = "HTTP server log", tip = "Show log from internal HTTP server",
                                                        slot = functools.partial(self.openApplication,'programHTTPLog'), enabled = 1)
        self.actionDefinitions['editScript'] = dict(text = "Edit/run script", tip = "Edit/run Python script within CCP4i2",
                                                    slot = functools.partial(self.openApplication,'editScript'), enabled = 1)
        self.actionDefinitions['sendReport'] = dict(text = "Send error report", tip = "Send error report to CCP4",
                                                    slot = self.openSendReport,enabled = self.isProjectViewer)
        self.actionDefinitions['update_core'] = dict(text = "Manage CCP4 updates", tip = "Examine, apply or remove CCP4 updates",
                                                     slot = self.um.manage, enabled = self.um.is_unlocked)
        self.actionDefinitions['import_task'] = dict(text = "Import task code", tip = "Load new task from compressed file",
                                                     slot = self.openImportTask, enabled = 1)
        self.actionDefinitions['export_task'] = dict(text = "Export task", tip = "Save task to compressed file",
                                                     slot = self.openExportTask, enabled = 1)
        self.actionDefinitions['list_tasks'] = dict(text = "List tasks", tip = "List tasks currently available in CCP4i2",
                                                    slot = self.listTasks, enabled = 1)
        self.actionDefinitions['grab_task'] = dict(text = "Grab task widget/report", tip = "Make screenshot of current task widget or report",
                                                   slot = self.grabWidget, enabled = self.isProjectViewer)
        self.setActionDefinition('compress_demo_data', dict(text = "Compress demo data", tip = "Create zip file to place on demo data download site",
                                                            slot = self.makeDemoData, enabled = 1))
        self.setActionDefinition('general_project', dict(text = "Create a General project", tip = "Create a project for odds and ends",
                                                         slot = self.makeGeneralProject, enabled = 1))
        self.setActionDefinition('new_project', dict(text = "New project", tip = "Create new CCP4 project",
                                                     slot = functools.partial(self.handleProjectMenu,'new_project'), enabled = 1))
        self.setActionDefinition('more_projects', dict(text = "More projects..",tip = "Show all projects",
                                                       slot = functools.partial(self.handleProjectMenu,'more_projects'), enabled = 1))
        self.setActionDefinition('manage_projects',dict(text = "Manage/open projects", tip = "Review, manage and open CCP4 projects",
                                                        slot = functools.partial(self.handleProjectMenu, 'manage_projects'), enabled = 1))
        self.setActionDefinition('export_project',dict(text = "Export project", tip = "Save project to a compressed file",
                                                       slot = self.handleProjectMenuExport, enabled = self.isProjectViewer, icon = 'export_arrow_new'))
        self.setActionDefinition('import_project',dict(text = "Import project", tip = "Load project from a compressed file",
                                                       slot = functools.partial(self.handleProjectMenu, 'import_project'), enabled = 1, icon = 'import_arrow_new'))
        self.setActionDefinition('serverSetup',dict(text="Configure servers for 'remote' run jobs", tip = "Specify host and mechanism to run remote jobs",
                                                    slot = self.openServerSetup, checked = self.isServerSetupOpen))
        self.setActionDefinition('pdb_redo_setup',dict(text="Set login tokens for PDB_REDO", tip = "This must be done at least once per year to enable PDB_REDO jobs to work",
                                                    slot = self.openConfigurePDBREDOTokens))
        self.setActionDefinition('redo_report',dict(text="Remake report", tip="Remake the task report",
                                                    slot=functools.partial(self.handleDeveloperTools, 'redo_report'),
                                                    enabled=self.reportAvailable))
        self.setActionDefinition('view_report_source', dict(text = "View report source", tip = "Display the report in a text viewer",
                                                            slot = functools.partial(self.handleDeveloperTools, 'view_report_source'),
                                                            enabled = self.reportAvailable))
        self.setActionDefinition('remake_cached_lookups', dict(text = "Remake cached lookups", tip = "Remake lookups, e.g. when a new module or datatype is added",
                                                            slot = functools.partial(self.handleDeveloperTools, 'remake_cached_lookups'),
                                                            enabled = self.isProjectViewer))
        self.setActionDefinition('view_test_report', dict(text = "Show project test report", tip = "Display the report from re-running the project",
                                                          slot = functools.partial(self.handleDeveloperTools, 'view_test_report'),
                                                          enabled = self.testReportAvailable))
        self.setActionDefinition('make_test_suite', dict(text="Convert project to test suite", tip="Add unittest template files",
                                                         slot=functools.partial(self.handleDeveloperTools, 'make_test_suite'),
                                                         enabled=self.isProjectViewer))
        self.setActionDefinition('demo_data_info', dict(text="Tutorials", tip="Description of downloadable demo data",
                                                        slot=functools.partial(self.showDemoDataInfo)))
        self.setActionDefinition('download_demo_data', dict(text="Download data", tip="Download data from internet",
                                                            slot=functools.partial(self.downloadDemoData)))
        self.setActionDefinition('auto_task_docs', dict(text="Recreate task docs index", tip="Auto generate docs/tasks/index.html",
                                                        slot=functools.partial(self.handleDeveloperTools, 'auto_task_docs')))

        self.actionDefinitions['resetZoom'] = dict(text = "Reset zoom",tip = "Reset zoom",
                                                slot = self.resetZoom, shortcut = "Ctrl+0", enabled = 1)

        self.actionDefinitions['zoomIn'] = dict(text = "Zoom in",tip = "Zoom in",
                                                slot = self.zoomIn, shortcut = QtGui.QKeySequence.ZoomIn, enabled = 1)

        self.actionDefinitions['zoomOut'] = dict(text = "Zoom out",tip = "Zoom out",
                                                slot = self.zoomOut, shortcut = QtGui.QKeySequence.ZoomOut, enabled = 1)


        self.setActionDefinition('clearStatusFiles', dict(text = "Clear shutdown/restart status", tip="Remove possibly corrupt status files",
                                                          slot = functools.partial(purgeStatusFiles,0)))

        otherBackShortCut = QtWidgets.QShortcut(QtGui.QKeySequence(self.tr("Ctrl+Left", "Back")),self)
        otherBackShortCut.activated.connect(self.historyBack)
        otherForwardShortCut = QtWidgets.QShortcut(QtGui.QKeySequence(self.tr("Ctrl+Right", "Forward")),self)
        otherForwardShortCut.activated.connect(self.historyForward)

    def resetZoom(self):
        if hasattr(self,"tab"):
            print(self.tab().currentWidget())
            self.tab().currentWidget().setZoomFactor(1.0)
            CCP4Modules.PREFERENCES().BROWSER_ZOOM_FACTOR.set(self.tab().currentWidget().zoomFactor())

    def zoomIn(self):
        if hasattr(self,"tab"):
            self.tab().currentWidget().setZoomFactor(self.tab().currentWidget().zoomFactor()*1.2)
            CCP4Modules.PREFERENCES().BROWSER_ZOOM_FACTOR.set(self.tab().currentWidget().zoomFactor())

    def zoomOut(self):
        if hasattr(self,"tab"):
            self.tab().currentWidget().setZoomFactor(self.tab().currentWidget().zoomFactor()/1.2)
            CCP4Modules.PREFERENCES().BROWSER_ZOOM_FACTOR.set(self.tab().currentWidget().zoomFactor())

    def updateActionEnabled(self, dummy1=None, dummy2=None):
        for actionName in ['print', 'run', 'find', 'save']:
            ifEnabled = self.actionDefinitions[actionName]['enabled']()
            action = self.findChild(QtWidgets.QAction,actionName)
            #print 'updateActionEnabled',actionName,ifEnabled,action
            if action is not None:
                action.setEnabled(ifEnabled)

    def zoomFactorChanged(self,val):
        CCP4Modules.PREFERENCES().BROWSER_ZOOM_FACTOR.set(self.tab().currentWidget().zoomFactor())

    def createWindow(self,type=None):
        view = CCP4WebView.CWebView()
        view.searchFound.connect(self.findFrame.searchCallback)
        view.zoomFactorChanged.connect(self.zoomFactorChanged)
        view.setZoomFactor(CCP4Modules.PREFERENCES().BROWSER_ZOOM_FACTOR)
        view.CustomMimeTypeRequested.connect(self.CustomMimeTypeRequested)
        #view.NewWindowRequested.connect(self.NewWindowRequested)
        view.StatusBarMessage.connect(self.StatusBarMessage)
        view.titleChanged.connect(self.setTabTitles)
        view.IconReady.connect(self.IconReady)
        self.newTab(view, "New tab")
        view.loadFinished.connect(self.loadFinished)
        return view

    def actionDefinition(self, actionName=''):
        return self.actionDefinitions.get(actionName,{})

    def setActionDefinition(self,actionName='',actionDefinition={}):
        self.actionDefinitions[actionName] = actionDefinition

    def getActionDef(self,name):
        return self.actionDefinitions.get(name,dict(text=name))

    def openBrowserWindow(self):
        b = CBrowserWindow()   # KJS : Broken ? 
        if hasattr(self,"handleFindNext"):
            b.loadFinished.connect(self.handleFindNext)
        #print 'openBrowserWindow',b,len(CBrowserWindow.Instances)
        b.show()
        #print 'openBrowserWindow raise_'
        b.raise_()

    def openI1Projects(self, mode='default'):
        if mode == 'default':
            self.openI1Projects1()
        else:
            fb = QtWidgets.QFileDialog()
#Not possible with native browser as far as I know. SJM 22/11/2018.
            if not CCP4Modules.PREFERENCES().NATIVEFILEBROWSER:
                fb.setOption(QtWidgets.QFileDialog.DontUseNativeDialog)
            fb.setFileMode(QtWidgets.QFileDialog.ExistingFile)
            fb.setFilter(QtCore.QDir.AllEntries | QtCore.QDir.Hidden | QtCore.QDir.NoDotAndDotDot)
            
            if fb.exec_():
                print(fb.selectedFiles())
                self.openI1Projects1(fb.selectedFiles()[0])
            
    @QtCore.Slot(str)
    def openI1Projects1(self, fileName=None):
        from qtgui import CCP4I1Projects
        from core import CCP4Utils
        if fileName is None:
            if sys.platform == "win32":
                fileName = os.path.normpath(os.path.join(CCP4Utils.getHOME(), 'CCP4', 'windows', 'directories.def'))
            else:
                fileName = os.path.normpath(os.path.join(CCP4Utils.getHOME(), '.CCP4', 'unix', 'directories.def'))
        for pV in CCP4I1Projects.CI1ProjectViewer.Instances:
            if pV.model().sourceFile == fileName:
                pV.show()
                pV.raise_()
                return
        pV = CCP4I1Projects.CI1ProjectViewer(fileName=fileName)
        CCP4I1Projects.CI1ProjectViewer.Instances[-1].show()
        CCP4I1Projects.CI1ProjectViewer.Instances[-1].raise_()
    
    def toolBar(self,name):
        return self.findChild(QtWidgets.QToolBar,name)

    def openProject(self,projectId,projectName=None):
        #print 'CMainWindow.openProject',projectId,projectName
        from qtgui import CCP4ProjectViewer
        for window in CCP4ProjectViewer.CProjectViewer.Instances:
            #print 'currently open',window,window.taskFrame.openJob.projectId,type(window.taskFrame.openJob.projectId)
            try:
                if window.taskFrame.openJob.projectId == projectId:
                    window.show()
                    window.raise_()
                    return
            except Exception as e:
                #print 'CMainWindow.openProject error',str(e)
                pass
        projectDir = CCP4Modules.PROJECTSMANAGER().db().getProjectInfo(projectId=projectId,mode='projectdirectory')
        if not os.path.exists(projectDir):
            QtWidgets.QMessageBox.warning(self,'','Error in opening project '+str(projectName)+'\nProject directory does not exist:\n'+str(projectDir))
            return
        if not os.path.exists(os.path.join(projectDir,'CCP4_JOBS')):
            QtWidgets.QMessageBox.warning(self,'','Error in opening project '+str(projectName)+'\nProject directory does not exist:\n'+os.path.join(projectDir,'CCP4_JOBS'))
            return
        p = CCP4ProjectViewer.CProjectViewer(projectId=projectId)
        p.resize(CCP4ProjectViewer.DEFAULT_WINDOW_SIZE[0],CCP4ProjectViewer.DEFAULT_WINDOW_SIZE[1])
        p.show()

    def openWorkflow(self):
        from core import CCP4WorkflowManagerGui
        CCP4WorkflowManagerGui.openWorkflowManagerGui()

    def openPatchComFile(self):
        from core import CCP4ComFilePatchManagerGui
        CCP4ComFilePatchManagerGui.openGui()

    def openCustomTasks(self):
        from qtgui import CCP4CustomTaskManagerGui
        CCP4CustomTaskManagerGui.openGui()

    def openImportJobs(self):
        from qtgui import CCP4ImportedJobManagerGui
        CCP4ImportedJobManagerGui.openGui()

    def handleProjectMenu(self, mode=''):
        from qtgui import CCP4ProjectManagerGui
        #print 'handleprojectMenu',mode
        if mode == 'new_project':
            if CCP4ProjectManagerGui.CNewProjectGui.insts is None:
                CCP4ProjectManagerGui.CNewProjectGui.insts = CCP4ProjectManagerGui.CNewProjectGui()
            CCP4ProjectManagerGui.CNewProjectGui.insts.clear()
            CCP4ProjectManagerGui.CNewProjectGui.insts.show()     # KJS : PyDev indicating no clear/show. Investigate this.
        elif mode == 'manage_projects':
            if CMainWindow.projectManagerDialog is None:
                CMainWindow.projectManagerDialog = CCP4ProjectManagerGui.CProjectManagerDialog()
            CMainWindow.projectManagerDialog.setMode('all')
            CMainWindow.projectManagerDialog.show()
            CMainWindow.projectManagerDialog.raise_()
        elif mode == 'more_projects':
            if CMainWindow.projectManagerDialog is None:
                CMainWindow.projectManagerDialog = CCP4ProjectManagerGui.CProjectManagerDialog()
            CMainWindow.projectManagerDialog.setMode('open')
            CMainWindow.projectManagerDialog.show()
            CMainWindow.projectManagerDialog.raise_()
        elif mode == 'import_project':
            if CMainWindow.projectManagerDialog is None:
                CMainWindow.projectManagerDialog = CCP4ProjectManagerGui.CProjectManagerDialog()
                CMainWindow.projectManagerDialog.hide()
            CMainWindow.projectManagerDialog.handleImportProject()
        elif mode == 'general_project':
            self.makeGeneralProject()

    def makeGeneralProject(self):
        from core import CCP4Utils
        name = 'General'
        status = CCP4Modules.PROJECTSMANAGER().projectStatus(name)
        if status == 0:
            QtWidgets.QMessageBox.information(self,'General project exists','General project already exists - use Project menu to open it')
            return
        elif status != 1:
            # Its in the database but broken in some way..
            QtWidgets.QMessageBox.warning(self,'General project exists','General project exists but may be corrupted - use Manage Projects tool from the Projects menu')
            return
        else:
            directory = os.path.join(CCP4Utils.getHOME(),'CCP4_General_Project')
            print('Making General Project in directory:',directory)
            try:
                projectId = CCP4Modules.PROJECTSMANAGER().createProject(projectName=name,projectPath=directory)
            except CException as e:
                e.warningMessage()
                return
            except:
                QtWidgets.QMessageBox.warning(self, 'Error creating project', 'Unknown error creating project')
                return
            else:
                QtWidgets.QMessageBox.information(self, 'General project created', 'General project created in directory\n' + directory)
                self.openProject(projectName='General')

    def openPreferences(self):
        if self.preferenceWindow is None:
            from qtgui import CCP4PreferencesGui
            self.preferenceWindow = CCP4PreferencesGui.CPreferencesWindow(self)
        #print 'CMainWindow.openPreferences',self.preferenceWindow
        self.preferenceWindow.show()
        self.preferenceWindow.raise_()

    def isPreferencesOpen(self):
        if self.preferenceWindow is not None and self.preferenceWindow.isVisible():
            return True
        else:
            return False

    def openConfigurePDBREDOTokens(self):
        dialog = QtWidgets.QDialog()
        dialog.setWindowTitle("Generate PDB_REDO tokens")
        layout = QtWidgets.QGridLayout()
        dialog.setLayout(layout)
        layout.addWidget(QtWidgets.QLabel("Please enter your PDB-REDO web services token and secret required by the PDB_REDO task"),0,0,1,2)
        layout.addWidget(QtWidgets.QLabel("PDB-REDO token id:"))
        tokenIdEdit = QtWidgets.QLineEdit()
        layout.addWidget(tokenIdEdit,layout.rowCount()-1,1)
        layout.addWidget(QtWidgets.QLabel("PDB-REDO token secret:"))
        secretEdit = QtWidgets.QLineEdit()
        layout.addWidget(secretEdit,layout.rowCount()-1,1)
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok|QtWidgets.QDialogButtonBox.Cancel)
        layout.addWidget(buttons,layout.rowCount(),0,1,2)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        ret = dialog.exec_()
        if ret == QtWidgets.QDialog.Accepted:
            if len(tokenIdEdit.text()) > 0 and len(secretEdit.text()) > 0:
                    CCP4Modules.PREFERENCES().PDB_REDO_TOKEN_ID.set(tokenIdEdit.text())
                    CCP4Modules.PREFERENCES().PDB_REDO_TOKEN_SECRET.set(secretEdit.text())
                    QtWidgets.QMessageBox.information(self,'Successfully set PDB_REDO token','Successfully set PDB_REDO token<br/>You should now be able use the PDB_REDO task.'.format(payload.get('expires')),QtWidgets.QMessageBox.Ok)


    def openServerSetup(self):
        if getattr(self,'serverSetupWindow',None) is None:
            from qtgui import CCP4JobControlGui
            self.serverSetupWindow = CCP4JobControlGui.CServerSetupWindow(self)
        #print 'CMainWindow.openPreferences',self.preferenceWindow
        self.serverSetupWindow.show()
        self.serverSetupWindow.raise_()

    def isServerSetupOpen(self):
        if getattr(self,'serverSetupWindow',None) is not None and self.serverSetupWindow.isVisible():
            return True
        else:
            return False

    def openConfig(self):
        if self.configWindow is None:
            from core import CCP4ConfigGui
            self.configWindow = CCP4ConfigGui.CConfigWindow(self)
        #print 'CMainWindow.openConfig',self.configWindow
        self.configWindow.show()
        self.configWindow.raise_()

    def isConfigOpen(self):
        if self.configWindow is not None and self.configWindow.isVisible():
            return True
        else:
            return False

    def openFileDialog(self):
        from qtgui import CCP4FileBrowser
        filter_list = []
        mimeTypeHandler = CCP4Modules.MIMETYPESHANDLER()
        if self.fileDialog is None:
            filter_list.append('All files (*.*)')
            self.fileDialog = CCP4FileBrowser.CFileDialog(self, title='Open file', filters=filter_list, addAll=False)
            self.fileDialog.selectFile.connect(self.handleOpenFileSelection)
        self.fileDialog.show()

    @QtCore.Slot(str)
    def handleOpenFileSelection(self, fileName=''):
        #print 'CCP4BrowserWindow.handleOpenFileSelection',fileName
        if isinstance(self, CBrowserWindow):
            self.openFile(fileName, internal=True)
        else:
            OPENFILE(fileName)

    def handleDeveloperTools(self,mode):
        print("########################################")
        print("########################################")
        print("CWebBrowser thread",QtCore.QThread.currentThread())
        print("########################################")
        print("########################################")
        #print 'handleDeveloperTools',mode
        if mode == 'redo_report':
            self.redoReport()
        elif mode == 'view_report_source':
            try:
                openJob = self.taskFrame.openJob
            except:
                return
            fileName = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=openJob.jobId,mode='REPORT')
            if os.path.exists(fileName):
                CCP4Modules.WEBBROWSER().openFile(fileName=fileName,format='text/plain')
        elif mode == 'remake_cached_lookups':
            try:
                from core.CCP4Modules import TASKMANAGER, PIXMAPMANAGER
                from core.CCP4DataManager import DATAMANAGER
                print("Off to remake CTaskManager cache...")
                TASKMANAGER().buildLookupFromScratch()
                print("...Completed")
                print("Off to remake CDataManager cache...")
                DATAMANAGER().buildClassLookupFromScratch()
                print("...Completed")
                print("Off to remake CPixmapManager cache...")
                PIXMAPMANAGER().buildCacheFromScratch()
                print("...Completed")
            except:
                print("Failed to remake caches")
        elif mode == 'view_test_report':
            testReportList = glob.glob(os.path.join(CCP4Modules.PROJECTSMANAGER().db().getProjectInfo(projectId=self.openJob.projectId,mode='projectdirectory'),'CCP4_test*.log'))
            if len(testReportList) > 0:
                CCP4Modules.WEBBROWSER().openFile(fileName=testReportList[0])
            return
        elif mode == 'make_test_suite':
            from core import CCP4ProjectBasedTesting
            testBuilder = CCP4ProjectBasedTesting.BuildTestSuite(projectId=self.getProject())  
            rv = testBuilder.run()
            if rv.maxSeverity() > SEVERITY_WARNING:
                rv.warningMessage(parent=self,windowTitle=self.windowTitle(),message='Error converting project to test suite')
            else:
                builder = CCP4ProjectBasedTesting.BuildTestSuite(projectId=self.getProject())
                try:
                    err = builder.run()
                except:
                    err.warningMessage()
                else:
                    projDir = CCP4Modules.PROJECTSMANAGER().db().getProjectDirectory(projectId=self.getProject())
                    self.testSuiteDialog = QtWidgets.QDialog(self)
                    self.testSuiteDialog.setWindowTitle('Creating a test suite')
                    self.testSuiteDialog.setLayout(QtWidgets.QVBoxLayout())
                    label = QtWidgets.QLabel("""A directory for test definitions has been created at\n""" +
                                         str(projDir)+'/CCP4_TEST_SYSTEM\n' +
                                         """You should edit in your tests and then 'Export' this project from 'Manage Projects'\n""" +
                                         """You can then run the test with 'testi2sys exported_file_path'""", self)
                    self.testSuiteDialog.layout().addWidget(label)
                    self.testSuiteDialog.show()
                    self.testSuiteDialog.raise_()
        elif mode == 'auto_task_docs':
            from core import CCP4TaskManager
            rv = CCP4TaskManager.CMakeDocsIndex().run()
            if rv.maxSeverity() > SEVERITY_WARNING:
                rv.warningMessage(parent=self,windowTitle=self.windowTitle(),message='Error recreating task documentation index page')

    def listTasks(self):
        import tempfile
        from core import CCP4TaskManager
        from core import CCP4Utils
        text = CCP4TaskManager.LISTTASKS(ifPrint=False)
        fileName = tempfile.mktemp(suffix='.txt')
        CCP4Utils.saveFile(fileName,text)
        widget = CCP4Modules.WEBBROWSER().openFile(fileName)
        widget.setFont(style='fixed_width')

    def openExportTask(self):
        from core import CCP4Utils
        import zipfile
        dirPath = QtWidgets.QFileDialog.getExistingDirectory(self,'Select pipeline directory to save as task compressed file',os.path.join(CCP4Utils.getCCP4I2Dir(),'pipelines')).__str__()
        if not os.path.isdir(dirPath):
            return
        if not os.path.relpath(dirPath,CCP4Utils.getCCP4I2Dir()).startswith('pipelines'):
            err = CException(self.__class__,209,dirPath)
            err.warningMessage(parent=self,windowTitle='Error creating compressed task file',message='Selected directory is not a pipeline in the currently running cccp4i2')
            return
        zipPath = dirPath + '.ccp4i2task.zip'
        if os.path.exists(zipPath):
            query = QtWidgets.QMessageBox.question(self,'Overwrite task compressed file?','Overwrite existing'+zipPath+'?',QtWidgets.QMessageBox.Cancel|QtWidgets.QMessageBox.Yes)
            if query == QtWidgets.QMessageBox.Cancel:
                return
            else:
                os.remove(zipPath)
        try:
            zip = zipfile.ZipFile(zipPath, mode='w')
        except:
            err = CException(self.__class__, 201, zipPath)
            err.warningMessage(parent=self,windowTitle='Error creating compressed task file',message='Creating'+str(zipPath))
            return
        try:
            CCP4Utils.zipDirectory(zip,dirPath,rootRelPath=CCP4Utils.getCCP4I2Dir())
        except:
            err = CException(self.__class__, 202, 'Saving', dirPath, 'to', zipPath)
            err.warningMessage(parent=self,windowTitle='Error creating compressed task file',message='Saving'+str(dirPath))
            return
        try:
            zip.close()
        except:
            return CErrorReport(self.__class__,203,zipPath)
            err.warningMessage(parent=self,windowTitle='Error creating compressed task file',message='Closing'+str(zipPath))
            return # ? Extra return here.
        info = QtWidgets.QMessageBox.information(self,'Saved compressed task file','Task saved to'+str(zipPath))

    def openImportTask(self):
        import shutil
        import zipfile
        from core import CCP4Utils
        filePath,selectedFilter = QtWidgets.QFileDialog.getOpenFileName(self,'Select compressed file containing task','',"Compressed task (*.ccp4i2task.zip)").__str__()
        if not os.path.isfile(filePath): return
        try:
            zip = zipfile.ZipFile(filePath,mode='r')
        except:
            err = CException(self.__class__,204,filePath)
            err.warningMessage(parent=self,windowTitle='Error reading compressed task file',message='Reading '+str(filePath))
            return
        targetPath = CCP4Utils.getCCP4I2Dir()
        err= CException()
        for zinfo in zip.infolist():
            if not zinfo.filename.startswith('pipelines'):
                err.append(self.__class__,206,str(zinfo.filename))
            else:
                if os.path.exists(os.path.join(targetPath,zinfo.filename)):
                    overwrite = QtWidgets.QMessageBox.question(self,'Overwrite existing task?','Overwrite existing '+str(zinfo.filename)+'?',QtWidgets.QMessageBox.Cancel|QtWidgets.QMessageBox.Yes)
                    if overwrite == QtWidgets.QMessageBox.Cancel:
                        err.append(self.__class__,207,zinfo.filename)
                        break
                    else:
                        try:
                            shutil.rmtree(os.path.join(targetPath,zinfo.filename))
                        except:
                            err.append(self.__class__,208,os.path.join(targetPath,zinfo.filename))
                            break
                try:
                    zip.extract(zinfo,targetPath)
                except:
                    err.append(self.__class__,205,filePath)
        zip.close()
        if len(err) > 0:
            err.warningMessage(parent=self,windowTitle='Error reading compressed task file',message='Extracting from '+str(filePath)+' to '+str(targetPath))
            return
        info = QtWidgets.QMessageBox.information(self,'New task installed','New '+ zinfo.filename+ ' task installed - Please restart CCP4i2')

    def showAbout(self):
        if not hasattr(self,'aboutDialog'):
            from core import CCP4Utils
            self.aboutDialog = QtWidgets.QDialog(self)
            self.aboutDialog.setLayout(QtWidgets.QVBoxLayout())
            self.aboutDialog.layout().setContentsMargins(1,1,1,1)
            #self.aboutDialog.layout().setSpacing(1)
            self.aboutDialog.setWindowTitle('About CCP4i2')
            topWidget = QtWidgets.QWidget()
            topWidget.setLayout(QtWidgets.QHBoxLayout())
            topWidget.setStyleSheet("QWidget { background-color:white; }")
            label =  QtWidgets.QLabel(self)
            fileName = os.path.normpath(os.path.join(CCP4Utils.getCCP4I2Dir(),'qticons','ccp4.png'))
            label.setPixmap(QtGui.QPixmap(fileName))
            topWidget.layout().addWidget(label)
            topRightWidget = QtWidgets.QWidget()
            topRightWidget.setLayout(QtWidgets.QVBoxLayout())
            label = QtWidgets.QLabel('CCP4i2',self)
            label.setStyleSheet("QLabel { font-size: 24px; font-weight: bold; }")
            topRightWidget.layout().addWidget(label)
            topWidget.layout().addWidget(topRightWidget)
            self.aboutDialog.layout().addWidget(topWidget)
            version = CCP4Utils.getProgramVersion('ccp4i2')
            date = CCP4Utils.getProgramVersion('ccp4i2',mode='date')
            label = QtWidgets.QLabel('Version '+version+' built on '+date+"\n"+'User interface to CCP4 Program Suite version '+ CCP4Utils.getProgramVersion('ccp4'),self)
            label.setStyleSheet("QLabel { font-size: 14px; font-style: italic; font-weight: bold; }")
            topRightWidget.layout().addWidget(label)
            label = QtWidgets.QLabel(self)
            label.setText('Copyright (C) 2001-2008 University of York, CCLRC\n' +
                          'Copyright (C) 2009-2010 University of York\n' +
                          'Copyright (C) 2011-2018 Science & Technology Facilities Council')
            bottomWidget = QtWidgets.QWidget()
            bottomWidget.setLayout(QtWidgets.QVBoxLayout())
            bottomWidget.layout().addWidget(label)
            licenseLabel = QtWidgets.QLabel(self)
            licenseLabel.setText('<a href="http://www.ccp4.ac.uk/download/licence.php">CCP4 License</a>')
            licenseLabel.setTextFormat(QtCore.Qt.RichText);
            licenseLabel.setTextInteractionFlags(QtCore.Qt.TextBrowserInteraction);
            licenseLabel.setOpenExternalLinks(True);
            bottomWidget.layout().addWidget(licenseLabel)
            self.aboutDialog.layout().addWidget(bottomWidget)

        self.aboutDialog.show()
        self.aboutDialog.raise_()

    def showDemoDataInfo(self):
        from core import CCP4Utils
        CCP4Modules.WEBBROWSER().loadWebPage(CCP4Modules.DEMODATAMANAGER().getOverviewPage())

    def downloadDemoData(self):
        CCP4Modules.DEMODATAMANAGER().downloadDemoData(self)

    def getProject(self):
        return None

    def isProjectViewer(self):
        return (self.getProject() is not None)

    def isNotProjectViewer(self):
        return True

    def reportAvailable(self):
        return False

    def testReportAvailable(self):
        return False

    def grabWidget(self):
        # Dummy method reimplemented in CProjectViewer
        pass
    
    def makeDemoData(self):
        CCP4Modules.DEMODATAMANAGER().makeDemoData(parentWidget=self)
    
    def listProcesses(self):
        CCP4Modules.JOBCONTROLLER().listLocalProcesses()
        from qtgui import CCP4JobControlGui
        if getattr(self,'listProcessesWindow',None) is None:
            self.listProcessesWindow = CCP4JobControlGui.CListProcesses(self)
        self.listProcessesWindow.load()
        self.listProcessesWindow.show()
        self.listProcessesWindow.raise_()

class TwoFingerFlickGesture(QtWidgets.QGesture):

    def __init__(self,parent=None):
        QtWidgets.QGesture.__init__(self,parent)

class TwoFingerFlickRecognizer(QtWidgets.QGestureRecognizer):

    def create(self,target):
        return TwoFingerFlickGesture()

    def recognize(self,state,watched,event):
        import time
        if event.type() == QtCore.QEvent.TouchEnd and self.inGesture:
            #print self.diff,(time.time() - self.startTime), self.diff/(time.time() - self.startTime)
            if abs(self.diff/(time.time() - self.startTime))>300:
                watched.flick.emit(self.diff)
            self.inGesture = False
            self.diff = 0.0
            return QtWidgets.QGestureRecognizer.FinishGesture
        if event.type() == QtCore.QEvent.TouchBegin:
            self.diff = 0.0
            if len(event.touchPoints()) == 2:
                self.inGesture = True
            else:
                self.inGesture = False
            self.startTime = time.time()
            return QtWidgets.QGestureRecognizer.TriggerGesture
        if event.type() == QtCore.QEvent.TouchUpdate:
            if len(event.touchPoints()) == 2:
                self.diff += event.touchPoints()[0].pos().x() -  event.touchPoints()[0].lastPos().x()
                self.inGesture = True
                return QtWidgets.QGestureRecognizer.TriggerGesture
        return QtWidgets.QGestureRecognizer.Ignore

class Flicker:

    def __init__(self):
        fftRecognizer = TwoFingerFlickRecognizer()
        self.tffType = QtWidgets.QGestureRecognizer.registerRecognizer(fftRecognizer)
        self.grabGesture(QtCore.Qt.PanGesture);
        self.grabGesture(QtCore.Qt.SwipeGesture);
        self.grabGesture(self.getTFFType())

    def getTFFType(self):
        return self.tffType

    def gestureEvent(self,event):
        tff = event.gesture(self.tffType)
        if tff:
            event.accept(tff);

    def event(self,e):
        if e.type() == QtCore.QEvent.Gesture:
            self.gestureEvent(e)
        return QtWidgets.QWidget.event(self,e)


#class CBrowserWindow(CMainWindow,Flicker):
class CBrowserWindow(CMainWindow):

    ERROR_CODES = {100 : {'description' : 'Invalid file name'}}
    Instances = set()
    Dummy = None

    flick = QtCore.Signal(float)
    #windowMoved = QtCore.Signal('QPos')
    loadFinished = QtCore.Signal()

    def __init__(self, parent=None, welcome=True):
        self.stack = None
        CMainWindow.__init__(self,parent)
        #Flicker.__init__(self)
        self.setWindowTitle(self.version)
        CBrowserWindow.Instances.add(self)
        self.stack = QtWidgets.QStackedWidget(self)
        #print 'CBrowserWindow.__init__ parent',self,parent,CBrowserWindow.Instances
        self.setCentralWidget(self.stack)
        self.destroyed.connect(CBrowserWindow.updateInstances)
        self.factory = None
        #QtWebKit.QWebSettings.globalSettings().setAttribute(QtWebKit.QWebSettings.PluginsEnabled, True)
        # Create menu and tool bar and add the 'general' actions
        self.mainToolBar = CToolBar(self,'main','Main toolbar')
        self.fileToolBar = CToolBar(self,'file','File')
        self.addToolBar(self.mainToolBar)
        self.addToolBarBreak()
        self.addToolBar(self.fileToolBar)
        if self.mainToolBar is not None:
            self.mainToolBar.extend(['back','forward','reload','find','sep','help_ccp4i2','open_tab'])
        self.mainToolBar.setMovable(False)
        self.fileToolBar.setMovable(False)
        self.mini = False

        self.editSplitter = CEditSplitter(None)
        self.editSplitter.addressEdit.returnPressed.connect(self.addressEdited)
        self.editSplitter.searchEdit.returnPressed.connect(self.searchEdited)

        if sys.platform == "darwin":
            self.mainToolBar.hide()
            self.webviewToolBar = CCP4WebToolBarButtons.CCP4WebToolBarButtons(fileName=os.path.join(os.path.dirname(__file__),"html_i2_browser_buttons.html"))
            dockWidget = QtWidgets.QDockWidget()
            dockWidget.setTitleBarWidget(QtWidgets.QWidget());

            dockLayoutWidget = QtWidgets.QWidget()
            dockLayout = QtWidgets.QVBoxLayout()
            dockLayout.setContentsMargins(0,0,0,0)
            dockLayoutWidget.setLayout(dockLayout)
            dockLayout.addWidget(self.webviewToolBar)
            dockLayout.addWidget(self.editSplitter)
            
            dockWidget.setWidget(dockLayoutWidget)

            self.addDockWidget(QtCore.Qt.TopDockWidgetArea,dockWidget)
            dockWidget.setFixedHeight(103)    # Don't want any scroolbars or resize decorations.
            def handleToolBarClick(buttonName):
                if self.findChild(QtWidgets.QAction, buttonName):
                    self.findChild(QtWidgets.QAction, buttonName).trigger()
            self.webviewToolBar.buttonClicked.connect(handleToolBarClick)
            self.webviewToolBar.loadFinished.connect(functools.partial(self.updateActionEnabled,None))
            self.webviewToolBar.setContextMenuPolicy(QtCore.Qt.NoContextMenu)

            self.toolBar('file').hide()
        else:
            self.toolBar('file').addWidget(self.editSplitter)

        self.setTab('nocontext')
        self.setHistoryActionAvailability()
        self.tab().currentChanged[int].connect(self.setTabTitles)
        sb = self.statusBar()
        self.searchEngineString = "http://www.google.com/search?q="
        self.drawFindTool()
        self.defaultSizePolicy = [self.sizePolicy().horizontalPolicy(),self.sizePolicy().verticalPolicy()]
        if welcome or self.tab().count() == 0:
            docstart = os.path.join('sphinx', 'build', 'html', 'index.html')
            self.loadWebPage(helpFileName=docstart, newTab=True)
        self.resize(960, 800)

        def handleFlick(flick):
           if flick < 0.0:
               self.historyBack()
           else:
               self.historyForward()

        self.flick.connect(handleFlick)

    @QtCore.Slot()
    @staticmethod
    def updateInstances(qobj):
        CBrowserWindow.Instances = set([window for window in CBrowserWindow.Instances if isAlive(window)])
        print('webbrowser.updateInstances',CBrowserWindow.Instances)

    def closeEvent(self,e):
        CBrowserWindow.Instances = []

    def close(self):
        if self.tab().count()==1:
            return CMainWindow.close(self)
        self.tab().closeTab(self.tab().currentIndex())
        if self.tab().count()>0:
            return False
        return CMainWindow.close(self)

    def setMini(self,mode=True):
        self.mini = mode
        tabs = self.findChildren(CTabWidget)
        if self.mini:
            self.mainToolBar.hide()
            self.fileToolBar.hide()
            self.resize(500,300)
            self.layout().setContentsMargins(0,0,0,0)
            self.layout().setSpacing(0)
            for tab in tabs:
                tab.cornerWidget().setText('Expand window')
        else:
            self.mainToolBar.show()
            self.fileToolBar.show()
            self.resize(800,800)
            for tab in tabs:
                tab.cornerWidget().setText('Close tab')

    def positionOver(self,widget):
        geom = widget.geometry()
        globalPos = widget.mapToGlobal(QtCore.QPoint(geom.x(),geom.y()))
        self.move(globalPos)
    
    def downloadTestData(self,**kw):
        pass
    
    def drawFindTool(self):
        self.findFrame = CFindFrame(self)
        self.findDock = QtWidgets.QDockWidget('Find',self)
        self.findDock.setWidget(self.findFrame)
        self.findDock.setAllowedAreas(QtCore.Qt.TopDockWidgetArea|QtCore.Qt.BottomDockWidgetArea)
        self.findDock.setTitleBarWidget(None)
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea,self.findDock,QtCore.Qt.Horizontal)
        self.findDock.close()
        self.findFrame.findNext.connect(self.handleFindNext)
        self.findFrame.findPrevious.connect(self.handleFindPrevious)
        #self.findFrame.highlightAll.connect(self.handleFindHighlightAll)

    def StatusBarMessage(self,message=None):
        sb = self.statusBar()
        sb.showMessage(message)

    def NewWindowRequested(self,view=None):
        if view:
            self.tab().addTab(view,view.title())
            view.CustomMimeTypeRequested.connect(self.CustomMimeTypeRequested)
            view.NewWindowRequested.connect(self.NewWindowRequested)
            view.StatusBarMessage.connect(self.StatusBarMessage)
            view.titleChanged.connect(self.setTabTitles)
            view.IconReady.connect(self.IconReady)
            self.tab().setCurrentWidget(view)

    def setTabTitles(self,title=None):
        #print 'setTabTitles',title
        for i in range(self.tab().count()):
            widget = self.tab().widget(i)
            if hasattr(widget,"title") and isinstance(widget.title, Callable):
                self.tab().setTabText(i,widget.title())
            self.setWindowTitle(self.version+self.tab().tabText(self.tab().currentIndex()))
            widget = self.tab().widget(self.tab().currentIndex())
            if hasattr(widget,"url") and isinstance(widget.url, Callable):
                self.editSplitter.addressEdit.setText(widget.url().toString())
        self.setHistoryActionAvailability()

    def tab(self):
        # Returns the CTabWidget for the current context (ie project)
        if self.stack is None:
            return None
        return self.stack.currentWidget()

    def contextTab(self,context):
        for indx in range(self.stack.count()):
            if str(self.stack.widget(indx).objectName()) == context:
                return self.stack.widget(indx)
        return None

    def setTab(self,context):
        tab = self.stack.currentWidget()
        if tab is not None and str(tab.objectName()) == context:
            return
        for indx in range(self.stack.count()):
            if str(self.stack.widget(indx).objectName()) == context:
                self.stack.setCurrentIndex(indx)
                return
        tab = CTabWidget(self,name=context,mini=self.mini)
        tab.currentChanged.connect(self.loadFinished)
        indx = self.stack.addWidget(tab)
        self.stack.setCurrentIndex(indx)
        #print 'setTab new',indx,tab,context
      
#-------------------------------------------------------------------
    def newTab(self, widget=None, title=None, context='nocontext', toolTip=None, icon=None, copen=True, url=None):
#-------------------------------------------------------------------
        if title is None:
            title = widget.title()
        #print 'CWebBrowser.newTab title',widget,title
        if icon is not None:
            idx = self.tab().addTab(widget,icon,title)
        else:
            idx = self.tab().addTab(widget,title)
        if toolTip is not None:
            self.tab().setTabToolTip(idx,toolTip)
        #if colour: self.tab().setTabTextColor(QtGui.QColor(colour))
        if copen:
            self.setTabIndex(idx)
            if url is not None:
                self.editSplitter.addressEdit.setText(url.toString())
            else:
                self.editSplitter.addressEdit.setText(title)
        self.tab().repaint()
        return idx

    def handleTabCornerWidget(self,clickBool = None):
        if self.mini:
            self.setMini(False)

    def tabWidgets(self):
        widgets = []
        for idx in range(0,self.tab().count()):
            widgets.append(self.tab().widget(idx))
        return widgets

#-------------------------------------------------------------------
    def setTabIndex(self, index=-1, fileName=''):
#-------------------------------------------------------------------
        if fileName and index < 0:
            index = self.fileOpenInTab(fileName)
        if index >= 0 and index <= self.tab().count():
            self.tab().setCurrentIndex(index)
            self.setWindowTitle(self.version+self.tab().widget(index).objectName())

#-------------------------------------------------------------------
    def fileOpenInTab(self, fileName):
#-------------------------------------------------------------------
        for idx in range(0,self.tab().count()):
            if str(self.tab().tabText(idx)) == fileName:
                return idx
        fileName = os.path.abspath(str(fileName))
        for idx in range(0,self.tab().count()):
            if self.tab().widget(idx).fileName == fileName:
                return idx
        return -1

    @QtCore.Slot('QUrl')
    def CustomMimeTypeRequested(self, url=None):
        #print 'CWebBrowser.CustomMimeTypeRequested',url.path(),url.scheme()
        if url is not None and url.scheme() == 'file':
            path = str(url.path())
            if os.path.exists(path):
                self.openFile(path, internal=True)
            elif os.path.splitext(path)[1] == '.i2com':
                com = os.path.split(os.path.splitext(path)[0])[1]
                #print 'CustomMimeTypeRequested i2com',com
                action = self.findChild(QtWidgets.QAction,com)
                if action is not None:
                    action.trigger()
                elif com in self.actionDefinitions and 'slot' in self.actionDefinitions[com]:
                    #print 'CustomMimeTypeRequested try',self.actionDefinitions[com]['slot']
                    try:
                        self.actionDefinitions[com]['slot']()
                    except:
                        pass

#-------------------------------------------------------------------
    def openApplication(self, application=None):
#-------------------------------------------------------------------
        from qtgui import CCP4TextViewer
        from qtgui import CCP4ErrorReportViewer
        widget = None
        application = str(application)
        idx = self.fileOpenInTab(application)
        #print 'openApplication',application,idx
        if idx >= 0:
            self.tab().setCurrentIndex(idx)
            #self.tab().widget(idx).open(self.tab().widget(idx).filename)
            return self.tab().widget(idx)
        if application == 'programPrintLog':
            widget = CCP4ErrorReportViewer.CPrintLogViewer(self)
            #print 'CWebBrowser.openApplication',widget
            widget.openThread(thread='main_thread')
            title = 'Program print log'
        elif application == 'programHTTPLog':
            widget = CCP4ErrorReportViewer.CPrintLogViewer(self)
            #print 'CWebBrowser.openApplication',widget
            widget.openThread(thread='HTTPServer')
            title = 'Program print log'
        elif  application == 'editScript':
            widget = CCP4TextViewer.CScriptViewer(self)
            widget.open()
            title= 'Edit script'
        if widget is not None:
            widget.setObjectName(title)
            self.newTab(widget,title)
            widget.show()
        return widget

#-------------------------------------------------------------------
    def openFile(self, fileName=None, format=None, title=None, toFront=False, internal=False):
#-------------------------------------------------------------------
        from core import CCP4Config
        # Is the file already displayed - force a redraw
        if fileName is None or fileName=='None':
            return None
        fileName = str(fileName)
        idx = self.fileOpenInTab(fileName)
        if idx >= 0:
            if self.tab().widget(idx).isFileModified():
                self.tab().widget(idx).reload()
            self.tab().setCurrentIndex(idx)
            if toFront: self.raise_()
            return self.tab().widget(idx)
        mimeTypeHandler = CCP4Modules.MIMETYPESHANDLER()
        #print "CCP4WebBrowser.open fileName",fileName,format,mimeTypeHandler,title
        if not os.path.exists(fileName):
            QtWidgets.QMessageBox.warning(self,self.windowTitle(),'File does not exist \n'+fileName)
            return None
        if os.path.isdir(fileName):
            format = 'dir'
        # If format not recognised try if file extension is recognised
        if mimeTypeHandler is None and format is None:
            # No handler so just throw it at web browser and hope
            format = 'text/html'
        else:
            if format is not None and not mimeTypeHandler.isSupportedFormat(format):
                format = None
            if format is None:
                format = mimeTypeHandler.formatFromFileExt(fileName)
            if format is None:
                format = 'text/html'
        #print 'CCP4WebBrowser.openFile format',format
        # If its a plugin format launch the plugin
        if ['text/html'].count(format):
            widget = self.loadWebPage(fileName,newTab=True)
        elif mimeTypeHandler.useDesktopServices(format):
            abs_fileName = os.path.abspath(fileName)
            #print 'calling QDesktopServices',abs_fileName
            url = QtCore.QUrl.fromLocalFile(abs_fileName)
            rv = QtGui.QDesktopServices.openUrl(url)
            if not rv:
                QtWidgets.QMessageBox.warning(None,self.windowTitle(),'Attempting to display file '+os.path.split(abs_fileName)[-1]+'\nusing desktop services failed')
            #.. and close this window if nothing else is displayed
            if not internal and self.stack.count() == 1:
                self.close()
            return None
        else:
            widgetClassList = mimeTypeHandler.getViewers(format)
            widgetClass = widgetClassList[0]
            if len(widgetClassList) == 0:
                QtWidgets.QMessageBox.warning(None,self.windowTitle(),'Sorry no way to display file type: '+format)
                return None
            else:
                #widgetClass = widgetClassList[0]
                #  PhilE hack to display non-MTZ files as text
                root, ext = os.path.splitext(str(fileName))
                if ext != '.mtz':
                    widgetClass = mimeTypeHandler.getViewers("text/plain")[0]
                         
            if isinstance(widgetClass,str):
                # This is a keyword for the launcher
                CCP4Modules.LAUNCHER().openInViewer(viewer=widgetClass,fileName=fileName)
                return None
            if CCP4Config.DEVELOPER():
                widget = widgetClass(self)
            else:
                try:
                    widget = widgetClass(self)
                except:
                    QtWidgets.QMessageBox.warning(None,self.windowTitle(),'Error opening display for file type: '+format)
                    return None
            #print 'CWebBrowser.openFile',widget
            if CCP4Config.DEVELOPER():
                rv = widget.open(fileName=fileName)
            else:
                try:
                    rv = widget.open(fileName=fileName)
                except CException as e:
                    widget.close()
                    e.warningMessage(windowTitle=self.windowTitle())
                except:
                    widget.close()
                    QtWidgets.QMessageBox.warning(None,self.windowTitle(),'Attempting to display file '+os.path.split(fileName)[-1]+'\n as '+format+' failed.\n')
                    return None
            if title is None:
                title = widget.title()
            self.newTab(widget,title,url=QtCore.QUrl.fromLocalFile(fileName),toolTip=fileName)
            widget.show()
            #print 'openFile size',widget.height(),widget.width()
            if toFront: self.raise_()
            return widget

    def reloadFile(self,fileName):
        # Beware this is used as slot by CCP4ProjectManagerGui so arguments are fixed
        #print 'WEBBROWSER.reloadFile',fileName
        if fileName is None or fileName=='None':
            return None
        fileName = str(fileName)
        idx = self.fileOpenInTab(fileName)
        if idx >= 0:
            if self.tab().widget(idx).isFileModified():
                self.tab().widget(idx).reload()
            self.tab().setCurrentIndex(idx)
            #if toFront: self.raise_()
            return self.tab().widget(idx)
        else:
            return None

    def showHelp(self,mode='ccp4i2',newTab=True):
        print('showHelp',mode)
        from core import CCP4Utils
        sph_pth = os.path.join(CCP4Utils.getCCP4I2Dir(), 'docs', 'sphinx', 'build', 'html')
        if mode == 'ccp4i2':
            url = QtCore.QUrl.fromLocalFile(os.path.join(sph_pth, 'index.html'))
            self.loadPage(url, newTab=newTab)
        elif mode == 'quickstart':
            url = QtCore.QUrl.fromLocalFile(os.path.join(sph_pth, 'quickstart', 'index.html'))
            print("Loading",url)
            self.loadPage(url, newTab=newTab)
        elif mode == 'quickexpert':
            url = QtCore.QUrl.fromLocalFile(os.path.join(sph_pth, 'quickexpert', 'index.html'))
            self.loadPage(url, newTab=newTab)
        elif mode == 'youtube':
            rv = QtGui.QDesktopServices.openUrl(QtCore.QUrl('https://www.youtube.com/watch?v=cxp0SGKkmG8'))
        elif mode == 'task_docs':
            url = QtCore.QUrl.fromLocalFile(os.path.join(sph_pth, 'tasks', 'index.html'))
            self.loadPage(url, newTab=newTab)
        elif mode == 'cloud_docs':
            url = QtGui.QDesktopServices.openUrl(QtCore.QUrl('http://ccp4serv6.rc-harwell.ac.uk/jscofe-dev/manuals/html-taskref/index.html'))
        elif mode == 'tutorial':
            url = QtCore.QUrl.fromLocalFile(os.path.join(sph_pth, 'general','tutorial.html'))
            self.loadPage(url, newTab=newTab)
        elif mode == 'ccp4_home':
            url =  QtCore.QUrl("http://www.ccp4.ac.uk")
            self.loadPage(url, newTab=newTab)
        elif mode == 'updates':
            url = QtCore.QUrl("http://www.ccp4.ac.uk/updates") 
            self.loadPage(url, newTab=newTab)
        elif mode == 'license':
            self.loadWebPage(helpFileName='general/license.html', newTab=newTab)

    def loadWebPage(self, fileName=None, helpFileName=None, target=None, newTab=0):
        from core import CCP4Utils
        #print 'loadWebPage',fileName,helpFileName,target
        if fileName is not None:
            if target is None and fileName.count('#'):
                fileName, target = fileName.split('#', 1)
        elif helpFileName is not None:
            if target is None and helpFileName.count('#'):
                helpFileName, target = helpFileName.split('#', 1)
            if os.path.splitext(helpFileName)[1] == '':
                helpFileName = helpFileName + '.html'
            if os.path.exists(helpFileName):
                fileName = helpFileName
            else:
                docDir = os.path.join(CCP4Utils.getCCP4I2Dir(), 'docs', 'sphinx', 'build' ,'html')
                fileName = os.path.join(docDir, helpFileName)
                idx = 0
                subDirList =  ['general', 'developers', 'tasks']
                while not os.path.exists(fileName) and idx < len(subDirList):
                    fileName = os.path.join(docDir, subDirList[idx], helpFileName)
                    idx += 1
                if idx > len(subDirList):
                    return None
        if target is not None:
            if target[0] != '#':
                target = '#' + target
        else:
            target = ''
        if os.path.exists(fileName):
            url = QtCore.QUrl.fromLocalFile(os.path.abspath(fileName))
            if len(target) > 1:
                url.setFragment(target[1:])
            #url = QtCore.QUrl(fileName+target)
        elif not fileName[0:4] == 'http':
            url = QtCore.QUrl("http://" + fileName)
        else:
            #url = QtCore.QUrl(file)
            url = file
        return self.loadPage(url, newTab=newTab)

    def loadPage(self, url, newTab=0):
        target = ''
        if url.scheme() == 'file':
            if url.path().find('#') > 0:
                urlpath, target = str(url.path().toAscii()).split('#', 1)
                url = QtCore.QUrl.fromLocalFile(urlpath)
                #print 'loadPage file',target
        if not url.isValid():
            err = CException(self.__class__, 100, url.path().toAscii())
            err.warningMessage()
            return
            #print "CWebBrowser.loadPage trying", url.path(), target
        if self.tab().count() > 0:
            widget = self.tab().currentWidget()
        else:
            widget= None
        if not newTab and isinstance(widget, CCP4WebView.CWebView):
            try:
                widget.setTarget(target)
                widget.load(url)
                #print 'CWebBrowser.loadPage widget.loaded'
            except:
                exc_type, exc_value, exc_tb = sys.exc_info()[:3]
                print(exc_type)
                print(exc_value)
                import traceback
                traceback.print_tb(exc_tb)
                return None
            return widget
        else:
            print(self.tab().currentWidget())
            view = CCP4WebView.CWebView()
            view.searchFound.connect(self.findFrame.searchCallback)
            view.zoomFactorChanged.connect(self.zoomFactorChanged)
            view.setZoomFactor(CCP4Modules.PREFERENCES().BROWSER_ZOOM_FACTOR)
            view.CustomMimeTypeRequested.connect(self.CustomMimeTypeRequested)
            #view.NewWindowRequested.connect(self.NewWindowRequested)
            view.StatusBarMessage.connect(self.StatusBarMessage)
            view.titleChanged.connect(self.setTabTitles)
            view.IconReady.connect(self.IconReady)
            view.load(url)
            view.setTarget(target)
            view.loadFinished.connect(self.loadFinished)

            #FIXME - Do not know why I have to do things differently. Why do new tabs get triggered anyway when address changes.
            #Added a "True" as this test was clearly broken. Now confused. SJM 9/5/2017
            if newTab or self.tab().count()==0:
                self.newTab(view, str(view.title()))
            return view

    def loadViaHTTP(self, fileName=None, target=None):
        from core import CCP4Utils
        relPath = os.path.relpath(fileName, os.path.join(CCP4Utils.getCCP4I2Dir(), 'docs', 'sphinx', 'build', 'html'))
        #print 'loadViaHTTP',relPath
        url =  QtCore.QUrl('http://127.0.0.1:43434/' + relPath + '#' + target)
        #print 'loadViaHTTP fragment',url.fragment()
        self.loadPage(url, True)

    def handleLinkClick(self, url):
        #print 'CCP4BrowserWindow.handleLinkClick',url
        pass

    @QtCore.Slot(tuple)
    def IconReady(self, args=None):
        icon,view = args
        for i in range(self.tab().count()):
            widget = self.tab().widget(i)
            if widget is view:
                self.tab().setTabIcon(i,icon)

    @QtCore.Slot()
    def searchEdited(self):
#FIXME PYQT - potential hairiness
        text =  re.sub(r'\s+',' ',self.editSplitter.searchEdit.text())
        #print 'QWebBrowser.searchEdited', text
        self.loadPage(QtCore.QUrl(self.searchEngineString + str(text)))

    @QtCore.Slot()
    def addressEdited(self):
        text = str(self.editSplitter.addressEdit.text())
        if os.path.exists(text.split('#')[0]):
            self.openFile(fileName=text, internal=True)
        else:
            self.loadPage(QtCore.QUrl(self.editSplitter.addressEdit.text()),newTab=False)

    def setHistoryActionAvailability(self):
        widget = self.tab().currentWidget()
        forwardAction = self.findChild(QtWidgets.QAction, 'forward')
        backAction = self.findChild(QtWidgets.QAction, 'back')
        if hasattr(widget,"history") and isinstance(widget.history, Callable):
            history = widget.history()
            #print 'setHistoryActionAvailability', backAction, history.canGoBack()
            if forwardAction:
                if history.canGoForward():
                    forwardAction.setEnabled(True)
                else:
                    forwardAction.setEnabled(False)
            if backAction:
                if history.canGoBack():
                    backAction.setEnabled(True)
                else:
                    backAction.setEnabled(False)
        else:
            # No notion of history so set actions disabled
            if forwardAction:
                forwardAction.setEnabled(False)
            if backAction:
                backAction.setEnabled(False)

    @QtCore.Slot()
    def historyForward(self):
        widget = self.tab().currentWidget()
        if widget and hasattr(widget, "history") and isinstance(widget.history, Callable):
            history = widget.history()
            if history.canGoForward():
                history.goToItem(history.forwardItem())
        self.setHistoryActionAvailability()

    @QtCore.Slot()
    def historyBack(self):
        widget = self.tab().currentWidget()
        if widget and hasattr(widget, "history") and isinstance(widget.history, Callable):
            history = widget.history()
            if history.canGoBack():
                history.goToItem(history.backItem())
        self.setHistoryActionAvailability()

    @QtCore.Slot()
    def reloadPage(self):
        widget = self.tab().currentWidget()
        if widget is None:
            return
        widget.page().triggerAction(QtWebEngineWidgets.QWebEnginePage.ReloadAndBypassCache)

#-------------------------------------------------------------------
    def currentWidget(self):
#-------------------------------------------------------------------
        if self.tab() is None:
            return None
        idx = self.tab().currentIndex()
        if idx < 0:
            return None
        return self.tab().widget(idx)

#-------------------------------------------------------------------
    def widgetIsRunable(self):
#-------------------------------------------------------------------
        w = self.currentWidget()
        if w is not None:
            return w.isRunable()
        else:
            return False

#-------------------------------------------------------------------
    def widgetIsSaveable(self):
#-------------------------------------------------------------------
        w = self.currentWidget()
        if w is not None:
            return w.isSaveable()
        else:
            return False
#-------------------------------------------------------------------
    def widgetIsPrintable(self):
#-------------------------------------------------------------------
        w = self.currentWidget()
        if w is not None:
            return w.isPrintable()
        else:
            return False
#-------------------------------------------------------------------
    def widgetIsSearchable(self):
#-------------------------------------------------------------------
        w = self.currentWidget()
        if w is not None:
            return w.isSearchable()
        else:
            return False
#-------------------------------------------------------------------
    def widgetIsScaleable(self):
#-------------------------------------------------------------------
        w = self.currentWidget()
        if w is not None:
            return w.isScaleable()
        else:
            return None

#-------------------------------------------------------------------
    def isCurrentWidget(self, mode='text'):
#-------------------------------------------------------------------
        widget = self.currentWidget()
        if widget is None:
            return 0
        if mode=='text':
            from qtgui import CCP4TextViewer
            if isinstance(widget, CCP4TextViewer.CTextViewer):
                return 1
            else:
                return 0
        if mode=='image':
            from qtgui import CCP4ImageViewer
            if isinstance(widget, CCP4ImageViewer.CImageViewer):
                return 1
            else:
                return 0
        else:
            return 0

#-------------------------------------------------------------------
    def handlePrint(self):
#-------------------------------------------------------------------
        idx = self.tab().currentIndex()
        if idx < 0:
            return
        printer = QtPrintSupport.QPrinter()
        pd = QtPrintSupport.QPrintDialog(printer)
        pdret = pd.exec_()
        if pdret == QtWidgets.QDialog.Accepted:
            painter = QtGui.QPainter()
            painter.begin(printer)
            self.tab().widget(idx).Print(painter)
            painter.end()

#-------------------------------------------------------------------
    def handleSave(self):
#-------------------------------------------------------------------
        from qtgui import CCP4FileBrowser
        widget = self.currentWidget()
        if widget is None:
            return
        defaultSuffix = ''
        filter_list = []
        ext_list =  widget.getFileExt()
        if ext_list:
            defaultSuffix = ext_list[0]
        ext_text = ''
        for ext in widget.getFileExt():
            ext_text = ext_text + '*.' + ext + ' '
        if ext_text:
            filter_list.append(widget.getLabel() + ' (' + ext_text + ')')
        self.saveWidget = CCP4FileBrowser.CFileDialog(self, title='Save ' + widget.objectName(), filters=filter_list, defaultSuffix=defaultSuffix, fileMode=QtWidgets.QFileDialog.AnyFile)
        self.saveWidget.selectFile.connect(self.currentWidget().Save)
        self.saveWidget.show()

    def handleRun(self):
        widget = self.currentWidget()
        if widget is None:
            return
        widget.run()

#-------------------------------------------------------------------
    def openFind(self):
#-------------------------------------------------------------------
        if self.findDock.isVisible():
            self.findDock.setVisible(0)
        else:
            self.findDock.setVisible(1)
            self.findDock.widget().findTextWidget.setFocus(QtCore.Qt.OtherFocusReason)

#-------------------------------------------------------------------
    def isFindFrameOpen(self):
#-------------------------------------------------------------------
        if self.findFrame is not None and self.findFrame.isVisible():
            return 1
        return 0

#-------------------------------------------------------------------
    def handleFind(self, direction=1):
#-------------------------------------------------------------------
        widget=self.currentWidget()
        if widget is None or not widget.isSearchable():
            return
        text = self.findFrame.text()
        #print 'CCP4WebBrowser.handleFind',text,direction
        widget.findText(subString=text, direction=direction, caseSensitive=self.findFrame.caseSensitive())

#-------------------------------------------------------------------
    def handleFindNext(self):
#-------------------------------------------------------------------
        self.handleFind(1)

#-------------------------------------------------------------------
    def handleFindPrevious(self):
#-------------------------------------------------------------------
        self.handleFind(-1)

#-------------------------------------------------------------------
    def handleFindHighlightAll(self):
#-------------------------------------------------------------------
        widget=self.currentWidget()
        if widget is None or not widget.isSearchable():
            return
        text = self.findFrame.text()
        rv = widget.findText(subString=text, highlightAll=1)

    def setFixedWidth(self, width=0):
        if width <= 0:
            self.setMaximumWidth(16777215)
            self.setMinimumWidth(0)
        else:
            self.setMaximumWidth(width + 10)
            self.setMinimumWidth(width + 10)
      
    def moveEvent(self, event):
        delta = event.pos()-event.oldPos()
        #print 'CWebBrowser.moveEvent', event.pos().x(), event.pos().y(), event.oldPos().x(), event.oldPos().y(), '*', delta.x(), delta.y()
        #self.windowMoved.emit(delta)
        event.ignore()

    def openSendReport(self):
        from qtgui import CCP4ErrorReportViewer
        widget = CCP4ErrorReportViewer.CSendJobError(self, projectId=self.getProject())
        widget.show()

    def handleProjectMenuExport(self):
        pass

    def openManageImportFiles(self):
        pass


class CFindFrame(QtWidgets.QWidget):

    findNext = QtCore.Signal()
    findPrevious = QtCore.Signal()

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        find_layout = QtWidgets.QHBoxLayout()
        margin = 2
        find_layout.setContentsMargins(margin, margin, margin, margin)
        find_layout.setSpacing(margin)
        #widget = QtWidgets.QPushButton('Close', self)
        #find_layout.addWidget(widget)
        find_layout.addWidget(QtWidgets.QLabel('Find:'))
        self.findTextWidget = QtWidgets.QLineEdit(self)
        self.findTextWidget.setClearButtonEnabled(True)
        find_layout.addWidget(self.findTextWidget)
        self.findTextWidget.returnPressed.connect(self.findNext.emit)
        @QtCore.Slot(str)
        def tChanged(s):
            self.findNext.emit()
        self.findTextWidget.textChanged.connect(tChanged)
        self.next = QtWidgets.QPushButton('Next', self)
        find_layout.addWidget(self.next)
        self.next.clicked.connect(self.findNext.emit)
        self.previous= QtWidgets.QPushButton('Previous', self)
        find_layout.addWidget(self.previous)
        self.previous.clicked.connect(self.findPrevious.emit)
        self.findCase = QtWidgets.QCheckBox('Match case', self)
        find_layout.addWidget(self.findCase)
        doneButton = QtWidgets.QPushButton("Done")
        find_layout.addWidget(doneButton)
        self.searchFailMessageBox = QtWidgets.QLabel('<span style="color:red;">Search string not found</span>')
        self.searchFailMessageBox.hide()
        find_layout.addWidget(doneButton)
        @QtCore.Slot()
        def hideFindWidget():
           self.findTextWidget.clear()
           self.hide()
        doneButton.clicked.connect(hideFindWidget)
        find_layout.addWidget(self.searchFailMessageBox)
        find_layout.addStretch(5)
        find_layout.addStretch(5)
        self.setLayout(find_layout)

    def text(self):
        return str(self.findTextWidget.text())

    def caseSensitive(self):
        return self.findCase.isChecked()
 
    @QtCore.Slot(bool)
    def searchCallback(self,t):
        if t or self.text()=="":
            self.searchFailMessageBox.hide()
        else:
            self.searchFailMessageBox.show()
