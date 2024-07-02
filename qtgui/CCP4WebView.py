from __future__ import print_function

"""
     CCP4WebView.py: CCP4 GUI Project
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
     Liz Potterton Jan 2010 - Copied from CCP4mg
"""
##@package CCP4WebView (QtWebKit) Web browser 'plugin' to view web pages
import os
import re
import traceback

from PySide6 import QtGui, QtWidgets, QtCore, QtWebEngineCore, QtWebEngineWidgets
from core.CCP4ErrorHandling import *
from core import CCP4Modules, CCP4Config

def isAlive(qobj):
    import shiboken6
    return shiboken6.isValid(qobj)

def setGlobalSettings():
    """
    webSettings = QtWebKit.QWebSettings.globalSettings()
    #print 'setGlobalSetting webSettings',webSettings,key,value
    preferences = CCP4Modules.PREFERENCES()
    webSettings.setFontSize(QtWebKit.QWebSettings.DefaultFontSize, preferences.REPORT_FONT_SIZE)
    webSettings.setFontSize(QtWebKit.QWebSettings.DefaultFixedFontSize, preferences.REPORT_FONT_SIZE)
    webSettings.setAttribute(QtWebKit.QWebSettings.DeveloperExtrasEnabled, True)
    """
    for window in CWebView.Instances:
        if isAlive(window):
            window.setLoggraphFont()

class CWebPage(QtWebEngineCore.QWebEnginePage):

    NavigationRequest = QtCore.Signal('QUrl')
    CustomMimeTypeRequested = QtCore.Signal('QUrl')

    def __init__(self, parent=None):
        QtWebEngineCore.QWebEnginePage.__init__(self, parent)
        #self.setLinkDelegationPolicy(QtWebEngineWidgets.QWebEnginePage.DelegateAllLinks)
        #self.setForwardUnsupportedContent(True)
        #self.settings().setAttribute(QtWebKit.QWebSettings.JavascriptEnabled ,1)
        #self.settings().setAttribute(QtWebKit.QWebSettings.PluginsEnabled ,1)

    def shouldInterruptJavaScript(self):
        # Reimplemented to prevent accasianally seen message box (by Phil) when creating DR reports
        # Implied from http://forum.qt.io/topic/21387/maximum-duration-for-a-script/3
        #http://stackoverflow.com/questions/9284511/reimplement-the-shouldinterruptjavascript-in-qt-c
        print('CWebPage.shouldInterruptJavaScript called')
        traceback.print_stack()
        #QApplication::processEvents(QEventLoop::AllEvents, 42)
        CCP4Modules.QTAPPLICATION().processEvents(QtCore.QEventLoop.AllEvents, 42)
        return False

    def javaScriptConsoleMessage(self, level, message, lineNumber, sourceID):
        return
        print(" ############ JAVASCRIPT MESSAGE ########### ")
        print("lineNumber: ", lineNumber, " sourceID: ", sourceID, " message: ", message, " level",level)
        print(" ########################################### ")

    def acceptNavigationRequest(self, url, navType, isMainFrame):
        print ('CWebPage.acceptNavigationRequest',url)
        # if CCP4Modules.PREFERENCES().EXTERNAL_FILES_IN_EXTERNAL_BROWSER and not request.url().isLocalFile():
        isXiaReport = (re.search('xia2.html', url.path()) is not None) or (re.search('xia2.html', url.query()) is not None)
        isPairefRep = (re.search('PAIREF_project.html', url.path()) is not None) or (re.search('PAIREF_project.html', url.query()) is not None)
        isArcimboldoReport = (re.search('arcimboldo.html', url.path()) is not None) or (re.search('arcimboldo.html', url.query()) is not None)
        isDuiReport = (re.search('dui_files.*?\d_report\.html', url.path()) is not None) or (re.search('dui_files.*?\d_report\.html', url.query()) is not None)
        isMRParseRe = (re.search('mrparse_i2.html', url.path()) is not None) or (re.search('mrparse_i2.html', url.query()) is not None)
        isDocx = (re.search('tabtest.docx', url.path()) is not None) or (re.search('tabtest.docx', url.query()) is not None)
        isReportFile = False
        if url.path() == "/database/" and url.host() == "127.0.0.1" and url.query().find("fileName=report.html") > -1:
            isReportFile = True
        isIntfReport = isXiaReport or isDuiReport
        if (CCP4Modules.PREFERENCES().EXTERNAL_FILES_IN_EXTERNAL_BROWSER and not url.isLocalFile() and not isReportFile and not isPairefRep and not isMRParseRe) or  isDocx or isIntfReport:
            # if not (url.path() == "/database/" and url.host() == "127.0.0.1"):
            #     url = QtCore.QUrl(url.path()) # Removing this seems to sort out the link problem.
            #print("OPEN ", url, url.path(),isMRParseRe)
            rv = QtGui.QDesktopServices.openUrl(url)
            if not rv:
                QtWidgets.QMessageBox.warning(None, 'Display url', 'Attempting to display external URL using non-CCP4i2 web browser failed')
            return False
        # path = url.path() # Redundant line of code.
#FIXME
        """
        format = CCP4Modules.MIMETYPESHANDLER().formatFromFileExt(path)
        if format and not CCP4Modules.MIMETYPESHANDLER().useWebBrowser(format):
            self.CustomMimeTypeRequested.emit(url)
            return False    # CCP4WebBrowser is not taking responsibility.
        """
#FIXME _ Bummer, this is called before entering load function in WebEngine. Thus, _blockLoading is always True
        """
        if self.parent().blockLoading():
            self.NavigationRequest.emit(url)
            return False
        """
        # This is broken if setHtml() has been used to set the page text and give
        # a local base href - the base href is returned by request.url()
        #self.view().extractCCP4Data(request.url())
        return True

    def clickOn(self,pos):
        # send a mouse click event to the web page # print 'CWebpage.clickOn', pos.x(),pos.y()
        event = QtGui.QMouseEvent(QtCore.QEvent.MouseButtonPress, pos, QtCore.Qt.LeftButton, QtCore.Qt.LeftButton, QtCore.Qt.NoModifier)
        CCP4Modules.QTAPPLICATION().sendEvent(self, event)
        event = QtGui.QMouseEvent(QtCore.QEvent.MouseButtonRelease, pos, QtCore.Qt.LeftButton, QtCore.Qt.LeftButton, QtCore.Qt.NoModifier)
        CCP4Modules.QTAPPLICATION().sendEvent(self, event)


class CWebView(QtWebEngineWidgets.QWebEngineView):

    IconReady = QtCore.Signal(tuple)
    NavigationRequest = QtCore.Signal('QUrl')
    CustomMimeTypeRequested = QtCore.Signal('QUrl')
    StatusBarMessage = QtCore.Signal(str)
    csvDownloadRequest = QtCore.Signal(str,str)
    nonReportNavigated = QtCore.Signal()
    reportNavigated = QtCore.Signal()
    searchFound = QtCore.Signal(bool)

    ERROR_CODES= {101 : {'description' : 'Error creating CCP4i2 Web Plugins'}}

    Instances = set()

    zoomFactorChanged = QtCore.Signal(float)

    def wheelEvent(self,e):
        if (e.modifiers() & QtCore.Qt.ControlModifier) or (e.modifiers() & QtCore.Qt.MetaModifier):
            if e.delta() > 0:
                self.setZoomFactor(self.zoomFactor()*1.2)
                self.zoomFactorChanged.emit(self.zoomFactor())
            else:
                self.setZoomFactor(self.zoomFactor()/1.2)
                self.zoomFactorChanged.emit(self.zoomFactor())
        else:
            QtWebEngineWidgets.QWebEngineView.wheelEvent(self,e)

    def __init__(self, parent=None, blockLoading=False):
        QtWebEngineWidgets.QWebEngineView.__init__(self,parent)
        CWebView.Instances.add(self)
        self.destroyed.connect(CWebView.updateInstances)
        self._blockLoading = blockLoading
        print("blockLoading",blockLoading)
        self.fileName = None
        self.report = None
        self.resetScroll = None
        self.mgpage = CWebPage(self)
        self.customContextMenu = None
        self.setPage(self.mgpage)
        self.mgpage.CustomMimeTypeRequested.connect(self.CustomMimeTypeRequested.emit)
        #self.mgpage.statusBarMessage.connect(self.StatusBarMessage.emit)
        self.mgpage.linkHovered.connect(self.StatusBarMessage.emit)
        #self.mgpage.linkClicked.connect(self.load)
        self.loadFinished.connect(self.LoadFinished)
        self.loadStarted.connect(self.LoadStarted)
        shortcut = QtGui.QShortcut(QtGui.QKeySequence(self.tr("Ctrl+C", "Copy")), self)
        shortcut.activated.connect(self.copyHighlighted)
        self.setTarget('')
        #self.mgpage.unsupportedContent.connect(self.handleUnsupportedContent)
        #Create Qt web plugin factory
        """
        try:
            from qtgui import CCP4WebPluginFactory
            factory = CCP4WebPluginFactory.CWebPluginFactory(self)
            self.page().setPluginFactory(factory)
        except CException as e:
            print('CCP4WebView Error creating web plugin factory')
        except Exception as e:
            print('CCP4WebView Error creating web plugin factory')
        """
        fontDataBase = QtGui.QFontDatabase()

        profile = QtWebEngineCore.QWebEngineProfile.defaultProfile()
        defaultSettings = profile.settings()
        standardFont = fontDataBase.font("Courier","",14)
        defaultSettings.setFontFamily(QtWebEngineCore.QWebEngineSettings.StandardFont, standardFont.family())
        if CCP4Modules.PREFERENCES().DISABLE_WEBGL:
            defaultSettings.setAttribute(QtWebEngineCore.QWebEngineSettings.WebGLEnabled,False)

    @staticmethod
    def updateInstances(qobj):
        CWebView.Instances = set([window for window in CWebView.Instances if isAlive(window)])

    def objectPath(self):
        return '' if self.mgpage is None else str(self.mgpage.objectName())

    def blockLoading(self):
        print("blockLoading",self._blockLoading); sys.stdout.flush()
        return self._blockLoading

    @QtCore.Slot()
    def LoadStarted(self):
        icon = QtGui.QIcon()
        self.IconReady.emit((icon, self))

    def setTarget(self, target=None):
        self.target=target

    @QtCore.Slot(bool)
    def LoadFinished(self, ok=None):
        nonReport = False
        for q in self.url().query().split("?"):
            if q.startswith("fileName="):
                if q.split("=")[1] != "report.html" and q.split("=")[1] !=  "report_tmp.html":
                    nonReport = True

        if self.url().path().startswith("/database/projectId"):
                nonReport = True
        if self.url().scheme() == "file" and os.path.basename(self.url().path()) != "report.html" and os.path.basename(self.url().path()) != "report_tmp.html":
            nonReport = True
        if nonReport:
            if self.history().canGoBack():
                print("########################################")
                print("Non-report file:")
                print(self.url())
                print("########################################")
                self.nonReportNavigated.emit()
        else:
            self.reportNavigated.emit()
        icon = self.icon()
        #print "CWebView.LoadFinished", self.target, ok,self.url().toLocalFile()
        if ok:
            return
        if self.target:
            #print 'currentFrame', self.page().mainFrame(), self.target.strip('#')
            self.page().mainFrame().scrollToAnchor(self.target.strip('#'))
            self.setTarget('')
            #self.load(myurl)
        self.IconReady.emit((icon, self))
        self.subJobReport = CSubJobReport(self)
#FIXME - More QWebChannel stuff?
        #self.page().currentFrame().addToJavaScriptWindowObject('SubJobReport', self.subJobReport)
        self.setLoggraphFont()
        if self.resetScroll is not None:
            self.page().mainFrame().setScrollBarValue(QtCore.Qt.Vertical, self.resetScroll)
            self.resetScroll = None

    def contextMenuEvent(self, event):
        """
        #print 'CWebView.contextMenuEvent', event.pos().x(), event.pos().y()
        hitTestResult = None
        alt_words = []
        frame = self.page().frameAt(event.pos())
        if frame:
            hitTestResult = frame.hitTestContent(event.pos())
        #print 'CWebView.contextMenuEvent',frame,hitTestResult
        #print 'CWebView.contextMenuEvent hitTestResult alt text',str(hitTestResult.alternateText())
        alt_text = hitTestResult.alternateText()
        if len(alt_text)>0:
            alt_words = (str(alt_text)).split()
        #print 'CWebView.contextMenuEvent alt text',alt_words
        if hitTestResult:
            webElement = hitTestResult.element()
        if len(alt_words) >= 2 and alt_words[0] == 'x-ccp4-widget/CDataFileView':
            if not self.customContextMenu:
                from qtgui import CCP4GuiUtils
                self.customContextMenu = QtWidgets.QMenu(self)
                CCP4GuiUtils.populateMenu(self, self.customContextMenu,['copy', 'view', 'help'], default_icon='')
            self.customContextMenu.popup(event.globalPos())
            event.accept()
        else:
            QtWebEngineWidgets.QWebEngineView.contextMenuEvent(self, event)
        """
        menu = QtWidgets.QMenu(self)
        reloadAct = menu.addAction("Reload")
        reloadAct.triggered.connect(self.reload)
        print(self.url())
        if self.url().scheme() == "file" and os.path.basename(self.url().path()) != "report.html":
            if self.history().canGoBack():
                backAct = menu.addAction("Back")
                backAct.triggered.connect(self.history().back)
        menu.popup(event.globalPos());

    def getActionDef(self, name=''):
        if name == 'copy':
            return dict(text=self.tr("Copy"), tip=self.tr('Copy data'), slot=self.dummyHandler)
        if name == 'help':
            return dict(text=self.tr("Help"), tip = self.tr('Help'), slot = self.dummyHandler)
        if name == 'view':
            return dict(text=self.tr("View"), tip=self.tr('View data'), slot=self.dummyHandler)

    def dummyHandler(self):
            pass

    def copyHighlighted(self):
        CCP4Modules.QTAPPLICATION().clipboard().setText(self.selectedText())

    def handleUnsupportedContent(self, qnr):
        #print 'handleUnsupportedContent',qnr,qnr.url(),'scheme',str(qnr.url().scheme()),'path',str(qnr.url().path())
        pass

    def isPrintable(self):
        return 1
  
    def isSaveable(self):
        return 1
  
    def isSearchable(self):
        return 1

    def isRunable(self):
        return 0

    def findText(self, subString='',direction=1, caseSensitive=0, wrapAroundDocument=1, highlightAll=0):
        def searchCallback(t):
            self.searchFound.emit(t)

        flags = QtWebEngineCore.QWebEnginePage.FindFlags(0)
        if direction < 0:
            flags = flags | QtWebEngineCore.QWebEnginePage.FindBackward 
        if caseSensitive:
            flags = flags | QtWebEngineCore.QWebEnginePage.FindCaseSensitively
        QtWebEngineWidgets.QWebEngineView.findText(self, subString, flags, searchCallback)

    def getFileExt(self):
        return None
        from qtcore import CCP4CustomMimeTypes
        return CCP4CustomMimeTypes.MimeTypesHandler().getMimeTypeInfo(name='text/html', info='fileExtensions')

    def getLabel(self):
        return None
        from qtcore import CCP4CustomMimeTypes
        return CCP4CustomMimeTypes.MimeTypesHandler().getMimeTypeInfo(name='text/html', info='description')

    def Save(self,fileName):
        from core import CCP4Utils
        CCP4Utils.saveFile(fileName, str(self.mgpage.currentFrame().toHtml()))

    def title(self):
        t = QtWebEngineWidgets.QWebEngineView.title(self)
        return str(t)

    def load(self, url, reportErrors=False):
        self.extractCCP4Data(url, reportErrors=reportErrors)
        #print 'WebView.load', str(url.toLocalFile()),self.report,self.report.resetBaseHref
        self.editHTTPPort(self.fileName)
        block = self._blockLoading
        self._blockLoading = False
        if self.report is None or self.report.resetBaseHref is None:
            QtWebEngineWidgets.QWebEngineView.load(self, url)
        else:
            url = QtCore.QUrl.fromLocalFile(self.report.resetBaseHref)
            QtWebEngineWidgets.QWebEngineView.load(self, url)
            #print 'CWebView.load fragment',url.fragment()
        if url.isLocalFile():
            self.fileName = str(url.toLocalFile())
        else:
            self.fileName = str(url)
        self._blockLoading = block

    def attemptImageLoad(self):
        self.page().runJavaScript('''
var imgElements = document.getElementsByTagName('img');
for (var i=0; i<imgElements.length;i++) {
   var imgElement = imgElements[i];
   var imgPath = imgElement.getAttribute('src');
   var d=new Date();
   $(imgElement).attr('src',imgPath+'?'+d.getTime());}''')

    """
    def reload(self):
        self.resetScroll = self.page().mainFrame().scrollBarValue(QtCore.Qt.Vertical)
        block = self._blockLoading
        self._blockLoading = False
        #if self.report is None or self.report.resetBaseHref is None:
        #  QtWebEngineWidgets.QWebEngineView.reload(self)
        #else:
        #  QtWebEngineWidgets.QWebEngineView.load(self,QtCore.QUrl.fromLocalFile(+self.report.resetBaseHref))
        self.load(self.url())
        self._blockLoading = block
    """

    def loadHtmlText(self):
        # http://stackoverflow.com/questions/2727080/how-to-get-qwebkit-to-display-image
        # answer from badcat that security feature does not allow local disk baseRef in setHtml
        print('CWebView loading text from', self.fileName)

    def editHTTPPort(self, fileName, port=None):
        from core import CCP4Utils
        if len(fileName) < 0:
            return
        try:
            text = CCP4Utils.readFile(fileName)
        except:
            return False
        if port is None: 
            port = CCP4Modules.HTTPSERVER().port
        changed = False
        #print 'CWebView.editHTTPPort',re.search('../../report_files',text)
        if re.search('../../report_files', text) is not None:
            changed = True
            text = re.sub('../../report_files', 'http://127.0.0.1:' + str(port) + '/report_files', text)
            #print 'CWebView.editHTTPPort after first sub',text[0:500]
        else:
            matchObj = re.search('(.*)127.0.0.1:(.*)/report_files(.*)', text)
            #print 'CWebView.editHTTPPort matchObj',matchObj
            if matchObj is not None:
                oldPort = matchObj.groups()[1]
                if oldPort != str(port):
                    changed = True
                    text = re.sub('127.0.0.1:'+oldPort+'/report_files','127.0.0.1:'+str(port)+'/report_files',text)
        if changed:
            CCP4Utils.backupFile(fileName,delete=True)
            CCP4Utils.saveFile(fileName,text)
        return changed

    def extractCCP4Data(self, url, reportErrors=False):
        #html_exts = CCP4Modules.MIMETYPESHANDLER().getMimeTypeInfo('text/html','fileExtensions')
        html_exts = ['html','htm']
        self.report = None
        self.fileName = os.path.normpath(str(url.toLocalFile()))
        ext = os.path.splitext(self.fileName)[1]
        if len(ext)>1:
            ext=ext[1:]
        if html_exts.count(ext) and os.path.exists(self.fileName):
            self.lastModTime = os.path.getmtime(self.fileName)
            from qtcore import CCP4Report
            self.report = CCP4Report.CReport()
            try:
                self.report.loadFromXmlFile(self.fileName)
                #print 'CWebView.extractCCP4Data: Loaded CCP4 data:',self.report.containsCCP4Data()
            except CException as e:
                if reportErrors: e.warningMessage(windowTitle=self.windowTitle())
        else:
            self.fileName = ''
            self.lastModTime = None

    def isFileModified(self):
        if len(self.fileName) == 0: 
            return False
        modTime = os.path.getmtime(self.fileName)
        if modTime > self.lastModTime:
            return True
        else:
            return False

    def close(self):
        pass

    def clear(self):
        self.mgpage = CWebPage(self)
        self.setPage(self.mgpage)

    def handleTabbedOpen(self):
        pass

    def handleTabbedClosed(self):
        pass

    def setLoggraphFont(self, font={}):
        #print 'setLoggraphFont',CCP4Modules.PREFERENCES().REPORT_FONT_SIZE-2, font
        fontSel = {'family' : 'Lucida Sans Unicode', 'size' : int(CCP4Modules.PREFERENCES().REPORT_FONT_SIZE - 2)}
        fontSel.update(font)
        from pimple import MGQTmatplotlib
        loggraphList = self.findChildren(MGQTmatplotlib.LogGraph)
        for loggraph in loggraphList:
            loggraph. applyPreferences(titleFontSel=fontSel, axesFontSel=fontSel, axesLabelFontSel=fontSel, legendFontSel=fontSel)

    def emitDownloadRequest(self, jobId, dataName):
        print('CWebView.emitDownloadRequest', jobId, dataName)
        self.csvDownloadRequest.emit(jobId, dataName)

class CSubJobReport(QtCore.QObject):

    def __init__(self,parent):
        QtCore.QObject.__init__(self,parent)

    @QtCore.Slot(str)
    def load(self,idText):
        ''' Invoked from xreport.js to merge sub-job report into paret report '''
        reportErr = True
        if idText[0:5] != 'jobId':
            print('ERROR in load sub-job report - CSubJobReport.load did not receive jobId')
            return
        jobId = int(idText[5:])
        from dbapi import CCP4DbUtils
        from report import CCP4ReportGenerator
        openJob = CCP4DbUtils.COpenJob(jobId=jobId)
        generator = CCP4ReportGenerator.CReportGenerator(jobId=openJob.jobId,jobStatus=openJob.status)
        if CCP4Config.DEVELOPER():
            reportFile, newPageOrNewData = generator.makeReportFile()
        try:
            reportFile, newPageOrNewData = generator.makeReportFile()
        except CException as e:
            if reportErr and e.maxSeverity()>SEVERITY_WARNING:
                e.warningMessage()
        except Exception as e:
            if reportErr:
                QtWidgets.QMessageBox.warning(self,self.parent().windowTitle(),'Unknown error creating report file for job number '+str(openJob.jobNumber))
        if os.path.exists(reportFile):
            err = self.merge(reportFile,idText)

    def merge(self,reportFile,idText):
        from lxml import etree
        from core import CCP4Utils
        '''
        parentTree = etree.parse(parentFile)
        print 'mergeIntoParent jobId',self.jobId
        # Find the span-link created by CCP4ReportParser.foldLinkLine
        path = 'body/div[@class="sub-job-list"]/span[@class="folder_link"]/a[@id="jobId'+str(self.jobId)+'"]'
        aEle = parentTree.xpath(path)
        print 'mergeIntoParent aEle',aEle
        if len(aEle) == 0: return CErrorReport(self.__class__,5,parentFile)
        aEle = aEle[0]
        label = aEle.text
        print 'mergeIntoParent',label
        insertEle = aEle.xpath('../..')[0]
        insertIndex = insertEle.index(aEle.xpath('..')[0])
        '''
        #myTree =  etree.parse(reportFile)
        myTree = CCP4Utils.openFileToEtree(reportFile)
        # Load the ccp4_data while we have the file as xml
        self.parent().report.loadFromEtree(myTree)
        # Convert body of sub-job report html to a div.subjob
        body = myTree.xpath('./body')[0]
        body.tag ='div'
        body.set('id',str(idText))
        body.set('class','subjob')
        subJobText = etree.tostring(body)
        #print 'CSubJobReport.merge subJobText(cut)',subJobText[0:200]
        ele = self.parent().page().currentFrame().findFirstElement('body div#'+idText)
        if ele:
            ele.replace(subJobText)
            ele.setStyleProperty('display','block')
            text = self.parent().page().currentFrame().toHtml().__str__()
            CCP4Utils.saveFile(self.parent().fileName+'.tmp',text,overwrite=True)

