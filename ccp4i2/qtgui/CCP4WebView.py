"""
Copyright (C) 2009-2010 University of York
Liz Potterton Jan 2010 - Copied from CCP4mg
"""

##@package CCP4WebView (QtWebKit) Web browser 'plugin' to view web pages

import os
import re
import sys
import traceback
import xml.etree.ElementTree as ET

from lxml import etree
from PySide2 import QtCore, QtGui, QtWebEngineWidgets, QtWidgets

from ..core import CCP4Utils
from ..core.CCP4ErrorHandling import CException, Severity
from ..core.CCP4Modules import HTTPSERVER
from ..core.CCP4Modules import PREFERENCES
from ..core.CCP4WarningMessage import warningMessage
from ..utils.QApp import QTAPPLICATION


def setGlobalSettings():
    for window in CWebView.Instances:
        if CCP4Utils.isAlive(window):
            window.setLoggraphFont()


class CWebPage(QtWebEngineWidgets.QWebEnginePage):

    NavigationRequest = QtCore.Signal('QUrl')
    CustomMimeTypeRequested = QtCore.Signal('QUrl')

    def __init__(self, parent=None):
        QtWebEngineWidgets.QWebEnginePage.__init__(self, parent)

    def shouldInterruptJavaScript(self):
        # Reimplemented to prevent accasianally seen message box (by Phil) when creating DR reports
        # Implied from http://forum.qt.io/topic/21387/maximum-duration-for-a-script/3
        #http://stackoverflow.com/questions/9284511/reimplement-the-shouldinterruptjavascript-in-qt-c
        print('CWebPage.shouldInterruptJavaScript called')
        traceback.print_stack()
        QTAPPLICATION().processEvents(QtCore.QEventLoop.AllEvents, 42)
        return False

    def acceptNavigationRequest(self, url, navType, isMainFrame):
        print ('CWebPage.acceptNavigationRequest',url)
        # if PREFERENCES().EXTERNAL_FILES_IN_EXTERNAL_BROWSER and not request.url().isLocalFile():
        isXiaReport = (re.search('xia2.html', url.path()) is not None) or (re.search('xia2.html', url.query()) is not None)
        isPairefRep = (re.search('PAIREF_project.html', url.path()) is not None) or (re.search('PAIREF_project.html', url.query()) is not None)
        isArcimboldoReport = (re.search('arcimboldo.html', url.path()) is not None) or (re.search('arcimboldo.html', url.query()) is not None)
        isDuiReport = (re.search(r'dui_files.*?\d_report\.html', url.path()) is not None) or (re.search(r'dui_files.*?\d_report\.html', url.query()) is not None)
        isMRParseRe = (re.search('mrparse_i2.html', url.path()) is not None) or (re.search('mrparse_i2.html', url.query()) is not None)
        isDocx = (re.search('tabtest.docx', url.path()) is not None) or (re.search('tabtest.docx', url.query()) is not None)
        isReportFile = False
        if url.path() == "/database/" and url.host() == "127.0.0.1" and url.query().find("fileName=report.html") > -1:
            isReportFile = True
        isIntfReport = isXiaReport or isDuiReport
        if (PREFERENCES().EXTERNAL_FILES_IN_EXTERNAL_BROWSER and not url.isLocalFile() and not isReportFile and not isPairefRep and not isMRParseRe) or  isDocx or isIntfReport:
            # if not (url.path() == "/database/" and url.host() == "127.0.0.1"):
            #     url = QtCore.QUrl(url.path()) # Removing this seems to sort out the link problem.
            #print("OPEN ", url, url.path(),isMRParseRe)
            rv = QtGui.QDesktopServices.openUrl(url)
            if not rv:
                QtWidgets.QMessageBox.warning(None, 'Display url', 'Attempting to display external URL using non-CCP4i2 web browser failed')
            return False

        # This is broken if setHtml() has been used to set the page text and give
        # a local base href - the base href is returned by request.url()
        #self.view().extractCCP4Data(request.url())
        return True

    def clickOn(self,pos):
        # send a mouse click event to the web page # print 'CWebpage.clickOn', pos.x(),pos.y()
        event = QtGui.QMouseEvent(QtCore.QEvent.MouseButtonPress, pos, QtCore.Qt.LeftButton, QtCore.Qt.LeftButton, QtCore.Qt.NoModifier)
        QTAPPLICATION().sendEvent(self, event)
        event = QtGui.QMouseEvent(QtCore.QEvent.MouseButtonRelease, pos, QtCore.Qt.LeftButton, QtCore.Qt.LeftButton, QtCore.Qt.NoModifier)
        QTAPPLICATION().sendEvent(self, event)


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
        shortcut = QtWidgets.QShortcut(QtGui.QKeySequence(self.tr("Ctrl+C", "Copy")), self)
        shortcut.activated.connect(self.copyHighlighted)
        self.setTarget('')
        #self.mgpage.unsupportedContent.connect(self.handleUnsupportedContent)
        #Create Qt web plugin factory
        fontDataBase = QtGui.QFontDatabase()
        defaultSettings = QtWebEngineWidgets.QWebEngineSettings.globalSettings()
        standardFont = fontDataBase.font("Courier","",14)
        defaultSettings.setFontFamily(QtWebEngineWidgets.QWebEngineSettings.StandardFont, standardFont.family())
        if PREFERENCES().DISABLE_WEBGL:
            defaultSettings.setAttribute(QtWebEngineWidgets.QWebEngineSettings.WebGLEnabled,False)
        defaultSettings.setAttribute(QtWebEngineWidgets.QWebEngineSettings.JavascriptCanAccessClipboard, True)

    @staticmethod
    def updateInstances(qobj):
        CWebView.Instances = set([window for window in CWebView.Instances if CCP4Utils.isAlive(window)])

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

        baseName = os.path.basename(self.url().path())
        if self.url().path().startswith("/database/projectId"):
            nonReport = True
        if (
            self.url().scheme() == "file"
            and baseName not in {"report.html", "report_tmp.html"}
        ):
            nonReport = True
        if (
            self.url().scheme() == "http"
            and "/database/projectid/" in self.url().path()
            and "/jobnumber/" in self.url().path()
            and "/file/" in self.url().path()
            and baseName not in {"blank.html", "report.html", "report_tmp.html"}
        ):
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
        QTAPPLICATION().clipboard().setText(self.selectedText())

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

        flags = QtWebEngineWidgets.QWebEnginePage.FindFlags(0)
        if direction < 0:
            flags = flags | QtWebEngineWidgets.QWebEnginePage.FindBackward 
        if caseSensitive:
            flags = flags | QtWebEngineWidgets.QWebEnginePage.FindCaseSensitively
        QtWebEngineWidgets.QWebEngineView.findText(self, subString, flags, searchCallback)

    def getFileExt(self):
        return None

    def getLabel(self):
        return None

    def Save(self,fileName):
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

    def editHTTPPort(self, fileName, port=None):
        if len(fileName) < 0:
            return
        try:
            text = CCP4Utils.readFile(fileName)
        except:
            return False
        if port is None: 
            port = HTTPSERVER().port
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
        #html_exts = MIMETYPESHANDLER().getMimeTypeInfo('text/html','fileExtensions')
        html_exts = ['html','htm']
        self.report = None
        self.fileName = os.path.normpath(str(url.toLocalFile()))
        ext = os.path.splitext(self.fileName)[1]
        if len(ext)>1:
            ext=ext[1:]
        if html_exts.count(ext) and os.path.exists(self.fileName):
            self.lastModTime = os.path.getmtime(self.fileName)
            from ..qtcore import CCP4Report
            self.report = CCP4Report.CReport()
            try:
                self.report.loadFromXmlFile(self.fileName)
                #print 'CWebView.extractCCP4Data: Loaded CCP4 data:',self.report.containsCCP4Data()
            except CException as e:
                if reportErrors: warningMessage(e, windowTitle=self.windowTitle())
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
        #print 'setLoggraphFont',PREFERENCES().REPORT_FONT_SIZE-2, font
        fontSel = {'family' : 'Lucida Sans Unicode', 'size' : int(PREFERENCES().REPORT_FONT_SIZE - 2)}
        fontSel.update(font)
        from ..pimple import MGQTmatplotlib
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
        from ..dbapi import CCP4DbUtils
        from ..report import CCP4ReportGenerator
        openJob = CCP4DbUtils.COpenJob(jobId=jobId)
        generator = CCP4ReportGenerator.CReportGenerator(jobId=openJob.jobId,jobStatus=openJob.status)
        try:
            reportFile, newPageOrNewData = generator.makeReportFile()
        except CException as e:
            if reportErr and e.maxSeverity()>Severity.WARNING:
                warningMessage(e)
        except Exception as e:
            if reportErr:
                QtWidgets.QMessageBox.warning(self,self.parent().windowTitle(),'Unknown error creating report file for job number '+str(openJob.jobNumber))
        if os.path.exists(reportFile):
            err = self.merge(reportFile,idText)

    def merge(self,reportFile,idText):
        #myTree =  etree.parse(reportFile)
        myTree = ET.parse(reportFile).getroot()
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
