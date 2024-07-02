from __future__ import print_function

"""
     CCP4DemoData.py: CCP4 GUI Project
     Copyright (C) 2015STFC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.sstac
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

"""
   Liz Potterton Sept 2015 - Create, download demo data
"""

import sys
import os
import glob
import time
import shutil
from PySide6 import QtGui, QtWidgets,QtCore
from core.CCP4ErrorHandling import *
from core import CCP4Utils
from core import CCP4Modules

class CDemoData:

    insts = None

    def __init__(self):
        self.builtInfo = False
        self.testDatasets = None
        self.restoreTestDatasets()

    def getOverviewPage(self):
        if not self.builtInfo: self.makeDemoDataInfo()
        #return os.path.join(CCP4Utils.getDotDirectory(),'demo_data','README.html')
        return os.path.join(CCP4Utils.getCCP4I2Dir(),'docs','sphinx','build','html','tutorials','index.html')
    
    def getTestDatasets(self):
        return self.testDatasets
    
    def loadTestDatasets(self):
        if self.testDatasets is None:
            #print 'CDemoData.loadTestDatasets'
            dirList = glob.glob(os.path.join(CCP4Utils.getCCP4I2Dir(),'demo_data','*'))
            dirList.extend(glob.glob(os.path.join(CCP4Utils.getDotDirectory(),'demo_data','*')))
            self.testDatasets = []
            for dr in dirList:
                #print 'loadTestDatasets',dr,'*',os.path.split(dr)[1]
                if os.path.isdir(dr) and os.path.split(dr)[1] != 'resources':
                    try:
                        infoText = CCP4Utils.readFile(os.path.join(dr,'INFO.txt'))
                        label = infoText.split('\n')[0]
                    except:
                        label = os.path.split(dr)[1]
                    self.testDatasets.append([os.path.split(dr)[1],label])
        return self.testDatasets

    def saveTestDatasets(self):
        from lxml import etree
        if self.testDatasets is None: self.loadTestDatasets()
        fileName = os.path.join(CCP4Utils.getCCP4I2Dir(),'demo_data','datasets_list.xml')
        #print 'saveTestDatasets',len(self.testDatasets)
        root = etree.Element('datasetList')
        for name,label in self.testDatasets:
            ele = etree.SubElement(root,'dataset')
            e = etree.SubElement(ele,'name')
            e.text = name
            e = etree.SubElement(ele,'label')
            e.text = label
        CCP4Utils.saveEtreeToFile(root,fileName)

    def restoreTestDatasets(self):
        fileName = os.path.join(CCP4Utils.getCCP4I2Dir(),'demo_data','datasets_list.xml')
        #print 'restoreTestDatasets',os.path.exists(fileName)
        if not os.path.exists(fileName):
            self.saveTestDatasets()
            return
        root = CCP4Utils.openFileToEtree(fileName)
        self.testDatasets = []
        for ele in root:
            self.testDatasets.append([ele.find('name').text,ele.find('label').text])
        #print 'restoreTestDatasets',self.testDatasets

    def copyDemoDataToProject(self,parentWidget=None,projectId=None,dataset=None):
        source = os.path.join(CCP4Utils.getCCP4I2Dir(),'demo_data',dataset)
        if not os.path.exists(source): source = os.path.join(CCP4Utils.getDotDirectory(),'demo_data',dataset)
        dest = os.path.join(CCP4Modules.PROJECTSMANAGER().db().getProjectInfo(projectId=projectId,mode='projectdirectory'),dataset)
        if os.path.exists(dest):
            QtWidgets.QMessageBox.information(parentWidget,'Download demo data','Test data was already copied to the project directory '+dest)
            return
        try:
            shutil.copytree(source,dest)
        except:
            QtWidgets.QMessageBox.warning(parentWidget,'Download demo data','Failed to copy test data to the project directory '+dest)
        else:
            QtWidgets.QMessageBox.information(parentWidget,'Download demo data','Test data copied to the project directory '+dest)

    def makeDemoDataInfo(self):
        sectionList =  ['Introduction','Molecular replacement','Ligands']
        infoFileList = glob.glob(os.path.join(CCP4Utils.getCCP4I2Dir(),'demo_data','*','INFO.txt'))
        infoFileList.extend(glob.glob(os.path.join(CCP4Utils.getDotDirectory(),'demo_data','*','INFO.txt')))
        for ii in range(len(infoFileList)-1,-1,-1):
            if infoFileList[ii].count('resources'): del infoFileList[ii]        
        info = []
        for fl in infoFileList:
            try:
                lines = CCP4Utils.readFile(fl).split('\n')
                if len(lines) < 3:
                    lines.append('Others')
                relPath = os.path.relpath(os.path.split(fl)[0],os.path.join(CCP4Utils.getDotDirectory(),'demo_data'))
                info.append([relPath, lines[0], lines[1], lines[2]])
                if not lines[2] in sectionList:
                    sectionList.append(lines[2])
            except:
                pass

        text = '''<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "DTD/xhtml1-transitional.dtd" >
<html>
<head>
   <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1"/>
   <meta name="AUTHOR" content="Liz"/>
   <title>Demo data overview</title>
<link rel="stylesheet" type="text/css" href="../docs/ccp4i2.css" title="CCP4i2" />
<link rev="made" href="mailto:ccp4gui@ccp4.ac.uk" />
</head>
<body>

<h1>CCP4i2 demo data</h1>

<p>The demo data sets provide a set of examples that you can copy to your project directory and that can be solved quickly as an introduction to CCP4i2.</p>

'''
        for section in sectionList:
            text = text + '<h4>' + section + '</h4>\n<ul>\n'
            for relPath,title,desc,sec in  info:
                if sec == section:         
                    text = text + '''<li><a href="./''' + relPath + '''/README.html">''' + title + '''</a>''' + desc + '''</li>\n\n'''
            text = text + '</ul>\n'
        text = text + '''
<address></address>
<!-- hhmts start -->Automatically created by CCP4i2: ''' + time.strftime("%a, %d %b %Y %H:%M:%S") + '''<!-- hhmts end -->
</body>
</html>
'''
        CCP4Utils.saveFile(os.path.join(CCP4Utils.getDotDirectory(),'demo_data','README.html'),text,overwrite=True)
        self.builtInfo=True


    def makeDemoData(self,parentWidget):
        dirPath = QtWidgets.QFileDialog.getExistingDirectory(parentWidget,'Select directory containing demo data').__str__()
        if not os.path.isdir(dirPath): return
        warningText = ''
        if not os.path.exists(os.path.join(dirPath,'README.html')):
            CCP4Utils.saveFile(os.path.join(dirPath,'README.html'),'''<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "DTD/xhtml1-transitional.dtd" >
<html>
<head>
   <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1"/>
   <meta name="AUTHOR" content="Martin"/>
   <title>Demo data: ***dataset name****</title>
<link rel="stylesheet" type="text/css" href="../../docs/ccp4i2.css" title="CCP4i2" />
<link rev="made" href="mailto:ccp4gui@ccp4.ac.uk" />
</head>
<body>

<h1>CCP4i2 demo data: ***Short title***</h1>

<p>***Short description***</p>

<h2>Data files</h2>
Data files for this tutorial are located in the ccp4i2/test/data/***dataset name**** directory.  When browsing to find files, the browser widget has a pull-down labelled "look in" which is populated by the projects known to your instance of ccp4i2 and (at the bottom) by CCP4I2_TOP: a shortcut to the ccp4i2 code base containing the demo data.
<ul>
<li><b>filename</b>Description</li>
</ul>

<address></address>
<!-- hhmts start -->Last modified: Thu Feb 26 14:37:05 GMT 2015 <!-- hhmts end -->
</body>
</html>
''')
            warningText = warningText + '''There was no README.html file in the directory.
A template README.html file has been added - please edit in details of the demo.\n\n'''
        if not os.path.exists(os.path.join(dirPath,'INFO.txt')):
            CCP4Utils.saveFile(os.path.join(dirPath,'INFO.txt'),'''First line is name of dataset
Second line is brief description for the overview page $CCP4I2/demo_data/README.html\n''')
            warningText = warningText + '''There was no INFO.txt file in the directory.
A template INFO.txt file has been added - please edit in details of the demo.\n\n'''
        if len(warningText) > 0:
            QtWidgets.QMessageBox.warning(parentWidget,'Compress demo data','Data not saved - some files missing:\n\n'+warningText)
            return
        try:
            import zipfile
            zip = zipfile.ZipFile(dirPath+'.ccp4i2_demo.zip',mode='w')
            CCP4Utils.zipDirectory(zip,dirPath,rootRelPath=os.path.split(dirPath)[0])
            zip.close()
        except:
            QtWidgets.QMessageBox.warning(parentWidget,'Compress demo data','Failed writing file:'+dirPath+'.ccp4i2_demo.zip')
        else:      
            QtWidgets.QMessageBox.information(parentWidget,'Compress demo data','Compressed demo data saved to'+dirPath+'.ccp4i2_demo.zip')


    def downloadDemoData(self,parentWidget=None,projectId=None):
        if getattr(self,'downloadDialog',None) is None:
            self.downloadDialog = CDownloadDemoDataDialog(parentWidget)
        self.downloadDialog.show()
        self.downloadDialog.raise_()
    

class CDownloadDemoDataDialog(QtWidgets.QDialog):
    ERROR_CODES = { 101 : {'description' : 'Failed to open compressed demo data file'},
                    102 : {'description' : 'Failed to extract data from compressed demo data file'},
                    103 : {'description' : 'Unknown error extracting data from compressed demo data file'},
                    104 : {'description' : 'Error downloading compressed demo data file'},
                    110 : {'description' : 'Failed reading URL'}}

    def __init__(self,parent):
        import functools
        from qtgui import CCP4Widgets
        self.path = ''
        self.makeTargetList()
        self.zipList = []
        QtWidgets.QDialog.__init__(self,parent)
        self.setWindowTitle('Download demo data')
        self.setModal(True)
        self.setLayout(QtWidgets.QVBoxLayout())
        self.layout().setContentsMargins(4,4,4,4)
        self.layout().setSpacing(4)
          
        line = QtWidgets.QHBoxLayout()
        self.layout().addLayout(line)
        line.addWidget(QtWidgets.QLabel('Enter URL for web page that has links to CCP4i2 demo data',self))
        line.itemAt(0).widget().setObjectName('italic')
        line = QtWidgets.QHBoxLayout()
        self.layout().addLayout(line)
        line.addWidget(QtWidgets.QLabel('Web page',self))
        self.source = QtWidgets.QComboBox(self)
        self.source.setMinimumWidth(400)
        self.source.setEditable(True)
        #for source in [ 'http://www.ysbl.york.ac.uk/~eap5/downloadPage.html' ]:
        for source in [ ]:
            self.source.addItem(source)
        line.addWidget(self.source)
        self.source.editTextChanged.connect(self.handleSourceChanged)
        self.search = QtWidgets.QPushButton('Search page',self)
        self.search.setToolTip('Search the page for downloadable compressed demo data files')
        self.search.clicked.connect(self.handleHtmlSearch)
        line.addWidget(self.search)
        line = QtWidgets.QHBoxLayout()
        line.addWidget(QtWidgets.QLabel('Download to:',self))
        self.targetButGroup = QtWidgets.QButtonGroup(self)
        id = -1
        for label,path,useable in self.targetList:
            id = id + 1
            if useable:
                rb = QtWidgets.QRadioButton(label,self)
                rb.setToolTip(path)
                self.targetButGroup.addButton(rb,id)
                line.addWidget(rb)
        self.targetButGroup.button(0).setChecked(True)
        self.layout().addLayout(line)
    
        line = QtWidgets.QHBoxLayout()
        line.addWidget(QtWidgets.QLabel('Downloadable demos..',self))
        line.itemAt(0).widget().setObjectName('italic')
        self.layout().addLayout(line)
          
        frame = QtWidgets.QFrame(self)
        frame.setFrameShape(QtWidgets.QFrame.Box)
        frame.setLineWidth(2)
        frame.setLayout(QtWidgets.QVBoxLayout())
        self.frame = QtWidgets.QFrame()
        self.frame.setLayout(QtWidgets.QVBoxLayout())
        frame.layout().addWidget(self.frame)
        self.seleLine = QtWidgets.QFrame()
        self.seleLine.setLayout(QtWidgets.QHBoxLayout())
        for butText,mode in [['Select all',True],['Unselect all',False]]:
            pb = QtWidgets.QPushButton(butText,self)
            pb.clicked.connect(functools.partial(self.selectAll,mode))
            self.seleLine.layout().addWidget(pb)
        frame.layout().addWidget(self.seleLine)    
        self.layout().addWidget(frame)
        self.setSelectVis()
        butBox = QtWidgets.QDialogButtonBox(self)
        but = butBox.addButton('Download',QtWidgets.QDialogButtonBox.ApplyRole)
        but.setDefault(False)
        but.clicked.connect(self.handleDownload)
        but = butBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
        but.clicked.connect(self.close)
        line = QtWidgets.QHBoxLayout()
        line.addWidget(butBox)
        self.layout().addLayout(line)

    def makeTargetList(self):
        # Do we have write access to the ccp4i2 source - otherwise use dot dir
        self.targetList = []
        fi = QtCore.QFileInfo(os.path.join(CCP4Utils.getCCP4I2Dir(),'demo_data'))
        self.targetList.append(['CCP4i2 installation',os.path.join(CCP4Utils.getCCP4I2Dir(),'demo_data'),fi.isDir() and fi.isWritable()])
        target = os.path.join(CCP4Utils.getDotDirectory(),'demo_data')
        if not os.path.exists(target): os.mkdir(target)
        self.targetList.append(['My CCP4i2 area',target,True])

    def handleFtpSearch(self):
        from ftplib import FTP
        self.host = str(self.source.currentText())
        self.path = ''
        if self.host.count('/'): self.host,self.path = self.host.split('/',1)
        if len(self.path)>0 and not self.path.endswith('/'): self.path = self.path + '/'
        self.path = self.path + '*.ccp4i2_demo.zip'
        #print 'CDownloadDemoData.handleSearch',self.host,self.path
        self.ftp = FTP(self.host)   # connect to host, default port
        self.ftp.login()               # user anonymous, passwd anonymous@
        self.retrlines = []
        self.ftp.retrlines('NLST '+self.path,self.handleRetrlines)
        #print 'CDownloadDemoData handleSearch retrlines',self.retrlines

    def handleRetrlines(self,line):
        self.retrlines.append(line)

    @QtCore.Slot()
    def handleHtmlSearch(self):
        if sys.version_info >= (3,0):
            import urllib.request, urllib.error, urllib.parse
        else:
            import urllib2
        from lxml import etree
        self.path = str(self.source.currentText())
        try:
            if sys.version_info >= (3,0):
                urlpath =urllib.request.urlopen(self.path)
            else:
                urlpath =urllib2.urlopen(self.path)
        except Exception as e0:
            if not self.path.startswith('http'):
                self.path = 'http://'+self.path
                try:
                    if sys.version_info >= (3,0):
                        urlpath =urllib.request.urlopen(self.path)
                    else:
                        urlpath =urllib2.urlopen(self.path)
                except Exception as e1:
                    raise CException(self.__class__,110,str(e0))
                else:
                    self.source.setEditText(self.path)
        pageText = urlpath.read().decode('utf-8')
        urlpath.close()
        pageNode = etree.fromstring(pageText)
        self.zipList = []
        for node in pageNode.iterfind('.//a'):
            if node.get('href') is not None and 'ccp4i2_demo.zip' in node.get('href'):
                self.zipList.append([node.get('href'),str(node.text)])
        #print 'handleHtmlSearch',self.zipList
        self.loadFrame()

    def loadFrame(self):
        for w in self.frame.findChildren(QtWidgets.QCheckBox):
            w.hide()
            w.deleteLater()
        for link, desc in self.zipList:
            cb = QtWidgets.QCheckBox(desc,self)
            self.frame.layout().addWidget(cb)
        self.setSelectVis()

    def setSelectVis(self):
        vis = len(self.zipList)>0
        for w in self.seleLine.findChildren(QtWidgets.QPushButton):
            w.setVisible(vis)

    @QtCore.Slot(bool)
    def selectAll(self,mode):
        for w in self.frame.findChildren(QtWidgets.QCheckBox):
            w.setChecked(mode)

    @QtCore.Slot()
    def handleDownload(self):
        if sys.version_info >= (3,0):
            import urllib.parse
        else:
            import urlparse
        downloadList = []
        err = CErrorReport()
        ii = -1
        for w in self.frame.findChildren(QtWidgets.QCheckBox):
            ii = ii + 1
            if w.isChecked() and not 'DOWNLOADED' in str(w.text()):
                if sys.version_info >= (3,0):
                    downloadList.append([w,urllib.parse.urljoin(os.path.split(self.path)[0]+'/',self.zipList[ii][0])])
                else:
                    downloadList.append([w,urlparse.urljoin(os.path.split(self.path)[0]+'/',self.zipList[ii][0])])
        #print 'handleDownload',self.path,downloadList
        if len(downloadList)==0:
            QtWidgets.QMessageBox.information(self,'Download demo data','No demo sets selected for download')
            return
        target = self.targetList[self.targetButGroup.checkedId()][1]
        for w,url in downloadList:
            try:
                self.downloadFile(url,target)
            except CException as e:
                err.extend(e)
            except Exception as e:
                err.append(self.__class__,103,'Extracting from '+url+' to '+target+'\n'+str(e))
            else:
                w.setText(str(w.text())+' DOWNLOADED')
                w.repaint()
        CCP4Modules.DEMODATAMANAGER().makeDemoDataInfo()
        self.testDatasets = None
        self.saveTestDatasets()
        if len(err) > 0:
            err.warningMessage('Download demo data','Error downloaing demo data',parent=self)
        else:
            self.close()

    def downloadFile(self,url,target=None,replace=False):
        if sys.version_info >= (3,0):
            import urllib.request, urllib.error, urllib.parse
        else:
            import urllib2
        import tempfile
        import zipfile
        tmpFile =tempfile.mktemp(suffix='ccp4i2_demo.zip')
        demoDir =   os.path.normpath(os.path.join(target,os.path.splitext(url.split('/')[-1])[0]))
        print('Downloading demo data from',url,'to',tmpFile,'unpack to',demoDir)
        if os.path.exists(demoDir):
            if replace:
                try:
                    shutil.move(demoDir,fileName+'.bak')
                except:
                    pass
            else:
                return
        #url = 'http://www.ysbl.york.ac.uk/~eap5/foo.ccp4i2_demo.zip'
        try:
            if sys.version_info >= (3,0):
                req = urllib.request.urlopen(url)
            else:
                req = urllib2.urlopen(url)
            fp = open(tmpFile,'wb')
            shutil.copyfileobj(req, fp)
            req.close()
            fp.close()
        except:
            raise CException(self.__class__,104,tmpFile)
        try:
            czip = zipfile.ZipFile(tmpFile,mode='r')
        except:
            raise CException(self.__class__,101,tmpFile)
        #for zinfo in czip.infolist():
        #  print 'zinfo',zinfo.filename,zinfo.date_time     
        try:
            czip.extractall(path=target)
        except:
            czip.close()
            raise CException(self.__class__,102,'Extracting from '+tmpFile+' to '+target)
        czip.close()

    @QtCore.Slot('QComboBox')
    def handleSourceChanged(self,source):
        if str(self.source.currentText()) != self.path:
            self.path = ''
            for w in self.frame.findChildren(QtWidgets.QCheckBox):
                w.hide()
                w.deleteLater()
        self.setSelectVis()
