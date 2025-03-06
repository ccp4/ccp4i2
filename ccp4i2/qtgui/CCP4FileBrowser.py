from __future__ import print_function

"""
     qtgui/CCP4FileBrowser.py: CCP4 Gui Project
     Copyright (C) 2001-2008 University of York, CCLRC
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
'''
Liz Potterton Mar 10 - Copied from ccp4mg/qtgui/MGWidgets.py
Liz Potterton Apr 12 - Rewrite from basics to fix the failure of select by double clicking
'''

##@package CCP4FileBrowser (QtGui) File browser - based on QtFileDialog but aware of ccp4 projects

import os
import re
import sys
import types
import functools
import requests
import datetime
from PySide2 import QtGui, QtWidgets,QtCore
from core import CCP4Modules
from core import CCP4File

class GenericWorker(QtCore.QObject):

    start = QtCore.Signal()
    finished = QtCore.Signal()

    def __init__(self, function, *args):
        super(GenericWorker, self).__init__()

        self.function = function
        self.args = args
        self.start.connect(self.run)

    @QtCore.Slot()
    def run(self):
        self.function(*self.args)
        self.finished.emit()

class CFileDialog(QtWidgets.QDialog):

    selectDownloadedFile = QtCore.Signal(str,dict)
    selectFile = QtCore.Signal(str)
    selectFiles = QtCore.Signal(tuple)

    NewDirectory = 101
    NewFolder = 102
    
    def __init__(self, parent=None, title='Open file', filters=[], defaultSuffix=None, defaultFileName='',
                 fileMode=QtWidgets.QFileDialog.ExistingFile, fileLabel=None, **kw):
        QtWidgets.QDialog.__init__(self, parent)
        if CCP4Modules.PREFERENCES().NATIVEFILEBROWSER:
            self.input = {'title' : title, 'filters' : filters, 'defaultSuffix':defaultSuffix,
                          'defaultFileName' : defaultFileName, 'fileMode' : fileMode,
                          'fileLabel' : fileLabel, 'kw' : kw}
            return
        self.setWindowTitle(title)
        self.useGetExistingDirectory = False
        layout = QtWidgets.QGridLayout()
        self.widget = CFileDialog1(self, fileLabel=fileLabel, projectCombo=kw.get('projectCombo', True))
        layout.addWidget(self.widget, 1, 0)
        self.setLayout(layout)
        if fileMode in [QtWidgets.QFileDialog.ExistingFile, QtWidgets.QFileDialog.ExistingFiles, QtWidgets.QFileDialog.Directory]:
            self.widget.fileDialog.setFileMode(fileMode)
        elif fileMode == CFileDialog.NewDirectory:
            self.widget.fileDialog.setFileMode(QtWidgets.QFileDialog.AnyFile)
            saveButton = self.widget.fileDialog.saveButton()
            if saveButton is not None:
                saveButton.setText('Create')
        elif fileMode == CFileDialog.NewFolder:
            self.widget.fileDialog.setFileMode(QtWidgets.QFileDialog.Directory)
            self.useGetExistingDirectory = True
        elif fileMode == QtWidgets.QFileDialog.AnyFile:
            self.widget.fileDialog.selectFile(defaultFileName)
            self.widget.fileDialog.setFileMode(QtWidgets.QFileDialog.AnyFile)
            self.widget.fileDialog.setAcceptMode(QtWidgets.QFileDialog.AcceptSave)
            self.widget.fileDialog.setConfirmOverwrite(False)
            if defaultSuffix is not None:
                self.widget.fileDialog.setDefaultSuffix(defaultSuffix)
        if fileMode in [QtWidgets.QFileDialog.Directory, CFileDialog.NewDirectory, CFileDialog.NewFolder]:
            self.widget.fileDialog.setOptions(QtWidgets.QFileDialog.ShowDirsOnly|QtWidgets.QFileDialog.DontUseNativeDialog)
        else:
            self.widget.fileDialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog)
        if 'saveButtonText' in kw:
            saveButton = self.widget.fileDialog.saveButton()
            if saveButton is not None:
                saveButton.setText(kw['saveButtonText'])
        self.suffix_lookup = {}
        filters = self.makeSuffixLookup(filters, addAll=kw.get('addAll', True))
        filters.append('All files ( *.* )')
        self.widget.fileDialog.setNameFilters(filters)
        if self.suffix_lookup:
            self.setSuffix(filters[0][0])
            self.filterSelected.connect(self.setSuffix)
        self.widget.fileDialog.filesSelected.connect(self.handleFilesSelected)
        self.widget.fileDialog.rejected.connect(self.handleRejected)

    @QtCore.Slot()
    def handleRejected(self):
        self.rejected.emit()
        self.close()
    
    def setFileMode(self,arg):
        self.widget.fileDialog.setFileMode(arg)

    def setAcceptMode(self,arg):
        self.widget.fileDialog.setAcceptMode(arg)
    
    def setConfirmOverwrite(self,arg):
        self.widget.fileDialog.setConfirmOverwrite(arg)

    @QtCore.Slot(str)
    def setSuffix(self, cfilter):
        cfilter = str(cfilter)
        #print "setSuffix",cfilter
        if cfilter in self.suffix_lookup:
            self.widget.fileDialog.setDefaultSuffix(self.suffix_lookup[cfilter])

    @QtCore.Slot(list)
    def handleFilesSelected(self,selectedFiles):
        files = self.filesToPy(selectedFiles)
        if files[0] != self.widget.download.get('fileName',''):
            self.widget.download = {}
        self.selectDownloadedFile.emit(files[0], self.widget.download)
        self.selectFile.emit(files[0])
        self.selectFiles.emit(files)
        self.close()

    def handleSaveButton(self):
        selectedFiles = self.widget.fileDialog.selectedFiles()
        files = self.filesToPy(selectedFiles)
        print('CFileDialog.handleSaveButton',files)
        if files[0] != self.widget.download.get('fileName', ''):
            self.widget.download = {}
        self.selectDownloadedFile.emit(files[0], self.widget.download)
        self.selectFile.emit(files[0])
        self.selectFiles.emit(files)
        self.close()

    def filesToPy(self,selectedFiles):
        files = []
        for ff in selectedFiles:
            f = str(ff)
            # QFileDialog should have added suffix for us if the defaultSuffix has been set
            files.append(f)
        return files

    def isVisible(self):
        return self.widget.fileDialog.isVisible()

    def show(self):
        if CCP4Modules.PREFERENCES().NATIVEFILEBROWSER:
            # OpenFileName with double extension messing up see: https://bugreports.qt.io/browse/QTBUG-44227
            filterText, secondExt = self.makeFilterText(self.input['filters'], truncateExt=(self.input['fileMode'] == QtWidgets.QFileDialog.AnyFile))
            if self.input['fileMode'] == QtWidgets.QFileDialog.ExistingFile:
                fileName,selectedFilter = QtWidgets.QFileDialog.getOpenFileName(self, self.input['title'], '', filterText)
            elif self.input['fileMode'] == QtWidgets.QFileDialog.AnyFile:
                #This should set the default output filename, but it is broken in Qt 4.8 on Mac (in C++). It does work in Qt5.7.
                fileName,selectedFilter = QtWidgets.QFileDialog.getSaveFileName(self, self.input['title'], self.input['defaultFileName'], filterText)
            else:  # This is the 'existingDirectory' option
                fileName = QtWidgets.QFileDialog.getExistingDirectory(self, self.input['title'], '')
            if fileName is not None and len(str(fileName)) > 0:
                fileName = str(fileName)
                if self.input['fileMode'] == QtWidgets.QFileDialog.AnyFile and secondExt is not None and not fileName.endswith(secondExt):
                    fileName = fileName + '.' + secondExt
                print('CFileDialog.show', fileName)
                self.selectFile.emit(fileName)
                self.selectFiles.emit(fileName)
                self.close()
            return
        if self.useGetExistingDirectory:
            rv = self.widget.fileDialog.getExistingDirectory(self)
            fileName = str(rv)
            self.selectFile.emit(fileName)
            self.selectFiles.emit(fileName)
            self.close()
        else:
            if sys.platform == "darwin":
                for t in self.widget.fileDialog.findChildren(QtWidgets.QAbstractItemView):
                    t.setIconSize(QtCore.QSize(16,16))
            
            self.widget.fileDialog.show()
            QtWidgets.QDialog.show(self)

    def makeFilterText(self, filterList, truncateExt=False):
        filterText = ''
        secondExt = None
        for ff in filterList:
            if not truncateExt:
                filterText = filterText + ff + ';;'
            else:
                try:
                    m = re.match(r'(.*)\((.*)\)', ff)
                    lab,exTxt = m.groups()
                    lab = lab + '('
                    exList = exTxt.split()
                    if exList[0].count('.') > 1:
                        secondExt = exList[0].split('.')[-1]
                    for ex in exList:
                        lab = lab + '*.'+ ex.split('.')[1] + ' '
                    filterText = filterText + lab[0:-1] + ');;'
                except Exception as e:
                    print('makeFilterText error', e)
                    filterText = filterText + ff + ';;'
        print('makeFilterText', filterList, filterText[0:-2], 'secondExt', secondExt)
        return filterText[0:-2], secondExt

    def makeSuffixLookup(self, filters, addAll=True):
        if len(filters) == 0:
            return []
        all_filters = "All files ("
        new_filters = []
        for filt in filters:
            #print "setNameFilters filt", filt
            if isinstance(filt,list):
                self.suffix_lookup[filt[0]] = filt[1]
                cfilter = str(filt[0])
            else:
                cfilter = str(filt)
            if len(cfilter) > 11 and cfilter[:11] == 'All files (' and (cfilter).rfind(')') > 0:
                pass
            else:
                new_filters.append(cfilter)
                left_par = cfilter.rfind('(')
                right_par = cfilter.rfind(')')
                if left_par > -1 and right_par > 0 and right_par > left_par:
                    if (cfilter[left_par + 1:right_par]).find('*') > -1:
                        all_filters = all_filters + " " + cfilter[left_par + 1:right_par]
                    else:
                        all_filters = all_filters + " " + cfilter
                else:
                    all_filters = all_filters + " " + cfilter
        all_filters = all_filters + " )"
        if addAll and len(new_filters) > 1:
            new_filters.append(all_filters)
        return new_filters

    def setDownloadMode(self, modeList=['ebiPdb'], projectId=None):
        self.widget.drawDownload(modeList=modeList, projectId=projectId)
    
    def addWidget(self, widget):
        self.layout().addWidget(widget, 0, 0)


class CFileDialog1(QtWidgets.QWidget):

    selectDownloadedFile = QtCore.Signal(str,dict)

    def __init__(self, parent=None, projectCombo=True, downloadOnlyMode=False, fileLabel=None):
        QtWidgets.QWidget.__init__(self,parent)
        self.DOWNLOAD_DEFINITIONS = {'ebiPdb' : {'label' : 'Download PDB from EBI-PDBe. 4-letter code:',
                                                 'url' : 'https://www.ebi.ac.uk/pdbe/entry-files/download/pdbCODE.ent',
                                                 'page' : 'www.ebi.ac.uk/pdbe', 'rename' : 'CODE.pdb', 'ext' : '.pdb'},
                                     'ebiSFs' : {'label' : 'Download structure factors from EBI-PDBe. 4-letter code:',
                                                 'url' : 'https://www.ebi.ac.uk/pdbe/entry-files/download/rCODEsf.ent',
                                                 'page' : 'www.ebi.ac.uk/pdbe', 'ext' : '.cif' },
                                     'Uppsala-EDS' : {'label' : 'Download map coefficients from Uppsala Electron Density Server.  4-letter code:' ,
                                                      'url' : 'http://eds.bmc.uu.se/eds/dfs/CODE13/CODE/CODE_sigmaa.mtz',
                                                      'page' : 'http://eds.bmc.uu.se/eds/', 'ext' : '_sigmaa.mtz' },
                                     'rcsbPdb' : {'label' : 'Download PDB from RCSB. 4-letter code:',
                                                  'url' :  'https://files.rcsb.org/download/CODE.pdb',
                                                  'page' : 'www.rcsb.org'},
                                     'uniprotFasta' : {'label' : 'Download sequence from uniprot database',
                                                       'url':'http://www.uniprot.org/uniprot/CODE.fasta',
                                                       'page' : 'www.uniprot.org', 'ext' : '.fasta'},
                                     'uniprotAFPdb' : {'label' : 'Download AlphaFold PDB File from Protein Structure db',
                                                       'url' : 'https://alphafold.ebi.ac.uk/files/AF-CODE-F1-model_VERSI.pdb',
                                                       'page' : 'alphafold.com', 'ext': '.pdb'} }
        self.projectId = None
        self.fileLabel = fileLabel
        self.downloadOnlyMode = downloadOnlyMode
        self.downloader = None
        self.download = {}
        self.downloadProgressFrame= None
        self.downloadProgress = None
        layout = QtWidgets.QVBoxLayout()
        layout.setSpacing(2)
        layout.setContentsMargins(2, 2, 2, 2)
        self.setLayout(layout)
        # This is the weird bit - CFileDialog0 must have the CFileDialog as parent
        # for it to appear in the CFileDialog diaolog box
        if not self.downloadOnlyMode:
            if projectCombo:
                self.drawProjectCombo()
            self.fileDialog = CFileDialog0(parent)
            layout.addWidget(self.fileDialog)
        self.ccp4i_combo = None

    def drawProjectCombo(self):
        ccp4_dirs, ccp4_aliases = CCP4Modules.PROJECTSMANAGER().getProjectsList()
        #print 'CFileDialog.__init__', ccp4_dirs,ccp4_aliases
        #if len(ccp4_dirs) > 0 or len(ccp4_aliases) > 0:
        layout0 = QtWidgets.QHBoxLayout()
        self.projectCombo = QtWidgets.QComboBox()
        self.projectCombo.setEditable(False)
        self.loadProjectCombo()
        ccp4i_label1 = QtWidgets.QLabel()
        ccp4i_label1.setText("Look in ccp4 project:")
        ccp4i_label2 = QtWidgets.QLabel(" or ...")
        layout0.addWidget(ccp4i_label1)
        layout0.addWidget(self.projectCombo)
        layout0.addWidget(ccp4i_label2)
        layout0.addStretch()
        self.layout().addLayout(layout0)
        self.projectCombo.currentIndexChanged[str].connect(self.projectComboChanged)
        CCP4Modules.PROJECTSMANAGER().projectsListChanged.connect(self.loadProjectCombo)

    @QtCore.Slot()
    def loadProjectCombo(self):
        self.projectCombo.clear()
        self.projectCombo.addItem('Full path..')
        ccp4_dirs, ccp4_aliases = CCP4Modules.PROJECTSMANAGER().getProjectsList()
        for project in ccp4_dirs:
            self.projectCombo.addItem(project)
    
    @QtCore.Slot(str)
    def projectComboChanged(self, alias):
        #print('projectComboChanged',alias)
        path = CCP4Modules.PROJECTSMANAGER().getProjectDirectory(projectName=str(alias))
        if path is not None:
            self.fileDialog.setDirectory(str(path))

    def drawDownload(self, modeList=[], projectId=None):
        #print 'CFileDialog1.drawDownload',modeList,projectId
        if modeList is None:
            return
        if not isinstance(modeList, list):
            modeList = [modeList]
        self.projectId = projectId
        self.downloadFrame = QtWidgets.QFrame(self)
        layout0 = QtWidgets.QHBoxLayout()
        self.downloadFrame.setLayout(layout0)
        layout0.setSpacing(1)
        layout0.setContentsMargins(0, 0, 0, 0)
        self.downloadCombo =  QtWidgets.QComboBox()
        self.downloadCombo.setEditable(False)
        for mode in modeList:
            try:
                self.downloadCombo.addItem(self.DOWNLOAD_DEFINITIONS[mode]['label'], mode)
            except:
                pass
        if self.downloadCombo.count() <= 0:
            self.downloadCombo.deleteLater()
            del self.downloadCombo
            return
        layout0.addWidget(self.downloadCombo)
        self.downloadCode =  QtWidgets.QLineEdit(self)
        self.downloadCode.setMinimumWidth(40)
        self.downloadCode.returnPressed.connect(self.handleDownload)
        layout0.addWidget(self.downloadCode)
        downloadButton = QtWidgets.QPushButton('Download', self)
        downloadButton.clicked.connect(self.handleDownload)
        viewButton = QtWidgets.QPushButton('View website', self)
        viewButton.clicked.connect(self.viewWebSite)
        if self.downloadOnlyMode:
            self.butFrame =  QtWidgets.QFrame(self)
            self.butFrame.setLayout(QtWidgets.QHBoxLayout())
            self.butFrame.layout().addWidget(downloadButton)
            self.butFrame.layout().addWidget(viewButton)
            self.layout().addWidget(self.butFrame)
        else:
            layout0.setStretchFactor(self.downloadCombo, 3)
            layout0.addWidget(downloadButton)
            layout0.addWidget(viewButton)
        #self.layout().insertLayout(0,layout0)
        self.layout().insertWidget(0, self.downloadFrame)
        self.diagDone=False

    def drawDownloadProgress(self):
        if self.downloadProgressFrame is None:
            self.downloadProgressFrame = QtWidgets.QFrame(self)
            layout0 = QtWidgets.QHBoxLayout()
            self.downloadProgressFrame.setLayout(layout0)
            layout0.setSpacing(1)
            layout0.setContentsMargins(0, 0, 0, 0)
            layout0.addWidget(QtWidgets.QLabel('Downloading:', self))
            self.downloadProgress = QtWidgets.QProgressBar(self.parent())
            self.downloadProgress.setMinimum(0)
            self.downloadProgress.setMaximum(100)
            layout0.addWidget(self.downloadProgress)
            interrupt = QtWidgets.QPushButton('Stop download', self)
            layout0.addWidget(interrupt)
            interrupt.clicked.connect(self.handleInterruptDownload)
            self.layout().insertWidget(0, self.downloadProgressFrame)
        else:
            self.downloadProgressFrame.show()
            self.downloadProgress.setValue(0)
        self.downloadFrame.hide()
        if self.downloadOnlyMode:
            self.butFrame.hide()

    def removeDownloadProgress(self):
        self.downloadProgressFrame.hide()
        self.downloadFrame.show()

    def downloadFileName(self, urlname, code=None, rename=None):
        from core import CCP4Utils
        if self.projectId is None:
            tmpDir = CCP4Utils.getTMP()
        else:
            tmpDir = os.path.join(CCP4Modules.PROJECTSMANAGER().getProjectDirectory(projectId=self.projectId),'CCP4_DOWNLOADED_FILES')
            if not os.path.exists(tmpDir):
                try:
                    os.mkdir(tmpDir)
                except:
                    tmpDir = CCP4Utils.getTMP()
        if code is not None and rename is not None:
            fileName = os.path.join(tmpDir,re.sub('CODE',code.lower(),rename))
        else:
            fileName = os.path.join(tmpDir,os.path.split(urlname)[1])
        #print 'downloadFileName urlname',urlname,fileName,os.path.exists(fileName)
        return fileName

    @QtCore.Slot()
    def handleDownload(self):
        nowis = datetime.datetime.now()
        cyear = nowis.year
        mode = self.downloadCombo.itemData(self.downloadCombo.currentIndex()).__str__()
        code = self.downloadCode.text().__str__()
        if len(code) < 4:
            warning = QtWidgets.QMessageBox.warning(self, 'Downloading data', 'Please enter code of at least 4 characters')
            return
        urlname = re.sub('CODE13', code[1:3].lower(), self.DOWNLOAD_DEFINITIONS[mode]['url'])
        if mode == "uniprotAFPdb":
            # Get the latest version possible. Limit loop by year, unclear what vrs will exist in future & how many at a time.
            # This code will handle max v5 in 2022, v9 is 2023, v14 afterwards.
            maxvrs = 15
            if cyear == 2022:
                maxvrs = 6
            elif cyear == 2023:
                maxvrs = 10
            urlname = re.sub('CODE', code.upper(), urlname)
            for i in reversed(range(3, maxvrs)):
                ckvrs = "v" + str(i)
                linea = re.sub('VERSI', ckvrs, urlname)
                req = requests.get(linea)
                if req.status_code == 200:
                    urlname = linea
                    break
        else:
            urlname = re.sub('CODE', code.lower(), urlname)
        rename = self.DOWNLOAD_DEFINITIONS[mode].get('rename',None)
        targetFile = self.downloadFileName(urlname, code=code, rename=rename)
        if os.path.exists(targetFile):
            question = QtWidgets.QMessageBox(self)
            question.setWindowTitle('Downloading data')
            question.setText('A file for ' + code + ' already exists.\nDo you want to use it or overwrite it?')
            question.addButton('Use existing file', QtWidgets.QMessageBox.YesRole)
            question.addButton('Delete file and re-download', QtWidgets.QMessageBox.NoRole)
            question.buttonClicked.connect(functools.partial(self.handleQueryExistingFile, code, mode, targetFile))
            question.show()
            question.raise_()
            return
        if self.downloader is None:
            self.drawDownloadProgress()

            self.downloader =  CDownloader()
            self.downloader.Finished.connect(functools.partial(self.handleDownloadFinished, code, mode, targetFile))
            self.downloader.Error.connect(functools.partial(self.handleDownloadError, code))
            self.downloader.ProgressChanged.connect(self.handleProgress)

            if hasattr(self,"my_thread"):
                self.my_thread.quit()

            self.my_thread = QtCore.QThread()
            self.my_thread.start()

            self.my_worker = GenericWorker(functools.partial(self.downloader.download, urlname))
            self.my_worker.moveToThread(self.my_thread)
            self.my_worker.finished.connect(self.my_thread.quit)
            self.my_worker.start.emit()

    @QtCore.Slot()
    def handleInterruptDownload(self):
        #print 'handleInterruptDownload'
        self.downloader.interrupt_dl = True
        self.downloader = None
        self.downloadThread = None
        self.removeDownloadProgress()

    @QtCore.Slot(str,str,str,'QAbstractButton')
    def handleQueryExistingFile(self, code, mode, targetFile, button):
        #print 'handleQueryExistingFile',button.text()
        if button.text().count('Delete'): 
            os.remove(targetFile)
            self.handleDownload()
        else:
            self.fileDialog.setDirectory(os.path.split(targetFile)[0])
            self.fileDialog.selectFile(targetFile)
            self.download = { 'code' : code, 'source' : mode, 'fileName' : targetFile }

    @QtCore.Slot(str,str)
    def handleDownloadError(self,code,message):
        self.downloader = None
        self.downloadThread = None
        self.removeDownloadProgress()
        #print 'handleDownloadError',code,message
        warning = QtWidgets.QMessageBox.warning(self,'Downloading data for ' + code,message)

    @QtCore.Slot(str,str,str,str)
    def handleDownloadFinished(self,code,mode,targetFile,tempFile):
        #print('handleDownloadFinished',code,mode,targetFile,tempFile)
        import shutil
        shutil.copyfile(tempFile,targetFile)
        self.downloader = None
        self.downloadThread = None
        self.removeDownloadProgress()
        self.download = {'code' : code, 'source' : mode, 'fileName' : targetFile}
        if not self.downloadOnlyMode:
            self.fileDialog.setDirectory(os.path.split(targetFile)[0])
            self.fileDialog.selectFile(targetFile)
        else:
            # Beware signal from CFileDialog1 class when used in downloadOnlyMode 
            self.selectDownloadedFile.emit(targetFile, self.download)

    @QtCore.Slot(int)
    def handleProgress(self,frac):
        try:
            self.downloadProgress.setValue(frac)
        except:
            pass

    @QtCore.Slot()
    def viewWebSite(self):
        mode = self.downloadCombo.itemData(self.downloadCombo.currentIndex()).__str__()
        CCP4Modules.WEBBROWSER().loadPage(QtCore.QUrl('http://' + self.DOWNLOAD_DEFINITIONS[mode]['page']), newTab=True)


class IconProvider(QtWidgets.QFileIconProvider):
    def __init__(self):
        QtWidgets.QFileIconProvider.__init__(self)
        self._filters = {}
        self._icon_cache = {}

    def setIconForFilter(self, icon, cfilter):
        self._filters[icon] = cfilter

    def setIconsForFilters(self, filters):
        self._filters = filters

    def icon(self, fileInfo):
        if isinstance(fileInfo, QtCore.QFileInfo):
            suffix = str(fileInfo.completeSuffix())
            if suffix == 'mtz':
                content = fileInfo.baseName().__str__().split(CCP4File.CDataFile.SEPARATOR)[-1]
                mimeType=CCP4Modules.MIMETYPESHANDLER().formatFromFileExt(ext=suffix, contentLabel=content)
            else:
                mimeType=CCP4Modules.MIMETYPESHANDLER().formatFromFileExt(ext=suffix)
            if mimeType is not None:
                icon = CCP4Modules.MIMETYPESHANDLER().icon(mimeType)
                if icon is not None:
                    return icon
        return QtWidgets.QFileIconProvider.icon(self, fileInfo)

    '''
    def icon(self,info):
        #print 'IconProvider.icon',info.filePath(),self._filters
        try:
            if info and hasattr(info,'filePath'):
                for k in self._filters:
                    for f in self._filters[k]:
                        if str(info.filePath()).endswith(f):
                            if not self._icon_cache.has_key(k):
                                self._icon_cache[k] = QtGui.QIcon(k)
                            return self._icon_cache[k]
        except:
            pass
        return QtWidgets.QFileIconProvider.icon(self,info)
    '''


class ProxyModel(QtCore.QSortFilterProxyModel):

    def __init__(self, parent=None):
        QtCore.QSortFilterProxyModel.__init__(self, parent)

    def filterAcceptsRow(self, sourceRow, sourceParent):
        if sourceParent.row() == -1 or sourceParent.column():
            return True
        index0 = self.sourceModel().index(sourceRow, 0, sourceParent)
        if self.sourceModel().isDir(index0):
            return True
        for f in self.sourceModel().nameFilters():
            if f == '*.*':
                return True
            try:
                if str(self.sourceModel().fileName(index0)).endswith(str(f).lstrip('*')):
                    return True
            except:
                pass
        return False


class CFileDialog0(QtWidgets.QFileDialog):

    def __init__(self, parent):
        QtWidgets.QFileDialog.__init__(self, parent)
        self.setSizeGripEnabled(0)
        iconProvider = IconProvider()
        #iconProvider.setIconsForFilters(CCP4Modules.MIMETYPESHANDLER().getIconsForFileFilters())
        self.setIconProvider(iconProvider)
        self.cleanupSidebar()

    def cleanupSidebar(self):
        '''Remove invalid directories from left side-bar'''
        # These can cause whole application to stall
        try:
            urlList = self.sidebarUrls()
        except:
            print('Error cleaning up side bar - no urlList')
            return
        newUrlList = []
        for url in urlList:
            try:
                if os.path.exists(str(url.path())):
                    newUrlList.append(url)
                else:
                    print('Removing non-existant directory from sidebar:', str(url.path()))
            except:
                try:
                    print('Removing problem directory from side bar:', str(url.path()))
                except:
                    print('Removing problem directory from side bar - could not print name')
        if len(newUrlList) != urlList:
            try:
                self.setSidebarUrls(newUrlList)
                print('Done resetting side bar')
            except:
                print('Error cleaning up side bar in setSidebarUrls')

    def saveButton(self):
        children = self.findChildren(QtWidgets.QPushButton, "")
        for child in children:
            if child.text() in [self.tr("&Save"), self.tr("&Choose"), self.tr("&Open")]:
                return child
        return None


class CDownloader(QtCore.QObject):

    Error = QtCore.Signal(str)
    Finished = QtCore.Signal(str)
    ProgressChanged = QtCore.Signal(int)
    ProgressBarRangeChanged = QtCore.Signal(int,int)

    def __init__(self):
        QtCore.QObject.__init__(self)
        self.interrupt_dl = False

    def download(self, url=None):
        if sys.version_info >= (3,0):
            import urllib.request, urllib.error, urllib.parse
        else:
            import urllib2
        import tempfile
        if not url:
            self.Error.emit('No file to download' )
            self.ProgressChanged.emit(0)
            return
        # Beware localFile is tuple of file handle and path name
        self.localFile = tempfile.mkstemp()
        try:
            # setup_proxies currently always returns None
            opener = setup_proxies()
            if opener:
                dl = opener.open(url)
            else:
                if False and str(url).startswith("https://"): #I get certificate problems with https and ccp4-python
                    import ssl
                    context = ssl._create_unverified_context()
                    if sys.version_info >= (3,0):
                        dl = urllib.request.urlopen(url, context=context)
                    else:
                        dl = urllib2.urlopen(url, context=context)
                else:
                    if sys.version_info >= (3,0):
                        dl = urllib.request.urlopen(url)
                    else:
                        dl = urllib2.urlopen(url)
        except:
            exc_type, exc_value, exc_tb = sys.exc_info()[:3]
            sys.stderr.write(str(exc_type) + '\n')
            sys.stderr.write(str(exc_value) + '\n')
            self.Error.emit('Error opening url ' + str(url) + '\nPlease check that you have internet connection.\nand that the file id is valid.\nAlso the data may not exist on the server.')
            return
        try:
            info = dl.info()
            contentLength = -1
            if sys.version_info > (3,0):
                for h in info.items():
                    if h[0].lower() == "content-length":
                        try:
                            cl = h[1]
                            contentLength = int(cl)
                            self.ProgressBarRangeChanged.emit(0, 100)
                        except:
                            self.ProgressBarRangeChanged.emit(0, 0)
                            self.ProgressChanged.emit(0)
                            contentLength = -1
                        break
            else:
                for h in info.headers:
                    if h.find('Content-Length:') > -1:
                        cl = h.split(':')[1].strip()
                        try:
                            contentLength = int(cl)
                            self.ProgressBarRangeChanged.emit( 0, 100)
                        except:
                            self.ProgressBarRangeChanged.emit(0, 0)
                            self.ProgressChanged.emit(0)
                            contentLength = -1
                        break
            readBytes = 0
#FIXME - Bytes might be fine with Python2 too.
            if sys.version_info > (3,0):
                cbuffer = b""
            else:
                cbuffer = ""
            CHUNKSIZE = 2097152
            read = dl.read(CHUNKSIZE)
            while len(read) > 0 and not self.interrupt_dl:
                os.write(self.localFile[0], read)
                cbuffer = cbuffer + read
                readBytes = readBytes + len(read)
                pc = int(float(readBytes)/contentLength * 100)
                if contentLength > -1:
                    self.ProgressChanged.emit(pc)
                read = dl.read(CHUNKSIZE)
            dl.close()
            if not self.interrupt_dl:
                fname = self.localFile[1]
                os.close(self.localFile[0])
                self.ProgressBarRangeChanged.emit(0, 100)
                self.ProgressChanged.emit(0)
                self.Finished.emit(fname)
                # Just in case
                self.interrupt_dl = False
            self.interrupt_dl = False
        except:
            self.interrupt_dl = False
            print("CCP4FileBrowser.CDownloader Error getting", url)
            raise
            exc_type, exc_value = sys.exc_info()[:2]
            print(exc_type)
            print(exc_value)
            os.close(self.localFile[0])
            self.Error.emit('Error getting file' + str(exc_value))
            self.ProgressChanged.emit(0)

def setup_proxies():
    return None
    if sys.version_info >= (3,0):
        import urllib.request, urllib.error, urllib.parse
    else:
        import urllib2
    from global_definitions import PM
    proxy_uri = PM('download_preferences').get('http_proxy')
    user = PM('download_preferences').get('http_user')
    passwd =  PM('download_preferences').get('http_password')
    if proxy_uri == "":
        return None
    if sys.version_info >= (3,0):
        proxy = urllib.request.ProxyHandler({'http': proxy_uri})
        proxy_auth_handler = urllib.request.HTTPBasicAuthHandler()
    else:
        proxy = urllib2.ProxyHandler({'http': proxy_uri})
        proxy_auth_handler = urllib2.HTTPBasicAuthHandler()
    # No idea if this is correct
    if user != "":
        proxy_auth_handler.add_password(None, uri=proxy_uri, user=user, passwd=passwd)
    if sys.version_info >= (3,0):
        opener = urllib.request.build_opener(proxy, proxy_auth_handler)
    else:
        opener = urllib2.build_opener(proxy, proxy_auth_handler)
    return opener
