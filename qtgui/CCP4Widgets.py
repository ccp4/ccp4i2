from __future__ import print_function

"""
     qtgui/CCP4Widgets.py: CCP4 Gui Project
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

##@package CCP4Widgets (QtGui) Collection of widgets for simple data types
from PySide6 import QtGui, QtWidgets,QtCore,QtSvg
from core import CCP4Data
from core import CCP4ModelData
from core import CCP4XtalData
from core import CCP4File
from qtgui import CCP4SequenceList
from qtgui import CCP4RefmacMultiAtomSelection
from core.CCP4Modules import PROJECTSMANAGER, PIXMAPMANAGER, QTAPPLICATION, TASKMANAGER, LAUNCHER
from core.CCP4Modules import WEBBROWSER, PREFERENCES, MIMETYPESHANDLER, COMFILEPATCHMANAGER

from core.CCP4ErrorHandling import *
import os
import re
import sys
import shutil
import glob
import functools
import traceback
import json

# Jon Agirre - 1 Feb 2019: Increased DRAGICONSIZE from 16 to 20
# as it was impossible to see the dragged icon on some systems
ICONBUTTONSIZE=24
DRAGICONSIZE=20
CONTRASTCOLOUR = '#9BFFFF'
ERRORCOLOUR = '#FF0000'

class MyProxyStyle(QtWidgets.QCommonStyle):
    def __init__(self):
        """Initialize all functions we're not overriding.
        This simply calls the corresponding function in self._style.
        """
        self._style = QtWidgets.QStyleFactory.create('macintosh')
        for method in ['drawComplexControl', 'drawItemPixmap',
                       'generatedIconPixmap', 'hitTestComplexControl',
                       'itemPixmapRect', 'itemTextRect', 'polish', 'styleHint',
                       'subControlRect', 'unpolish', 'drawItemText',
                       'sizeFromContents', 'drawPrimitive','drawControl','subElementRect',]:
            target = getattr(self._style, method)
            setattr(self, method, functools.partial(target))
        QtWidgets.QCommonStyle.__init__(self)

    def pixelMetric(self, metric, option, widget):
        if metric == QtWidgets.QStyle.PM_SmallIconSize or metric == QtWidgets.QStyle.PM_LargeIconSize:
            return 8
        else:
            return self._style.pixelMetric(metric, option, widget)

def bool_kw_decorator(f):
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        f(*args, **kwargs)
    return wrapper

class CBaseWidget:

    CHARSIZE = None
    STRETCH = -1
    TABLE_EDIT_MODE = None

    def __init__(self, dragType=None):
        self._dragType = dragType
        if self.dragType() is not None:
            self.setAcceptDrops(1)
        if CBaseWidget.CHARSIZE is None:
            someText = "/99/Z/999.2/HG51"
            dummy = QtWidgets.QLineEdit()
            fn = dummy.font()
            fm = QtGui.QFontMetrics(fn)
            CBaseWidget.CHARSIZE = fm.horizontalAdvance(someText) / len(someText)

    def setCharWidth(self, charWidth=None, mode='fixed'):
        #print 'CBaseWidget.setCharWidth',self,charWidth
        #Add small extra (2), widgets never seem wide enough
        if charWidth is not None and charWidth > 0:
            setting = self.CHARSIZE * (charWidth + 2)
            if mode in ('fixed', 'minimum'):
                self.setMinimumWidth(setting)
            if mode in ('fixed', 'maximum'):
                self.setMaximumWidth(setting)

    def dragType(self):
        #print 'CBaseWidget.dragType', self.objectName(), self._dragType
        if self._dragType is None:
            return None
        elif not isinstance(self._dragType, list):
            return self._dragType
        else:
            return self._dragType[0]

    def dragTypeList(self):
        #print 'CBaseWidget.dragTypeList', self.objectName(), self._dragType
        if self._dragType is None:
            return []
        elif not isinstance(self._dragType, list):
            return [self._dragType]
        else:
            return self._dragType

    def dropTypes(self):
        if getattr(self,'model', None) is not None:
            return [str(self.model.__class__.__name__)[1:]]
        else:
            return []

    def getFilesFromMimeData(self, mimeData):
        fileList = []
        if mimeData.hasUrls():
            urlList = mimeData.urls()
            for url in urlList:
                if url.isLocalFile():
                    path = url.toLocalFile()
                    mimeTypeList = MIMETYPESHANDLER().classListFromFileName(path)
                    for mimeType in mimeTypeList:
                        if sys.version_info > (3,0):
                            fileList.append([path, mimeType])
                        else:
                            fileList.append([path, mimeType.encode('ascii', 'ignore')])
        return fileList

    def makeDragEtree(self, path, mimeType):
        from lxml import etree
        rel, base = os.path.split(path)
        fileEle = etree.Element(str(mimeType).encode('ascii', 'ignore'))
        ele = etree.Element('relPath')
        ele.text = rel
        fileEle.append(ele)
        ele = etree.Element('baseName')
        ele.text = base
        fileEle.append(ele)
        text = etree.tostring(fileEle)
        return text

    def acceptMimeData(self, mimeData):
        #print('CBaseData.acceptMimeData', repr(self), self.dragTypeList(), self.dropTypes())
        # for item in mimeData.formats(): print 'CBaseData.acceptMimeData data:',str(item),mimeData.data(item)
        dragTypeList = self.dragTypeList()
        if dragTypeList is not None:
            for dragType in dragTypeList:
                if mimeData.hasFormat(dragType):
                    return mimeData.data(dragType).data()
            if len(self.dropTypes()) > 0:
                for item in self.dropTypes():
                    if mimeData.hasFormat(item):
                        return mimeData.data(item).data()
            elif isinstance(self, CIconButton) and len(self.parent().dropTypes()) > 0:
                for item in self.parent().dropTypes():
                    if mimeData.hasFormat(item):
                        return mimeData.data(item).data()
            elif mimeData.hasFormat('jobId'):
                dropText = mimeData.data('jobId').data()
                from lxml import etree
                tree = etree.fromstring(dropText)
                for groupName in ['outputFiles','outputData']:
                    groupEle = tree.find(groupName)
                    if groupEle is not None:
                        for ele in groupEle:
                            if ele.tag in dragTypeList:
                                return etree.tostring(ele)
        fileList = self.getFilesFromMimeData(mimeData)
        #print('acceptMimeData fileList',fileList, self.dropTypes())
        for path,mimeType in fileList:
            if mimeType in self.dragTypeList():
                return self.makeDragEtree(path,mimeType)
            elif len(self.dropTypes()) > 0:
                for item in self.dropTypes():
                    if mimeType == item:
                        return self.makeDragEtree(path, mimeType)
        return None

    def dragEnterEvent(self,event):
        print("dragEnter")
        if self.acceptMimeData(event.mimeData()) is not None:
            event.accept()
        else:
            event.ignore()

    def dragMoveEvent(self, event):
        if self.acceptMimeData(event.mimeData()) is not None:
            event.setDropAction(QtCore.Qt.CopyAction)
            event.accept()
        else:
            event.ignore()

    def dropEvent(self,event):
        print("Drop!!!!!!!")
        dropData = self.acceptMimeData(event.mimeData())
        if dropData is not None:
            #print 'CBaseWidget.dropEvent',dropData,type(dropData),type(self)
            self.acceptDropDataSignal.emit(dropData)
            event.setDropAction(QtCore.Qt.CopyAction)
            event.accept()
        else:
            event.ignore()

    def makeToolTip(self,qualis):
        if 'toolTip' in qualis and qualis['toolTip'] is not None and qualis['toolTip'] is not NotImplemented:
            tip = qualis['toolTip'] + '.'
        else:
            tip = ''
        for key in ['default', 'min', 'max']:
            if key in qualis and qualis[key] is not None and qualis[key] is not NotImplemented:
                tip = tip + ' '+key+':'+str(qualis[key])
        return tip

    def parentTaskWidget(self):
        from qtgui import CCP4TaskWidget
        from qtgui import CCP4ContainerView
        w = self.parent()
        while w is not None and isinstance(w, QtWidgets.QWidget):
            if isinstance(w, (CCP4TaskWidget.CTaskWidget, CCP4ContainerView.CContainerView)):
                return w
            else:
                w = w.parent()
        return None

    def parentTaskFrame(self):
        from qtgui import CCP4ProjectViewer
        w = self.parent()
        while w is not None and isinstance(w, QtWidgets.QWidget):
            if isinstance(w, (CCP4ProjectViewer.CTaskFrame)):
                return w
            else:
                w = w.parent()
        return None

    def reset(self):
        pass

    @QtCore.Slot()
    @bool_kw_decorator
    def validate(self, isValid=None, **kw):
        pass


class CPixmapManager:
    '''Load all icon source files from qticons and return QPixmap when requested'''
    insts = None
    SIZE = 32

    def __init__(self):
        CPixmapManager.insts = self
        self.pixmaps = {}
        self.pngFiles = []
        self.svgFiles = []
        self.loadCache()
        self.pixmaps = {}
        #self.loadPixmaps()

    def loadCache(self, dir='', subDirectories = ['']):
        if not dir:
            dir = os.path.join(os.environ['CCP4I2_TOP'], 'qticons')
        try:
            with open(os.path.join(dir, "CachedPixmapPaths.json"), "r") as pixmapCacheFile:
                self.cachedFilenames=json.loads(pixmapCacheFile.read())
        except:
            self.buildCacheFromScratch(dir, subDirectories)

    def buildCacheFromScratch(self, dir='', subDirectories = []):
        if not dir:
            dir = os.path.join(os.environ['CCP4I2_TOP'], 'qticons')
            subDirectories = ['']
        self.cachedFilenames = {}
        for subd in subDirectories:
            pngPaths = glob.glob(os.path.join(dir, subd, '*.png'))
            for fileName in pngPaths:
                name, extension = os.path.splitext(os.path.split(fileName)[-1])
                self.cachedFilenames[name] = os.path.join(subd, name+extension)
            svgPaths = glob.glob(os.path.join(dir, subd, '*.svg'))
            for fileName in svgPaths:
                name, extension = os.path.splitext(os.path.split(fileName)[-1])
                self.cachedFilenames[name] = os.path.join(subd, name+extension)
        cacheFilename = os.path.join(dir, "CachedPixmapPaths.json")
        with open(cacheFilename,"w") as pixmapCacheFile:
            pixmapCacheFile.write(json.dumps(self.cachedFilenames))

    def loadPixmaps(self, dir='', width=32, height=0):
        if not dir:
            dir = os.path.join(os.environ['CCP4I2_TOP'], 'qticons')
            subDirectories = ['']
        self.pixmaps = {}
        if height<=0:
            height = width
        for name in self.cachedFilenames:
            filename = self.cachedFilenames[name]
            if filename.endswith("png"):
                self.pixmaps[name] = os.path.join(dir, filename)
            elif filename.endswith("svg"):
                self.pixmaps[name] = self.loadSvg(os.path.join(dir, filename))

    def loadSvg(self,fileName):
        svg = QtSvg.QSvgRenderer()
        svg.load(fileName)
        pixmap = QtGui.QPixmap(CPixmapManager.SIZE, CPixmapManager.SIZE)
        pixmap.fill(QtGui.QColor(0, 0, 0, 0))
        painter = QtGui.QPainter(pixmap)
        svg.render(painter)
        painter.end()
        return pixmap

    def getPixmap(self,name):
        dir = os.path.join(os.environ['CCP4I2_TOP'], 'qticons')
        #print 'getPixmap',name,self.pixmaps.has_key(name)
        if name in self.cachedFilenames:
            filename = self.cachedFilenames[name]
            if not name in self.pixmaps:
                if filename.endswith("png"):
                    self.pixmaps[name] = os.path.join(dir, filename)
                elif filename.endswith("svg"):
                    self.pixmaps[name] = self.loadSvg(os.path.join(dir, filename))
            return self.pixmaps[name]
        elif 'blank' in self.pixmaps:
            return self.pixmaps['blank']
        else:
            return None

class CIntValidator(QtGui.QIntValidator):
    '''Reimplementation of Qt validator'''

    def __init__(self, bottom=None, top=None, parent=None):
        QtGui.QIntValidator.__init__(self, bottom, top, parent)

    def validate(self, input, pos):
        if len(input) == 0:
            return (QtGui.QValidator.Acceptable, pos)
        return QtGui.QIntValidator.validate(self, input, pos)


class CDoubleValidator(QtGui.QDoubleValidator):
    '''Reimplementation of Qt validator'''

    def __init__(self, bottom=None, top=None, decimals=3, parent=None):
        if bottom is None:
            bottom = float('-inf')
        if top is None:
            top = float('inf')
        QtGui.QDoubleValidator.__init__(self, bottom, top, decimals, parent)
        self.setNotation(QtGui.QDoubleValidator.StandardNotation)

    def validate(self, input, pos):
        if len(input) == 0:
            return (QtGui.QValidator.Acceptable, pos)
        try:
            float(input)
        except:
            return QtGui.QDoubleValidator.validate(self, input, pos)
        return QtGui.QDoubleValidator.validate(self, input, pos)


class CClickableFrame(QtWidgets.QFrame):

    contextMenu = QtCore.Signal(int,int)

    def __init__(self, parent):
        QtWidgets.QFrame.__init__(self, parent)

    def mouseReleaseEvent(self,event):
        # Reimplementing QWidget.mouseReleaseEvent because QFrame does not have required signal
        if event.button() ==  QtCore.Qt.RightButton:
            self.contextMenu.emit(event.globalX(), event.globalY())
        QtWidgets.QFrame.mouseReleaseEvent(self, event)


class CLabel(QtWidgets.QLabel,CBaseWidget):

    acceptDropDataSignal = QtCore.Signal(object)
    dataChanged = QtCore.Signal()

    rightMousePress = QtCore.Signal('QMouseEvent')
    rightMouseRelease = QtCore.Signal('QMouseEvent')
#FIXME - I put this one in to satisfy connect statements using CLabel.editingFinished. But they are probably bogus and probably ought to be removed.
    editingFinished = QtCore.Signal()

    def __init__(self, parent=None, qualifiers={}, **kw):
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(kw)
        QtWidgets.QLabel.__init__(self, parent)
        CBaseWidget.__init__(self, dragType=qualis.get('dragType', None))
        if 'text' in qualis:
            if qualis['text'] is not None:
                self.setText(qualis['text'])
        self.setCharWidth(qualis.get('charWidth', None))
        if not qualis.get('editable', True):
            styleSheet = qualis.get('styleSheet', None)
            outlineFrame = True
        else:
            styleSheet = qualis.get('styleSheet', None)
            outlineFrame = qualis.get('outlineFrame', None)
        if styleSheet is not None:
            self.setStyleSheet(styleSheet)
        if outlineFrame is not None:
            self.setFrameStyle(QtWidgets.QFrame.StyledPanel | QtWidgets.QFrame.Plain)

    def mousePressEvent(self, event):
        if event.button() == QtCore.Qt.RightButton:
            self.rightMousePress.emit(event)
            event.accept()
        else:
            QtWidgets.QLabel.mousePressEvent(self, event)

    def mouseReleaseEvent(self,event):
        if event.button() == QtCore.Qt.RightButton:
            self.rightMouseRelease.emit(event)
        QtWidgets.QLabel.mouseReleaseEvent(self,event)

    def setValue(self, text=None):
        if text is None:
            self.setText('')
        else:
            self.setText(str(text))

    def getValue(self):
        return self.text().__str__()

class CItalicLabel(CLabel):
    def __init__(self, text=None, parent=None):
        CLabel.__init__(self, parent)
        if text is not None:
            self.setText(text)
        # use style sheet to set style for QLabel#italic
        self.setObjectName('italic')

class CBoldLabel(CLabel):
    def __init__(self, text=None, parent=None):
        CLabel.__init__(self, parent)
        if text is not None:
            self.setText(text)
        # use style sheet to set style for QLabel#bold
        self.setObjectName('bold')

class CPushButton(QtWidgets.QPushButton, CBaseWidget):

    acceptDropDataSignal = QtCore.Signal(object)
    dataChanged = QtCore.Signal()

    rightMousePress = QtCore.Signal('QMouseEvent')

    def __init__(self, parent=None, qualifiers={}, **kw):
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(kw)
        QtWidgets.QPushButton. __init__(self,parent)
        if 'text' in qualis and qualis['text'] is not None:
            self.setText(qualis['text'])
        tip = self.makeToolTip(qualis)
        if len(tip) > 0:
            self.setToolTip(tip)

    def mousePressEvent(self,event):
        if event.button() == QtCore.Qt.RightButton:
            self.rightMousePress.emit(event)
            event.accept()
        else:
            QtWidgets.QPushButton.mousePressEvent(self, event)

class CMultiAtomSelection(CCP4RefmacMultiAtomSelection.MultiAtomSelection, CBaseWidget):
    def __init__(self, parent=None, qualifiers={},**kw):
        qualis= {}
        qualis.update(qualifiers)
        qualis.update(kw)
        CCP4RefmacMultiAtomSelection.MultiAtomSelection.__init__(self, parent)
        CBaseWidget.__init__(self, dragType=qualis.get('dragType', None))
        self.setCharWidth(qualis.get('charWidth', None))
        self.customContextMenu = qualis.get('customContextMenu', False)
        tip = self.makeToolTip(qualis)
        if len(tip) > 0:
            self.setToolTip(tip)

class CSequenceTable(CCP4SequenceList.SequenceTable, CBaseWidget):

    acceptDropDataSignal = QtCore.Signal(object)
    dataChanged = QtCore.Signal()

    def __init__(self, parent=None, qualifiers={},**kw):
        qualis= {}
        qualis.update(qualifiers)
        qualis.update(kw)
        CCP4SequenceList.SequenceTable.__init__(self, parent)
        CBaseWidget.__init__(self, dragType=qualis.get('dragType', None))
        self.setCharWidth(qualis.get('charWidth', None))
        self.customContextMenu = qualis.get('customContextMenu', False)
        tip = self.makeToolTip(qualis)
        if len(tip) > 0:
            self.setToolTip(tip)

class CLineEdit(QtWidgets.QLineEdit, CBaseWidget):

    acceptDropDataSignal = QtCore.Signal(object)
    dataChanged = QtCore.Signal()

    rightMousePress = QtCore.Signal('QMouseEvent')
    contextMenuRequest = QtCore.Signal(int,int)
    editSignal = QtCore.Signal()

    def __init__(self, parent=None, qualifiers={},**kw):
        qualis= {}
        qualis.update(qualifiers)
        qualis.update(kw)
        QtWidgets.QLineEdit.__init__(self, parent)
        CBaseWidget.__init__(self, dragType=qualis.get('dragType', None))
        self.setCharWidth(qualis.get('charWidth', None))
        self.customContextMenu = qualis.get('customContextMenu', False)
        tip = self.makeToolTip(qualis)
        if len(tip) > 0:
            self.setToolTip(tip)

        self.editingFinished.connect(self.editSignal.emit)

    def mousePressEvent(self,event):
        if event.button() == QtCore.Qt.RightButton:
            self.rightMousePress.emit(event)
            event.accept()
        else:
            QtWidgets.QLineEdit.mousePressEvent(self,event)

    def contextMenuEvent(self,event):
        if self.customContextMenu:
            self.contextMenuRequest.emit(event.globalX(), event.globalY())
        else:
            QtWidgets.QLineEdit.contextMenuEvent(self, event)

    def setValue(self,value=None):
        #traceback.print_stack(limit=5)
        if value is None or value == 'None':
            self.setText('')
        else:
            self.setText(str(value))

    def getValue(self):
        return str(self.text())

class CCheckBox(QtWidgets.QCheckBox, CBaseWidget):

    acceptDropDataSignal = QtCore.Signal(object)
    dataChanged = QtCore.Signal()

    rightMousePress = QtCore.Signal('QMouseEvent')
    contextMenuRequest = QtCore.Signal(int,int)
    editSignal = QtCore.Signal()

    def __init__(self, parent=None, qualifiers={},**kw):
        QtWidgets.QCheckBox.__init__(self, parent)
        qualis={}
        qualis.update(qualifiers)
        qualis.update(kw)
        self.onValue = qualis.get('onValue', True)
        self.offValue = qualis.get('offValue', False)
        tip = self.makeToolTip(qualis)
        if len(tip) > 0:
            self.setToolTip(tip)

        self.clicked.connect(self.editSignal.emit)

    def mousePressEvent(self,event):
        if event.button() == QtCore.Qt.RightButton:
            self.rightMousePress.emit(event)
            event.accept()
        else:
            QtWidgets.QCheckBox.mousePressEvent(self, event)


    def contextMenuEvent(self, event):
        self.contextMenuRequest.emit(event.globalX(), event.globalY())

    def setValue(self, value=None):
        if value is None:
            pass
        elif value == self.onValue:
            self.setChecked(1)
        elif value == self.offValue:
            self.setChecked(0)

    def getValue(self):
        i = int(self.isChecked())
        if i:
            return self.onValue
        else:
            return self.offValue


class CCheckBoxUneditable(CCheckBox):

    def mousePressEvent(self, event):
        event.accept()

    def mouseReleaseEvent(self, event):
        event.accept()


class CComboBox(QtWidgets.QComboBox, CBaseWidget):

    acceptDropDataSignal = QtCore.Signal(object)
    dataChanged = QtCore.Signal()

    rightMousePress = QtCore.Signal('QMouseEvent')
    editSignal = QtCore.Signal()

    FONTHEIGHT = None

    def __init__(self, parent=None, qualifiers={}, **kw):
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(kw)
        self.onlyEnumerators = qualis.get('onlyEnumerators', False)
        QtWidgets.QComboBox.__init__(self, parent)
        CBaseWidget.__init__(self,dragType=qualis.get('dragType', None))
        values = qualis.get('enumerators', [])
        if len(values) > 0:
            self.populate(values, qualis.get('menuText', []))
        self.setEditable(not(self.onlyEnumerators))
        #print 'CComboBox editable', parent, qualis.get('onlyEnumerators',None),self.isEditable(),self.insertPolicy()
        default = qualis.get('default', None)
        if default is not None:
            indx = self.findData(default)
            if indx >= 0:
                self.setCurrentIndex(indx)
        self.currentIndexChanged.connect(self.handleMenuChange)
        tip = self.makeToolTip(qualis)
        if len(tip) > 0:
            self.setToolTip(tip)

        @QtCore.Slot()
        def currentChangedInter(int):
            self.editSignal.emit()
        self.currentIndexChanged[int].connect(currentChangedInter)

    def showPopup(self):
        '''This is attempt to fix the CDataFileView popup listview being drawn too short after several different
        CDataFileViews have been clicked on. Even after we have explicitly set the height of the listview it can appear with
        not all items displayed and the scroll arrows at top/bottom. This seems to be fixed by setting wrapping True but
        why I dont know. This same method is copied in CFinishedJobsCombo'''
        QtWidgets.QComboBox.showPopup(self)
        if sys.platform != 'darwin':
            return
        listViewList = self.findChildren(QtWidgets.QListView)
        if len(listViewList) > 0:
            listViewList[0].setWrapping(True)
            if CComboBox.FONTHEIGHT is None:
                CComboBox.FONTHEIGHT = listViewList[0].fontMetrics().height() + 2
            listViewList[0].setMinimumHeight(self.count() * CComboBox.FONTHEIGHT + 5)
            # print 'CComboBox.showPopup height', self.count(), CComboBox.FONTHEIGHT, listViewList[0].height()
            listViewList[0].setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
            listViewList[0].horizontalScrollBar().setEnabled(False)
            listViewList[0].horizontalScrollBar().setValue(0)

    def populate(self, menuItems=[], menuText=[]):
        if menuText is None:
            menuText = []
        self.menuText = menuText
        self.blockSignals(True)
        self.clear()
        for ii in range(len(menuItems)):
            if ii < len(menuText):
                self.addItem(menuText[ii], menuItems[ii])
            else:
                self.addItem(str(menuItems[ii]), menuItems[ii])
        self.blockSignals(False)

    def mousePressEvent(self, event):
        if event.button() == QtCore.Qt.RightButton:
            self.rightMousePress.emit(event)
            event.accept()
        else:
            QtWidgets.QComboBox.mousePressEvent(self, event)

    def keyReleaseEvent(self, event):
        #print 'CComboBox.keyPressEvent', event.key()
        self.dataChanged.emit()
        QtWidgets.QComboBox.keyReleaseEvent(self, event)

    def focusOutEvent(self, event):
        self.dataChanged.emit()
        QtWidgets.QComboBox.focusOutEvent(self, event)

    def contextMenuEvent(self, event):
        self.contextMenuRequest.emit(event.globalX(), event.globalY())

    def setValue(self, value=None):
        self.blockSignals(True)
        if value is None:
            self.setEditText('')
            self.setCurrentIndex(-1)
        else:
            ic = self.findText(str(value))
            if ic >= 0:
                text = value
            else:
                ic = self.findData(value)
                if ic >= 0 and ic < len(self.menuText):
                    text = self.menuText[ic]
                else:
                    text = value
            self.setCurrentIndex(ic)
            self.setEditText(str(text))
        self.blockSignals(False)

    def getValue(self):
        return str(self.currentText())

    def wheelEvent(self, event):
        # Intercept wheel event 'cos mostly gets input intended to scroll a frame 'CComboBox.wheelEvent'
        event.ignore()

    '''
    def mouseReleaseEvent(self,event):
      #print 'CCombo.mouseReleaseEvent',event.x(),event.y()
      QtWidgets.QComboBox.mouseReleaseEvent(self,event)
    '''

    @QtCore.Slot(int)
    def handleMenuChange(self, item):
        self.dataChanged.emit()


class CTextEdit(QtWidgets.QTextEdit, CBaseWidget):

    acceptDropDataSignal = QtCore.Signal(object)
    dataChanged = QtCore.Signal()

    editSignal = QtCore.Signal()

    def placeholderText(self):
        return self._placeholderText

    def setPlaceholderText(self,t):
        self._placeholderText = t

    def paintEvent(self,e):
        if not self.toPlainText() and self._placeholderText and not self.hasFocus():
            p = QtGui.QPainter(self.viewport())
            pen = p.pen()
            pen.setColor(QtGui.QColor(128,128,128))
            p.setPen(pen)
            p.drawText(self.geometry(), QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop, self._placeholderText)
        else:
            QtWidgets.QTextEdit.paintEvent(self,e);

    def __init__(self, parent=None, qualifiers={}, **kw):
        self._placeholderText = ""
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(kw)
        QtWidgets.QTextEdit.__init__(self, parent)
        CBaseWidget.__init__(self, dragType=qualis.get('dragType', None))
        if not qualis.get('editable', True):
            self.setReadOnly(True)
        tip = self.makeToolTip(qualis)
        if len(tip) > 0:
            self.setToolTip(tip)

        self.textChanged.connect(self.editSignal.emit)

    def setValue(self, value=None):
        self.blockSignals(True)
        if value is None:
            self.document().clear()
        else:
            self.document().setPlainText(value)
        self.blockSignals(False)

    def getValue(self):
        text = str(self.toPlainText())
        if len(text) == 0:
            return None
        else:
            return text

    def focusOutEvent(self, event):
        self.textChanged.emit()
        QtWidgets.QTextEdit.focusOutEvent(self, event)


class CRadioButtonGroup(QtWidgets.QButtonGroup, CBaseWidget):

    acceptDropDataSignal = QtCore.Signal(object)
    dataChanged = QtCore.Signal()

    editSignal = QtCore.Signal()

    def __init__(self, parent=None, qualifiers={}, **kw):
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(kw)
        QtWidgets.QButtonGroup.__init__(self, parent)
        self.setExclusive(True)
        CBaseWidget.__init__(self, dragType=qualis.get('dragType',None))
        self.toolTip = self.makeToolTip(qualis)
        self.lastButtonId = -1

        self.buttonReleased.connect(self.editSignal.emit)

    def addRadioButton(self, value=None, text=None):
        if text is not None and text is not NotImplemented:
            but = QtWidgets.QRadioButton(text,self.parent())
        else:
            but = QtWidgets.QRadioButton(self.parent())
        but.setToolTip(self.toolTip)
        but.setObjectName(str(value))
        self.lastButtonId = self.lastButtonId + 1
        QtWidgets.QButtonGroup.addButton(self,but,self.lastButtonId)
        return but

    def removeRadioButton(self, value=None):
        for but in self.buttons():
            if str(but.objectName()) == value:
                QtWidgets.QButtonGroup.removeButton(self, but)
                return

    def setValue(self, value=None):
        #print 'CRadioButtonGroup.setValue',value,self.buttons()
        value = str(value)
        for but in self.buttons():
            if str(but.objectName()) == value:
                self.blockSignals(True)
                but.setChecked(True)
                self.blockSignals(False)
                return

    def getValue(self):
        indx = int(self.checkedId())
        v = str(self.button(indx).objectName())
        return v

    def getButton(self,value):
        if isinstance(value,CCP4Data.CData): value = value.get()
        for but in self.buttons():
            if str(but.objectName()) == value:
                return but
        return None

    def setToolTip(self,tip):
        for but in self.buttons(): but.setToolTip(tip)


class CBooleanRadioButtonGroup(CRadioButtonGroup):

    def __init__(self,parent=None,qualifiers={},**kw):
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(kw)
        CRadioButtonGroup.__init__(self,parent=parent,qualifiers=qualis)
        menuText = qualis.get('menuText',[NotImplemented,NotImplemented])
        self.addRadioButton(True,menuText[0])
        self.addRadioButton(False,menuText[1])


class CFolder(QtWidgets.QFrame):

    folderToggled = QtCore.Signal()

    MARGIN = 4

    def __init__(self,parent=None,title='Folder',open=1,titleBar=1,toggle=[],toggleFunction=[],icon=None):
        QtWidgets.QFrame.__init__(self,parent)
        self.setObjectName(title)
        layout = QtWidgets.QVBoxLayout()
        layout.setSpacing(CFolder.MARGIN)
        layout.setContentsMargins(CFolder.MARGIN,CFolder.MARGIN,CFolder.MARGIN,CFolder.MARGIN)
        if titleBar:
            self.titleBar = CClickableFrame(self)
            self.titleBar.setFrameStyle(QtWidgets.QFrame.Box|QtWidgets.QFrame.Plain)
            layout.addWidget(self.titleBar)
            title_layout = QtWidgets.QHBoxLayout()
            title_layout.setSpacing(CFolder.MARGIN)
            title_layout.setContentsMargins(CFolder.MARGIN,CFolder.MARGIN,CFolder.MARGIN,CFolder.MARGIN)
            if icon is not None:
                pix = PIXMAPMANAGER().getPixmap(icon)
                if pix is not None:
                    ico = QtGui.QIcon(pix)
                    self.titleBar.layout().addWidget(ico,0,QtCore.Qt.AlignRight)
            self.titleIcon = QtWidgets.QLabel(self)
            title_layout.addWidget(self.titleIcon)
            self.titleLabel = QtWidgets.QLabel(title,self)
            title_layout.addWidget(self.titleLabel)
            title_layout.addStretch(5)
            self.titleBar.setLayout(title_layout)
        else:
            self.titleBar = None
            self.titleIcon = None
        self.contents = QtWidgets.QFrame(self)
        self.contentsLayout = None
        layout.addWidget(self.contents)
        #layout.addStretch(5)
        if not open:
            open = 0
        self.open = 1-open
        self.toggleFolder()
        self.setLayout(layout)
        if self.titleBar is not None:
            self.titleBar.mouseRelease.connect(self.toggleFolder)
        if len(toggle)>0:
            if len(toggle) == 1:
                self.parent().setToggle(self, toggle[0], 'open', [True])
            elif len(toggle) == 2:
                self.parent().setToggle(self, toggle[0], toggle[1], [True])
            else:
                self.parent().setToggle(self, toggle[0], toggle[1], toggle[2])
        if len(toggleFunction):
            self.parent().setToggle(self, toggleFunction[0], toggleFunction=toggleFunction[1], parameterList=toggleFunction[2])

    def title(self):
        return self.titleLabel.text().__str__()

    def setTitleColour(self,colour):
        if self.titleBar:
            self.titleBar.setStyleSheet("QFrame { background-color:"+colour +"; }")

    def setContentsLayout(self,layout):
        self.contentsLayout = layout
        self.contents.setLayout(layout)

    @QtCore.Slot()
    def toggleFolder(self):
        self.open = 1 - self.open
        if self.open:
            self.contents.show()
            self.resetIcon('toc-minus')
        else:
            self.contents.hide()
            self.resetIcon('toc-plus')
        self.folderToggled.emit()

    def resetIcon(self, feature):
        if not self.titleIcon:
            return
        pix = QtGui.QPixmap(PIXMAPMANAGER().getPixmap(feature))
        self.titleIcon.setPixmap(pix)

    def openFolder(self):
        if not self.open:
            self.toggleFolder()

    def closeFolder(self):
        if self.open:
            self.toggleFolder()

    def isOpen(self):
        return self.open


class CIconButton(QtWidgets.QToolButton,CBaseWidget):

    acceptDropDataSignal = QtCore.Signal(object)
    dataChanged = QtCore.Signal()

    rightMousePress = QtCore.Signal('QMouseEvent')
    rightMouseRelease = QtCore.Signal('QMouseEvent')
    leftMousePress = QtCore.Signal('QMouseEvent')
    leftMouseRelease = QtCore.Signal('QMouseEvent')

    '''Sub-class QToolButton for drag'n'drop'''

    def __init__(self, parent=None, icon=None, dragType=None):
        QtWidgets.QToolButton.__init__(self, parent)
        CBaseWidget.__init__(self, dragType=dragType)
        self.setPopupMode(QtWidgets.QToolButton.InstantPopup)
        # The menu indicator on Linux is making the icon small. The flollowing ought to remove menu indicator but does
        # not due to bug that is worked round by the reimplementation of paintEvent() below
        self.setStyleSheet("QToolButton::menu-indicator { image: none; }")
        self.setToolTip('Left mouse to drag, right mouse for menu')
        if icon is not None:
            self.setIcon(icon)
        self.dragStartPosition = None

    # From https://bugreports.qt-project.org/browse/QTBUG-2036
    def paintEvent(self, event):
        p = QtWidgets.QStylePainter(self)
        opt = QtWidgets.QStyleOptionToolButton()
        self.initStyleOption(opt)
        p.drawComplexControl(QtWidgets.QStyle.CC_ToolButton,opt)

    def sizeHint(self):
        return QtCore.QSize(ICONBUTTONSIZE, ICONBUTTONSIZE)

    def mousePressEvent(self, event):
        #print 'CIconButton.mousePressEvent', event.button(), self.dragType()
        if event.button() == QtCore.Qt.LeftButton:
            self.leftMousePress.emit(event)
            if self.dragType() is not None:
                self.dragStartPosition = event.pos()
                event.accept()
        elif event.button() == QtCore.Qt.RightButton:
            self.rightMousePress.emit(event)
            QtWidgets.QToolButton.mousePressEvent(self, event)
        return

    def mouseReleaseEvent(self, event):
        if event.button() == QtCore.Qt.RightButton:
            self.rightMouseRelease.emit(event)
        if event.button() == QtCore.Qt.LeftButton:
            self.leftMouseRelease.emit(event)
        QtWidgets.QToolButton.mouseReleaseEvent(self, event)

    def mouseMoveEvent(self,event):
        if self.dragType() is None:
            return
        if self.dragStartPosition is not None and \
                            (event.pos() - self.dragStartPosition).manhattanLength() < QTAPPLICATION().startDragDistance():
            return
        self.parent().startDrag()


class CViewWidget(QtWidgets.QFrame, CBaseWidget):

    acceptDropDataSignal = QtCore.Signal(object)
    dataChanged = QtCore.Signal()

    editSignal = QtCore.Signal()

    MODEL_CLASS = CCP4Data.CData
    MARGIN = 4
    STRETCH = -1

    def __init__(self, parent=None, qualifiers={}, **kw):
        QtWidgets.QFrame.__init__(self, parent)
        CBaseWidget.__init__(self)
        self.editable = qualifiers.get('editable', True)
        self._presetWidth = qualifiers.get('width', -1)
        self.model=None
        self.widgets = {}
        self.lastWarning = None
        self.isValid = None
        self.blockUpdateView = False
        self.validityMessage = None
        #self.clicked.connect(self.editSignal.emit)

    def enterEvent(self, event):
        if self.model is None or not self.editable:
            return
        # print 'enterEvent', self.model.objectName()
        try:
            self.parentTaskWidget().setMessage(self.validityMessage,self.model.objectPath())
        except:
            pass

    def leaveEvent(self, event):
        if self.model is None or not self.editable:
            return
        #print 'leaveEvent', self.model.objectName()
        try:
            self.parentTaskWidget().unsetMessage(self.model.objectPath())
        except:
            pass

    def setValue(self, value=None):
        if self.blockUpdateView:
            return
        if hasattr(self, 'widget'):
            self.widget.blockSignals(True)
            self.widget.setValue(value)
            self.widget.blockSignals(False)
        else:
            if value is None:
                return
            for item in list(value.keys()):
                w = self.widgets.get(item,None)
                if w is not None:
                    w.blockSignals(True)
                    w.setValue(value[item])
                    w.blockSignals(False)
        self.lastWarning = None

    def getValue(self):
        if  isinstance(self.model, CCP4Data.CBaseData):
            data = self.widget.getValue()
        else:
            data = {}
            for item in self.model.dataOrder():
                w = self.widgets.get(item, None)
                if w is not None:
                    data[item] = w.getValue()
        return data

    def getModel(self):
        return self.model

    def setModel(self, model):
        if isinstance(model, self.MODEL_CLASS) or model is None:
            if self.model is not None:
                self.connectUpdateViewFromModel(False)
            self.model = model
            if model is not None:
                self.connectUpdateViewFromModel(True)
            self.validate()

    def modelObjectPath(self):
        if self.model is None:
            return ''
        else:
            return self.model.objectPath()

    def reset(self):
        if not isinstance(self.model, CCP4Data.CBaseData):
            for w in list(self.widgets.keys()):
                self.widgets[w].reset()


    def connectUpdateViewFromModel(self, mode=True, propagate=False):
        if self.model is None:
            return
        self.blockUpdateView = not(mode)
        if mode:
            self.model.dataChanged.connect(self.updateViewFromModel)
        else:
            pass#self.model.dataChanged.disconnect(self.updateViewFromModel)
        if propagate:
            p = self.parent()
            while p is not None and not isinstance(p,(QtWidgets.QDialog,QtWidgets.QMainWindow)):
                if isinstance(p,CViewWidget):
                    try:
                        p.connectUpdateViewFromModel(mode)
                    except:
                        p = None
                    else:
                        p = p.parent()
                else:
                    p = p.parent()

    def unsetModel(self):
        if self.model is not None:
            # ARGHHH - how to avoid disconnect (besides not changing the model)?
            self.connectUpdateViewFromModel(False)
            self.model = None

    @QtCore.Slot()
    def updateViewFromModel(self):
        if self.model is None:
            return
        if self.blockUpdateView:
            return
        #print 'CViewWidget.updateViewFromModel',self.model.objectName(),self.model.get0()
        #  traceback.print_stack(limit=6)
        self.blockSignals(True)
        self.setValue(self.model.get0())
        self.blockSignals(False)
        self.lastWarning = None

    @QtCore.Slot()
    def updateModelFromView(self):
        if self.model is None:
            return
        #print 'CViewWidget.updateModelFromView',self.model.objectName(),self.getValue()
        self.connectUpdateViewFromModel(False)
        isValid = None
        try:
            self.model.set(self.getValue())
        except:
            isValid = False
        self.connectUpdateViewFromModel(True)
        #print 'CViewWidget.updateModelFromView',self.model.objectName(),self.model,self.getValue(),isValid
        self.validate(isValid=isValid)

    def updateModelFromText(self):
        pass

    def stretchFactor(self):
        return self.STRETCH

    def validate(self,isValid=None,excludeWidgets=[],report=None,reportMessage=True):
        if self.model is None:
            return False
        #print 'CViewWidget.validate', self.model.objectName(),isValid,report
        if not getattr(self,'editable', True):
            return True
        #print 'CViewWidget.validate', self.model.objectName(),type(self.model),isValid,self.model.get(),self.widgets
        if isValid is None:
            v = self.model.validity(self.model.get())
            update = False
            #print 'CViewWidget.validate maxSeverity', v.maxSeverity(), v.report(ifStack=False,mode=2)
            if v.maxSeverity() > SEVERITY_WARNING:
                if self.isValid is None or self.isValid:
                    update = True
                self.setProperty("isValid", False)
                self.setProperty("hasWarning", False)
                self.isValid = False
                self.validityMessage = v.report(ifStack=False, user=True, mode=2, minSeverity=SEVERITY_UNDEFINED_ERROR)
            elif v.maxSeverity()==SEVERITY_WARNING:
                # For warning message do not highlight widget but do set the warning text
                if self.isValid is None or not self.isValid:
                    update = True
                self.setProperty("isValid", True)
                self.setProperty("hasWarning", True)
                self.isValid = True
                self.validityMessage = v.report(ifStack=False, user=True, mode=2, minSeverity=SEVERITY_UNDEFINED_ERROR)
            else:
                if self.isValid is None or not self.isValid: update = True
                self.setProperty("isValid",True)
                self.setProperty("hasWarning",False)
                self.isValid = True
                self.validityMessage = None
        else:
            update = True
            self.setProperty("isValid",isValid)
            self.isValid = isValid
            self.validityMessage = report
        #print 'CViewWidget.validate', self.model.objectName(),self.isValid,self.validityMessage,reportMessage
        # It does not update correctly unless we do this --
        # See http://stackoverflow.com/questions/1595476/are-qts-stylesheets-really-handling-dynamic-properties
        if update:
            self.updateValidityIndicator()
            try:
                self.parentTaskWidget().updateMessage(report,self.model.objectPath())
            except:
                pass
        if isinstance(self.widgets,dict):
            for key,w in list(self.widgets.items()):
                if key not in excludeWidgets:
                    w.validate(isValid,reportMessage=False)
        else:
            for w in self.widgets:
                w.validate(isValid,reportMessage=False)
        return self.isValid

    def updateValidityIndicator(self):
        self.style().unpolish(self)
        self.style().polish(self)
        self.update()

    def setFocus1(self,reason,indx):
        #print 'CViewWidget.setFocus1',self,self.model.objectName(),reason,indx
        self.setFocus(reason)

    def setToolTip(self,tip):
        #print 'CViewWidget.setToolTip',self.model.objectName(),tip,self.toolTip()
        if getattr(self,'widget',None) is not None:
            try:
                self.widget.setToolTip(tip+'\n'+self.widget.toolTip())
            except:
                pass
        try:
            for w in list(self.widgets.values()):
                w.setToolTip(tip+'\n'+w.toolTip())
            for l in self.findChildren(QtWidgets.QLabel):
                l.setToolTip(tip)
        except:
            pass


class CSeparator(CViewWidget):
    '''Hold container widgets in auto-generated gui'''
    from core import CCP4Container

    MODEL_CLASS = CCP4Container.CContainer

    def __init__(self, parent=None, qualifiers={}, **kw):
        CViewWidget.__init__(self,parent,)
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(kw)
        layout = QtWidgets.QHBoxLayout()
        layout.setSpacing(CViewWidget.MARGIN)
        layout.setContentsMargins(CViewWidget.MARGIN,CViewWidget.MARGIN,CViewWidget.MARGIN,CViewWidget.MARGIN)
        self.setLayout(layout)
        self.label = CLabel(self)
        self.label.setStyleSheet("QFrame { background-color:"+CONTRASTCOLOUR +"; }")
        layout.addWidget(self.label)

    @QtCore.Slot()
    def updateViewFromModel(self):
        if self.model is not None:
            self.label.setValue(str(self.model.objectName()))


class CComplexLineWidget(CViewWidget):
  MARGIN = 4
  MODEL_CLASS = CCP4Data.CData
  TABLE_EDIT_MODE = 'window'

  def dragEnterEvent(self,event):
        if self.acceptMimeData(event.mimeData()) is not None:
            event.accept()
        else:
            event.ignore()


  def dragMoveEvent(self, event):
        if self.acceptMimeData(event.mimeData()) is not None:
            event.setDropAction(QtCore.Qt.CopyAction)
            event.accept()
        else:
            event.ignore()

  def dropEvent(self,event):
        dropData = self.acceptMimeData(event.mimeData())
        if dropData is not None:
            self.acceptDropDataSignal.emit(dropData)
            event.setDropAction(QtCore.Qt.CopyAction)
            event.accept()
        else:
            event.ignore()


  def __init__(self,parent=None,qualifiers={}):
    qualis = {}
    qualis.update(qualifiers)
    self.iconButton = None
    self.iconMenu = None
    self.widgets = {}
    CViewWidget.__init__(self,parent,qualis)
    #print 'CComplexLineWidget.__init__',self.MODEL_CLASS, qualis.get('stack',False)
    self._stacked = qualis.get('stacked',False)
    self._dragType = qualis.get('dragType',None)
    self._dropTypes = qualis.get('dropTypes',[])
    if self._dragType is None:
      ty = (self.MODEL_CLASS.__name__).split('.')[-1]
      if ty[0]=='C':
        self._dragType= ty[1:]
      else:
        self._dragType= ty
    self.iconName = qualis.get('iconName',None)


    if qualis.get('gridLayout',False):
      layout = QtWidgets.QGridLayout()
    elif  qualis.get('vboxLayout',False):
      layout = QtWidgets.QVBoxLayout()
    else:
      layout = QtWidgets.QHBoxLayout()
    layout.setSpacing(CComplexLineWidget.MARGIN)
    layout.setContentsMargins(CComplexLineWidget.MARGIN,CComplexLineWidget.MARGIN,CComplexLineWidget.MARGIN,CComplexLineWidget.MARGIN)
    self.setLayout(layout)


    iconButton = qualis.get('iconButton',None)
    if iconButton is None: iconButton = (parent is not None and (not isinstance(parent,CComplexLineWidget) or isinstance(parent,CListView)))
    if iconButton:
      self.createIconButton()
      '''
      if qualis.get('guiLabel',NotImplemented) is not NotImplemented:
        layout.addWidget(QtWidgets.QLabel(qualis['guiLabel'],self))
      '''
      self.setFrameShape(QtWidgets.QFrame.StyledPanel)

      self.iconButton.rightMousePress.connect(self.updateMenu)
      if self.editable:
        self.iconButton.acceptDropDataSignal.connect(self.acceptDropData)
        self.acceptDropDataSignal.connect(self.acceptDropData)


  def createIconButton(self):
    if self.iconName is None and self.dragType() is not None:
      self.iconName = self.dragType()
    #print 'CComplexLineWidget.createIconButton',self,type(self),self.dragType(),'iconName',self.iconName
    from qtgui import CCP4GuiUtils
    icon =CCP4GuiUtils.createIcon(name=self.iconName)
    self.iconButton = CIconButton(self,icon=icon,dragType= self._dragType)
    self.iconButton.setIcon(icon)
    layout = self.layout()
    if isinstance(layout,QtWidgets.QGridLayout):
      layout.addWidget(self.iconButton,0,0)
    else:
      layout.addWidget(self.iconButton)

    self.iconMenu = QtWidgets.QMenu()
    self.iconButton.setMenu(self.iconMenu)
    height = self.iconButton.height()
    self.iconButton.setIconSize(QtCore.QSize(height,height))
    #print 'createIconButton',self.iconButton.iconSize().height(),self.iconButton.height()

  def createListButtons(self,editor= False):
    buttonBox = QtWidgets.QFrame()
    if CListView.SIDE_BOX:
      buttonBox.setLayout(QtWidgets.QVBoxLayout())
    else:
      buttonBox.setLayout(QtWidgets.QHBoxLayout())
    buttonBox.layout().setSpacing(0)
    buttonBox.layout().setContentsMargins(0,0,0,0)
    #buttonBox.layout().addStretch(1)
    menu =  [['list_add_grey','Add item','append'],['list_delete_grey','Remove item','delete']]
    if editor: menu.append(['bullet_arrow_right','Open editor','editor'])
    for item in menu:
      icon = QtGui.QIcon(PIXMAPMANAGER().getPixmap(item[0]))
      button = QtWidgets.QToolButton(self)
      button.setIcon(icon)
      button.setMaximumHeight(ICONBUTTONSIZE)
      button.setMaximumWidth(ICONBUTTONSIZE)
      #button.setText(item[1])
      button.setToolTip(item[1])
      buttonBox.layout().addWidget(button)
      button.clicked.connect(functools.partial(self.handleButtonClick,item[2]))
    if not CListView.SIDE_BOX: buttonBox.layout().addStretch(1)
    if editor: buttonBox.layout().itemAt(2).widget().setCheckable(True)

    return buttonBox

  def handleButtonClick(self,item):
    pass

  def startDrag(self):
    mimeData = self.createMimeData()
    #print 'CComplexLineWidget.startDrag',mimeData,mimeData is None
    if mimeData is None: return
    drag = QtGui.QDrag(self)
    drag.setMimeData(mimeData)
    drag.setPixmap(self.iconButton.icon().pixmap(DRAGICONSIZE,DRAGICONSIZE))
    drag.targetChanged.connect(self.changedTarget)
    drag.exec_(QtCore.Qt.CopyAction)

  @QtCore.Slot('QWidget')
  def changedTarget(self,widget):
    #print 'CComplexLineWidget.changedTarget',widget
    pass

  def createMimeData(self):
    dragText = self.dragData()
    #print 'CComplexLineWidget.createMimeData',repr(self),self.dragType(),dragText
    if dragText is None: return None
    data = QtCore.QByteArray()
    data.append(dragText)
    mimeData = QtCore.QMimeData()
    # With mime type as text the data can be dropped on desktop
    # but the type of the data is lost
    #mimeData.setData('text/plain',data)
    mimeData.setData(self.dragType(),data)
    return mimeData

  def isWidgetOpen(self):
    return True

  def hideWidgetLabel(self):
    if self.isWidgetOpen():
      return 'Hide'
    else:
      return 'Show'

  def handleHideWidget(self):
    pass

  def getActionDef(self,name):

    def e(): return self.dragType() is not None
    if name == 'copy':
      return dict (
        text = self.tr("Copy"),
        tip = self.tr('Copy data'),
        enabled = e,
        slot = self.copy
      )
    def e(): return (self.editable and self.clipboardLoaded())
    if name == 'paste':
      return dict (
        text = self.tr("Paste"),
        tip = self.tr('Paste data'),
        enabled = e,
        slot = self.paste
      )
    def e(): return (self.model.qualifiers('helpFile') is not NotImplemented )
    if name == 'help':
      return dict (
        text = self.tr("Help"),
        tip = self.tr('Help'),
        slot = self.helpMethod,
        enabled = e,
      )
    if name == 'clear':
      return dict (
        text = self.tr("Clear"),
        tip = self.tr('Clear all data'),
        slot = self.clear,
        enabled = self.editable
      )
    def e(): return (self.model is not None and hasattr(self.model,'exists') and self.model.exists())
    if name[0:4] == 'view':
      if name == 'view':
        text = 'View'
      elif name == 'view_text':
        text = 'View as text'
      else:
        text = "View in "+name[5:]
      return dict (
        text = self.tr(text),
        tip = self.tr('View data'),
        slot = functools.partial(self.openViewer,name),
        enabled = e
      )
    if name == 'handleStack':
      return dict (
        text = self.tr("Alternate view"),
        tip = self.tr('View data'),
        slot = self.parent().handleStack,
        checkable = True,
        checked = self.parent().isAltViewOpen
      )
    if name == 'hide':
      return dict (
        text = self.hideWidgetLabel,
        tip = self.tr('Open/close the widget'),
        slot = self.handleHideWidget
      )
    if name == 'show_list':
      return dict (
        text = self.tr("Show list"),
        tip = self.tr('Enter multiple data items'),
        slot = self.setListVisible,
        enabled = self.showListEnabled,
        checkable = True,
        checked = self.isListVisible
      )

  def getMenuDef(self):
    if self.editable:
      menu = ['clear','copy','paste','help']
    else:
      menu = ['copy','help']
    if self._stacked: menu.insert(0,'handleStack')
    return menu

  @QtCore.Slot('QMouseEvent')
  def updateMenu(self,event=None):
    #print 'CDataFileView.updateMenu',self,type(self),self.getMenuDef()
    from qtgui import CCP4GuiUtils
    self.iconMenu.clear()
    CCP4GuiUtils.populateMenu(self,menuWidget=self.iconMenu,definition=self.getMenuDef(),getActionDef=self.getActionDef)
    if sys.platform == "darwin":
        proxyStyle = MyProxyStyle()
        self.iconMenu.setStyle(proxyStyle)
    self.iconMenu.popup(event.globalPos())


  def copy(self):
    #print 'CComplexLineWidget.copy',self.dragType()
    if self.dragType() is None: return
    mimeData = self.createMimeData()
    if mimeData is not None:
      cb = QTAPPLICATION().clipboard()
      cb.setMimeData(mimeData)


  def paste(self):
    #print 'CComplexLineWidget.paste'
    if self.dragType() is None: return 0
    mimeData = QTAPPLICATION().clipboard().mimeData()
    textData = self.acceptMimeData(mimeData)
    #print('CComplexLineWidget.paste textData',textData)
    if textData is not None:
      self.acceptDropData(textData)

  def clipboardLoaded(self):
    if self.dragType() is None: return 0
    try:
       mimeData = QTAPPLICATION().clipboard().mimeData()
       textData = self.acceptMimeData(mimeData)
    except:
       return 0
    if textData is not None:
      return 1
    else:
      return 0

  @QtCore.Slot()
  def helpMethod(self):
    helpFile = self.model.qualifiers('helpFile')
    if helpFile is not NotImplemented:
      WEBBROWSER().loadWebPage(helpFileName=helpFile)
    else:
      print('No help file defined for ',self.model.__class__)

  def dragData(self):
    #print 'CComplexLineWidget.dragData',self.model.xmlText(pretty_print=False)
    return self.model.xmlText(pretty_print=False)

  @QtCore.Slot(object)
  def acceptDropData(self,textData):
    #print 'CComplexLineWidget.acceptDropData',textData
    from lxml import etree
    tree = etree.fromstring(textData)
    tree.tag = self.dragType()
    #print 'CComplexLineWidget.acceptDropData',textData,tree.tag
    self.connectUpdateViewFromModel(False)
    self.model.unSet()
    self.model.setEtree(tree)
    #print 'CComplexLineWidget.acceptDropData model',self.model
    self.connectUpdateViewFromModel(True)
    self.updateViewFromModel()
    #print 'CComplexLineWidget.acceptDropData validity',self.model.validity(self.model.get()).report()
    self.validate()

  def openViewer(self,**kw):
    print('CComplexLineWidget.openViewer should by reimplemented in subclass')

  def clear(self):
    self.model.blockSignals(True)
    self.model.unSet()
    self.model.blockSignals(False)
    #print 'CComplexLineWidget.clear',self.model.objectName(),self.model
    self.updateViewFromModel()
    self.validate()


  def setModel(self,model):
    #print 'CComplexLineWidget.setModel',model,self.widgets.keys()
    #if model is not None: print model.objectName()
    if model is None or isinstance(model,self.MODEL_CLASS):
      CViewWidget.setModel(self,model)

      if model is not None:
        model.dataChanged.connect(self.validate)
        toolTip = model.qualifiers('toolTip')
        #SJM 07/07/2017 Commented this out because it triggers https://fg.oisin.rc-harwell.ac.uk/tracker/?func=detail&group_id=92&aid=1266&atid=118
        """
        if toolTip is not NotImplemented and toolTip is not None  and self.iconButton is not None:
          self.iconButton.setToolTip(toolTip+'\n'+self.iconButton.toolTip())
        """
        #print 'CComplexLineWidget.setModel', self.widgets.items()
        for key,w in list(self.widgets.items()):
          #print 'CComplexLineWidget.setModel',key,w
          if isinstance(w,CViewWidget): w.setModel(model.get(key))
      else:
        for key,w in list(self.widgets.items()):
           if isinstance(w,CViewWidget): w.setModel(None)

  def getModel(self):
    return self.model

  def setToolTip(self,toolTip):
    try:
      self.iconButton.setToolTip(toolTip+'\n'+self.iconButton.toolTip())
    except:
      pass
    CViewWidget.setToolTip(self,toolTip)

  def taskProjectId(self):
    # Widget may be in a task interface (in which case it has a parent CTaskWidget
    p =  self.parentTaskWidget()
    if p is not None:
      return p.projectId()
    else:
      #or in a report, in which case it has a parent CTaskFrame
      p=self.parentTaskFrame()
      if p is not None:
        return p.projectId()
      return None


class CDataFileView(CComplexLineWidget):

  editSignal = QtCore.Signal()

  STRETCH = 5
  MODEL_CLASS = CCP4File.CDataFile
  ERROR_CODES = { 100 : { 'description' : 'Error copying file' },
                 101 : { 'description' : 'Unable to access ProjectManager to set project combo box' }
                 }
  MAXCHARS = 9999

  def __init__(self,parent=None,model=None,qualifiers={},**kw):

    qualis = qualifiers

    CComplexLineWidget.__init__(self,parent,qualifiers=qualis)
    self.setAcceptDrops(True)
    self.ifJobCombo = qualis.get('jobCombo',True)
    self.ifInfo = qualis.get('ifInfo',False)
    self.autoInfoOnFileImport=qualis.get('autoInfoOnFileImport',True)
    toolTip = qualis.get('toolTip',NotImplemented)
    self.browser = None
    self.databaseBrowser= None
    self.jobChooser = None
    self.model = None
    self.role = qualifiers.get('role',None)

    if qualis.get('vboxLayout',False):
      layout = QtWidgets.QHBoxLayout()
      iconWidgetItem = self.layout().takeAt(0)
      if iconWidgetItem is not None:
        layout.addWidget(iconWidgetItem.widget())
        self.layout().addLayout(layout)
    else:
      layout = self.layout()

    #print 'CDataFileView.init qualifiers',qualifiers
    if qualifiers.get('stackButton',None) is not None: layout.addWidget(qualifiers['stackButton'])

    if qualis.get('guiLabel',None) not in (None,NotImplemented):
      label = QtWidgets.QLabel(qualis['guiLabel'],self)
      label.setMinimumWidth(100)
      layout.addWidget(label)

    if self.editable:
      if self.ifJobCombo:
        self.nComboFiles = 0
        self.jobCombo = CComboBox(self,qualis,onlyEnumerators=True)
        self.jobCombo.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.jobCombo.setToolTip('Select a file output by previous job from menu')
        #self.jobCombo.setStyleSheet("QComboBox { combobox-popup: 0; }");
        #self.jobCombo.setMaxVisibleItems(20)
        self.jobCombo.setEditable(False)
        # http://stackoverflow.com/questions/11252299/pyqt-qcombobox-setting-number-of-visible-items-in-dropdown
        #self.jobCombo.setStyleSheet("QComboBox { combobox-popup: 0; }")
        self.jobCombo.setMaxVisibleItems(10)
        #if toolTip is not NotImplemented: self.jobCombo.setToolTip(toolTip)
        layout.addWidget(self.jobCombo,5)
        self.setFocusProxy(self.jobCombo)
      else:
        self.fileLineEdit = CLineEdit(self,qualifiers=qualis)
        self.fileLineEdit.setToolTip('Name of file or full path for file not associated with project')
        if toolTip is not NotImplemented: self.fileLineEdit.setToolTip(toolTip)
        layout.addWidget(self.fileLineEdit,5)
        self.setFocusProxy(self.fileLineEdit)
    else:
      if self.ifJobCombo:
        self.jobLabel = CLabel(self,qualifiers=qualis,editable=False)
        self.jobLabel.setToolTip('Using the best file available at finish of this job')
        layout.addWidget(self.jobLabel,5)
        PROJECTSMANAGER().db().fileUpdated.connect(self.handleDbFileUpdated)
      else:
        self.fileLineEdit = CLabel(self,qualifiers=qualis,editable=False)
        layout.addWidget(self.fileLineEdit,5)


    if self.editable:
      #print  'browseDb',qualis.get('browseDb',None)
      if qualis.get('browseDb',None) is not None:
        browseDb = qualis['browseDb']
      elif model is None:
        browseDb=True
      else:
        browseDb = model.qualifiers('fromPreviousJob')

      if qualis.get('browseFiles', None) is not None:
          browseFiles = qualis['browseFiles']
      else:
          browseFiles = True

      ifDownload = PREFERENCES().NATIVEFILEBROWSER and model is not None and model.qualifiers('downloadModes') is not NotImplemented and browseFiles
      ifDemoData = PREFERENCES().NATIVEFILEBROWSER

      if ifDemoData:
          browseFilesText = 'Browse files, press and hold to browse for demo data'
      else:
          browseFilesText = 'Browse files'

      iconDefn =  [ [ ifDownload , 'download', 'Download file' , self.downloadGui ],
                    [ browseFiles , 'file_manager', browseFilesText, self.openBrowser],
                    [ browseDb, 'database' ,'Browse database', self.openDatabaseBrowser ] ]


      for draw,iconName,toolTip,action in iconDefn:
        if draw:
          icon = QtGui.QIcon(PIXMAPMANAGER().getPixmap(iconName))
          browserButton = QtWidgets.QToolButton(self)
          browserButton.setIcon(icon)
          browserButton.setMaximumHeight(ICONBUTTONSIZE)
          browserButton.setMaximumWidth(ICONBUTTONSIZE)
          browserButton.setToolTip(toolTip)
          layout.addWidget(browserButton)
          if iconName == "file_manager" and ifDemoData:
              menu = QtWidgets.QMenu();
              testAction = QtGui.QAction("Browse for demo data",browserButton);
              @QtCore.Slot(bool)
              def handleOpenTrigger():
                  self.openBrowser(os.path.join(os.environ["CCP4"],"share/ccp4i2/demo_data"))
              testAction.triggered.connect(handleOpenTrigger)
              menu.addAction(testAction);
              browserButton.setMenu(menu);
              browserButton.setPopupMode(QtWidgets.QToolButton.DelayedPopup);
          if iconName == "file_manager":
              @QtCore.Slot(object)
              def handleOpenTrigger2(act):
                  act()
              browserButton.clicked.connect(functools.partial(handleOpenTrigger2,action))
          else:
              browserButton.clicked.connect(action)

      if not self.ifJobCombo:
        self.fileLineEdit.acceptDropDataSignal.connect(self.acceptDropData)
        self.fileLineEdit.editingFinished.connect(self.updateModelFromView)
        self.fileLineEdit.textEdited.connect(self.updateModelFromView1)

    self.setModel(model)
    self.validate()

    if self.editable and self.ifJobCombo:
      PROJECTSMANAGER().db().jobFinished.connect(self.handleJobFinished)
      self.jobCombo.acceptDropDataSignal.connect(self.acceptDropData)
      self.jobCombo.currentIndexChanged[int].connect(self.handleJobComboChange)

    PROJECTSMANAGER().db().fileUpdated.connect(self.handleFileUpdated)
    if self.parentTaskWidget() is not None and hasattr(self.parentTaskWidget(),"followFromJobUpdated"):
      self.parentTaskWidget().followFromJobUpdated.connect(self.handleFollowFrom)

    #self.fileChanged.connect(self.editSignal.emit)

  @QtCore.Slot()
  def downloadGui(self):
    '''Display a window to download from web server'''
    from qtgui import CCP4FileBrowser
    projectId = self.parentTaskWidget().projectId()
    downloadDir = os.path.join(PROJECTSMANAGER().db().getProjectDirectory(projectId),'CCP4_DOWNLOAD')
    self.dowloadDialog = QtWidgets.QDialog(self)
    self.dowloadDialog.setLayout(QtWidgets.QVBoxLayout())
    self.dowloadDialog.setWindowTitle('Download file from database')
    self.dowloadDialog.layout().addWidget(QtWidgets.QLabel('Downloaded file will be saved to '+downloadDir+' directory',self))
    dialog1 = CCP4FileBrowser.CFileDialog1(self.dowloadDialog,False,True)
    self.dowloadDialog.layout().addWidget(dialog1)
    dialog1.drawDownload(self.model.qualifiers('downloadModes'),projectId=projectId)
    self.dowloadDialog.show()
    dialog1.selectDownloadedFile.connect(self.handleDownloadFinished)

  @QtCore.Slot(str,str)
  def handleDownloadFinished(self,fileName,downloadInfo):
    #print 'CDataFileView.handleDownloadFinished',fileName,downloadInfo
    self.dowloadDialog.close()
    self.dowloadDialog.deleteLater()
    del self.dowloadDialog
    self.handleBrowserOpenFile(fileName,downloadInfo)


  def setToolTip(self,toolTip):
    try:
      self.iconButton.setToolTip(toolTip+'\n'+self.iconButton.toolTip())
    except:
      pass
    if getattr(self,'jobCombo',None) is not None:
      self.jobCombo.setToolTip(toolTip+'\n'+self.jobCombo.toolTip())
    elif getattr(self,'fileLineEdit',None) is not None:
      self.fileLineEdit.setToolTip(toolTip+'\n'+self.fileLineEdit.toolTip())
    elif getattr(self,'jobLabel',None) is not None:
      self.jobLabel.setToolTip(toolTip+'\n'+self.jobLabel.toolTip())
    for l in self.findChildren(QtWidgets.QLabel):
        l.setToolTip(toolTip)

  @QtCore.Slot(str,str)
  def handleFollowFrom(self,contextJobId,projectId):
    '''Handle change to the followFrom job'''
    print('CDataFileView.handleFollowFrom',self.model.objectName(),contextJobId,self.model.qualifiers('fromPreviousJob'))
    from dbapi import CCP4DbApi
    if self.model.qualifiers('fromPreviousJob'):
      if contextJobId is None:
        self.model.unSet()
      else:
        mimeType = self.model.qualifiers('mimeTypeName')
        if CCP4DbApi.FILETYPES_TEXT.count(mimeType):
          fileType = CCP4DbApi.FILETYPES_TEXT.index(mimeType)
        else:
          fileType = CCP4DbApi.FILETYPES_TEXT[0]
        subType =  self.model.qualifiers('requiredSubType')
        contentFlag =  self.model.requiredContent()
        #print 'CDataFileView.handleFollowFrom',self.model.objectName(),fileType,subType,type(subType),contentFlag,type(contentFlag)

        fileIdList = PROJECTSMANAGER().db().getFileByJobContext(contextJobId=contextJobId,fileType=fileType,subType=subType,contentFlag=contentFlag,projectId=projectId)
        #print 'CDataFileView.handleFollowFrom',self.model.objectName(),mimeType,fileType,fileIdList
        if len(fileIdList)>0:
          fileInfo = PROJECTSMANAGER().db().getFileInfo(fileId=fileIdList[0],mode=['jobid','filename','relpath','projectid','annotation','filecontent','filesubtype'])
          self.model.set( { 'baseName' : fileInfo['filename'],
                          'relPath' :  fileInfo['relpath'],
                          'project' : fileInfo['projectid'],
                          'annotation' : fileInfo['annotation'],
                          'dbFileId' :fileIdList[0],
                          'contentFlag' :fileInfo['filecontent'],
                          'subType' :fileInfo['filesubtype'] } )
        else:
          self.model.unSet()
    self.updateViewFromModel()
    self.validate()

  def reset(self):
    if self.editable:
      if self.ifJobCombo:
        self.loadJobCombo()
        #self.jobCombo.setCurrentIndex(0)
      else:
        self.fileLineEdit.clear()
    self.validate()



  def getMenuDef(self):
    '''Return text list of actions required on menu'''
    if self.editable:
      menu = ['clear','view','sep','copy','paste','help']
      #menu = ['clear','view','sep','copy','paste','help']
    else:
      if self.role is not None and self.role == 'output':
        menu = ['view','sep','copy','editLabel','export','help']
      else:
        menu = ['view','sep','copy','export','help']
    if self._stacked: menu.insert(0,'handleStack')
    if self.ifInfo: menu.insert(menu.index('sep'),'annotate')
    return menu

  def getActionDef(self,name):
    '''Return definition of action to appear on menu'''
    if name == 'annotate':
      return dict (
        text = self.tr("Annotate "),
        tip = self.tr('Add information on file, particularly its provenance'),
        slot = self.openInfo,
        enabled = self.isImportedFile
      )
    elif name == 'export':
      return dict (
        text = self.tr("Export"),
        tip = self.tr('Copy and rename file'),
        slot = self.handleExport
      )
    elif name == 'editLabel':
      return dict (
        text = self.tr("Edit label"),
        tip = self.tr('Edit file label'),
        slot = self.handleEditLabel
      )
    else:
      return CComplexLineWidget.getActionDef(self,name)


  def filterText(self):
    # make the filters text for QFileDialog
    if self.model.qualifiers('mimeTypeDescription') is not NotImplemented:
      text = self.model.qualifiers('mimeTypeDescription') +' ('
    else:
      text = '('
    if self.model.qualifiers('fileExtensions') is not NotImplemented and len(self.model.qualifiers('fileExtensions'))>0:
      for ext in self.model.qualifiers('fileExtensions'):
        text = text + '*.'+ext+' '
      text = text[0:-1]+')'
    else:
      text = text + '*)'
    #print 'filterText',text
    return [ text ]


  def clear(self):
    '''Clear any selection'''
    if self.editable:
      if self.ifJobCombo:
        self.jobCombo.setCurrentIndex(0)
      else:
        self.fileLineEdit.setValue(None)
    self.updateModelFromView()
    #print 'CDataFileView.clear',self.model

  def isImportedFile(self):
    #print 'CDataFile.isImportedFile',self.model.objectName(),self.model.dbFileId.__str__()
    if not self.model.dbFileId.isSet():
      return True
    else:
      return False

  def handleExport(self):
    '''Handle menu option to export file'''
    from qtgui import CCP4FileBrowser
    fileType = PROJECTSMANAGER().db().getFileInfo(fileId=self.model.dbFileId,mode='fileType')
    filters = MIMETYPESHANDLER().getMimeTypeInfo(fileType,'filter')
    defaultSuffix =  MIMETYPESHANDLER().getMimeTypeInfo(fileType,'fileExtensions')[0]
    #print 'CDataFileView.handleExport',fileType,filters,defaultSuffix
    fileBrowser = CCP4FileBrowser.CFileDialog(parent=self,
                                      title='Export '+self.model.baseName.__str__(),
                                      filters = [filters],
                                      defaultSuffix = defaultSuffix,
                                      fileMode = QtWidgets.QFileDialog.AnyFile)
    fileBrowser.selectFile.connect(functools.partial(self.exportData,self.model.__str__(),str(self.model.dbFileId)))
    fileBrowser.show()

  @QtCore.Slot(str,str,str)
  def exportData(self,myFileName,fileId,exportFileName):
    # This is copy of code in CCP4ProjectWidget
    #print 'CDataFileView.exportData',myFileName,exportFileName,fileId
    if os.path.splitext(exportFileName) !=  os.path.splitext(myFileName):
      exportFileName = os.path.splitext(exportFileName)[0] +  os.path.splitext(myFileName)[1]
    try:
      shutil.copyfile(myFileName,exportFileName)
    except:
      e = CException(self.__class__,100,'From: '+str(myFileName)+' to: '+str(exportFileName),name=self.modelObjectPath())
      e.warningMessage('Copying file',parent=self)
    else:
      PROJECTSMANAGER().db().createExportFile(fileId=fileId,exportFilename=exportFileName)
      fileInfo = PROJECTSMANAGER().db().getFileInfo(fileId=fileId,mode=['jobid','projectname'])
      from dbapi import CCP4DbUtils
      CCP4DbUtils.makeJobBackup(jobId=fileInfo['jobid'],projectName=fileInfo['projectname'])

  def handleEditLabel(self):
    from qtgui import CCP4Widgets
    d = CCP4Widgets.CEditFileLabel(parent=self,fileId=self.model.dbFileId.__str__(),fileLabel=self.jobLabel.getValue())

  @QtCore.Slot(dict)
  def handleDbFileUpdated(self,args):
    '''handle db signal that file attributes have changed - possible change to annotation'''
    if args['key'] != 'annotation': return
    if self.model is None or args['fileId'] != self.model.dbFileId.__str__(): return
    #print 'CDataFileView.handleDbFileUpdated updating'
    self.updateJobLabel()


  def setModel(self,model):
    #if model is not None: print 'CDataFileView.setModel',model.objectName(),model,isinstance(model,self.MODEL_CLASS)
    if isinstance(model,self.MODEL_CLASS) or model is None:
      if self.model is not None:
          self.connectUpdateViewFromModel(False)
      self.model = model
      if self.editable and self.ifJobCombo: self.loadJobCombo()
      self.updateViewFromModel()
      if model is not None: self.connectUpdateViewFromModel(True)


  @QtCore.Slot()
  def updateViewFromModel(self):
    if self.model is None:
      if self.ifJobCombo:
        pass
      else:
        self.fileLineEdit.setValue('')
      return

    #print 'CDataFileView.updateViewFromModel',self.model.objectName(),self.model
    self.blockSignals(True)
    if self.ifJobCombo:
      if self.editable:
        self.updateJobCombo()
      else:
        self.updateJobLabel()
    else:
      if self.model.project.isSet():
        projectName = str(self.model.project)
      else:
        projectName = None
      #print 'updateViewFromModel',projectName,type(projectName),'self.taskProject',self.taskProjectId()

      #print 'CDataFileView.updateViewFromModel',self.model.objectName(),type(self.model),self.model.annotaton.isSet(),self.model.annotaton.__str__()
      if self.model.annotation.isSet():
        fileName = self.model.annotation.__str__()
      else:
        fileName = self.model.__str__()
      self.fileLineEdit.setValue(fileName)
    self.blockSignals(False)

  def updateJobLabel(self):
    #print 'updateJobLabel',self.model.objectName(),self.model.isSet(), self.model,'*', self.role

    if not self.model.isSet():
      self.jobLabel.setValue('')
      return

    if self.role is None:
      container = self.model.parentContainer()
      if container is not None and container.objectName() == 'outputData':
        self.role = 'output'
      else:
       self.role = 'input'


    jobId = None
    jobInfo = {}
    if self.model.dbFileId.isSet():
      try:
        fileInfo = PROJECTSMANAGER().db().getFileInfo(fileId=str(self.model.dbFileId),mode=['taskname','jobnumber','annotation','filename'])
      except CException as e:
        print('updateJobLabel', e.report())
        fileInfo = {}
      #print 'updateJobLabel',self.model.objectName(),fileInfo
      if fileInfo.get('annotation',None) is not None:
        if self.role == 'output':
          self.jobLabel.setValue(str(fileInfo['annotation']))
        else:
          self.jobLabel.setValue(fileInfo.get('jobnumber','')+' '+str(fileInfo['annotation']))
      elif fileInfo.get('taskname',None) is not None :
        if self.role == 'output':
          self.jobLabel.setValue(self.model.qualifiers('guiLabel'))
        else:
          self.jobLabel.setValue(fileInfo.get('jobnumber','')+' '+TASKMANAGER().getTitle(fileInfo['taskname']))
      else:
        try:
          fileInfo = PROJECTSMANAGER().db().getFileInfo(fileId=str(self.model.dbFileId), mode=['relpath', 'filename'])
        except CException as e:
          print('updateJobLabel', e.report())
          fileInfo = {}
        if fileInfo.get('relpath',None) is None: fileInfo['relpath']=''
        if fileInfo.get('filename',None) is None: fileInfo['filename']=''
        self.jobLabel.setValue(os.path.join(fileInfo['relpath'], fileInfo.get('filename','')))
    elif self.model.annotation.isSet():
      self.jobLabel.setValue(str(self.model.annotation))
    elif not self.model.isSet():
      self.jobLabel.setValue('.. is not used')
    else:
      self.jobLabel.setValue(str(self.model))

  def updateJobCombo(self):
    '''Called from updateViewFromModel() to update the jobCombo when the model is changed'''
    if self.model is None: return
    if getattr(self,'jobCombo',None) is None: return
    #print 'updateJobCombo',self.model.objectName(),self.model,self.jobCombo.currentIndex()

    self.jobCombo.blockSignals(True)
    if not self.model.isSet():
      #print 'CDataFileView.updateJobCombo model not set',self.model.objectName()
      self.jobCombo.setCurrentIndex(0)
      self.validate()
      self.jobCombo.blockSignals(False)
      return

    fileId = self.model.dbFileId.get()
    #print 'CDataFileView.updateJobCombo',self.model.objectName(),fileId
    if fileId is not None:
      try:
        if CCP4Data.varToUUID(self.jobCombo.itemData(self.jobCombo.currentIndex())) == fileId:
          self.jobCombo.blockSignals(False)
          self.validate()
          return
        #print 'skip update in updateJobCombo'
      except:
        pass
      indx = self.jobCombo.findData(fileId)

      #print 'CDataFileView.updateJobCombo getDbFileId',self.model.objectName(),self.model,fileId,indx
      if indx>=0:
        self.jobCombo.setCurrentIndex(indx)
      else:
        # This file not yet in the combo list
        try:
          fileInfo = PROJECTSMANAGER().db().getFileInfo(fileId=fileId,mode=['taskname','jobnumber','annotation','relpath','filename','projectid'])
        except:
          print('Failed retrieving file data from Db data',fileId,self.model.objectPath(),self.model)
          self.jobCombo.setCurrentIndex(0)
        else:
          #print 'CDataFileView.updateJobCombo fileInfo',fileInfo
          if fileInfo['projectid'] != self.taskProjectId():
            title = PROJECTSMANAGER().db().getProjectInfo(projectId=fileInfo['projectid'],mode='projectname')+' '
          else:
            title = ''
          if fileInfo['annotation'] is not None:
            title = title + str(fileInfo['jobnumber'])+' '+fileInfo['annotation']
          elif fileInfo['taskname'] is not None :
            title = title + str(fileInfo['jobnumber'])+' '+TASKMANAGER().getTitle(fileInfo['taskname'])
          else:
            if fileInfo['relpath'] is None: fileInfo['relpath']=''
            title = title + os.path.join(fileInfo['relpath'],fileInfo['filename'])
          self.jobCombo.insertItem(1,title,fileId)
          self.nComboFiles = self.nComboFiles + 1
          self.jobCombo.setCurrentIndex(1)
    else:
      if self.model.annotation.isSet():
        text = str(self.model.annotation)
      else:
        text = str(self.model)
      shortText = text[-CDataFileView.MAXCHARS:]
      # Beware menu text may cwbwgin with include the job number
      indx = self.jobCombo.findText(shortText,QtCore.Qt.MatchEndsWith)
      if indx<0: indx = self.jobCombo.findData(text)
      #print 'CDataFileView.updateJobCombo indx',self.model.objectName(),text,indx
      if indx>=0:
        self.jobCombo.setCurrentIndex(indx)
        self.jobCombo.blockSignals(False)
        self.validate()
        return
      else:
        self.nComboFiles = self.nComboFiles + 1
        self.jobCombo.insertItem(1,shortText)
        self.jobCombo.setItemData(1,str(self.model),QtCore.Qt.UserRole+1)
        if self.model.__dict__.get('sourceFileName',None) is not None:
          #print 'updateJobCombo setting sourceFileName',self.model.__dict__['sourceFileName']
          self.jobCombo.setItemData(1,str(self.model.__dict__['sourceFileName']),QtCore.Qt.UserRole+2)
          if self.model.__dict__.get('sourceFileAnnotation',None) is not None:
            self.jobCombo.setItemData(1,str(self.model.__dict__['sourceFileAnnotation']),QtCore.Qt.UserRole+3)
        self.jobCombo.setCurrentIndex(1)
        self.validate()
    self.jobCombo.blockSignals(False)

  def setMaxChars(self):
    self.show()
    someText = 'foo/bar'
    fm = QtGui.QFontMetrics(self.jobCombo.font())
    #print 'setMaxChars',
    layout = self.layout()
    #for i in range(layout.count()):
    #  print layout.itemAt(i).widget(),layout.itemAt(i).widget().width(),
    #print ' self',self.width()
    CDataFileView.MAXCHARS = (((self.width()-40) * len(someText)) / fm.horizontalAdvance(someText))-8
    #print 'setMaxChars',self.jobCombo.width(),self.width(),fm.width(someText),CDataFileView.MAXCHARS

  @QtCore.Slot(str)
  def updateModelFromView1(self,text):
    self.updateModelFromView()

  @QtCore.Slot()
  def updateModelFromView(self):
    if self.model is None or not self.editable: return
    self.connectUpdateViewFromModel(mode=False)
    if hasattr(self,'jobCombo'):
      self.handleJobComboChange()
    else:
      fullPath = self.fileLineEdit.getValue()
      #print 'CDataFileView.updateModelFromView',fullPath
      self.model.setFullPath(fullPath)
    self.connectUpdateViewFromModel(mode=True)
    self.validate()
    #print 'CDataFileView.updateModelFromView',PREFERENCES().AUTO_INFO_ON_FILE_IMPORT

  def openViewer(self,mode=None):
    '''Open file viewer'''
    if not self.model.isSet(): return
    fileName = self.model.fullPath.get()
    #print 'CDataFileView.openViewer',fileName,mode,self.taskProjectId()
    if not os.path.exists(fileName): return
    #  PhilE hack to display MTZ files in viewhkl
    root, ext = os.path.splitext(str(fileName))
    if ext == '.mtz' and mode[0:4] == 'view':
      mode = 'view_viewhkl'
    if  mode  in ['view_text']:
      WEBBROWSER().openFile(fileName=fileName,format='text/plain')
    elif mode is not None and mode[0:4] == 'view' and len(mode)>4:
      # Open in external viewer
      LAUNCHER().openInViewer(mode[5:],fileName=fileName,projectId=self.taskProjectId(),guiParent=self)
    else:
      from qtgui import CCP4WebBrowser
      CCP4WebBrowser.OPENFILE(fileName,toFront=True)


  @QtCore.Slot(str)
  def openBrowser(self, theDir=''):
    '''Open file browser - dependent on PREFERENCES().NATIVEFILEBROWSER'''
    if self.model.qualifiers('mimeTypeDescription') is not NotImplemented:
      title = 'Open '+self.model.qualifiers('mimeTypeDescription')
    else:
      title = 'Open file'
    ifInput = self.model.qualifiers('mustExist')
    try:
      if not ifInput:
        ifInput = (str(self.model.parentContainer().objectName()) == 'inputData')
    except:
      pass
    if PREFERENCES().NATIVEFILEBROWSER:
      if self.model.qualifiers('isDirectory') and ifInput:
        fileName = QtWidgets.QFileDialog.getExistingDirectory ( self, title, '' )
      elif ifInput:
        fileName,selectedFilter = QtWidgets.QFileDialog.getOpenFileName ( self, title, theDir, self.filterText()[0] )
      else:
        fileName,selectedFilter = QtWidgets.QFileDialog.getSaveFileName ( self, title, '', self.filterText()[0])
      if fileName is not None and len(str(fileName)) > 0: self.handleBrowserOpenFile(str(fileName),None)
      return
    if self.browser is not None:
      #self.browser.exec_()
      self.browser.show()
      self.browser.raise_()
      return
    from qtgui import CCP4FileBrowser
    try:
      if not ifInput:
        ifInput = (str(self.model.parentContainer().objectName()) == 'inputData')
    except:
      pass
    #print 'CDataFileView.openBrowser',self.model.parent().objectName(),ifInput, self.model.qualifiers('isDirectory'),self.model.qualifiers('mustExist')
    if self.model.qualifiers('isDirectory'):
      if self.model.qualifiers('mustExist'):
        fileMode = QtWidgets.QFileDialog.Directory
      else:
        fileMode = CCP4FileBrowser.CFileDialog.NewDirectory
      #print 'CDataFileView.openBrowser fileMode',fileMode
      self.browser = CCP4FileBrowser.CFileDialog(self,
         title=title,
         filters= self.filterText(),
         fileMode = fileMode,
         defaultSuffix='',defaultFileName='')
      self.browser.selectDownloadedFile.connect(self.handleSelectDownloadedFile)
    elif ifInput:
      self.browser = CCP4FileBrowser.CFileDialog(self,
         title=title,
         filters= self.filterText(),
         defaultSuffix='',defaultFileName='')
      #   fileLabel = self.model.qualifiers('fileLabel'),
      self.browser.setStyleSheet("")
      self.browser.selectDownloadedFile.connect(self.handleSelectDownloadedFile)
      downloadModes = self.model.qualifiers('downloadModes')
      if downloadModes is not None and downloadModes is not NotImplemented and len(downloadModes)>0:
        projectId = None
        tw = self.parentTaskWidget()
        if tw is not None: projectId = tw.projectId()

        self.browser.setDownloadMode(modeList=downloadModes,projectId=projectId)

        if sys.version_info > (3,0):
            demo_data_dir = os.path.normpath(os.path.join(os.environ["CCP4"],"share/ccp4i2/demo_data"))
        else:
            demo_data_dir = unicode(os.path.normpath(os.path.join(os.environ["CCP4"],"share/ccp4i2/demo_data")))
        urls = self.browser.widget.fileDialog.sidebarUrls()
        url_paths = []

        for url in urls:
            if sys.version_info > (3,0):
                py_path = url.path()
            else:
                py_path = unicode(url.path())
            url_paths.append(os.path.normpath(py_path))
        if not demo_data_dir in url_paths:
            urls.append(QtCore.QUrl.fromLocalFile(demo_data_dir))
            self.browser.widget.fileDialog.setSidebarUrls(urls)
    else:
      if len(self.model.qualifiers('fileExtensions'))>0:
        defaultSuffix= self.model.qualifiers('fileExtensions')[0]
      else:
        defaultSuffix = ''
      self.browser = CCP4FileBrowser.CFileDialog(self,
         title=title,
         filters= self.filterText(),
         defaultSuffix=defaultSuffix,
         fileMode = QtWidgets.QFileDialog.AnyFile )
      self.browser.selectDownloadedFile.connect(self.handleSelectDownloadedFile)
    #self.browser.exec_()
    self.browser.show()

  def closeBrowser(self):
    # Programmatically close the browser
    if self.browser is not None: self.browser.close()

  # The method called from connect as pass through to handleBrowserOpenFile
  # Is here to ensure the arguments from the connect are consistent
  @QtCore.Slot(str,str)
  def handleSelectDownloadedFile(self,fileName,downloadInfo):
    self.handleBrowserOpenFile(fileName,downloadInfo)

  def handleBrowserOpenFile(self,fileName,downloadInfo,autoInfoOnFileImport=True,validate=True,updateView=True):
    '''Handle importing a file seelcted by the file browser (or other route)'''
    #print('CDataFileView.handleBrowserOpenFile',fileName,downloadInfo,autoInfoOnFileImport,validate,updateView)
    if downloadInfo is not None and downloadInfo.get('code',None) is not None:
      sourceFileAnnotation = 'Downloaded from '+str(downloadInfo.get('source','unknown'))+' Id: '+str(downloadInfo.get('code','unknown'))
    elif self.model.qualifiers('isDirectory'):
      sourceFileAnnotation = ''
    else:
      guiLabel = 'Object'
      if self.model.qualifiers('guiLabel') is not NotImplemented: guiLabel = self.model.qualifiers('guiLabel')
      if self.parentTaskWidget() is not None and self.parentTaskWidget().jobNumber() is not None:
        sourceFileAnnotation = guiLabel + ' loaded from ' + os.path.split(fileName)[1] + ' by job ' + self.parentTaskWidget().jobNumber()
      else:
        sourceFileAnnotation = ''
    try:
      self.connectUpdateViewFromModel(False)
      self.model.setFullPath(fileName)
      if len(sourceFileAnnotation)>0:
          self.model.annotation = sourceFileAnnotation
      self.connectUpdateViewFromModel(True)
    except CException as e:
      if isinstance(self.model,CCP4XtalData.CMiniMtzDataFile):
        pass
      else:
        e.warningMessage('Attempting to select file',parent=self)
    #print 'CDataFileView.handleBrowserOpenFile fullPath',self.model.fullPath
    if updateView: self.updateViewFromModel()
    if validate: self.validate()
    #print 'CDataFileView.handleBrowserOpenFile',autoInfoOnFileImport,PREFERENCES().AUTO_INFO_ON_FILE_IMPORT,self.model.dbFileId.isSet()
    if (self.autoInfoOnFileImport and autoInfoOnFileImport and PREFERENCES().AUTO_INFO_ON_FILE_IMPORT) and not self.model.dbFileId.isSet():
      self.openInfo(sourceFileAnnotation=sourceFileAnnotation)
    else:
      # Save to be used by PROJECTSMANAGER.importFiles() to load to Db when user runs the task
      if len(sourceFileAnnotation)>0: self.model.__dict__['sourceFileAnnotation'] = sourceFileAnnotation
    print('CDataFileView.handleBrowserOpenFile DONE')
    if hasattr(self,'parent') and callable(self.parent) and hasattr(self.parent(),'widgets') and 'structure' in self.parent().widgets and hasattr(self.parent().widgets['structure'],'widgets') and 'selection' in self.parent().widgets['structure'].widgets and hasattr(self.parent().widgets['structure'].widgets['selection'],'applySelection') and callable(self.parent().widgets['structure'].widgets['selection'].applySelection):
        self.parent().widgets['structure'].widgets['selection'].applySelection()
    return

  @QtCore.Slot(dict)
  def handleFileUpdated(self,args):
    '''Update combo box if file annotation changed'''
    if not args['key'] == 'annotation': return
    if self.editable and self.ifJobCombo:
      indx = self.jobCombo.findData(args['fileId'])
      #print 'CDataFileView.handleFileChanged',fileInfo,indx
      if indx<0: return
      text = str(self.jobCombo.itemText(indx)).split()[0] + ' '+args['value']
      self.jobCombo.setItemText(indx,text)
    elif not self.editable:
      pass

  def loadJobCombo(self):
    '''Load list of recent data objects into combo box'''
    #print 'loadJobCombo',self.taskProjectId(),self.model
    #if self.model is not None: print self.model.objectName(), self.model.__dict__.get('sourceFileAnnotation','no sourceFileAnnotation')
    if not self.editable or not self.ifJobCombo: return
    self.jobCombo.blockSignals(True)
    if self.jobCombo.count()>0:
      currentText = str(self.jobCombo.currentText())
    else:
      currentText = ''
    resetIndex = 0
    self.jobCombo.clear()
    self.nComboFiles = 0
    if self.model is None or self.model.qualifiers('allowUndefined'):
      self.jobCombo.addItem('..is not used',-1)
    else:
      self.jobCombo.addItem('..must be selected',-1)
    if self.taskProjectId() is not None and self.model is not None:
      jobListInfo =PROJECTSMANAGER().db().getJobsWithOutputFiles(projectId=self.taskProjectId(),fileType = self.model.qualifiers('mimeTypeName'),subType=self.model.qualifiers('requiredSubType'),contentFlag=self.model.requiredContent(),importFiles=True)
      #print 'loadJobCombo',self.model.objectName(),jobListInfo
      for item in jobListInfo[0:10]:
        if item['annotation'] is not None:
          title = item['annotation']
        else:
#FIXME
          title = MIMETYPESHANDLER().getMimeTypeInfo(name=item['filetype'],info='description')+' from '+TASKMANAGER().getTitle(item['taskname'])
        self.jobCombo.addItem(str(item['jobnumber'])+' '+title,item['fileid'])
        self.nComboFiles = self.nComboFiles + 1
        if str(item['jobnumber'])+' '+title == currentText: resetIndex = self.nComboFiles
      if len(jobListInfo)>10:
        self.jobCombo.insertSeparator(1000)
        self.jobCombo.addItem('More data..')
    #if self.model is not None: print 'loadJobCombo',self.model.objectName(),resetIndex
    self.jobCombo.setCurrentIndex(resetIndex)
    self.jobCombo.blockSignals(False)

  @QtCore.Slot(int)
  def handleJobComboChange(self,indx0=None):
    '''Handle a user selection from combo box'''
    #print 'handleJobComboChange',indx0
    if indx0 is None:
        indx = int(self.jobCombo.currentIndex())
    else:
        indx = int(indx0)
    if 'sourceFileAnnotation' in self.model.__dict__:
        del self.model.__dict__['sourceFileAnnotation']
    if indx == 0:
        self.model.unSet()
    elif str(self.jobCombo.itemText(indx)) == 'More data..':
        self.jobCombo.blockSignals(True)
        self.jobCombo.setCurrentIndex(0)
        self.jobCombo.blockSignals(False)
        if self.jobChooser is None:
            self.jobChooser = CJobChooser(self)
            self.jobChooser.fileSelected.connect(self.loadFileFromDb)
        self.jobChooser.populate()
        self.jobChooser.show()
    else:
        #print 'handleJobComboChange < nComboFiles',indx
        self.connectUpdateViewFromModel(False)
        fileId = CCP4Data.varToUUID(self.jobCombo.itemData(indx))
        #print 'CDataFileView.handleJobComboChange fileId',fileId
        if fileId is not None and len(fileId) > 0:
            self.loadFileFromDb(fileId,annotation=str(self.jobCombo.currentText()))
        else:
            #quuid
            self.model.blockSignals(True)
            self.model.contentFlag.unSet()
            fullPath = CCP4Data.varToUUID(self.jobCombo.itemData(indx,QtCore.Qt.UserRole+1))
            #print 'CDataFileView.handleJobComboChange fullPath',fullPath
            self.model.setFullPath(fullPath)
            self.model.annotation.set(str(self.jobCombo.currentText()))
            self.model.blockSignals(False)
            self.model.dataChanged.emit()
            try:
                self.model.__dict__['sourceFileName'] = self.jobCombo.itemData(indx,QtCore.Qt.UserRole+2).__str__()
                #print 'handleJobComboChange sourceFileName',self.model.__dict__['sourceFileName']
                self.model.__dict__['sourceFileAnnotation'] = self.jobCombo.itemData(indx,QtCore.Qt.UserRole+3).__str__()
            except:
                pass
            self.model.blockSignals(False)
    self.connectUpdateViewFromModel(True)
    self.validate()

  @QtCore.Slot(str,str)
  def loadFileFromDb(self,fileId,annotation=None):
    '''Update model based on db fileId from eith the combo box or the database search tool'''
    #print 'CDataFileView.loadFileFromDb',fileId
    if sys.version_info >= (3,0) and type(fileId) == bytes:
        fileId = fileId.decode()
    else:
        fileId = str(fileId)
    if fileId is None:
      self.model.unSet()
    else:
      fileInfo = 'UNSET'
      try:
        fileInfo = PROJECTSMANAGER().db().getFileInfo(fileId,mode=['projectid','relpath','filename','annotation','filesubtype','filecontent','filetype'])
        #print 'CDataFileView.loadFileFromDb',fileInfo
        self.model.blockSignals(True)
        self.model.set(project=fileInfo['projectid'],relPath=fileInfo['relpath'],baseName=fileInfo['filename'])
        if annotation is not None:
          self.model.annotation.set(annotation)
        else:
          self.model.annotation.set(fileInfo['annotation'])
        #print 'CDataFileView.loadFileFromDb set filename'
        self.model.subType.set(fileInfo.get('filesubtype',None))
        self.model.contentFlag.set(fileInfo.get('filecontent',None))
        #print 'CDataFileView.loadFileFromDb set contentFlag',self.model.contentFlag
        self.model.dbFileId.set(fileId)
        #print 'CDataFileView.loadFileFromDb set dbFileId'
        self.model.blockSignals(False)
        self.model.dataChanged.emit()
      except:
        print('Error retrieving file info from Db (CDataFileView.loadFileFromDb) for',fileId,self.model.objectPath(),str(self))
        self.model.unSet()
        self.jobCombo.blockSignals(True)
        self.jobCombo.setCurrentIndex(0)
        self.jobCombo.blockSignals(False)

  @QtCore.Slot(dict)
  def handleJobFinished(self,args):
    #print 'CDataFileView.handleJobFinished',args
    '''Add new file from just finished job to the combo box'''
    try:
      if args.get('projectId','x') != self.parentTaskWidget().projectId(): return
    except:
      return
    if not self.editable or not self.ifJobCombo or args.get('parentJobId','XXXX') is not None: return
    resetIndex = self.jobCombo.currentIndex()
    jobListInfo = PROJECTSMANAGER().db().getJobsWithOutputFiles(projectId=args['projectId'],jobId=args['jobId'],
                fileType = self.model.qualifiers('mimeTypeName'),subType=self.model.qualifiers('requiredSubType'),
                contentFlag=self.model.requiredContent())
    #print 'handleJobFinished jobListInfo',self.model.objectName(),jobListInfo
    if len(jobListInfo)==0: return

    self.jobCombo.blockSignals(False)
    insertIndex = 1
    try:
      sourceFileName = self.jobCombo.itemData(1,QtCore.Qt.UserRole+3).__str__()
      if len(sourceFileName)>0: insertIndex = 2
    except:
      pass
    for item in jobListInfo:
        if item['annotation'] is not None:
          title = item['annotation']
        else:
          title = MIMETYPESHANDLER().getMimeTypeInfo(name=item['filetype'],info='description')+' from '+TASKMANAGER().getTitle(item['taskname'])
        self.jobCombo.insertItem(insertIndex,str(item['jobnumber'])+' '+title,item['fileid'])
        self.nComboFiles = self.nComboFiles + 1
        if resetIndex>=insertIndex: resetIndex+=1

    self.jobCombo.setCurrentIndex(resetIndex)
    self.jobCombo.blockSignals(False)

  def showEditAnnotation(self):
    pass

  def openInfo(self,label=None,sourceFileAnnotation=''):
    '''Open the file info window for user to enter provenance of newly imported file'''
    #print 'CDataFileView.openInfo',self.model,self.model.exists(),self.model.getSourceFileName()
    if not self.model.exists(): return
    if self.model.dbFileId.isSet():
      importId = PROJECTSMANAGER().db().getFileInfo(fileId=self.model.dbFileId,mode='importId')
    else:
      importId = None
    #print 'CDataFileView.openInfo importId',importId
    from qtgui import CCP4ImpExpWidgets
    self.infoWidget = CCP4ImpExpWidgets.CImportInfoDialog(self,self.model,importId=importId,label=label,sourceFileAnnotation=sourceFileAnnotation)
    self.infoWidget.setProjectId(self.taskProjectId())
    self.infoWidget.show()
    self.infoWidget.raise_()

  @QtCore.Slot(object)
  def acceptDropData(self,textData):
    '''Reimplement drop to handle drag from outside i2'''
    #print('CDataFileView.acceptDropData',textData)
    from lxml import etree
    tree = etree.fromstring(textData)
    if tree.find('dbFileId') is not None and len(tree.find('dbFileId').text)>0:
      CComplexLineWidget.acceptDropData(self,textData)
    else:
      try:
        path = os.path.join(tree.xpath('//relPath')[0].text.__str__(),tree.xpath('//baseName')[0].text.__str__())
        #print 'CDataFileView.acceptDropData path',path
      except:
        pass
      else:
        self.handleBrowserOpenFile(path,{})

  def createMimeData(self):
    '''Reimplement drag/drop to add url for drop outside i2'''
    # For CDataFiles add mimeType and url
    mimeData = CComplexLineWidget.createMimeData(self)
    url = QtCore.QUrl()
    url.setPath(self.model.__str__())
    mimeData.setUrls([url])
    return mimeData

  @QtCore.Slot()
  def openDatabaseBrowser(self):
    '''Open window to select file from another project'''
    if self.databaseBrowser is not None:
      self.databaseBrowser.show()
      self.databaseBrowser.raise_()
      return
    from qtgui import CCP4DatabaseBrowser
    self.databaseBrowser = CCP4DatabaseBrowser.CDatabaseBrowser(parent=self,title='Open file from another project',
                                   fileType=self.model.qualifiers('mimeTypeName'),projectId=self.taskProjectId())
    self.databaseBrowser.fileSelected.connect(self.handleDatabaseBrowserOpenFile)

  @QtCore.Slot(str)
  def handleDatabaseBrowserOpenFile(self,dbFileId):
    '''Handle file selected in database browser'''
    #print 'CDataFileView.handleDatabaseBrowserOpenFile',dbFileId
    self.databaseBrowser.hide()
    self.connectUpdateViewFromModel(False)
    self.loadFileFromDb(dbFileId)
    self.connectUpdateViewFromModel(True)
    self.updateViewFromModel()
    self.validate(isValid=True)


class CJobChooser(QtWidgets.QDialog):

  fileSelected = QtCore.Signal(str)

  '''Popout window to display all files when CDataFileView.jobCombo not showing enough'''

  def __init__(self,parent):
    QtWidgets.QDialog.__init__(self,parent)
    self.setLayout(QtWidgets.QHBoxLayout())
    self.list = QtWidgets.QListWidget(self)
    self.layout().addWidget(self.list)
    self.list.itemDoubleClicked.connect(self.handleSelection)

  def populate(self):
    self.list.clear()
    jobListInfo =PROJECTSMANAGER().db().getJobsWithOutputFiles(projectId=self.parent().taskProjectId(),fileType = self.parent().model.qualifiers('mimeTypeName'),contentFlag=self.parent().model.requiredContent(),importFiles=True)
    #print 'CJobChooser.populate',jobListInfo
    for item in jobListInfo:
      if item['annotation'] is not None:
        title = item['annotation']
      else:
        title = TASKMANAGER().getTitle(item['taskname'])
      listItem = QtWidgets.QListWidgetItem(str(item['jobnumber'])+' '+title)
      listItem.setData(QtCore.Qt.UserRole,item['fileid'])
      self.list.addItem(listItem)

  @QtCore.Slot('QListWidgetItem')
  def handleSelection(self,listItem):
    fileId = CCP4Data.varToUUID(listItem.data(QtCore.Qt.UserRole))
    if fileId is not None:
      self.fileSelected.emit(fileId)
      self.close()


class CStackedWidget(CViewWidget):
  '''Use QStackedLayout to enable alternative widgets in same space in task widget'''

  def __init__(self,parent):
    CViewWidget.__init__(self,parent)
    from qtgui import CCP4TaskWidget
    self.setLayout(QtWidgets.QStackedLayout())

  @QtCore.Slot(bool)
  def handleStack(self, state=None):
    if state is None:
        state = not self.isAltViewOpen()
    indx =self.layout().currentIndex()
    if self.layout().currentWidget() is not None:
        self.layout().currentWidget().setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)
    if state:
        self.layout().setCurrentIndex(1)
    else:
        self.layout().setCurrentIndex(0)
    self.layout().currentWidget().setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
    self.adjustSize()

  def createIconButton(self):
    icon = QtGui.QIcon(PIXMAPMANAGER().getPixmap('view_forward'))
    stackButton = QtWidgets.QToolButton(self)
    stackButton.setIcon(icon)
    stackButton.setToolTip('Alternate view')
    stackButton.setIconSize(QtCore.QSize(ICONBUTTONSIZE,ICONBUTTONSIZE))
    stackButton.setMaximumHeight(ICONBUTTONSIZE)
    stackButton.setMaximumWidth(ICONBUTTONSIZE)
    stackButton.clicked.connect(self.handleStack)
    return stackButton

  def isAltViewOpen(self):
    indx = self.layout().currentIndex()
    #print 'CDataFileView.isAltViewOpen',indx
    return self.layout().currentIndex() == 1

  def addWidget(self,widget):
    self.layout().addWidget(widget)

  def getValue(self):
    value = {}
    for i in range(self.layout().count()):
      w = self.layout().widget(i)
      value.update(w.getValue)
    #print 'CStackedWidget.getValue',value
    return value

  def setValue(self,value):
    #print 'CStackedWidget.setValue',value
    for i in range(self.layout().count()):
      w = self.layout().widget(i)
      w.setValue(value)



  @QtCore.Slot()
  def updateViewFromModel(self):
    for indx in range(self.layout().count()):
      self.layout().widget(indx).updateViewFromModel()
    self.validate()

  @QtCore.Slot()
  def updateModelFromView(self):
    for indx in range(self.layout().count()):
      self.layout().widget(indx).updateModelFromView()


class CListViewListWidget( QtWidgets.QListWidget):

  rightMousePress = QtCore.Signal('QMouseEvent')
  leftMousePress = QtCore.Signal('QMouseEvent')
  edit = QtCore.Signal(int,int)
  insert = QtCore.Signal(int)
  delete = QtCore.Signal(int)

  '''Component of CListView - use a QListWidget to show items in list'''
  def __init__(self,parent=None,qualifiers={}):
    QtWidgets.QListWidget.__init__(self,parent)
    self.setDragEnabled(True)
    if qualifiers.get('editable',True):
      self.setAcceptDrops(True)
      self.setSelectionMode( QtWidgets.QListWidget.SingleSelection)
      #self.setDragDropMode(QtWidgets.QAbstractItemView.DropOnly)
      #self.setDefaultDropAction(QtCore.Qt.CopyAction)
      self.popupMenu = None
      self.rightMousePress.connect(self.showContextMenu)
      self.leftMousePress.connect(self.openEdit)
    else:
      self.setSelectionMode( QtWidgets.QListWidget.NoSelection)

  def makeItem(self,indx):
    if indx>=0:
      self.insertItem(indx,'')
    else:
      self.addItem('')
    #print 'CListViewListWidget.makeItem',indx,self.count(),repr(self)

  def setCurrentRow(self,row):
    self.blockSignals(True)
    self.setCurrentItem(self.item(row))
    self.blockSignals(False)


  def mousePressEvent(self,event):
    #print 'CListViewListWidget.mousePressEvent',event.button()
    if event.button() == QtCore.Qt.RightButton:
      self.rightMousePress.emit(event)
      event.accept()

    elif event.button() == QtCore.Qt.LeftButton:
      self.leftMousePress.emit(event)
      event.accept()

    else:
      QtWidgets.QListWidget.mousePressEvent(self,event)

  @QtCore.Slot('QMouseEvent')
  def openEdit(self,event):
    indx = self.indexAt(event.pos())
    self.edit.emit(indx.row(),indx.row())

  @QtCore.Slot('QMouseEvent')
  def showContextMenu(self,event):
    row = self.indexAt(event.pos()).row()
    column = self.indexAt(event.pos()).column()
    if self.popupMenu is None:
      self.popupMenu = QtWidgets.QMenu(self)
    else:
      self.popupMenu.clear()
    a = self.popupMenu.addAction('Edit')
    a.triggered.connect(functools.partial(self.edit.emit,row,column))
    a = self.popupMenu.addAction('Insert above')
    a.triggered.connect(functools.partial(self.insert.emit,row))
    a = self.popupMenu.addAction('Insert below')
    a.triggered.connect(functools.partial(self.insert.emit,row+1))
    a = self.popupMenu.addAction('Delete')
    a.triggered.connect(functools.partial(self.delete.emit,row))

    if sys.platform == "darwin":
        proxyStyle = MyProxyStyle()
        self.popupMenu.setStyle(proxyStyle)
    self.popupMenu.popup(event.globalPos())


class CListViewTableWidget( QtWidgets.QTableWidget):

  edit = QtCore.Signal(int)
  insert = QtCore.Signal(int)
  delete = QtCore.Signal(int)
  rightMousePress = QtCore.Signal('QMouseEvent')
  currentRowChanged = QtCore.Signal(int)

  '''Component of CListView - use a QTableWidget to show items in list'''
  def __init__(self,parent=None,qualifiers={}):
    #print 'CListViewTableWidget.__init__',qualifiers
    QtWidgets.QTableWidget.__init__(self,parent)
    #self.setDragEnabled(True)
    #self.setAcceptDrops(True)
    self.setSelectionMode( QtWidgets.QAbstractItemView.SingleSelection)
    self.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
    self.setColumnCount(len(qualifiers['tableItems']))
    self.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
    self.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
    if 'columnHeaders' in qualifiers:
      labelList = []
      for item in  qualifiers['columnHeaders']: labelList.append(item)
      self.setHorizontalHeaderLabels(labelList)
    else:
      self.setHorizontalHeader(None)
    self.verticalHeader().hide()
    #self.setDragDropMode(QtWidgets.QAbstractItemView.DropOnly)
    #self.setDefaultDropAction(QtCore.Qt.CopyAction)
    self.popupMenu = None
    self.rightMousePress.connect(self.showContextMenu)
    self.currentItemChanged.connect(self.handleItemChanged)
    #print 'CListViewTableWidget.__init__ done'

  @QtCore.Slot('QListWidgetItem','QListWidgetItem')
  def handleItemChanged(self,current,previous):
    #print 'CListViewTableWidget.handleItemChanged',current,previous
    if current is not None:
      self.currentRowChanged[int].emit(int(current.row()))

  def makeItem(self,indx):
    self.setRowCount(int(self.rowCount())+1)

  def setCurrentRow(self,row):
    self.blockSignals(True)
    self.setCurrentCell(row,0)
    self.blockSignals(False)

  def keyReleaseEvent(self,event):
    #print 'CListViewTableWidget.keyReleaseEvent',event.key()
    if event.key() == QtCore.Qt.Key_Delete:
      row = self.itemAt(self.cursor().pos()).row()
      #print 'CListViewTableWidget.keyReleaseEvent row',row
      self.delete.emit(row)

  def mousePressEvent(self,event):
    if event.button() == QtCore.Qt.RightButton:
      self.rightMousePress.emit(event)
      event.accept()
    else:
      QtWidgets.QTableWidget.mousePressEvent(self,event)

  @QtCore.Slot('QMouseEvent')
  def showContextMenu(self,event):
    row = self.indexAt(event.pos()).row()
    if self.popupMenu is None:
      self.popupMenu = QtWidgets.QMenu(self)
    else:
      self.popupMenu.clear()
    a = self.popupMenu.addAction('Edit')
    a.triggered.connect(functools.partial(self.edit.emit,row))
    a = self.popupMenu.addAction('Insert above')
    a.triggered.connect(functools.partial(self.insert.emit,row))
    a = self.popupMenu.addAction('Insert below')
    a.triggered.connect(functools.partial(self.insert.emit,row+1))
    a = self.popupMenu.addAction('Delete')
    a.triggered.connect(functools.partial(self.delete.emit,row))

    if sys.platform == "darwin":
        proxyStyle = MyProxyStyle()
        self.popupMenu.setStyle(proxyStyle)
    self.popupMenu.popup(event.globalPos())


class CListView(CComplexLineWidget):
  LINE_HEIGHT = 32
  MODEL_CLASS = CCP4Data.CList
  SIDE_BOX = False

  def __init__(self,parent=None,model=None,qualifiers={},editorQualifiers={}):
    #print 'CListView qualifiers',qualifiers
    qualis = {'gridLayout' : True, 'iconName' : 'List' }
    qualis.update(qualifiers)
    #print 'CListView init',repr(model)
    CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualis)
    self.setModel(model)
    self.viewMode = qualifiers.get('mode','list')
    if self.viewMode in ['table','tree']: self.tableItems = qualifiers['tableItems']
    #print 'CListView.__init__ viewMode',self.viewMode,qualis
    self.layout().removeWidget(self.iconButton)
    title = qualis.get('title',None)
    if title is None: title = qualis.get('label',None)
    line = QtWidgets.QHBoxLayout()
    if self.iconButton is not None:
      line.addWidget(self.iconButton)
      self.iconButton.leftMousePress.connect(self.updateMenu)
    if title is not None: line.addWidget(CItalicLabel(title,self))
    self.showListButton = QtWidgets.QToolButton(self)
    self.showListButton.setText('Show list')
    self.showListButton.setToolTip('Show a list to enable entering more than one item')
    self.showListButton.clicked.connect(functools.partial(self.setListVisible,None))
    line.insertWidget(1,self.showListButton)
    line.addStretch(2)
    if CListView.SIDE_BOX:
      self.layout().addLayout(line,0,0,1,2)
    else:
      self.layout().addLayout(line,0,0)
    # Put menu on left mouse button (as well as right)
    if self.iconButton is not None:
      self.iconButton.leftMousePress.connect(self.updateMenu)
      self.iconButton.setToolTip('Click for menu')
    #print 'CListView',self.viewMode
    if self.viewMode =='table':
      self.listWidget = CListViewTableWidget(self,qualifiers=qualifiers)
    elif self.viewMode =='tree':
      self.listWidget = CListViewTreeWidget(self,qualifiers=qualifiers)
    else:
      self.listWidget = CListViewListWidget(self,qualifiers=qualifiers)
    if CListView.SIDE_BOX:
      self.layout().addWidget(self.listWidget,1,1)
    else:
      self.layout().addWidget(self.listWidget,1,0)
    self.populateListWidget(resize=True)
    self.itemAddedToModel = False
    self.editItemIndex = 0
    if self.editable:
      self.buttonBox = self.createListButtons()
      self.editButton = self.buttonBox.layout().itemAt(2).widget()
      self.editButton = None
      if  CListView.SIDE_BOX:
        self.layout().addWidget(self.buttonBox,1,0)
      else:
        self.layout().addWidget(self.buttonBox,2,0)
      #if len(self.model) == 0:
      #  self.model.addItem()
      #  self.itemAddedToModel = True
      self.editorStack = QtWidgets.QStackedWidget()
      if CListView.SIDE_BOX:
        self.layout().addWidget(self.editorStack,2,0,1,2)
      else:
        self.layout().addWidget(self.editorStack,3,0)
      from core import CCP4DataManager
      editorClassName = qualis.get('editorClassName',None)
      #print "\n\nCListView.init qualis['editorClassName']", editorClassName
      if editorClassName is not None:
        editorClass = CCP4DataManager.DATAMANAGER().getClass(className=editorClassName)
        #print "ECN",editorClassName, editorClass
      else:
        editorClass = None
      if editorClass is None:
        editorClass = CCP4DataManager.DATAMANAGER().getWidgetClass(modelClass=model.subItemClass())
        #print "ECN",model.subItemClass(), editorClass
      #print 'CListView.init',model.objectName(),editorClass
      if editorClass is not None:
        if model is not None:
          ediQualis = model.subItemQualifiers()
        else:
         ediQualis = {}
        #print 'CListView.init',model.objectName(),ediQualis,len(model),self.editItemIndex
        ediQualis.update(editorQualifiers)
        self.editor = editorClass(self,qualifiers=ediQualis)
        #print 'CListView.init editor',self.editor
        self.editorStack.addWidget(self.editor)
        self.editorStack.setVisible(True)

      self.setListVisible(qualifiers.get('listVisible',False))

      if self.model is not None and len( self.model)>0:
         self.handleRowChange(row=0,force=True)
         self.connectDataChanged()
      self.listWidget.currentRowChanged[int].connect(self.handleRowChange)
      self.listWidget.edit.connect(self.handleRowChange)
      self.listWidget.insert.connect(self.addLine)
      self.listWidget.delete.connect(self.deleteLine)
      self.setEditorVisible(True)



  def setModel(self,value):
    self.model = value
    if self.model is not None:
      if self.model.__len__() == 0:
        self.model.addItem()
        if self.model.qualifiers('listMinLength') is NotImplemented or self.model.qualifiers('listMinLength')==0:
          #print 'CListView.setModel setting allowUndefined',repr(self.model),self.model
          self.model[0].setQualifier('allowUndefined',True)


  def handleListClicked(self,index):
    #print 'CListView.handleListClicked',index
    #print 'CListView.handleListClicked',index.row()
    pass

  def updateShowListButton(self,visible=None,enabled=None):
    #print 'updateListButton',visible,self.listWidget.isVisible(),len(self.model)
    if visible is None:
      visible = self.listWidget.isVisible()
    if visible:
      self.showListButton.setText('Hide list')
      if enabled is None:
        if self.model is None or len(self.model)>1:
          enabled = False
        else:
          enabled = True
    else:
      self.showListButton.setText('Show list')
      enabled = True
    self.showListButton.setEnabled(enabled)

  def clear(self):
    self.model.blockSignals(True)
    self.model.unSet()
    try:
      self.editor.model.unSet()
    except:
      pass
    self.model.blockSignals(False)
    #print 'CComplexLineWidget.clear',self.model.objectName(),self.model
    self.updateViewFromModel()
    self.validate()

  @QtCore.Slot(bool)
  def setListVisible(self,visible=None):
    #print('setListVisible',visible,'isVisible',self.listWidget.isVisible(),'model',repr(self.model))
    #if self.model is None or len(self.model)==0:
    #  visible = False
    if self.model is not None and len(self.model)>1:
      visible = True
    elif visible is None:
      visible = not(self.listWidget.isVisible())
    if visible: self.populateListWidget()
    self.listWidget.setVisible(visible)
    if self.editable: self.buttonBox.setVisible(visible)
    self.updateShowListButton(visible=visible)


  def isListVisible(self):
    return self.listWidget.isVisible()

  def showListEnabled(self):
    if len(self.model)>1:
      return False
    else:
      return True

  def setEditorVisible(self,visible=None):
    #print 'setEditorVisible',visible,'isVisible',self.editorStack.isVisible(),self.model
    #if self.model is None or len(self.model)==0:
    #  visible = False
    if not hasattr(self,'editorStack'): return
    if visible is None:
      visible = not(self.editorStack.isVisible())
    self.editorStack.setVisible(visible)
    if self.editButton is not None:
      if visible:
        self.editButton.setIcon(QtGui.QIcon(PIXMAPMANAGER().getPixmap('bullet_arrow_down')))
      else:
        self.editButton.setIcon(QtGui.QIcon(PIXMAPMANAGER().getPixmap('bullet_arrow_right')))


  @QtCore.Slot('QModelIndex',int,bool)
  def handleRowChange(self,row,indx1=None,force=False):
    #print 'CListView.handleRowChange',self.editItemIndex,self.model.get()
    #print 'CListView.handleRowChange input',row,force
    #traceback.print_stack(limit=5)
    #self.editor.updateModelFromView()
    if (row == self.editItemIndex) and not force:
      #self.setEditorVisible(True)
      #self.editor.setFocus1(QtCore.Qt.PopupFocusReason,indx1)
      if hasattr(self,'editor'): self.editor.setFocus(QtCore.Qt.PopupFocusReason)
      return
    self.connectDataChanged(False)
    if self.editor.model is not None: self.editor.connectUpdateViewFromModel(False)
    #print 'CListView.handleRowChange calling setModel',repr(self.model[row])
    self.editor.setModel(self.model[row])
    self.editor.updateViewFromModel()
    self.editor.validate()
    #print 'CListView.handleRowChange done calling setModel',self.editor,row,repr(self.model[row]),self.model[row]
    self.setEditorVisible(True)
    self.editItemIndex = row
    self.listWidget.setCurrentRow(self.editItemIndex)
    #self.editor.setFocus1(QtCore.Qt.PopupFocusReason,indx1)
    #print 'handleRowChange setFocus'
    self.editor.setFocus(QtCore.Qt.PopupFocusReason)
    self.connectDataChanged(True)
    self.editor.connectUpdateViewFromModel(True)

  @QtCore.Slot(int)
  def addLine(self,indx=-1):
    maxLength = self.model.qualifiers('listMaxLength')
    if maxLength is not NotImplemented and maxLength is not None and len(self.model) >= maxLength : return

    if indx>=0:
      self.editItemIndex = indx
    else:
      self.editItemIndex = len(self.model)
    #print 'CListView.addLine indx',indx,'editItemIndex',self.editItemIndex,self.model
    self.connectDataChanged(False)

    # Add row to listWidget and make it selected
    self.listWidget.makeItem(indx)
    self.listWidget.setCurrentRow(self.editItemIndex)
    #print 'CListView.addLine count after setCurrentRow',self.listWidget.count()
    # Create item in the model and make it the editor model
    #self.connectDataChanged(False)
    self.model.blockSignals(True)
    obj = self.model.addItem(index=indx)
    #self.connectDataChanged(True)
    self.model.blockSignals(False)
    #print 'CListView.addLine indx,obj',indx,obj,type(obj),self.model,self.editItemIndex,self.listWidget.count()
    self.editor.setModel(obj)
    self.editor.updateViewFromModel()
    self.setEditorVisible(True)
    # Update the new row in the listWidget with value of new model item
    if self.viewMode =='table':
      self.populateListWidget()
    else:
      self.updateListItem()

    # Reset the connect to update the listWidget
    #print 'addLine setFocus'
    self.editor.setFocus(QtCore.Qt.PopupFocusReason)
    self.connectDataChanged(True)
    self.updateShowListButton()

  @QtCore.Slot(int)
  def deleteLine(self,indx):
    minLength = self.model.qualifiers('listMinLength')
    if minLength is NotImplemented or minLength is None: minLength=0
    if len(self.model) <= minLength : return
    self.connectDataChanged(False)
    if self.viewMode == 'list': obj = self.listWidget.takeItem(indx)
    # Delete the item in the model
    try:
      obj = self.model.pop(indx)
    except:
      # Likely failed because below listMin
      self.connectDataChanged(True)
      return
    obj.deleteLater()
    # Determine the new current row
    self.editItemIndex = min( indx,len(self.model)-1)
    # Set the model item on the editor
    if self.editItemIndex>=0:
      self.editor.setModel(self.model[self.editItemIndex])
      self.editor.updateViewFromModel()
    else:
      self.editor.setModel(None)
      self.setEditorVisible(False)
    # Just repopulate table (?is there a more efficient approach)
    if self.viewMode == 'table': self.populateListWidget()
    # Set the selected listWidget row
    self.listWidget.setCurrentRow(self.editItemIndex)
    self.updateShowListButton()
    self.validate()
    # Reset the connect to update the listWidget
    self.connectDataChanged(True)

  @QtCore.Slot(str)
  def handleButtonClick(self,button):
    if button in ['append','Insert before']:
      if button == 'Insert before':
        indx = self.listWidget.currentRow()
      else:
        indx = -1
      self.addLine(indx)
    elif button == 'delete':
      #if len(self.model) == 1: return
      row = self.listWidget.currentRow()
      self.deleteLine(row)
    elif button == 'editor':
      self.setEditorVisible()


  @QtCore.Slot()
  @QtCore.Slot(bool)
  def populateListWidget(self,resize=False):
    #print 'populateListWidget resize',resize
    if self.model is None: return
    if self.viewMode == 'table':
      self.listWidget.setRowCount(len(self.model))
      for row in range(len(self.model)):
        tableTextItems = self.model[row].getTableTextItems()
        if tableTextItems is not None:
          #print 'populateListWidget tableTextItems',tableTextItems
          for col in range(0,len(tableTextItems)):
            self.listWidget.setItem(row,col,QtWidgets.QTableWidgetItem(tableTextItems[col]))
      if resize:
        self.listWidget.horizontalHeader().resizeSections(QtWidgets.QHeaderView.ResizeToContents)
        self.listWidget.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Interactive)
    else:
      self.listWidget.blockSignals(True)
      self.listWidget.clear()
      #print 'populateListWidget',self.model.objectName(),len(self.model)
      for n in range(len(self.model)):
        #print 'populateListWidget',n
        try:
          value = self.model[n].getTextItem()
        except Exception as e:
          #print 'CListView.populateListWidget',e
          value = '--'
        #print 'CListView.populateListWidget',self.model.objectName(),n,value
        self.listWidget.addItem(value)
        self.setValidity(n)
      #if getattr(self,'editItemIndex',None) is not None: self.handleRowChange(0,True)
      self.listWidget.blockSignals(False)

  @QtCore.Slot()
  def updateViewFromModel(self):
    #for ii in range(len(self.model)): print repr(self.model[ii]),self.model[ii]
    self.populateListWidget()
    # Beware major change to list may have zapped the current self.editor.model
    # Beware no editor in non-editable mode
    try:
      if hasattr(self,'editor'):
        if not self.model.count(self.editor.model):
          self.handleRowChange(0)
        else:
          self.editor.updateViewFromModel()
    except Exception as e:
      print('ERROR updating editor in CListView.updateViewFromModel')
      print(e)

  @QtCore.Slot()
  def updateModelFromView(self):
    # The model should be kept in sync at the appropriate times (or
    # we are in poo)
    pass

  @QtCore.Slot()
  def updateListItem(self):
    #print 'CListView.updateListItem',self.editItemIndex,self.model[self.editItemIndex]
    if self.viewMode =='table':
      tableTextItems = self.model[self.editItemIndex].getTableTextItems()
      #print 'CListView.updateListItem tableTextItems',tableTextItems
      for ic in range(0,len(tableTextItems)):
        w = self.listWidget.item(self.editItemIndex,ic)
        #print 'CListView.updateListItem',ic,w,tableTextItems[ic]
        if w is not None: w.setText(tableTextItems[ic])
    else:
      value = self.model[self.editItemIndex].getTextItem()
      listItem = self.listWidget.item(self.editItemIndex)
      #print 'CListView.updateListItem value',self.editItemIndex,value,'listItem',listItem,self.listWidget.count(),repr(self.listWidget)
      if listItem is not None: listItem.setText(value)
    self.setValidity(self.editItemIndex)
    self.validate()

  def connectDataChanged(self,mode=True):
    #print 'CListView.connectDataChanged',mode,self.editItemIndex,repr(self.editor.model)
    if self.editor.model is None: return
    if mode:
      self.editor.model.dataChanged.connect(self.updateListItem)
    else:
      pass#self.editor.model.dataChanged.disconnect(self.updateListItem)


  def setValidity(self,row,isValid=None):
    if isValid is None:
      validity = self.model[row].validity(self.model[row].get())
      if validity.maxSeverity()>SEVERITY_WARNING:
        isValid = False
      else:
        isValid = True
    #print 'CListView.setValidity',self.model.objectName(),row,isValid
    if isinstance(self.listWidget,QtWidgets.QListWidget):
      widgetRow = self.listWidget.item(row)
    else:
      widgetRow = None
    if widgetRow is not None:
      if isValid:
#FIXME PYQT - or maybe None? This used to set QVariant.
        widgetRow.setData(QtCore.Qt.DecorationRole,"")
      else:
        widgetRow.setData(QtCore.Qt.DecorationRole,QtGui.QIcon(PIXMAPMANAGER().getPixmap('error_star')))

    # Check overall list validity
    validity = self.model.validity(self.model.__dict__['_value'])
    previous = self.listWidget.property('isValid')
    #print 'CListView.setValidity',previous,validity.maxSeverity()
    if validity.maxSeverity()>SEVERITY_WARNING:
      self.listWidget.setProperty('isValid',False)
    else:
      self.listWidget.setProperty('isValid',True)
    self.listWidget.style().unpolish(self)
    self.listWidget.style().polish(self)
    self.listWidget.update()


  def getNofLines(self):
    if self.viewMode=='table':
      return self.listWidget.rowCount()
    else:
      return self.listWidget.count()

  def setNofLines(self,nLines):
    if self.viewMode == 'table':
      self.listWidget.setRowCount(nLines)
    else:
      # Add or delete lines to get to NLines
      curNLines = self.listWidget.count()
      #print 'exFrame.setNofLines',curNLines,nLines,self.parent()
      if nLines == curNLines:
        return
      elif nLines > curNLines:
        for ii in range(curNLines,nLines):
          self.addLine()

      elif nLines < curNLines and nLines>=0:
        for ii in range(curNLines,nLines,-1):
          #print 'exFrame.setNofLines deleteLine',ii-1
          self.deleteLine(ii-1)

  def getValue(self,index=-1):
    value = []
    if index<0:
      first = 0
      last = self.body.layout().count()
    else:
      first = index
      last = index + 1
    if self.viewMode =='table':
      for row in range(first,last):
        v = {}
        for col in range(self.listWidget.coloumnCount()):
          v[self.tableItems[col]] = str(self.listWidget.item(row,col).text())
        value.append(v)
    else:
      for row in range(first,last):
        value.append(str(self.listWidget.item(row).text()))

    return value

  def setValue(self,value=[],index=-1):

     #print 'CListView.setValue',value
     if self.viewMode=='table':
       self.listWidget.setRowCount(len(value))
       for row in range(len(value)):
         tableItemValues = value[row].getTableTextItems()
         for col in range(len(tableItemValues)):
           self.listWidget.setItem(row,col,QtWidgets.QTableWidgetItem(tableItemValues[col]))
     else:
       self.listWidget.clear()
       for rowValue in value:
         self.listWidget.addItem(rowValue)

  def getMenuDef(self):
    if self.editable:
      menu = ['show_list','clear','copy','paste','help']
    else:
      menu = ['show_list','copy','help']
    if self._stacked: menu.insert(0,'handleStack')
    return menu

  def isWidgetOpen(self):
    return self.listWidget.isVisible()

  def handleHideWidget(self):
    if self.isWidgetOpen():
      self.listWidget.hide()
      self.editorStack.hide()
      for ii in range(self.buttonBox.layout().count()):
        self.buttonBox.layout().itemAt(ii).widget().hide()
    else:
      self.listWidget.show()
      self.editorStack.show()
      for ii in range(self.buttonBox.layout().count()):
        self.buttonBox.layout().itemAt(ii).widget().show()

  @QtCore.Slot(bool,bool)
  def validate(self,isValid=None,reportMessage=True):
    if not getattr(self,'editable',True): return True
    CComplexLineWidget.validate(self,isValid=isValid,reportMessage=reportMessage)
    editorStack = getattr(self,'editorStack',None)
    if editorStack is not None:  editorStack.currentWidget().validate()


class CFileListView(CListView):

  def __init__(self,parent=None,model=None,qualifiers={},editorQualifiers={}):
    qualis = { 'editorClassName' : CDataFileView }
    qualis.update(qualifiers)
    CListView.__init__(self,parent,model,qualis,editorQualifiers=editorQualifiers)

class CTreeViewAbstractItemModel(QtCore.QAbstractItemModel):

  def __init__(self,parent,model):
     QtCore.QAbstractItemModel.__init__(self,parent)
     self.rootItem = model
     self.rootItem.setAbstractModelParent(self)
     self._headerData = []
     self.unsetCurrentItem()
     self.itemMap = {}

  def setHeaderData(self,value):
    self._headerData = value

  def columnCount(self,parent):
    return len(self._headerData)

  def headerData(self,section,orientation,role):
#FIXME PYQT - or maybe None? This used to return QVariant.
    if orientation != QtCore.Qt.Horizontal: return None
    value = self._headerData[section].get(role,None)
    if value is None:
#FIXME PYQT - or maybe None? This used to return QVariant.
      return None
    else:
      return value

  def childCount(self):
    return self.rootItem.__len__()

  def child(self,row):
    return self.rootItem.__getitem__(row)

  def data(self, index, role):
    if not index.isValid():
      return None
    item = index.internalPointer()
    return item.data(index.column(),role)

  def index(self, row, column, parent):
    if row < 0 or column < 0 or row >= self.rowCount(parent) or column >= self.columnCount(parent):
      return QtCore.QModelIndex()
    if not parent.isValid():
      parentItem = self.rootItem
    else:
      parentItem = parent.internalPointer()
    childItem = parentItem.child(row)
    #print 'CProjectModel.index',row, column, parent,childItem.objectPath()
    if childItem:
      return self.createIndex(row, column, childItem)
    else:
      return QtCore.QModelIndex()

  def parent(self, index):
    if not index.isValid():
      return QtCore.QModelIndex()
    childItem = index.internalPointer()
    parentItem = childItem.abstractModelParent()
    if parentItem == self.rootItem:
      return QtCore.QModelIndex()
    return self.createIndex(parentItem.row(), 0, parentItem)


  #MN Is something like this needed to avoid identical model items
  #appearing to be different entities
  #def createIndex(self, row, column, item):
  #  itemStr = str(row) + '_' + str(column) + '_' + str(item)
  #  if not itemStr in self.itemMap:
  #    self.itemMap[itemStr] = super(CTreeViewAbstractItemModel,self).createIndex(row,column,item)
  #  return self.itemMap[itemStr]


  def rowCount(self, parent):
    #if parent.column() > 0:
    #  return 0
    if not parent.isValid():
      parentItem = self.rootItem
    else:
      parentItem = parent.internalPointer()
    return parentItem.childCount()

  def row(self):
    return 0

  def unsetCurrentItem(self):
    self.unsetCurrentItem0(self.rootItem)

  def unsetCurrentItem0(self,node):
    for n in range(node.childCount()):
      if 'currentItem' in node.child(n).__dict__:
        del node.child(n).__dict__['currentItem']
      self.unsetCurrentItem0(node.child(n))

  def setCurrentItem(self,obj):
    if isinstance(obj,QtCore.QModelIndex): obj = obj.internalPointer()
    if obj is None: return
    obj.__dict__['currentItem'] = True

  def currentItem(self):
    rv = self.currentItem0(self.rootItem)
    #print 'CTreeViewAbstractItemModel.currentItem',rv
    return rv

  def currentItem0(self,node):
    if node.childCount() == 0: return None
    for n in range(node.childCount()):
      if 'currentItem' in node.child(n).__dict__ and node.child(n).__dict__['currentItem']:
        return [node.child(n),n]
    for n in range(node.childCount()):
      rv = self.currentItem0(node.child(n))
      if rv is not None:
        rv.append(n)
        return rv
    return None


  def onlyOneItem(self):
    node = self
    while node.childCount() > 0:
      if node.childCount() >1: return False
      node = node.child(0)
    return True

  def lastItem(self):
    # Find the last item in the tree to be set as the default object in
    # the editor
    if self.childCount() == 0:
      return QtCore.QModelIndex()
    else:
      indx = self.child(self.childCount()-1)
      if indx.childCount()>0:
        return indx.child(indx.childCount()-1)
      else:
        return indx

  def dataObj2ModelIndex(self,dataObj):
    qVar = str(repr(dataObj))
    indexList = self.match(self.index(0,0,QtCore.QModelIndex()),QtCore.Qt.UserRole,qVar,1)
    #print 'dataObj2ModelIndex',repr(dataObj),indexList
    if len(indexList)>0: return indexList[0]
    #top = self.index(0,0,QtCore.QModelIndex())
    for ir in range(self.rootItem.childCount()):
      node = self.index(ir,0,QtCore.QModelIndex())
      indexList = self.match( self.index(0,0,node),QtCore.Qt.UserRole,qVar,1)
      #print 'dataObj2ModelIndex',ir,indexList
      if len(indexList)>0: return indexList[0]
    return None

  '''
  def removeRows(self, position=0, count=1, parent=QtCore.QModelIndex()):
    node = self.nodeFromIndex(parent)
    self.beginRemoveRows(parent, position, position + count - 1)
    node.childItems.pop(position)
    self.endRemoveRows()
  '''

class CTreeViewTreeView(QtWidgets.QTreeView):

  rightMousePress = QtCore.Signal('QMouseEvent')
  leftMouseRelease = QtCore.Signal('QModelIndex')
  edit = QtCore.Signal('QModelIndex')
  insert = QtCore.Signal('QModelIndex',int)
  delete = QtCore.Signal('QModelIndex')

  def __init__(self,parent):
    QtWidgets.QTreeView.__init__(self,parent)
    self.setSelectionMode(QtWidgets.QAbstractItemView.NoSelection)
    self.popupMenu = None
    self.rightMousePress.connect(self.showContextMenu)

  def keyReleaseEvent(self,event):
    #print 'CListViewTableWidget.keyReleaseEvent',event.key(),QtCore.Qt.Key_Delete
    if event.key() == QtCore.Qt.Key_Delete:
      pos = self.mapFromGlobal(self.cursor().pos())-QtCore.QPoint(0,self.header().height())
      modelIndex = self.indexAt(pos)
      #print 'CListViewTableWidget.keyReleaseEvent item',pos,modelIndex
      if modelIndex is None: return
      #print 'CListViewTableWidget.keyReleaseEvent',indices
      self.delete.emit(modelIndex)
      event.accept()
    else:
      event.ignore()

  def mousePressEvent(self,event):
    if event.button() == QtCore.Qt.RightButton:
        self.rightMousePress.emit(event)
    QtWidgets.QTreeView.mousePressEvent(self,event)

  def mouseReleaseEvent(self,event):
    if event.button() == QtCore.Qt.LeftButton:
      modelIndex = self.indexAt(event.pos())
      #print 'CTreeViewTreeView emiting leftMouseRelease'
      self.leftMouseRelease.emit(modelIndex)
    QtWidgets.QTreeView.mouseReleaseEvent(self,event)


  @QtCore.Slot('QMouseEvent')
  def showContextMenu(self,event):
    if not self.parent().editable: return
    modelIndex = self.indexAt(event.pos())
    dataObject = modelIndex.internalPointer()
    iEditor = 0
    for editorDefn in self.parent().editorDefnList:
      if isinstance(dataObject,editorDefn.get('modelClass')): break
      iEditor += 1
    #print 'showContextMenu dataObject,iEditor',dataObject,iEditor
    if self.popupMenu is None:
      self.popupMenu = QtWidgets.QMenu(self)
    else:
      self.popupMenu.clear()
    a = self.popupMenu.addAction('Edit')
    a.triggered.connect(functools.partial(self.edit.emit,modelIndex))
    label = self.parent().editorDefnList[iEditor].get('label')
    a = self.popupMenu.addAction('Insert '+label+' above')
    a.triggered.connect(functools.partial(self.insert.emit,modelIndex, self.parent().editorDefnList[iEditor].get('name'),0))
    a = self.popupMenu.addAction('Insert '+label+' below')
    a.triggered.connect(functools.partial(self.insert.emit,modelIndex, self.parent().editorDefnList[iEditor].get('name'),1))
    for editorDefn in self.parent().editorDefnList[iEditor+1:]:
      a = self.popupMenu.addAction('Append '+editorDefn.get('label','?'))
      a.triggered.connect(functools.partial(self.insert.emit,modelIndex,editorDefn.get('name','?')),-1)
    a = self.popupMenu.addAction('Delete')
    a.triggered.connect(functools.partial(self.delete.emit,modelIndex))

    if sys.platform == "darwin":
        proxyStyle = MyProxyStyle()
        self.popupMenu.setStyle(proxyStyle)
    self.popupMenu.popup(event.globalPos())


class CTreeView(CComplexLineWidget):
  LINE_HEIGHT = 32
  SIDE_BOX = False

  def __init__(self,parent=None,model=None,qualifiers={}):
    from core import CCP4DataManager
    qualis = {'gridLayout' : True, 'iconName' : 'List' }
    qualis.update(qualifiers)
    #print 'CListView qualis',qualis
    CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualis)
    self.setModel(model)
    editorClosable = False
    self.editorDefnList = qualifiers.get('editors')
    self.columnHeaders = qualifiers.get('columnHeaders')
    # The title/icon line
    self.layout().removeWidget(self.iconButton)
    title = qualis.get('title',None)
    line = QtWidgets.QHBoxLayout()
    if self.iconButton is not None: line.addWidget(self.iconButton)
    #self.showListIcons = {}
    #self.showListIcons['enabled'] = CCP4GuiUtils.createIcon(name='ShowList')
    #self.showListIcons['disabled'] = CCP4GuiUtils.createIcon(name='List_grey')
    self.showListButton = QtWidgets.QToolButton(self)
    self.showListButton.setText('Show list')
    self.showListButton.setToolTip('Show a list to enable entering more than one item')
    #self.showListButton.setIconSize(QtCore.QSize(ICONBUTTONSIZE,ICONBUTTONSIZE))
    #self.showListButton.setMaximumSize(QtCore.QSize(ICONBUTTONSIZE,ICONBUTTONSIZE))
    self.showListButton.clicked.connect(functools.partial(self.setListVisible,None))
    line.insertWidget(1,self.showListButton)
    if title is not None: line.addWidget(CItalicLabel(title,self))
    line.addStretch(2)
    if CTreeView.SIDE_BOX:
      self.layout().addLayout(line,0,0,1,2)
    else:
      self.layout().addLayout(line,0,0)
    # Put menu on left mouse button (as well as right)
    self.iconButton.leftMousePress.connect(self.updateMenu)
    self.iconButton.setToolTip('Click for menu')
    # The tree widget
    self.abstractModel = CTreeViewAbstractItemModel(self,self.model)
    self.abstractModel.setHeaderData(qualifiers.get('columnHeaders',[]))
    self.listWidget = CTreeViewTreeView(self)
    self.listWidget.setModel(self.abstractModel)
    section=-1
    for column in qualifiers.get('columnHeaders',[]):
      section+=1
      width = column.get('width',None)
      if width is not None: self.listWidget.setColumnWidth(section,width)
    if CTreeView.SIDE_BOX:
      self.layout().addWidget(self.listWidget,1,1)
    else:
      self.layout().addWidget(self.listWidget,1,0)
    #self.abstractModel.setCurrentRow()
    self.listWidget.expandAll()
    self.editItemIndex = [None]
    # Edit buttons
    if self.editable:
      self.buttonBox = QtWidgets.QFrame(self)
      if CTreeView.SIDE_BOX:
        self.buttonBox.setLayout(QtWidgets.QVBoxLayout())
      else:
        self.buttonBox.setLayout(QtWidgets.QHBoxLayout())
      self.buttonBox.layout().setContentsMargins(0,0,0,0)
      editMenu = [['list_add_grey','Add item','append'],['list_delete_grey','Remove item','delete']]
      if editorClosable: editMenu.append(['bullet_arrow_down','Open editor','editor'])
      for item in editMenu:
        icon = QtGui.QIcon(PIXMAPMANAGER().getPixmap(item[0]))
        button = QtWidgets.QToolButton(self)
        button.setIcon(icon)
        button.setMaximumHeight(ICONBUTTONSIZE)
        button.setMaximumWidth(ICONBUTTONSIZE)
        button.setToolTip(item[1])
        self.buttonBox.layout().addWidget(button)
        if item[2] == 'append' and len(self.editorDefnList)>0:
          menu = QtWidgets.QMenu(self)
          for editorDefn in self.editorDefnList:
            action = menu.addAction('Add '+editorDefn.get('label'))
            action.setObjectName(editorDefn.get('name'))
            action.triggered.connect(functools.partial(self.handleButtonClick,'list_add_'+editorDefn.get('name')))
          button.setMenu(menu)
          button.setPopupMode(QtWidgets.QToolButton.InstantPopup)
        else:
          button.clicked.connect(functools.partial(self.handleButtonClick,item[2]))
      if not CListView.SIDE_BOX: self.buttonBox.layout().addStretch(1)
      if CTreeView.SIDE_BOX:
        self.layout().addWidget(self.buttonBox,1,0)
      else:
        self.layout().addWidget(self.buttonBox,2,0)
      if editorClosable:
        self.editButton = self.buttonBox.itemAt(2).widget()
        self.editButton.setCheckable(True)
      else:
        self.editButton = None
      # Set up the editors for individual items in a QStackedWidget
      self.editorStack = QtWidgets.QStackedWidget()
      if CTreeView.SIDE_BOX:
        self.layout().addWidget(self.editorStack,2,0,1,2)
      else:
        self.layout().addWidget(self.editorStack,3,0)
      from core import CCP4DataManager
      for editorDefn in self.editorDefnList:
        editor = CCP4DataManager.DATAMANAGER().widget(modelClass=editorDefn.get('modelClass'),parentWidget=self.editorStack)
        editor.setObjectName(editorDefn.get('modelClass').__name__)
        #print 'CTreeView.init editor',editorDefn,editor
        self.editorStack.addWidget(editor)
      self.editorStack.setVisible(True)
      self.handleRowChange(dataObj = self.abstractModel.lastItem())
      if self.model is None or self.abstractModel.onlyOneItem():
        self.setListVisible(False)
      else:
        self.setListVisible(True)
      self.connectDataChanged()
#FIXME - SJM 16/9/2020, I don't think this can ever be true.
      if hasattr(self.listWidget,"currentRowChanged"):
          self.listWidget.currentRowChanged.connect(self.handleRowChange)
      self.listWidget.edit.connect(self.handleRowChange)
      self.listWidget.insert.connect(self.addLine)
      self.listWidget.delete.connect(self.deleteLine)
      self.listWidget.leftMouseRelease.connect(self.handleRowChange)

  def resetAbstractModel(self):
    self.abstractModel = CTreeViewAbstractItemModel(self,self.model)
    self.abstractModel.setHeaderData(self.columnHeaders)
    self.listWidget.setModel(self.abstractModel)

  def setModel(self,value):
    self.model = value
    if self.model is not None and self.model.__len__() == 0:
      self.model.addItem()
    if self.model.qualifiers('listMinLength') is NotImplemented or self.model.qualifiers('listMinLength')==0:
      #print 'CTreeView.setModel setting allowUndefined',self.model.objectPath(),self.model
      self.model[0].setQualifier('allowUndefined',True)

  @QtCore.Slot(bool)
  def setListVisible(self,visible=None):
    #print('setListVisible',visible,'isVisible',self.listWidget.isVisible(),'model',self.model)
    #if self.model is None or len(self.model)==0:
    #  visible = False
    if not self.abstractModel.onlyOneItem():
      visible = True
    elif visible is None:
      visible = not(self.listWidget.isVisible())
    #if visible: self.populateListWidget()
    self.listWidget.setVisible(visible)
    self.buttonBox.setVisible(visible)
    self.editorStack.setCurrentIndex(1)
    '''
    if visible:
      self.editButton.setIcon(QtGui.QIcon(PIXMAPMANAGER().getPixmap('bullet_arrow_down')))
    else:
      self.editButton.setIcon(QtGui.QIcon(PIXMAPMANAGER().getPixmap('bullet_arrow_right')))
    '''
    self.update()
    self.updateShowListButton(visible=visible)

  def isListVisible(self):
    return self.listWidget.isVisible()

  def updateShowListButton(self,visible=None,enabled=None):
    #print 'updateListButton',visible,self.listWidget.isVisible(),len(self.model)
    if visible is None:
      visible = self.listWidget.isVisible()
    if visible:
      self.showListButton.setText('Hide list')
      if enabled is None:
        if not self.abstractModel.onlyOneItem():
          enabled = False
        else:
          enabled = True
    else:
      self.showListButton.setText('Show list')
      enabled = True
    self.showListButton.setEnabled(enabled)


  def showListEnabled(self):
    return self.abstractModel.onlyOneItem()


  def getMenuDef(self):
    if self.editable:
      menu = ['show_list','clear','copy','paste','help']
    else:
      menu = ['show_list','copy','help']
    if self._stacked: menu.insert(0,'handleStack')
    return menu

  def setEditorVisible(self,visible=None):
    if visible is None: visible = not(self.editorStack.isVisible())
    self.editorStack.setVisible(visible)
    if self.editButton is not None:
      if visible:
        self.editButton.setIcon(QtGui.QIcon(PIXMAPMANAGER().getPixmap('bullet_arrow_down')))
      else:
        self.editButton.setIcon(QtGui.QIcon(PIXMAPMANAGER().getPixmap('bullet_arrow_right')))


  def setEditor(self,mode,obj=None):
    w = self.editorStack.findChild(QtWidgets.QWidget,mode)
    if w is not None:
      self.editorStack.setCurrentWidget(w)
      if obj is not None:
        w.setModel(obj)
        w.updateViewFromModel()
      self.setEditorVisible(True)
      w.setFocus(QtCore.Qt.PopupFocusReason)
      w.validate()
    else:
      print('Error in setEditor for mode:',mode)
      print('traceback.print_stack:')
      traceback.print_stack()
    return w

  @QtCore.Slot('QModelIndex',object)
  def handleRowChange(self,modelIndex=None,dataObj=None):
    #traceback.print_stack(limit=8)
    if modelIndex is not None and modelIndex.isValid():
      dataObj = modelIndex.internalPointer()
    #print 'CTreeView.handleRowChange',modelIndex,dataObj
    self.connectDataChanged(False)
    self.setEditor(dataObj.__class__.__name__,dataObj)
    self.abstractModel.unsetCurrentItem()
    self.abstractModel.setCurrentItem(dataObj)
    self.connectDataChanged(True)

  @QtCore.Slot('QModelIndex',str,int,object)
  def addLine(self,modelIndex=None,mode=None,position=-1,currentItem=None):
    # mode is name value of definition in self.editorDefnList
    # modelIndex or currentItem (CTreeAbstractModel.currentItem) are a selected item
    # that a new item will be placed before (position=0) or after (position=1) or
    # new item will be appended if position=-1
    #print 'CTreeView.addLine',modelIndex,mode,position,currentItem
    if mode is None: mode = self.editorDefnList[0].get('name')
    self.connectDataChanged(False)
    #print 'CTreeView.addLine dataObj',modelIndex,mode,repr(dataObj)
    if mode == self.editorDefnList[0].get('name'):
      listObj = self.abstractModel.rootItem
      if position==-1 or (modelIndex is None and currentItem is None):
        insertIndex =listObj.__len__()
        listIndex = -1
      else:
        if modelIndex is not None and modelIndex.isValid():
          dataObj = modelIndex.internalPointer()
        elif currentItem is not None:
          dataObj = currentItem[0]
        insertIndex = dataObj.row()
        insertIndex += position
        if insertIndex >= listObj.__len__():
          listIndex = -1
        else:
          listIndex = insertIndex
      self.abstractModel.beginInsertRows(QtCore.QModelIndex(),insertIndex,insertIndex)
      newDataObj = listObj.addItem(index=listIndex)
      self.abstractModel.endInsertRows()
      # Try getting Qt to understand the newly inserted object has a child too
      newModelIndex = self.abstractModel.index(insertIndex,0,QtCore.QModelIndex())
      self.abstractModel.beginInsertRows(newModelIndex,0,0)
    else:
      dataObj = None
      modelClass = None
      for ediDef in self.editorDefnList:
        if ediDef.get('name') == mode: modelClass = ediDef.get('modelClass')
      # We are going to need modelIndex and dataObj when adding item below the top level
      if modelIndex is None and currentItem is None:
        print('ERROR in CTreeView.addLine - no  selected model index for insertion point')
        return
      elif modelIndex is not None and modelIndex.isValid():
        dataObj = modelIndex.internalPointer()
      elif currentItem is not None:
        dataObj = currentItem[0]
        modelIndex = self.abstractModel.dataObj2ModelIndex(dataObj)
      if dataObj is None or modelIndex is None:
        print('ERROR in CTreeView.addLine - no consistent selected model index/ data object for insertion point')
        return
      # If selected item is of the same type as item to be added then go
      # up the heirarchy one rung for the parent
      if isinstance(dataObj,modelClass):
        insertIndex = dataObj.row() + position
        modelIndex = self.abstractModel.parent(modelIndex)
        dataObj = modelIndex.internalPointer()
        listObj = dataObj.childListObject()
        if insertIndex>=len(listObj):
          listIndex = -1
        else:
          listIndex = insertIndex
      else:
        listObj = dataObj.childListObject()
        insertIndex = listObj.childCount()
        listIndex = -1
      self.abstractModel.beginInsertRows(modelIndex,insertIndex,insertIndex)
      newDataObj = listObj.addItem(index=listIndex)
    self.abstractModel.endInsertRows()
    self.listWidget.expandAll()

    self.setEditor(newDataObj.__class__.__name__,newDataObj)
    self.abstractModel.unsetCurrentItem()
    self.abstractModel.setCurrentItem(newDataObj)
    self.updateShowListButton()
    self.validate()
    self.connectDataChanged(True)

  @QtCore.Slot('QModelIndex',object)
  def deleteLine(self,modelIndex=None,dataObj=None):
    # Delete the item in the model
    if modelIndex is not None and modelIndex.isValid():
      dataObj = modelIndex.internalPointer()
    if dataObj is None: return
    if modelIndex is None:  modelIndex = self.abstractModel.dataObj2ModelIndex(dataObj)
    if modelIndex is None: return
    self.connectDataChanged(False)
    self.listWidget.blockSignals(True)
    listObj = dataObj.parent()
    listIndex = dataObj.row()
    #print 'CTreeView.deleteLine',repr(dataObj),modelIndex,repr(listObj),listIndex
    # If this is the model for the editor need to clear up before deleting it
    # Otherwise risk runtime error trying to access deleted C++ object
    # But if deleting dataObj fails (eg goes below listMinLength) then this is not necessary
    editorModel = None
    for idx in range(self.editorStack.count()):
      if self.editorStack.widget(idx).model == dataObj:
         editorModel = dataObj
         editorIndex = idx
         self.editorStack.widget(idx).setModel(None)
    #print 'deleteLine',repr(dataObj),repr(listObj),listIndex
    self.abstractModel.beginRemoveRows(modelIndex.parent(),listIndex,listIndex)
    try:
      listObj.__delitem__(listIndex)
    except CException as e:
      #print 'deleteLine __delitem__ failed'
      if e.count(code=107):
         QtWidgets.QMessageBox.warning(self,self.windowTitle(),'Can not delete item\nList must contain at least one item')
      else:
        e.warningMessage(message='Error attempting to delete item',windowTitle='Error editing',parent=self)
      self.abstractModel.endRemoveRows()
      if editorModel is not None: self.editorStack.widget(editorIndex).setModel(editorModel)
    except Exception as e:
       QtWidgets.QMessageBox.warning(self,self.windowTitle(),'Unknown error attempting to delete item')
    else:
      #print 'deleteLine __delitem__ ok'
      self.abstractModel.endRemoveRows()
      self.listWidget.update()
      if listIndex>=listObj.__len__(): listIndex = listObj.__len__()-1
      if listIndex>=0:
        self.handleRowChange(dataObj=listObj.__getitem__(listIndex))
      self.updateShowListButton()
    self.connectDataChanged(True)
    self.validate()
    self.listWidget.blockSignals(False)
    #self.handleRowChange()

  def handleButtonClick(self,button):
    rv = self.abstractModel.currentItem()
    #print 'handleButtonClick',button,'currentItem',rv
    if button == 'delete':
      if rv is None: return
      self.deleteLine(dataObj=rv[0])
    elif button == 'editor':
      self.setEditorVisible()
    elif button[0:9] == 'list_add_':
      self.addLine(currentItem=rv,mode=button[9:],position=-1)

  @QtCore.Slot()
  def updateViewFromModel(self):
    #self.populateListWidget()
    self.listWidget.update()
    self.editorStack.currentWidget().updateViewFromModel()

  @QtCore.Slot()
  def updateModelFromView(self):
    # The model should be kept in sync at the appropriate times
    pass

  @QtCore.Slot()
  def updateEditListItem(self):
    self.listWidget.update()
    self.validate()

  @QtCore.Slot()
  def updateListItem(self):
    self.listWidget.update()

  def connectDataChanged(self,mode=True):
    #print 'CListView.connectDataChanged',mode,self
    editor = self.editorStack.currentWidget()
    if editor.model is None: return
    if mode:
      editor.model.dataChanged.connect(self.updateEditListItem)
    else:
      pass#editor.model.dataChanged.disconnect(self.updateEditListItem)

  def setValidity(self,row,isValid=None):
    return
    if isValid is None:
      validity = self.model[row].validity(self.model[row].get())
      if validity.maxSeverity()>SEVERITY_WARNING:
        isValid = False
      else:
        isValid = True
    #print 'CListView.setValidity',self.model.objectName(),row,validity.report(),isValid
    widgetRow = self.listWidget.item(row)
    if widgetRow is not None:
      if isValid:
#FIXME PYQT - or maybe None? This used to set QVariant.
        widgetRow.setData(QtCore.Qt.DecorationRole,"")
      else:
        widgetRow.setData(QtCore.Qt.DecorationRole,QtGui.QIcon(PIXMAPMANAGER().getPixmap('error_star')))

  def validate(self,isValid=None,reportMessage=True):
    if not getattr(self,'editable',True): return True
    CComplexLineWidget.validate(self,isValid=isValid,reportMessage=reportMessage)
    editorStack = getattr(self,'editorStack',None)
    if editorStack is not None:  editorStack.currentWidget().validate()

  def getValue(self,index=-1):
    print('CTreeView.getValue doing nothing')
    return []

  def setValue(self,value=[],index=-1):
     print('CTreeView.setValue doing nothing',value)


class CStringView(CViewWidget):

    MODEL_CLASS = CCP4Data.CString

    def __init__(self,parent=None,model=None,qualifiers={}):
      qualis = {}
      qualis.update(qualifiers)
      CViewWidget.__init__(self,parent,qualifiers=qualis)

      layout = QtWidgets.QHBoxLayout()
      layout.setSpacing(CViewWidget.MARGIN)
      layout.setContentsMargins(CViewWidget.MARGIN,CViewWidget.MARGIN,CViewWidget.MARGIN,CViewWidget.MARGIN)
      #print 'CStringView qualis',qualis.get('guiMode',None)
      if qualis.get('editable',True) or qualis.get('guiMode','NotAMultiLineText') == 'multiLine':
        self.mode = 'edit'
        if qualifiers.get('enumerators',[]) or (model is not None and len(model.qualifiers('enumerators'))>0): self.mode = 'combo'
        if qualis.get('modelClass',None) is not None:
          enumerators = CCP4Data.classQualifier(qualis['modelClass'],'enumerators')
          if enumerators is not None and len(enumerators)>0: self.mode = 'combo'
        if qualis.get('guiMode',None) is not None: self.mode =  qualis.get('guiMode')
      else:
        self.mode = 'label'
      #print 'CStringView mode',qualis,self.mode
      if self.mode == 'combo':
        comboQualis = {}
        if 'modelClass' in qualifiers:
          #comboQualis.update(qualifiers['modelClass']['qualifiers'])
          comboQualis.update(CCP4Data.classQualifiers(qualis['modelClass']))
        else:
          comboQualis.update(qualifiers)
        self.widget = CComboBox(self,qualifiers=comboQualis)
        # Use the CComboBox signal that covers edit and menu change signals
        self.widget.dataChanged.connect(self.updateModelFromView)
        layout.addWidget(self.widget)
      elif ['radio','multiLineRadio'].count(self.mode):
        if self.mode == 'multiLineRadio':
          layout = QtWidgets.QVBoxLayout()
          layout.setSpacing(CViewWidget.MARGIN)
          layout.setContentsMargins(CViewWidget.MARGIN,CViewWidget.MARGIN,CViewWidget.MARGIN,CViewWidget.MARGIN)
        labelText = qualis.get('label',NotImplemented)
        if labelText is not NotImplemented and labelText is not None:
          label = QtWidgets.QLabel(qualis.get('label'),self)
          layout.addWidget(label)
        self.widget = CRadioButtonGroup(self)
        #self.widget.buttonReleased[int].connect(self.updateModelFromView1)
        self.widget.buttonReleased.connect(self.updateModelFromView1)
      elif self.mode == 'label':
        labelText = qualis.get('label')
        if labelText is not NotImplemented and labelText is not None:
          label = QtWidgets.QLabel(labelText,self)
          layout.addWidget(label)
        self.widget = CBoldLabel(parent=self)
        layout.addWidget(self.widget)
      elif self.mode == 'multiLine':
        self.widget =CTextEdit(self)
        layout.addWidget(self.widget)
        self.widget.textChanged.connect(self.updateModelFromView)
        self.setFocusProxy(self.widget)
      else:
        if model is not None:
          self.widget = CLineEdit(self,qualifiers=model.qualifiers())
        else:
          self.widget = CLineEdit(self )
        if qualis.get('charWidth',None) is not None: self.widget.setCharWidth(qualis['charWidth'])
        self.widget.textEdited.connect(self.updateModelFromView)
        layout.addWidget(self.widget)
        self.setFocusProxy(self.widget)
      self.setLayout(layout)
      if self.mode == 'combo':
        self.populateComboBox(model,modelClass=qualifiers.get('modelClass',None))
      elif ['radio','multiLineRadio'].count(self.mode):
        self.populateRadioButtonBox(model,modelClass=qualifiers.get('modelClass',None))
      elif self.mode == 'label' and model is not None:
        self.widget.setText(str(model.getMenuValue()))
      self.setModel(model)

    def setToolTip(self,tip):
      CViewWidget.setToolTip(self,tip)
      self.widget.setToolTip(tip)

    @QtCore.Slot()
    def updateViewFromModel(self):
      #print 'CStringView.updateViewFromModel',self.blockUpdateView
      if self.blockUpdateView: return
      #import traceback
      #traceback.print_stack(limit=10)
      if self.model is None: return
      self.widget.blockSignals(True)
      if self.model is not None and self.model.get() is not None:
        if self.mode == 'label':
          # Beware model may be CInt or CFloat
          text = str(self.model.getMenuValue())
        else:
          text = self.model.get()
        if self.mode == 'combo':
          try:
            indx = self.model.qualifiers('enumerators').index(text)
            #print 'CStringView.updateViewFromModel indx',indx
            self.widget.setValue(self.model.qualifiers('menuText')[indx])
          except:
            self.widget.setValue(text)
        else:
          self.widget.setValue(text)
      else:
        self.widget.setValue(None)
      self.widget.blockSignals(False)

    def getWidgetText(self):
      if self.mode == 'combo':
        if self.widget.isEditable():
          value = str(self.widget.currentText())
          try:
            value = self.model.qualifiers('enumerators')[self.model.qualifiers('menuText').index(value)]
          except:
            pass
        else:
          value = self.widget.itemData(self.widget.currentIndex()).__str__()
      elif ['radio','multiLineRadio','multiLine'].count(self.mode):
        value = self.widget.getValue()
      else:
        txt = self.widget.text()
        try:
          value = str(txt)
        except Exception as e:
          value = ''
          print('ERROR invalid character')
          print(e)
      return value

    def updateModelFromText(self):
      self.updateModelFromView()

    @QtCore.Slot(str)
    def updateModelFromView1(self,text):
      self.updateModelFromView()

    @QtCore.Slot()
    def updateModelFromView(self):
      if self.model is None: return
      #traceback.print_stack(limit=6)
      value = self.getWidgetText()
      #print 'CStringView.updateModelFromView',self.model.objectName(),value
      self.connectUpdateViewFromModel(False,propagate=True)
      try:
        if isinstance(value,str) and len(value) == 0:
          self.model.set(value=None)
        else:
          #self.model.set(value=self.model.PYTHONTYPE(value))
          self.model.set(value=value)
      except CException as e:
        report=e.report(ifStack=False,mode=2)
        #print 'CStringView updateModelFromView error',value,type(value),report
        self.validate(False,report=report)
      else:
        self.validate()
      self.connectUpdateViewFromModel(True,propagate=True)

    def populateComboBox(self,model=None,modelClass=None):
      #print 'populateComboBox',self.model.qualifiers('enumerators'),self.model.qualifiers('menuText')
      if model is not None:
        menuText = model.qualifiers('menuText')
        enumerators = model.qualifiers('enumerators')
      elif modelClass is not None:
        menuText = CCP4Data.classQualifier(modelClass,'menuText')
        enumerators = CCP4Data.classQualifier(modelClass,'enumerators')
      else:
        return
      #MN: Note that the "frozen" (ie. running/run) version of this widget will not have a populate method
      if hasattr(self.widget,'populate'): self.widget.populate(enumerators,menuText)

    def populateRadioButtonBox(self,model=None,modelClass=None):
      if model is not None:
        menuText = model.qualifiers('menuText')
        enumerators = model.qualifiers('enumerators')
      elif modelClass is not None:
        menuText = CCP4Data.classQualifier(modelClass,'menuText')
        enumerators = CCP4Data.classQualifier(modelClass,'enumerators')
      #print 'populateRadioButtonBox',menuText,enumerators
      if self.mode == 'radio':
        for ii in range(0,min(len(enumerators),len(menuText))):
          but = self.widget.addRadioButton(enumerators[ii],menuText[ii])
          #print 'populateRadioButtonBox',but
          self.layout().addWidget(but)
      else:
        for ii in range(0,min(len(enumerators),len(menuText))):
          line = QtWidgets.QHBoxLayout()
          but = self.widget.addRadioButton(enumerators[ii],menuText[ii])
          #print 'populateRadioButtonBox',but
          line.addWidget(but)
          self.layout().addLayout(line)


class CFloatView(CStringView):

    MODEL_CLASS = CCP4Data.CFloat

    def __init__(self,parent=None,model=None,qualifiers={}):
        qualis = qualifiers
        CStringView.__init__(self,parent,model=model,qualifiers=qualis)
        #if self.editable: self.setValidator(qualis)

    def setModel(self,model):
      CViewWidget.setModel(self,model)
      #if self.editable: self.setValidator()

    def setValidator(self,qualifiers={}):
        if self.model is None: return
        minValue=self.model.qualifiers('min')
        maxValue=self.model.qualifiers('max')
        if minValue is None and maxValue is None: return
        validator = CDoubleValidator(bottom=minValue,top=maxValue,parent=self)
        if self.mode in ['combo','edit']:
          self.widget.setValidator(validator)

    @QtCore.Slot()
    def updateModelFromView1(self):
      # This is called after each key input but updating a CFloat value
      # can lead to autocompleting eg '1' to '1.0' which sticks in a zero
      # when the user is typing.  So don't do the model update and trust
      # that the call to updateModelFromView when widget losses focus will
      # update the model
      return


class CIntView(CStringView):

    MODEL_CLASS = CCP4Data.CInt

    def __init__(self,parent=None,model=None,qualifiers={}):
        #print 'CIntView.__init__',parent,model,qualifiers
        qualis = qualifiers
        CStringView.__init__(self,parent,model=model,qualifiers=qualis)
        #if self.editable: self.setValidator(qualis)

    def setModel(self,model):
      CViewWidget.setModel(self,model)
      #if self.editable: self.setValidator()

    def setValidator(self,qualifiers={}):
      if self.model is None: return
      minValue=self.model.qualifiers('min')
      maxValue=self.model.qualifiers('max')
      validator = QtGui.QIntValidator(self)
      if minValue is not None: validator.setBottom(minValue)
      if maxValue is not None: validator.setTop(maxValue)
      if self.mode in ['combo','edit']:
        self.widget.setValidator(validator)


class CBooleanView(CViewWidget):

    MODEL_CLASS = CCP4Data.CBoolean

    def __init__(self,parent=None,model=None,qualifiers={}):
        #print 'CBooleanView.__init__',model.objectName()
        qualis = qualifiers
        CViewWidget.__init__(self,parent,qualis)
        layout = QtWidgets.QHBoxLayout()
        layout.setSpacing(CViewWidget.MARGIN)
        layout.setContentsMargins(CViewWidget.MARGIN,CViewWidget.MARGIN,CViewWidget.MARGIN,CViewWidget.MARGIN)
        if self.editable:
          self.widget = CCheckBox(self)
          self.widget.stateChanged.connect(self.updateModelFromView1)
        else:
          self.widget = CCheckBoxUneditable(self)
        layout.addWidget(self.widget)
        self.setLayout(layout)
        if model is not None:
          self.setModel(model)
          self.updateViewFromModel()

    @QtCore.Slot()
    def updateViewFromModel(self):
      #print 'CBooleanView.updateViewFromModel',self.model.objectName(),self.model.get()
      if self.model is not None and self.model.get() is not None:
        self.widget.blockSignals(True)
        self.widget.setChecked( bool(self.model.get()) )
        self.widget.blockSignals(False)

    @QtCore.Slot()
    def updateModelFromView1(self):
      self.updateModelFromView()

    @QtCore.Slot()
    def updateModelFromView(self):
      #print 'CBooleanView.updateModelFromView',self.model.objectName(),self.widget.isChecked()
      if self.model is None: return
      self.connectUpdateViewFromModel(False)
      self.model.set(value=self.widget.isChecked())
      self.connectUpdateViewFromModel(True)


class CRangeView(CViewWidget):

   MODEL_CLASS = CCP4Data.CRange

   def __init__(self,parent=None,model=None,qualifiers={}):
     qualis = qualifiers
     CViewWidget.__init__(self,parent,qualifiers=qualis)
     layout = QtWidgets.QHBoxLayout()
     layout.setSpacing(CViewWidget.MARGIN)
     layout.setContentsMargins(CViewWidget.MARGIN,CViewWidget.MARGIN,CViewWidget.MARGIN,CViewWidget.MARGIN)
     #print 'CRangeView.__init__',model

     if isinstance(model,CCP4Data.CIntRange):
       self.widgets['start'] = CIntView(self,model.start,qualifiers=qualis)
       self.widgets['end'] = CIntView(self,model.end,qualifiers=qualis)
     else:
       self.widgets['start'] = CFloatView(self,model.start,qualifiers=qualis)
       self.widgets['end'] = CFloatView(self,model.end,qualifiers=qualis)


     layout.addWidget(self.widgets['start'])
     layout.addWidget(QtWidgets.QLabel('to',self))
     layout.addWidget(self.widgets['end'])
     self.setLayout(layout)

     if model is not None: self.setModel(model)


class CGenericGridView(CViewWidget):

   MODEL_CLASS = None
   ITEMS = []

   def __init__(self,parent=None,model=None,qualifiers={}):
     qualis = qualifiers
     CViewWidget.__init__(self,parent,qualifiers=qualis)
     #print 'CGenericView.__init__',model,model.objectName()
     self.setLayout(QtWidgets.QGridLayout())
     self.layout().setSpacing(CViewWidget.MARGIN)
     self.layout().setContentsMargins(CViewWidget.MARGIN,CViewWidget.MARGIN,CViewWidget.MARGIN,CViewWidget.MARGIN)
     maxRowLength = 0
     for row in self.getWidgetItems():
       maxRowLength = max(maxRowLength,len(row))
     from core import CCP4DataManager
     iRow = -1
     for rowItems in self.getWidgetItems():
        iRow = iRow + 1
        iCol = -1
        for item in rowItems:
          iCol = iCol + 1
          l = QtWidgets.QLabel(item,self)
          self.layout().addWidget(l,iRow,iCol)
          iCol = iCol + 1
          #self.widgets[item] =CLineEdit(self)
          try:
            self.widgets[item] =CCP4DataManager.DATAMANAGER().widget(model=model.get(item),parentWidget=self)
          except:
            try:
              self.widgets[item] =CStringView(parent=self,model=model.get(item))
            except:
              pass
          if item in self.widgets:
            #print 'CGenericView.__init__',item,self.widgets[item]
            self.widgets[item].editingFinished.connect(self.updateModelFromView)
            #print 'CGenericGridView',len(rowItems),maxRowLength,item == rowItems[-1]
            if len(rowItems) < maxRowLength and item == rowItems[-1]:
              # Is last item on a short tine so give it the remaining columns
              colSpan = 1 + (maxRowLength-len(rowItems))*2
            else:
              colSpan = 1
            self.layout().addWidget(self.widgets[item],iRow,iCol,1,colSpan)
          else:
            print('CGenericView.__init__',item,'No widget')
     if model is not None: self.setModel(model)

   @QtCore.Slot()
   def updateViewFromModel(self):
     #print 'CGenericGridView.updateViewFromModel'
     if self.model is not None:
       for item in list(self.widgets.keys()):
         self.widgets[item].updateViewFromModel()

   @QtCore.Slot()
   def updateModelFromView(self):
     if self.model is not None:
       for item in list(self.widgets.keys()):
         self.widgets[item].updateModelFromView()

   def getWidgetItems(self):
     return self.__class__.ITEMS

class CI2XmlHeaderView(CGenericGridView):

   MODEL_CLASS = CCP4File.CI2XmlHeader
   ITEMS = [ ['function','ccp4iVersion'],
              ['pluginName','pluginVersion'],
              ['pluginTitle'],
              ['userId','creationTime'],
              ['comment'] ]


class CPatchSelectionView(CComplexLineWidget):

  MODEL_CLASS = CCP4Data.CPatchSelection

  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = {}
    qualis.update(qualifiers)
    qualis['vboxLayout'] = True
    CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualis)
    self.setFrameShape(QtWidgets.QFrame.Box)
    icon = self.layout().takeAt(0).widget()
    icon.deleteLater()
    line = QtWidgets.QHBoxLayout()
    lab = QtWidgets.QLabel('Use control parameters and/or apply command file patches from one of the following:',self)
    # Beware the 'italic' label is used in handlePatchListChanged
    lab.setObjectName('italic')
    line.addWidget(lab)
    self.layout().addLayout(line)
    self.patchLayout= QtWidgets.QVBoxLayout()
    self.layout().addLayout(self.patchLayout)
    if model is not None: self.setModel(model)
    self.viewWidgetDict = {}
    COMFILEPATCHMANAGER().listChanged.connect(self.handlePatchListChanged)

  def setModel(self,model):
    if model is not None:
      self.draw(model)
    CComplexLineWidget.setModel(self,model)

  def draw(self,model):
    for patchName,patchTitle in model.getPatches():
      line = QtWidgets.QHBoxLayout()
      if self.editable:
        cb = CCheckBox(self)
      else:
        cb = CCheckBoxUneditable(self)
      if patchTitle is not None:
        cb.setText(str(patchTitle))
      else:
        cb.setText(str(patchName))
      cb.setObjectName(str(patchName))
      cb.clicked.connect(functools.partial(self.updateModelFromView,str(patchName)))
      but = QtWidgets.QPushButton('View',self)
      but.clicked.connect(functools.partial(self.viewPatch,str(patchName)))
      line.addWidget(cb)
      line.addStretch(1)
      line.addWidget(but)
      self.patchLayout.addLayout(line)
    #self.layout().addStretch(5)

  @QtCore.Slot()
  def updateViewFromModel(self):
    if self.model is None: return
    #print 'CPatchSelectionView.updateViewFromModel',self.model.patch
    for cb in self.findChildren(QtWidgets.QCheckBox):
      cb.setChecked(False)
    cb = self.findChild(QtWidgets.QCheckBox,self.model.patch.__str__())
    if cb is not None: cb.setChecked(True)

  @QtCore.Slot(str)
  def updateModelFromView(self,patchName=None):
    if patchName is not None:
      cb0 = self.findChild(QtWidgets.QCheckBox,patchName)
      #print 'CPatchSelectionView.updateModelFromView',cb0
      if cb0 is not None and cb0.isChecked():
        for cb in self.findChildren(QtWidgets.QCheckBox):
          if cb != cb0: cb.setChecked(False)
    if self.model is None: return
    self.connectUpdateViewFromModel(False)
    self.model.patch.unSet()
    for cb in self.findChildren(QtWidgets.QCheckBox):
      if cb.isChecked(): self.model.patch.set(str(cb.objectName()))

    self.connectUpdateViewFromModel(True)

  @QtCore.Slot(str)
  def viewPatch(self,patchName):
    if self.viewWidgetDict.get(patchName,None) is None:
      from qtgui import CCP4ComFilePatchManagerGui
      self.viewWidgetDict[patchName] = CCP4ComFilePatchManagerGui.CCreatePatchDialog(self,new=False)
      self.viewWidgetDict[patchName].loadPatch(patchName)
    self.viewWidgetDict[patchName].show()
    self.viewWidgetDict[patchName].raise_()

  @QtCore.Slot()
  def handlePatchListChanged(self):
    # refresh the model CPatchSelection.patchsForTask list
    self.model.set(self.model.fix({'taskName' : self.model.taskName.__str__(), 'patch' : self.model.patch.__str__() }))

    widgetlist = self.findChildren(QtWidgets.QWidget)
    for widget in widgetlist:
      if widget.objectName() != 'italic': widget.deleteLater()
    for idx in range(self.patchLayout.count()-1,-1,-1):
      try:
        item = self.patchLayout.takeAt(idx)
        item.layout().deleteLater()
        #print 'handlePatchListChanged deleted',idx
      except:
        print('CPatchSelectionView.handlePatchListChanged failed',idx)

    # redraw
    self.draw(self.model)

class CFollowFromJobView(CComplexLineWidget):
  MODEL_CLASS = CCP4Data.CFollowFromJob
  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = {'iconName':'job'}
    qualis.update(qualifiers)
    CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualis)
    self.setFrameShape(QtWidgets.QFrame.Box)
    #self.setLineWidth(2)
    self.layout().addWidget(QtWidgets.QLabel('Use data from job',self))
    self.projectId = qualifiers.get('projectId',None)
    if self.projectId is None:
      self.projectId =self.parentTaskWidget().projectId()
    self.combo = CFinishedJobsCombo(self,self.projectId)
    self.listView=None
    self.force = False
    self.setModel(model)
    self.layout().addWidget(self.combo)
    self.layout().addWidget(QtWidgets.QLabel('as input below..',self))
    self.layout().addStretch(1)
    self.combo.load()
    self.combo.currentIndexChanged[int].connect(functools.partial(self.updateModelFromView,True))
    self.combo.blockSignals(True)


  def getActionDef(self,name=None):
    #print 'CFollowFromJobView.updateMenu',self,type(self),self.getMenuDef()
    if name == 'paste':
      pasteText = 'Paste'
      mimeData = QTAPPLICATION().clipboard().mimeData()
      for item in mimeData.formats():
        #print 'CFollowFromJobView.updateMenu data:',str(item),mimeData.data(item)
        if str(item).startswith('taskParameters') and str(item)[15:]==self.parentTaskWidget().taskName():
          from lxml import etree
          root = etree.fromstring(str(mimeData.data(item)))
          jobNo = root.find('jobNumber').text
          projectName = root.find('projectName').text
          if root.find('projectId').text == self.parentTaskWidget().projectId():
            pasteText = 'Paste parameters from job '+jobNo
          else:
            pasteText = 'Paste parameters from '+projectName+' '+jobNo
      def e(): return (self.editable and self.clipboardLoaded())
      return dict (
        text = self.tr(pasteText),
        tip = self.tr('Paste data'),
        enabled = e,
        slot = self.paste
      )
    else:
      return CComplexLineWidget.getActionDef(self,name)



  def dropTypes(self):
    return ['taskParameters_'+self.parentTaskWidget().taskName()]

  def currentJobId(self):
    return self.combo.currentJobId()

  def reset(self):
    followFrom = PROJECTSMANAGER().db().getProjectFollowFromJobId(projectId=self.combo.projectId)
    self.model.set(followFrom)

  def getValue(self):
    return self.combo.currentJobId()

  @QtCore.Slot()
  def updateViewFromModel(self):
    self.combo.blockSignals(True)
    jobId = self.model.get()
    #print 'CFollowFromJobView.updateViewFromModel',jobId,self.combo.signalsBlocked()
    if jobId is not None:
      self.combo.set(self.model.get())
    else:
      self.combo.set(-1)
    self.combo.blockSignals(False)

  def updateInputFiles(self,jobId,projectId=None,force=False):
    #if self.model.get() is None: return
    #print 'CFollowFromJobView.updateInputFiles',jobId,self.model.get()
    if self.parentTaskWidget() is not None and hasattr(self.parentTaskWidget(),"followFromJobUpdated"):
        self.parentTaskWidget().followFromJobUpdated.emit(jobId,projectId)

  @QtCore.Slot(bool,int)
  def updateModelFromView(self,updateInputFiles=False,dum=12):
    jobId = self.combo.currentJobId()
    #print 'CFollowFromJobView.updateModelFromView',jobId,'*',type(jobId)
    if jobId is None:
      self.model.unSet()
      if updateInputFiles: self.updateInputFiles(None,force=True)
      return
    elif jobId != 'more':
      self.model.set(jobId)
      if updateInputFiles: self.updateInputFiles(jobId,force=True)
    else:
      if self.listView is None:
        self.listDialog = QtWidgets.QDialog(self)
        self.listDialog.setWindowTitle('Use data from job..')
        #self.listDialog.setWindowFlags(QtCore.Qt.FramelessWindowHint)
        self.listDialog.setLayout(QtWidgets.QHBoxLayout())
        self.listDialog.layout().setSpacing(0)
        self.listDialog.layout().setContentsMargins(0,0,0,0)
        from dbapi import CCP4DbApi
        self.listView = CFinishedJobsListWidget(self.listDialog,self.projectId,jobStatus = [CCP4DbApi.JOB_STATUS_FINISHED,CCP4DbApi.JOB_STATUS_INTERRUPTED])
        self.listDialog.layout().addWidget(self.listView)
        #self.listDialog.setModal(True)
        #self.listView.set(self.combo.currentJobId())
        self.listView.currentItemChanged.connect(self.handleListDialog)
      self.combo.blockSignals(True)
      self.combo.set(jobId=self.model.get())
      self.combo.blockSignals(False)
      self.listDialog.show()
      self.listDialog.raise_()
      self.listView.load()
      self.listDialog.setFocus(QtCore.Qt.OtherFocusReason)
      self.listDialog.move(self.mapToGlobal(self.combo.pos()))

  @QtCore.Slot('QListWidgetItem','QListWidgetItem')
  def handleListDialog(self,currentItem,previousItem):
    #jobId = CCP4Data.varToUUID(currentItem.data(QtCore.Qt.UserRole))
    if currentItem is None:
      self.model.unSet()
      return
    text = currentItem.text().__str__()
    jobId = PROJECTSMANAGER().db().getJobId(projectId=self.projectId,jobNumber=text.split()[0])
    #print 'CFollowFromJobView.handleListDialog', jobId
    if jobId is None: return
    self.listDialog.close()
    self.combo.addJob(jobId)
    self.combo.set(jobId)
    self.model.set(jobId)
    self.updateInputFiles(jobId,force=True)

  @QtCore.Slot(object)
  def acceptDropData(self,textData):
    #print 'CFollowFromJob.acceptDropData',textData, textData.count('taskParameters')
    if textData.count('taskParameters') > 0:
      from lxml import etree
      root = etree.fromstring(textData)
      #print 'CFollowFromJob.acceptDropData',root.find('taskName').text,self.parentTaskWidget().taskName()
      try:
        if root.find('taskName').text == self.parentTaskWidget().taskName():
          self.parentTaskWidget().loadControlParameters(jobId = root.find('jobId').text)
      except:
        return
    else:
      CComplexLineWidget.acceptDropData(self,textData)
      self.updateInputFiles(jobId=self.model.get(),force=True)


class CFinishedJobsListWidget(QtWidgets.QListWidget):

  def __init__(self,parent,projectId,jobStatus=[]):
    QtWidgets.QListWidget.__init__(self,parent)
    self.projectId = projectId
    self.jobStatus = jobStatus
    self.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
    PROJECTSMANAGER().db().jobFinished.connect(self.handleJobFinished)

  def load(self):
    self.clear()
    jobInfoList = PROJECTSMANAGER().db().getProjectJobListInfo(projectId=self.projectId,
          jobStatus=self.jobStatus,topLevelOnly=True,
          mode=['jobid','jobnumber','taskname'],order='DESC')
    #print 'CFinishedJobsListWidget.load',len(jobInfoList)
    self.blockSignals(True)
    for jobInfo in jobInfoList:
      title = str(jobInfo['jobnumber'])+' '+TASKMANAGER().getTitle(jobInfo['taskname'])
      #print 'CFinishedJobsListWidget.load',title,jobInfo['jobid']
      item = QtWidgets.QListWidgetItem(title,self)
      item.setData(QtCore.Qt.UserRole,jobInfo['jobid'])
      self.addItem(item)
      #self.addItem(title)
    self.blockSignals(False)
    #print 'CFinishedJobsListWidget.load',self.count()

  @QtCore.Slot(dict)
  def handleJobFinished(self,args):
    if args['projectId'] != self.projectId or args['status'] not in self.jobStatus: return
    self.addJob(args['jobId'])

  def addJob(self,jobId):
    try:
      jobInfo = PROJECTSMANAGER().db().getJobInfo(jobId=jobId,mode=['jobnumber','taskname'])
      jobNumber = jobInfo['jobnumber']
      taskName = jobInfo['taskname']
    except:
      return
    self.blockSignals(True)
    item = QtWidgets.QListWidgetItem(str(jobNumber)+' '+TASKMANAGER().getTitle(taskName),self)
    self.insertItem(0,item)
    self.blockSignals(False)
    #self.extras.append({'jobid':jobId,'jobnumber':jobNumber,'taskname':taskName})

  def set(self,jobId):
    row = self.findData(jobId)
    #print 'CFinishedJobsListWidget.set',jobId,row
    if row>=0: self.setCurrentRow(row)


  def findData(self,jobId):
    for i in range(self.count()):
      var = self.item(i).data(QtCore.Qt.UserRole)
      itemJobId = CCP4Data.varToUUID(var)
      if itemJobId == jobId: return i
    return -1

  def selectedJobs(self):
    selectedIndices =  self.selectionModel().selectedIndexes()
    #print 'CFinishedJobsListWidget.selectedJobs',selectedIndices
    jobList = []
    for idx in selectedIndices:
      jobId = CCP4Data.varToUUID(self.model().data(idx,QtCore.Qt.UserRole))
      if jobId is not None: jobList.append(jobId)
    return jobList

class CJobTitleView(CStringView):
  STRETCH = 5
  MODEL_CLASS = CCP4Data.CJobTitle

class CFinishedJobsCombo(QtWidgets.QComboBox):

  def __init__(self,parent,projectId):
    QtWidgets.QComboBox.__init__(self,parent)
    self.setEditable(False)
    self.projectId = projectId
    #self.extras = []
    self.more = False
    PROJECTSMANAGER().db().jobFinished.connect(self.handleJobFinished)
    PROJECTSMANAGER().db().jobUpdated.connect(self.handleJobUpdated)


  def showPopup(self):
      '''This should replicate method in CComboBox to prevent list view being drawn too short'''
      QtWidgets.QComboBox.showPopup(self)
      if sys.platform != 'darwin' : return
      listViewList = self.findChildren(QtWidgets.QListView)
      if len(listViewList)>0:
        listViewList[0].setWrapping(True)
        if CComboBox.FONTHEIGHT is None: CComboBox.FONTHEIGHT = listViewList[0].fontMetrics().height()+2
        listViewList[0].setMinimumHeight(self.count()*CComboBox.FONTHEIGHT+5)
        #print 'CComboBox.showPopup height',self.count(),CComboBox.FONTHEIGHT,listViewList[0].height()

  def load(self):
    MAXJOBS = 10
    #print 'CFinishedJobCombo.load'
    from dbapi import CCP4DbApi
    jobInfoList = PROJECTSMANAGER().db().getProjectJobListInfo(projectId=self.projectId,
          jobStatus=[CCP4DbApi.JOB_STATUS_FINISHED],topLevelOnly=True,maxJobs=MAXJOBS+1,
          mode=['jobid','jobnumber','taskname','jobtitle'],order='DESC')
    #print 'CFinishedJobCombo.load',jobInfoList
    self.blockSignals(True)
    self.clear()
    self.addItem('No',-1)
    self.insertSeparator(1)
    for item in jobInfoList[0:10]:
      if item['jobtitle'] is not None:
        title = str(item['jobnumber'])+' '+item['jobtitle']
      else:
        title = str(item['jobnumber'])+' '+TASKMANAGER().getTitle(item['taskname'])
      self.addItem(title,item['jobid'])

    if len(jobInfoList)>MAXJOBS:
      self.insertSeparator(999)
      self.addItem('More jobs..')
      self.more = True
    else:
      self.more = False
    self.blockSignals(False)

  @QtCore.Slot(dict)
  def handleJobFinished(self,args):
    from dbapi import CCP4DbApi
    if args['projectId'] != self.projectId or args['status']!= CCP4DbApi.JOB_STATUS_FINISHED: return
    self.addJob(args['jobId'])

  @QtCore.Slot(dict)
  def handleJobUpdated(self,args):
    if  args['projectId'] != self.projectId or args.get('key') != 'jobtitle': return
    self.blockSignals(True)
    jobId = self.currentJobId()
    self.load()
    self.set(jobId)
    self.blockSignals(False)

  def addJob(self,jobId,jobNumber=None,taskName=None):
    indx = self.findData(jobId)
    #print 'CFinishedJobsCombo.addJob',jobId,indx
    if indx>=0: return
    if jobNumber is None or taskName is None:
      try:
        jobInfo = PROJECTSMANAGER().db().getJobInfo(jobId=jobId,mode=['jobnumber','taskname'])
        jobNumber = jobInfo['jobnumber']
        taskName = jobInfo['taskname']
      except:
        return
    title = jobNumber+' '+TASKMANAGER().getTitle(taskName)
    self.blockSignals(True)
    self.insertItem(1,title,jobId)
    self.blockSignals(False)
    #self.extras.append({'jobid':jobId,'jobnumber':jobNumber,'taskname':taskName})

  def set(self,jobId=None,jobNumber=None,taskName=None):
    #print 'CFinishedJobsCombo.set',jobId,jobNumber,taskName
    if jobId is None or jobId==-1:
      indx = 0
    else:
      indx = self.findData(jobId)
      if indx<0:
        self.addJob(jobId,jobNumber,taskName)
      indx = self.findData(jobId)

    if self.currentIndex() == indx: return
    self.blockSignals(True)
    self.setCurrentIndex(indx)
    self.blockSignals(False)
    #print 'CFinishedJobsCombo.set DONE'

  def currentJobId(self):
    if self.currentIndex() == 0:
      return None
    if self.itemData(self.currentIndex()) is None:
      return 'more'
    jobId = CCP4Data.varToUUID(self.itemData(self.currentIndex()))
    if len(jobId) ==0:
      #Do not htink we should get here any more with QVariant/QString api 2.
      return 'more'
    else:
      return jobId

class CExePathListView(CListView):
  MODEL_CLASS = CCP4File.CExePathList
  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = {
               'columnHeaders':['Program name','Use this version'],
               'tableItems' : ['exeName','exePath'],
               'editorClassName' : 'CExePathView'
               }
    qualis.update(qualifiers)
    #print 'CExePathListView.__init__ model',model
    CListView.__init__(self,parent,model=model,qualifiers=qualis)



class CExePathView(CComplexLineWidget):
  MODEL_CLASS = CCP4File.CExePath

  def __init__(self,parent=None,model=None,qualifiers={}):
    CComplexLineWidget.__init__(self,parent,qualifiers)
    self.widgets = {}
    #line = QtWidgets.QHBoxLayout()
    #iconWidgetItem = self.layout().takeAt(0)
    #line.addWidget(iconWidgetItem.widget())
    self.layout().addWidget(QtWidgets.QLabel('Name',self))
    self.widgets['exeName'] = CStringView(parent=self)
    self.widgets['exeName'].setToolTip('The name of a program')
    self.layout().addWidget(self.widgets['exeName'])
    self.widgets['exeName'].setMaximumWidth(150)
    qualis = { 'jobCombo' : False, 'browseDb' : False , 'toolTip' : 'Select the (non-standard) program to use','autoInfoOnFileImport': False }
    self.widgets['exePath'] =  CDataFileView(parent=self,qualifiers=qualis)
    #self.widgets['exePath'].layout().takeAt(0)
    self.layout().addWidget(self.widgets['exePath'])
    self.setModel(model)


class CSearchPathListView(CListView):
  MODEL_CLASS = CCP4File.CSearchPathList
  def __init__(self,parent=None,model=None,qualifiers={}):
    #print 'CSearchPathListView'
    qualis = {
               'columnHeaders':['Executable','Search path'],
               'tableItems' : ['name','path'],
               'editorClassName' : 'CSearchPathView'
               }
    qualis.update(qualifiers)
    CListView.__init__(self,parent,model=model,qualifiers=qualis)




class CSearchPathView(CComplexLineWidget):
  MODEL_CLASS = CCP4File.CSearchPath

  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = { 'vboxLayout' : True }
    qualis.update(qualifiers)
    CComplexLineWidget.__init__(self,parent,qualis)
    self.widgets = {}
    line = QtWidgets.QHBoxLayout()
    iconWidgetItem = self.layout().takeAt(0)
    line.addWidget(iconWidgetItem.widget())
    line.addWidget(QtWidgets.QLabel('Executable name',self))
    self.widgets['name'] = CLineEdit(parent=self)
    line.addWidget(self.widgets['name'])
    #print 'CSearchPathView.init',repr(self.layout())
    self.layout().addLayout(line)
    line = QtWidgets.QHBoxLayout()
    qualis = { 'jobCombo' : False }
    self.widgets['path'] =  CDataFileView(parent=self,qualifiers=qualis)
    line.addWidget(self.widgets['path'])
    self.layout().addLayout(line)

class CFramelessWarning(QtWidgets.QDialog):
  def __init__(self,parent=None,message=None,position=None):
    QtWidgets.QDialog.__init__(self,parent=parent)
    #print 'CFramelessWarning',message
    self.setWindowFlags(QtCore.Qt.FramelessWindowHint | QtCore.Qt.Dialog)
    self.setFocusPolicy(QtCore.Qt.NoFocus)
    self.setLayout(QtWidgets.QGridLayout())
    self.layout().setSpacing(3)
    self.layout().setContentsMargins(3,3,3,3)
    self.label = QtWidgets.QLabel(self)
    self.layout().addWidget(self.label,1,0)
    self.label.setText(message)
    self.setGeometry(position.x(),position.y(),150,10)
    self.show()
    self.timer = QtCore.QTimer()
    self.timer.setSingleShot(True)
    self.timer.timeout.connect(self.close)
    self.timer.start(3000)

class CJobSelectionCombo(QtWidgets.QComboBox):
  def __init__(self,parent,projectId=None,ifOneJob=False):
    QtWidgets.QComboBox.__init__(self,parent=parent)
    self.projectId = None
    if projectId is not None: self.setProjectId(projectId)

  def setProjectId(self,projectId=None):
    #print 'CJobSelectionCombo.setProjectId',projectId
    if projectId is None:
      projectId = self.projectId
      if projectId is None: return
    elif projectId == self.projectId:
      return
    else:
      self.projectId = projectId
    self.clear()

    jobInfoList = PROJECTSMANAGER().db().getProjectJobListInfo(projectId=projectId,topLevelOnly=True,order='ASC',mode=['jobnumber','taskname','jobid'])
    #print 'CJobSelectionCombo.setProjectId jobInfoList',jobInfoList
    for jobInfo in jobInfoList:
      self.addItem(jobInfo['jobnumber']+' '+jobInfo['taskname'],jobInfo['jobid'])

    #print 'CJobSelectionCombo.setProjectId',self.rootModelIndex(),
    #print 'child',self.rootModelIndex().child(0,0),
    #print 'model',self.rootModelIndex().child(0,0).model()

  def getSelection(self):
    qvar = self.itemData(self.currentIndex())
    jobId=str(qvar)
    #print 'CJobSelectionCombo.getSelection',jobId
    return jobId



class CJobSelectionLineEdit(QtWidgets.QLineEdit):
    def __init__(self,parent,projectId=None,ifOneJob=False):
      QtWidgets.QLineEdit.__init__(self,parent)
      self.projectId = projectId
      self.ifOneJob = ifOneJob
      if ifOneJob:
        self.setToolTip("Enter job number")
      else:
        self.setToolTip("Enter list of jobs e.g. '27-29,31'")

    def setProjectId(self,projectId):
      self.projectId = projectId

    def getSelection(self):
      errList = []
      seleList = []
      text = str(self.text())
      splitList = text.split(',')
      #print 'CExportJobSelection.getSelection',text,type(text)
      for item in splitList:
        #print 'CExportJobSelection.getSelection split',item
        rSplit = item.split('-')
        if len(rSplit)==1:
          try:
            jobId = PROJECTSMANAGER().db().getJobId(projectId=self.projectId,jobNumber=item.strip())
          except:
            errList.append(item)
          else:
            seleList.append(jobId)
        elif len(rSplit)==2:
          try:
            jobList = PROJECTSMANAGER().db().getJobsInRange(projectId=self.projectId,jobNumberRange=[rSplit[0].strip(),rSplit[1].strip()])
          except:
            errList.append(item)
          else:
            if len(jobList)==0:
              errList.append(item)
            else:
              seleList.extend(jobList)
        else:
           errList.append(item)
      #print 'CExportJobSelection.getSelection', seleList,errList
      return seleList,errList

class CGenericDataView(CViewWidget):

  def __init__(self,parent=None,qualifiers={},**kw):
    CViewWidget.__init__(self,parent=parent,qualifiers=qualifiers,**kw)
    self.setLayout(QtWidgets.QVBoxLayout())
    self.widgets={}
    line = QtWidgets.QHBoxLayout()
    self.widgets['dataType'] = QtWidgets.QComboBox(self)


class CEditFileLabel(QtWidgets.QDialog):
  # Allow user to edit the job title or file annotation
  # Usually called from CProjectWidget
  def __init__(self,parent=None,fileId=None,jobId=None,position=None,fileLabel=None):
    QtWidgets.QDialog.__init__(self,parent=parent)
    self.fileId = fileId
    self.jobId = jobId
    #self.edited = False
    self.setLayout(QtWidgets.QVBoxLayout())
    if fileId is not None:
      fileInfo = PROJECTSMANAGER().db().getFileInfo(fileId=fileId,mode=['annotation','jobid','fileclass'])
      from core import CCP4DataManager
      fileClass = CCP4DataManager.DATAMANAGER().getClass(fileInfo.get('fileclass',''))
      if fileClass is not None:
        fileTypeLabel = fileClass.QUALIFIERS['guiLabel']
      else:
        fileTypeLabel = 'file'
      jobInfo =  PROJECTSMANAGER().db().getJobInfo(jobId=fileInfo.get('jobid',None),mode=['jobnumber','jobtitle'])
      self.layout().addWidget(QtWidgets.QLabel('Enter label for '+fileTypeLabel+' from job '+ \
                                          str(jobInfo.get('jobnumber',''))+' '+ str(jobInfo.get('jobtitle','')),self ))
      label = fileInfo.get('annotation','')

    else:
      jobInfo = PROJECTSMANAGER().db().getJobInfo(jobId=jobId,mode=['jobtitle','jobnumber'])
      self.layout().addWidget(QtWidgets.QLabel('Enter label for job number '+str(jobInfo.get('jobnumber','')),self))
      label = jobInfo.get('jobtitle','')
    self.edit = QtWidgets.QLineEdit(self)
    self.edit.setMinimumWidth(300)
    self.layout().addWidget(self.edit)
    bb = QtWidgets.QDialogButtonBox(self)
    self.layout().addWidget(bb)
    b = bb.addButton(QtWidgets.QDialogButtonBox.Save)
    b.clicked.connect(self.save)
    b = bb.addButton(QtWidgets.QDialogButtonBox.Cancel)
    b.clicked.connect(self.reject)
    if label is not None:
      self.edit.setText(label)
      self.oldLabel = label
    else:
      self.oldLabel = ''
    self.show()

  #def handleEdit(self,text):
  #  self.edited = True

  @QtCore.Slot()
  def save(self):
    label = str(self.edit.text()).strip()
    label = re.sub('\n',' ',label)
    if label != self.oldLabel:
      #print 'CEditFileLabel.save',label
      if self.fileId is not None:
        PROJECTSMANAGER().db().updateFile(fileId=self.fileId,key='annotation',value=label)
      else:
        PROJECTSMANAGER().db().updateJob(jobId=self.jobId,key='jobtitle',value=label)
    self.accept()

class SearchBoxCompleter(QtCore.QObject):
  def __init__(self,parent=None):
    QtCore.QObject.__init__(self,parent)

    self.editor = parent
    self.popup = QtWidgets.QTreeWidget()
    self.popup.setWindowFlags(QtCore.Qt.Popup);
    self.popup.setFocusPolicy(QtCore.Qt.NoFocus);
    self.popup.setFocusProxy(parent);
    self.popup.setMouseTracking(True);

    self.popup.setColumnCount(1);
    self.popup.setUniformRowHeights(True);
    self.popup.setRootIsDecorated(False);
    self.popup.setEditTriggers(QtWidgets.QTreeWidget.NoEditTriggers);
    self.popup.setSelectionBehavior(QtWidgets.QTreeWidget.SelectRows);
    self.popup.setFrameStyle(QtWidgets.QFrame.Box | QtWidgets.QFrame.Plain);
    self.popup.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff);
    self.popup.header().hide();

    self.popup.installEventFilter(self);

    self.popup.itemClicked.connect(self.doneCompletion)
    self.hits = []

  @QtCore.Slot(bool)
  def doneCompletion(self,clicked=False):
    self.popup.hide()
    if not clicked and self.editor.text() in self.hits:
      #print "doneCompletion (editor)",self.editor.text()
      QtCore.QMetaObject.invokeMethod(self.editor,"returnPressed")
      self.editor.dataSelected.emit()
    else:
      item = self.popup.currentItem()
      if item:
        #print "doneCompletion",item.text(0)
        self.editor.setText(item.text(0))
        self.editor.changed.emit()
        self.popup.hide()
        QtCore.QMetaObject.invokeMethod(self.editor,"returnPressed")
        self.editor.dataSelected.emit()

  def eventFilter(self,obj, ev):
    if obj != self.popup:
      return False

    if ev.type() == QtCore.QEvent.MouseButtonPress:
        self.popup.hide()
        self.editor.setFocus()
        return True

    if ev.type() == QtCore.QEvent.KeyPress:

        consumed = False;
        key = ev.key();
        if key == QtCore.Qt.Key_Enter or key == QtCore.Qt.Key_Return:
            self.doneCompletion(True)
            consumed = True

        elif key == QtCore.Qt.Key_Escape:
            self.editor.setFocus()
            popup.hide()
            consumed = True

        elif key == QtCore.Qt.Key_Up:
            return False
        elif key == QtCore.Qt.Key_Down:
            return False
        elif key == QtCore.Qt.Key_Home:
            return False
        elif key == QtCore.Qt.Key_End:
            return False
        elif key == QtCore.Qt.Key_PageUp:
            return False
        elif key == QtCore.Qt.Key_PageDown:
            return False

        else:
            self.editor.setFocus()
            self.editor.event(ev)
            #self.popup.hide()
            return True

        return consumed

    return False

  def showCompletion(self, choices):

    if len(choices)==0:
        return

    pal = self.editor.palette();
    color = pal.color(QtGui.QPalette.Disabled, QtGui.QPalette.WindowText);

    self.popup.setUpdatesEnabled(False);
    self.popup.clear();
    for i in range(len(choices)):
        item = QtWidgets.QTreeWidgetItem();
        item.setText(0, choices[i]);
        self.popup.insertTopLevelItem(i,item)

    self.popup.resizeColumnToContents(0);
    self.popup.adjustSize();
    self.popup.setUpdatesEnabled(True);

    h = self.popup.sizeHintForRow(0) * min(7, len(choices)) + 3;
    self.popup.resize(self.editor.width(), h);

    self.popup.move(self.editor.mapToGlobal(QtCore.QPoint(0, self.editor.height())));
    self.popup.setFocus();
    self.popup.show();
    self.popup.setCurrentItem(self.popup.topLevelItem(0))

class CSearchBox(CLineEdit):

  changed = QtCore.Signal()
  dataSelected = QtCore.Signal()

  def __init__(self,parent=None):
    CLineEdit.__init__(self,parent)
    self.textChanged.connect(self.doSearch)
    self.completer = SearchBoxCompleter(self)
    self.textChanged.connect(self.doSearch)

  def setHitList(self,hitList):
    self.completer.hits = []
    self.completer.hits.extend(hitList)

  @QtCore.Slot(str)
  def doSearch(self,text):
    #print 'doSearch',text
    self.changed.emit()
    self.completer.popup.hide()
    if len(text)<1:
      return
    choicesl = []
    for ch in self.completer.hits:
      if str(text).lower() in ch.lower() or str(text).lower().replace(' ','-') in ch.lower():
        choicesl.append(ch)
    choices = choicesl
    self.completer.showCompletion(choices)

class CTreeItem:
  def __init__(self,parent=None,data={},children=[]):
    #if len(data)>0:
    #  print 'CTreeItem.__init__',self,parent,str(data[QtCore.Qt.DisplayRole])
    self.parent = parent
    self.myData = {}
    if isinstance(data,dict):
      self.myData.update(data)
    else:
      self.myData[QtCore.Qt.DisplayRole] = data
    self.children = []
    for c in children:
      if isinstance(c,(list,tuple)):
        self.appendChild(CTreeItem(self,c[0],c[1]))
      else:
        self.appendChild(CTreeItem(self,c))

  def childCount(self):
    return len(self.children)

  def child(self,row):
    if row>=0 and row<len(self.children):
      return self.children[row]
    else:
      return None

  def appendChild(self,item):
     self.children.append(item)

  def row(self):
    if self.parent is not None:
      return self.parent.children.index(self)
    else:
      return 0

  def data(self,role):
#FIXME PYQT - or maybe None? This used to return QVariant.
    return self.myData.get(role,"")

  def findChild(self,role,value):
    for child in self.children:
      if role in child.myData and child.myData[role] == value:
        return child
    return None

class CTasksModel(QtCore.QAbstractItemModel):

  def __init__(self,parent):
     QtCore.QAbstractItemModel.__init__(self,parent)
     from core import CCP4TaskManager
     modelAsList = CCP4TaskManager.TASKMANAGER().taskTree(shortTitles=True)
     model = []
     for module in modelAsList:
        taskList = []
        #print 'CTasksModel',module
        for name,title in module[2]:
          taskList.append( { QtCore.Qt.DisplayRole: title, QtCore.Qt.UserRole : name } )
        model.append( [ { QtCore.Qt.DisplayRole: module[1], QtCore.Qt.UserRole : module[0] }, taskList ] )
     self.rootItem = CTreeItem(self,'root',model)

  def columnCount(self,parent):
    return 1

  def childCount(self):
    return self.rootItem.childCount()

  def rowCount(self,index):
    if not index.isValid():
      return self.childCount()
    else:
      return index.internalPointer().childCount()

  def child(self,row):
    return self.rootItem.child(row)

  def data(self, index, role):
    if not index.isValid():
      return None
    item = index.internalPointer()
    #return item.data(index.column(),role)
    return item.data(role)

  def index(self, row, column, parent):
    if row < 0 or column < 0 or row >= self.rowCount(parent) or column >= self.columnCount(parent):
      return QtCore.QModelIndex()
    if not parent.isValid():
      parentItem = self.rootItem
    else:
      parentItem = parent.internalPointer()
    childItem = parentItem.child(row)
    #print 'CProjectModel.index',row, column, parent,childItem.getName()
    if childItem:
      return self.createIndex(row, column, childItem)
    else:
      return QtCore.QModelIndex()

  def parent(self,index):
    if not index.isValid():
      return QtCore.QModelIndex()
    childItem = index.internalPointer()
    parentItem = childItem.parent
    if parentItem == self.rootItem:
      return QtCore.QModelIndex()
    return self.createIndex(parentItem.row(), 0, parentItem)

class CTreeComboBox(QtWidgets.QComboBox):
  # Represent a tree in a combobox?
  #from  http://qt.shoutwiki.com/wiki/Implementing_QTreeView_in_QComboBox_using_Qt-_Part_2

  def __init__(self,parent):
    QtWidgets.QComboBox.__init__(self,parent)
    self.skipNextHide = False
    self.setView(QtWidgets.QTreeView(self))
    self.view().viewport().installEventFilter(self)
    QtWidgets.QComboBox.resize(self,200,30)


  def eventFilter(self,object,event):
    if event.type() == QtCore.QEvent.MouseButtonPress and object == self.view().viewport():
        mouseEvent = QtGui.QMouseEvent(event)
        index = self.view().indexAt(mouseEvent.pos())
        #print 'eventFilter',self.view().visualRect(index).contains(mouseEvent.pos())
        if  not self.view().visualRect(index).contains(mouseEvent.pos()):
           self.skipNextHide = True
    return False

  def showPopup(self):
    #self.setRootModelIndex(self.model().rootItem.index())
    QtWidgets.QComboBox.showPopup(self)

  def hidePopup(self):
    self.setRootModelIndex(self.view().currentIndex().parent())
    self.setCurrentIndex(self.view().currentIndex().row())
    if self.skipNextHide:
      self.skipNextHide = False
    else:
      QtWidgets.QComboBox.hidePopup(self)
