from __future__ import print_function

"""
     qtgui/guiUtils.py: CCP4 Gui Project
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009 University of York

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

##@package CCP4GuiUtils (QtGui) Assorted Gui utilities
import os
import types
import functools
import sys
if sys.version_info >= (3,7):
    from collections.abc import Callable
else:
    from collections import Callable
from PySide2 import QtGui, QtWidgets,QtCore,QtSvg
from qtgui import CCP4Widgets

LOADED_PIXMAPS = {}
ICON_EXTENSIONS= ['.svg','.png']

def checkState(inp):
    if inp:
        return QtCore.Qt.Checked
    else:
        return QtCore.Qt.Unchecked

def fromCheckState(inp):
    if inp == QtCore.Qt.Unchecked:
        return 0
    else:
        return 1

#----------------------------------------------------------------------
def fromMenu(text,menu_list,alias_list):
#----------------------------------------------------------------------
    # return the short alias for a given long menu text
    if menu_list.count(text):
        return alias_list[menu_list.index(text)]
    else:
        return text


#----------------------------------------------------------------------
def createAction(name='',parent=None,definition={},icon_path='',default_icon=''):
#----------------------------------------------------------------------
    if definition:
        adef = definition
    else:
        adef = parent.getActionDef(name)
    #print "createAction",name,adef
    # Is it a 'build-in' action
    if 'action' in adef:
        #print "createAction build-in",name,adef['action']
        return adef['action']
    text = ''
    if 'text' in adef:
        if isinstance(adef['text'], Callable):
            text = adef['text']()
        else:
            text = adef['text']
    # Is it grouped with other actions
    if 'group' in adef:
        group=parent.findChild(QtWidgets.QActionGroup,adef['group'])
        if not group:
            group = QtWidgets.QActionGroup(parent)
            group.setObjectName(adef['group'])
        a = group.addAction(text)
    else:
        a = QtWidgets.QAction(text,parent)
    a.setObjectName(name)
    if 'shortcut' in adef and adef['shortcut']:
        if isinstance(adef['shortcut'],int):
            a.setShortcut(adef['shortcut'])
        else:
            a.setShortcut(QtGui.QKeySequence(adef['shortcut']))
    icon = createIcon(name,adef,icon_path=icon_path,default_icon=default_icon)
    if icon:
        a.setIcon(icon)
    tip = adef.get('toolTip','')
    if tip:
        a.setToolTip(tip)
    if 'checkable' in adef and adef['checkable']:
        a.setCheckable(1)
        if 'checked' in adef and isinstance(adef["checked"], Callable):
            res = adef["checked"]()
            if type(res)==int or type(res)==bool:
                a.setChecked(res)
    if 'enabled' in adef:
        if isinstance(adef["enabled"], Callable):
            res = adef["enabled"]()
        else:
            res = adef["enabled"]
        if type(res)==int or type(res)==bool:
            a.setEnabled(res)
    if adef.get('deleteLater',0):
        a.deleteLater()
    # Sort out a signal
    if 'signal' in adef:
        signal = adef['signal']
    else:
        signal = 'triggered()'
#FIXME!!!
    if adef.get('slot',''):
        if isinstance(adef['slot'],list) and len(adef['slot'])>1:
            a.triggered.connect(functools.partial(adef['slot'][0],adef['slot'][1]))
        else:
            a.triggered.connect(adef['slot'])
    return a

#----------------------------------------------------------------------
def createIcon(name=None,adef={},icon_path='',default_icon='unknown'):
#----------------------------------------------------------------------
    #print 'guiUtils.createIcon',name
    if 'icon' in adef and adef['icon']:
        icon_name =  str(adef['icon'])
    elif name is not None:
        icon_name = str(name)
    else:
        icon_name = 'unknown'
    if not icon_path:
        for item in ['CCP4I2_TOP']:
            if item in os.environ:
                icon_path = os.path.join(os.path.abspath(os.environ[item]),'qticons')
    #print 'guiUtils.createIcon',icon_path,icon_name
    for ext in ICON_EXTENSIONS:
        filename = os.path.join(icon_path,icon_name+ext)
        if os.path.exists(filename):
            #print 'guiUtils.createIcon',filename
            if ext == '.svg':
                pix = loadSvg(filename)
            else:
                pix =  QtGui.QPixmap(filename)
            ico =  QtGui.QIcon(pix)
            return ico
    if default_icon:
        for item in ['CCP4I2_TOP']:
            if item in os.environ:
                file = os.path.join(os.path.abspath(os.environ[item]),'qticons',default_icon+'.png')
                return QtGui.QIcon(file)
    else:
        return None

def loadSvgWithSize(fileName,width,height):
    svg = QtSvg.QSvgRenderer()
    svg.load(fileName)
    pixmap = QtGui.QImage(width,height,QtGui.QImage.Format_ARGB32)
    pixmap.fill(QtGui.QColor(0, 0, 0, 0))
    painter = QtGui.QPainter(pixmap)
    painter.setRenderHints(QtGui.QPainter.Antialiasing);
    svg.render(painter)
    painter.end()
    #pixmap = pixmap.scaledToHeight(height,QtCore.Qt.SmoothTransformation)
    return QtGui.QPixmap.fromImage(pixmap)

def loadSvg(fileName):
    ICONSIZE = 24
    svg = QtSvg.QSvgRenderer()
    svg.load(fileName)
    pixmap = QtGui.QPixmap(ICONSIZE,ICONSIZE)
    pixmap.fill(QtGui.QColor(0,0,0,0))
    painter = QtGui.QPainter(pixmap)
    svg.render(painter)
    painter.end()
    return pixmap

#----------------------------------------------------------------------
def clearMenu(menu=None):
#----------------------------------------------------------------------
    submenus = menu.findChildren(QtWidgets.QMenu)
    menu.clear()
    for sm in submenus:
        sm.clear()

#----------------------------------------------------------------------
def populateMenu(parent=None,menuWidget=None,definition=[],
                 getActionDef=None,default_icon='',info={}):
#----------------------------------------------------------------------
    #print "populateMenu",definition
    #menuWidget.hovered.connect(handleMenuHover)
#FIXME!!!
    if not definition:
        menuWidget.addAction(createAction('dummy',parent,{ 'text' : '  --  '} ,default_icon=default_icon))
    elif not isinstance(definition,list):
        menuWidget.aboutToShow.connect(functools.partial(definition,menuWidget))
    else:
        for item in definition:
            a = None
            if item == 'sep':
                #print 'Not adding a menu separator',definition
                menuWidget.addSeparator()
            elif isinstance(item,list) and len(item)>0:         
                sub_menu=menuWidget.addMenu(item[0])
                if len(item) == 1:
                    # This is a sub without actions defined
                    sub_menu.setEnabled(0)
                elif isinstance(item[1],str):
                    populateMenu(parent=parent,menuWidget=sub_menu,definition = item[1:],getActionDef=getActionDef,default_icon=default_icon,info=info)
                elif isinstance(item[1],types.MethodType):
                    sub_menu.aboutToShow.connect(functools.partial(item[1],(sub_menu,item[0])))
                else:
                    for item2 in item[1:]:
                       sub_menu2=sub_menu.addMenu(item2[0])
                       if len(item2) == 1:
                           sub_menu2.setEnabled(0)
                       else:
                           populateMenu(parent=parent,menuWidget=sub_menu2,definition = item2[1:],getActionDef=getActionDef,default_icon=default_icon,info=info)
            else:
                if getActionDef:
                    adef = getActionDef(item)
                else:
                    adef = parent.getActionDef(item)
                #print 'populateMenu adef ',item,adef
                if adef:
                    a = adef.get('action',None)
                    if a is None:
                        a = createAction(item,parent,adef,default_icon=default_icon)
                        adef['action'] = a
                    else:
                        # The action exists but lets update the checked status
                        if 'text' in adef and isinstance(adef["text"], Callable):
                            a.setText(adef['text']())
                        if 'checked' in adef and isinstance(adef["checked"], Callable):
                            res = adef["checked"]()
                            #print "populateMenu checked",item,res
                            if type(res)==int or type(res)==bool:
                                a.setChecked(res)
                        if 'enabled' in adef:
                            #print 'reload adef enabled',item,adef["enabled"]
                            if isinstance(adef["enabled"], Callable):
                                res = adef["enabled"]()
                            else:
                                res = adef["enabled"]
                            #print 'reload adef enabled',item,res
                            if type(res)==int or type(res)==bool:
                                a.setEnabled(res)
                if a:
                    menuWidget.addAction(a)

def handleMenuHover(action):
    tip = action.toolTip()
    if not tip: return
    QtWidgets.QToolTip.showText(QtGui.QCursor.pos(),tip)

#----------------------------------------------------------------------
def populateToolBar(parent=None,toolBarWidget=None,definition=[],
                 getActionDef=None,default_icon='',info={}):
#----------------------------------------------------------------------
    #print 'populateToolBar',definition
    if not definition:
        toolBarWidget.addAction(createAction('dummy',parent,{ 'text' : '  --  '} ,default_icon=default_icon))
    elif not isinstance(definition,list):
        toolBarWidget.aboutToShow.connect(functools.partial(definition,toolBarWidget))   # KJS : Typos fixed
    else:
        for item in definition:
            a = None
            if item == 'sep':
                toolBarWidget.addSeparator()
            elif isinstance(item,list) and len(item)>0:         
                #  ?????????????  handling splitter
                pass
            else:
                if getActionDef:
                    adef = getActionDef(*[item], **info)
                else:
                    adef = parent.getActionDef(*[item], **info)
                if adef:
                    #print "createAction adef",adef
                    a = parent.findChild(QtWidgets.QAction,item)
                    #print "populateMenu find action",item,a
                    if not a:
                        a = createAction(item,parent,adef,default_icon=default_icon)
                    else:
                        # The action exists but lets update the checked status
                        if 'text' in adef and isinstance(adef["text"], Callable):
                            a.setText(adef['text']())
                        if 'checked' in adef and isinstance(adef["checked"], Callable):
                            res = adef["checked"]()
                            #print "populateMenu checked",item,res
                            if type(res)==int or type(res)==bool:
                                a.setChecked(res)
                        if 'enabled' in adef:
                            #print 'reload adef enabled',item,adef["enabled"]
                            if isinstance(adef["enabled"], Callable):
                                res = adef["enabled"]()
                            else:
                                res = adef["enabled"]
                            #print 'reload adef enabled',item,res
                            if type(res)==int or type(res)==bool:
                                a.setEnabled(res)
                if a:
                    toolBarWidget.addAction(a)

#-------------------------------------------------------------------
def loadPixmap(name='',group='actions',icon_path='',width=32,height=0):
#-------------------------------------------------------------------
    if not height: height = width
    size_key = str(width)+'_'+str(height)
    if not (name in LOADED_PIXMAPS and size_key in LOADED_PIXMAPS[name]):
        if not icon_path:
            icon_path = os.path.join(os.path.abspath(os.environ['CCP4I2_TOP']),'qticons')
        file = ''
        for ext in ICON_EXTENSIONS:
            f = os.path.join(icon_path,group,name+ext)   
            if os.path.exists(f):
                file = f
                break
        if not file:
            file = os.path.join(os.environ['CCP4I2_TOP'],'qticons','unknown.png')
            if  not os.path.exists(file):
                return None
        if name not in LOADED_PIXMAPS:
            LOADED_PIXMAPS[name] = {}
        LOADED_PIXMAPS[name][size_key] = QtGui.QPixmap(file).scaled(width,height)
    return LOADED_PIXMAPS[name][size_key]

#-------------------------------------------------------------------
def loadResource(self,url='',mode=QtGui.QTextDocument.StyleSheetResource):
#-------------------------------------------------------------------
    qu = QtCore.QUrl(url)
    self.style = self.textBrowser.loadResource(mode, qu)


#-------------------------------------------------------------------
def setWidgetValue(widget,value):
#-------------------------------------------------------------------
    #print 'setWidgetValue',widget,value
    if isinstance(widget,QtWidgets.QComboBox):
        ic = widget.findText(value)
        if ic < 0:
            ic = widget.findData(value)
        if ic >= 0:
            widget.setCurrentIndex(ic)
    elif isinstance(widget,QtWidgets.QLineEdit):
        widget.setText(str(value))
    elif isinstance(widget,QtWidgets.QLabel):
        widget.setText(str(value))
    elif isinstance(widget,QtWidgets.QAbstractButton):
        widget.setChecked(value)
    elif isinstance(widget,QtWidgets.QDoubleSpinBox):
        widget.setValue(float(value))
    elif isinstance(widget,QtWidgets.QSpinBox):
        widget.setValue(int(value))
    elif isinstance(widget,CCP4Widgets.CDataFileView):
        widget.setDataFile(value)
    '''
    elif isinstance(widget,mgSelectionCombo):
      widget.setText(value)
    elif isinstance(widget,mgColourCombo):
      widget.set_colour(value)
    elif isinstance(widget,mgAtomPicker):
      widget.setSelection(atomName=value)
    '''

#-------------------------------------------------------------------
def getWidgetValue(widget=None,default='',baseType='text'):
#-------------------------------------------------------------------
    value = default
    if isinstance(widget,QtWidgets.QComboBox):
        #value = str(widget.itemData(widget.currentIndex()))
        variant = widget.itemData(widget.currentIndex())
        if variant:
            if baseType == 'int':
                value = variant
            elif baseType == 'real':
                value = variant
            else:
                value = str(variant)
            #print 'getWidgetValue menu variant',variant,value
        if not value:
            value = str(widget.currentText())
    elif isinstance(widget,QtWidgets.QLabel):
        value = str(widget.text())
    elif isinstance(widget,QtWidgets.QLineEdit):
        value = str(widget.text())
    elif isinstance(widget,QtWidgets.QAbstractButton):
        value = int(widget.isChecked())
    elif isinstance(widget,QtWidgets.QDoubleSpinBox):
        value = float(widget.value())
    elif isinstance(widget,QtWidgets.QSpinBox):
        value = int(widget.value())
    elif isinstance(widget,CCP4Widgets.CDataFileView):
        value = widget.getDataFile()
    '''
    elif isinstance(widget,mgSelectionCombo):
        value = str(widget.text())
    elif isinstance(widget,mgColourCombo):
        value = str(widget.get_colour())
    elif isinstance(widget,mgAtomPicker):
        value = str(widget.selection())
    '''
    return value

def getWidgetSignal(widget=None,default='clicked()'):
    signal = default
    if isinstance(widget,QtWidgets.QComboBox):
        signal = 'currentIndexChanged(int)'
    elif isinstance(widget,QtWidgets.QLabel):
        signal = 'linkActivated(const QString)'
    elif isinstance(widget,QtWidgets.QLineEdit):
        signal = 'editingFinished()'
    elif isinstance(widget,QtWidgets.QAbstractButton):
        signal = 'clicked(bool)'
    elif isinstance(widget,QtWidgets.QDoubleSpinBox):
        signal = 'valueChanged(double)'
    elif isinstance(widget,QtWidgets.QSpinBox):
        signal = 'valueChanged(int)'
    return signal


