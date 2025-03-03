from __future__ import print_function


"""
     qtgui/CCP4AnnotationWidgets.py: CCP4 Gui Project
     Copyright (C) 2016 STFC

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
from PySide2 import QtCore,QtGui, QtWidgets
from core.CCP4Modules import *
from core.CCP4ErrorHandling import *
from core import CCP4Annotation
from qtgui import CCP4Widgets


class CAnnotationView(CCP4Widgets.CComplexLineWidget):

    textChanged = QtCore.Signal()

    MODEL_CLASS = CCP4Annotation.CAnnotation
    STRETCH = 5

    def __init__(self,parent=None,model=None,qualifiers={}):
        qualis = {'charWidth': -1}
        qualis.update(qualifiers)
        ifMultiLine = qualis.get('multiLine',False)
        if ifMultiLine: qualis['gridLayout'] = True
        CCP4Widgets.CComplexLineWidget.__init__(self,parent,qualis)
        if self.editable:
            if ifMultiLine:
                self.lineEdit = CCP4Widgets.CTextEdit(self)
                self.layout().addWidget(self.lineEdit,1,0,1,2)
                if qualis.get('title',None) is not None:
                    self.layout().addWidget(CCP4Widgets.CItalicLabel(qualis['title'],self),0,1)
                self.lineEdit.textChanged.connect(self.updateModelFromView)
            else:
                self.lineEdit = CCP4Widgets.CLineEdit(self,{'charWidth':qualis['charWidth'] } )
                self.layout().addWidget(self.lineEdit)
                self.layout().setStretchFactor(self.lineEdit,5)
                self.lineEdit.editingFinished.connect(self.updateModelFromView)
        else:
            if ifMultiLine:
                if qualis.get('title',None) is not None:
                    self.layout().addWidget(CCP4Widgets.CItalicLabel(qualis['title'],self),0,1)
                self.lineEdit = CCP4Widgets.CTextEdit(self) # KJS : Another bug. Seems to be multiple problems with this code. Be surprised if it's actually used anywhere.
                self.lineEdit.setReadOnly(True)
                self.layout().addWidget(self.lineEdit,1,0,1,2)
            else:
                self.lineEdit = CCP4Widgets.CLabel(self,{'charWidth':qualis['charWidth'] })
                self.layout().addWidget(self.lineEdit)
                self.layout().setStretchFactor(self.lineEdit,5)
        if model is not None:
            self.setModel(model)

    def sizePolicy(self):
        #print 'using CAnnotationView.sizePolicy'
        return QtWidgets.QSizePolicy.MinimumExpanding

    def handleEdit(self):
        self.textChanged.emit()

    def getMenuDef(self):
        return ['copy','paste','help']

    def updateViewFromModel(self):
        if self.model is None:
            return
        text = self.model.text.get()
        #print 'CAnnotationView.updateViewFromModel',text
        self.lineEdit.blockSignals(True)
        self.lineEdit.setValue(text)
        self.lineEdit.blockSignals(False)

    @QtCore.Slot()
    def updateModelFromView(self):
        if self.model is None:
            return
        value = {'text' : self.lineEdit.getValue()}
        #print 'CAnnotationView.updateModelFromView',value
        self.connectUpdateViewFromModel(False)
        self.model.set(value=value)
        self.connectUpdateViewFromModel(True)

    def updateModelFromText(self):
        self.updateModelFromView()

class CTimeView(CCP4Widgets.CViewWidget):
    MODEL_CLASS = CCP4Annotation.CTime

    def __init__(self,parent=None,model=None,qualifiers={}):
        qualis = qualifiers
        CCP4Widgets.CViewWidget.__init__(self,parent,qualifiers=qualis)
        layout = QtWidgets.QHBoxLayout()
        layout.setSpacing(CCP4Widgets.CViewWidget.MARGIN) # KJS : This one seems to have a problem. Investigate.
        layout.setContentsMargins(CCP4Widgets.CViewWidget.MARGIN, CCP4Widgets.CViewWidget.MARGIN, CCP4Widgets.CViewWidget.MARGIN, CCP4Widgets.CViewWidget.MARGIN)
        self.setLayout(layout)
        if self.editable:
            self.widget = QtWidgets.QDateEdit(self)
            self.widget.dateTimeChanged.connect(self.updateModelFromView1)
        else:
            self.widget = CCP4Widgets.CLabel(self)
        self.layout().addWidget(self.widget)
        #print 'CTimeView.__init__',self.layout().count()
        self.setModel(model)
    
    def setValue(self,value=None):
        t = QtCore.QDateTime()
        t.setTime_t(value)
        if self.editable:
            self.widget.setDateTime(t)
        else:
            self.widget.setValue(str(t))

    def getValue(self):
        return  self.widget.dateTime().toTime_t()

    def updateViewFromModel(self):
        self.widget.blockSignals(True)
        self.widget.setDateTime(self.model.QDateTime)
        self.widget.blockSignals(False)

    @QtCore.Slot('QDateTime')
    def updateModelFromView1(self,date):
        self.updateModelFromView()

    def updateModelFromView(self):
        self.connectUpdateViewFromModel(False)
        self.model.QDateTime = self.widget.dateTime()
        self.connectUpdateViewFromModel(True)


class CMetaDataTagView(CCP4Widgets.CComplexLineWidget):
    MODEL_CLASS = CCP4Annotation.CMetaDataTag

    def __init__(self,parent=None,model=None,qualifiers={}):
        CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualifiers)
        self.widgets['tag'] = CCP4Widgets.CComboBox(self)
        self.layout().addWidget(self.widgets['tag'])
        self.setModel(model)

    def setModel(self,model):
        #print 'CMetaDataTagView.setModel',model
        #try:
        #  print 'CMetaDataTagView.setModel',model.parent().__dict__['_value'],repr(model)
        #except:
        #  pass
        # These only need connecting on first call to setModel() - if this is a list editor then
        # likely to be called multiple times
        if self.model is not None:
            self.model.enumeratorsUpdated.disconnect(self.updateCombo)
        if self.model is None and model is not None:
            self.widgets['tag'].dataChanged.connect(self.updateTag)
            self.widgets['tag'].returnPressed.connect(self.handleReturn)
            self.widgets['tag'].focusOut.connect(self.handleReturn)
        CComplexLineWidget.setModel(self,model)   # KJS : Another one...
        if self.model is not None:
            self.updateCombo()
            self.model.enumeratorsUpdated.connect(self.updateCombo)

    @QtCore.Slot()
    def handleReturn(self):
        #print 'CMetaDataTagView.handleReturn',self.widgets['tag'].currentText()
        # addEnumerator() creates new tag if necessary and return index of the tag
        idx = self.model.addEnumerator(str(self.widgets['tag'].currentText()))
        from qtgui import CCP4I1Projects
        #print 'CMetaDataTagView.handleReturn',CCP4I1Projects.CI1PREFERENCES().tagList
        #print 'CMetaDataTagView.handleReturn idx',idx
        if idx >= 0:
            self.updateCombo()
            self.widgets['tag'].setCurrentIndex(idx)
        self.updateTag()

    @QtCore.Slot()
    def updateTag(self):
        #print 'updateTag ',self.widgets['tag'].getValue()
        self.model.tag.set( self.widgets['tag'].getValue(),checkValidity=False)

    @QtCore.Slot()
    def updateCombo(self):
        self.widgets['tag'].populate(self.model.getEnumerators())


class CMetaDataTagListView(CCP4Widgets.CListView):

    MODEL_CLASS = CCP4Annotation.CMetaDataTagList

    def __init__(self,parent=None,model=None,qualifiers={},editorQualifiers={}):
        qualis = { 'mode' : 'list' }
        qualis.update(qualifiers)
        CCP4Widgets.CListView.__init__(self,parent=parent,model=model,qualifiers=qualis,editorQualifiers=editorQualifiers)
        self.editor.dataChanged.connect(self.populateListWidget)

    @QtCore.Slot()
    def updateCombo(self):
        self.editor.widget.updateCombo()

   
class CDateRangeView(CCP4Widgets.CComplexLineWidget):

    MODEL_CLASS = CCP4Annotation.CDateRange

    def __init__(self,parent=None,model=None,qualifiers={}):
        qualis = {}
        qualis.update(qualifiers)
        CCP4Widgets.CComplexLineWidget.__init__(self,parent=parent,qualifiers=qualis)
        self.layout().setSpacing(1)
        for item in [ 'year', 'month', 'day']:
            self.widgets[item] = CCP4Widgets.CComboBox(parent=self,qualifiers=qualis)
            self.layout().addWidget(self.widgets[item])
        self.layout().addWidget(QtWidgets.QLabel('+/-',self))
        self.widgets['year'].setMaximumWidth(QtGui.QFontMetrics(self.widgets['year'].font()).width('2222')+30)
        self.widgets['month'].setMaximumWidth(QtGui.QFontMetrics(self.widgets['month'].font()).width('September')+30)
        self.widgets['day'].setMaximumWidth(QtGui.QFontMetrics(self.widgets['day'].font()).width('28')+30)
        self.widgets['month'].currentIndexChanged[int].connect(self.updateDayCombo)
        for item in [ 'yearRange', 'monthRange', 'dayRange']:
            self.widgets[item] = CCP4Widgets.CIntView(parent=self,qualifiers=qualis)
            self.widgets[item].setMaximumWidth(QtGui.QFontMetrics(self.widgets[item].font()).width('9999'))
            self.layout().addWidget(self.widgets[item])
            self.layout().addWidget(QtWidgets.QLabel(item[0:-5]+'s',self))
        self.setModel(model)

    def setModel(self,model=None):
        CCP4Widgets.CComplexLineWidget.setModel(self,model=model)
        if model is None: return
        self.widgets['year'].populate(model.year.qualifiers('enumerators'))
        self.widgets['month'].populate(['']+model.month.qualifiers('enumerators'))
        self.updateDayCombo(self.widgets['month'].currentIndex())
        self.widgets['day'].populate([' ']+model.day.qualifiers('enumerators'))
        self.updateViewFromModel()

    @QtCore.Slot(int)
    def updateDayCombo(self,monthIndex):
        if monthIndex<1:
            return
        currentIndex = self.widgets['day'].currentIndex()
        dayList = list(range(1,self.model.MONTHLENGTH[monthIndex]+1))
        self.model.day.setQualifier('enumerators', dayList)
        menu = ['']
        if monthIndex > 0:
            for item in range(1, self.model.MONTHLENGTH[monthIndex]+1):
                menu.append(item)
        self.widgets['day'].populate(menu)
        self.widgets['day'].setCurrentIndex(min(currentIndex,self.model.MONTHLENGTH[monthIndex]))

    #Reimplement getValue()/setValue() to convert a data value of None to zero length string

    def setValue(self,value=None):
        #print 'CDateRangeView.setValue',value
        if self.blockUpdateView: return
        if value is None: return
        for item,w in  list(self.widgets.items()):
            w.blockSignals(True)
            if value[item] is None:
                w.setValue('')
            else:
                w.setValue(value[item])
            w.blockSignals(False)
        self.lastWarning = None

    def getValue(self):
        data = {}
        for item,w in list(self.widgets.items()):
            w = self.widgets.get(item,None)
            data[item] = w.getValue()
            if len( data[item]) == 0:
                data[item]=None
        #print 'CDateRangeView.getValue',data
        return data
