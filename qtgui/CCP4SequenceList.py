from __future__ import print_function

import sys
import functools
import textwrap

from PySide2 import QtCore, QtGui, QtWidgets

"""
This is a simple implementation of an edit for a sequence data model (name,no of copies, description, sequence).

The view is name, no. of copies, description, and a button to view/edit the sequence.

Double-click starts editing a cell. For the type this brings up a combo box with "PROTEIN","RNA","DNA".
All other cells use default editors.

The view button brings up a text view of the sequence with editing turned off. You can turn on editing if you wish.
There are Save/Don't save button.

"""

class SequenceModel(QtCore.QAbstractTableModel):

    def canDropMimeData(self, data, action, row, column, parent):
        return True

    def dropMimeData(self, data, action, row, column, parent):
        return True

    def setEditable(self,val):
        self.beginResetModel()
        self._editable = val
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.endResetModel()
        self.dataChanged.emit(tl,br)

    def __init__(self,parent=None):
        QtCore.QAbstractTableModel.__init__(self,parent)
        self._data = []
        self._editable = True

    def headerData(self,section, orientation, role):
        if role == QtCore.Qt.DisplayRole and orientation == QtCore.Qt.Horizontal:
            if section == 0:
                return "Name"
            elif section == 1:
                return "Copies"
            elif section == 2:
                return "Description"
            elif section == 3:
                return "Type"
            elif section == 4:
                return "Sequence"
        else:
            return QtCore.QAbstractTableModel.headerData(self,section, orientation, role)

    def removeRows(self,idxs):

        self.beginResetModel()
        indices = [ x.row() for x in idxs ]
        dataNew = []
        for i,row in zip(range(len(self._data)),self._data):
            if not i in indices:
                dataNew.append(row)
        self._data = dataNew
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.endResetModel()
        self.dataChanged.emit(tl,br)

    def getRow(self,i):
        return self._data[i]

    def rowCount(self,parent=None):
        return len(self._data)

    def columnCount(self,parent):
        return 5

#FIXME - These two methods could be combined by use of UserRole
    def setSequenceData(self,i,val):
        self.beginResetModel()
        self._data[i][4] = val
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.endResetModel()
        self.dataChanged.emit(tl,br)

    def setData(self,index,val,role):
        if role == QtCore.Qt.EditRole:
            self._data[index.row()][index.column()] = val
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.dataChanged.emit(tl,br)
        return QtCore.QAbstractTableModel.setData(self,index,val,role)

#FIXME - These two methods could be combined by use of UserRole
    def sequenceData(self,i):
        return self._data[i][4]

    def data(self,index,role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.ToolTipRole and index.column()>3:
            tip = self._data[index.row()][2]+"\n\n"+"\n".join(textwrap.wrap(self._data[index.row()][4][:1600],80))
            if len(self._data[index.row()][4])>1600:
                tip = tip + "\n..."
            return tip
        if index.column() > 3:
            return None
        if role == QtCore.Qt.DisplayRole or role == QtCore.Qt.EditRole:
            return self._data[index.row()][index.column()]
        return None

    def addItem(self,item):
        self.beginResetModel()
        self._data.append(list(item))
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.endResetModel()
        self.dataChanged.emit(tl,br)

    def flags(self,index):
        if not self._editable:
            return QtCore.QAbstractTableModel.flags(self,index) &~ QtCore.Qt.ItemIsEditable
        if not index.isValid():
            return QtCore.Qt.ItemIsEnabled


        return QtCore.QAbstractItemModel.flags(self,index) | QtCore.Qt.ItemIsEditable

class SquareButton(QtWidgets.QToolButton):

    def __init__(self,text,parent=None):
        QtWidgets.QToolButton.__init__(self,parent)
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        self.setText(text)

    def paintEvent(self,_paintEvent ):
        styleOption = QtWidgets.QStyleOptionToolButton()
        QtWidgets.QToolButton.initStyleOption(self,styleOption)
        styleOption.features &= ~QtWidgets.QStyleOptionToolButton.HasMenu

        painter = QtGui.QPainter(self)
        self.style().drawComplexControl( QtWidgets.QStyle.CC_ToolButton, styleOption, painter, self )

class ComboDelegate(QtWidgets.QStyledItemDelegate):

    def __init__(self,parent=None):
        QtWidgets.QStyledItemDelegate.__init__(self,parent)

    def createEditor(self,parent, option, index):
        cb = QtWidgets.QComboBox(parent)
        row = index.row()
        cb.addItem("PROTEIN")
        cb.addItem("RNA")
        cb.addItem("DNA")
        return cb

    def setEditorData(self,editor, index):
        if sys.version_info > (3,0):
            currentText = index.data(QtCore.Qt.EditRole)
        else:
            currentText = str(index.data(QtCore.Qt.EditRole))
        cbIndex = editor.findText(currentText)
        if (cbIndex >= 0):
           editor.setCurrentIndex(cbIndex)

    def setModelData(self,editor, model, index):
        model.setData(index, editor.currentText(), QtCore.Qt.EditRole)

class NoSpaceDelegate(QtWidgets.QStyledItemDelegate):

    def __init__(self,parent=None):
        QtWidgets.QStyledItemDelegate.__init__(self,parent)

    def createEditor(self,parent, option, index):
        editor = QtWidgets.QLineEdit(parent)
        validator = QtGui.QRegExpValidator(QtCore.QRegExp("\\S+"), self )
        editor.setValidator(validator)
        editor.setFrame(False)
        return editor

    def setEditorData(self,editor, index):
        value = index.data(QtCore.Qt.EditRole)
        editor.setText(value)

    def setModelData(self,editor, model, index):
        model.setData(index, editor.text(), QtCore.Qt.EditRole)

class SpinDelegate(QtWidgets.QStyledItemDelegate):

    def __init__(self,parent=None):
        QtWidgets.QStyledItemDelegate.__init__(self,parent)

    def createEditor(self,parent, option, index):
        editor = QtWidgets.QSpinBox(parent)
        editor.setFrame(False)
        editor.setMinimum(0)
        editor.setMaximum((2<<30)-1)
        return editor

    def setEditorData(self,editor, index):
        value = index.data(QtCore.Qt.EditRole)
        editor.setValue(value)

    def setModelData(self,editor, model, index):

        spinBox = editor
        spinBox.interpretText()
        value = spinBox.value()
        model.setData(index, value, QtCore.Qt.EditRole)

class SequenceTableView(QtWidgets.QTableView):

    showSignal = QtCore.Signal(int)
    deleteSignal = QtCore.Signal()
    urlsDroppedSignal = QtCore.Signal(list)
    textDroppedSignal = QtCore.Signal(str)

    def paintEvent(self,e):
        if self.model().rowCount() == 0:
            p = QtGui.QPainter(self.viewport())
            fm = QtGui.QFontMetrics(self.font())
            h = fm.height()
            r = self.rect()
            r.moveTop(-h)
            p.drawText(r, QtCore.Qt.AlignCenter, "There are currently no sequences in this asymmetric unit description.\nClick '+' below to add sequence(s).\n\nOr drag and drop sequence/coordinate files into here.")
        else:
            QtWidgets.QTableView.paintEvent(self,e)

    def keyPressEvent(self,e):
        if ((e.key()==QtCore.Qt.Key_Delete) or (e.key()==QtCore.Qt.Key_Backspace)) and len(self.selectionModel().selectedRows())>0:
            self.deleteSignal.emit()

    @QtCore.Slot(str)
    def showSequence(self,val):
        win = QtWidgets.QDialog()
        layout = QtWidgets.QVBoxLayout()
        tb = QtWidgets.QTextEdit()
        tb.setStyleSheet("QTextEdit{font-family:courier;}")
        tb.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse|QtCore.Qt.TextSelectableByKeyboard)
        layout.addWidget(tb)
        dbb = QtWidgets.QDialogButtonBox()
        layout.addWidget(dbb)
#This is to mimic the idea of clicking a pen icon to edit on a phone/tablet/web page.
        if self.model()._editable:
            closeButton = dbb.addButton(QtWidgets.QDialogButtonBox.Save)
            discardButton = dbb.addButton(QtWidgets.QDialogButtonBox.Discard)
            editButton = dbb.addButton("Edit",QtWidgets.QDialogButtonBox.ActionRole)
            self.editing = False
            @QtCore.Slot()
            def toggleEditing():
                self.editing = not self.editing
                if self.editing:
                    editButton.setText("Stop editing")
                    tb.setTextInteractionFlags(QtCore.Qt.TextEditorInteraction)
                else:
                    editButton.setText("Edit")
                    tb.setTextInteractionFlags(QtCore.Qt.TextBrowserInteraction)
            editButton.clicked.connect(toggleEditing)
            closeButton.clicked.connect(win.accept)
            discardButton.clicked.connect(win.reject)
        else:
            closeButton = dbb.addButton(QtWidgets.QDialogButtonBox.Close)
            closeButton.clicked.connect(win.accept)
        tb.setText(self.model().sequenceData(val))
        tb.setLineWrapMode(QtWidgets.QTextEdit.FixedColumnWidth)
        tb.setLineWrapColumnOrWidth(80)
        win.setLayout(layout)
        rv = win.exec_()
        if rv == QtWidgets.QDialog.Accepted and self.model()._editable:
            self.model().setSequenceData(val,"".join(str(tb.toPlainText()).split()))

    def __init__(self,parent=None):
        QtWidgets.QTableView.__init__(self,parent)
        self.showSignal.connect(self.showSequence)
        self.setEditTriggers(QtWidgets.QAbstractItemView.DoubleClicked|QtWidgets.QAbstractItemView.SelectedClicked)
        textDelegate = NoSpaceDelegate(self)
        self.setItemDelegateForColumn(0,textDelegate)
        spinDelegate = SpinDelegate(self)
        self.setItemDelegateForColumn(1,spinDelegate)
        comboDelegate = ComboDelegate(self)
        self.setItemDelegateForColumn(3,comboDelegate)
        self.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.viewport().setAcceptDrops(True)
        #self.horizontalHeader().setStretchLastSection(True)
        self.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)

    def dragMoveEvent(self, e):
        if e.mimeData().hasUrls() or e.mimeData().hasText():
            e.acceptProposedAction()

    def dragEnterEvent(self, e):
        if e.mimeData().hasUrls() or e.mimeData().hasText():
            e.acceptProposedAction()

    def dropEvent(self, e):
        if e.mimeData().hasText():
            e.acceptProposedAction()
            self.textDroppedSignal.emit(e.mimeData().text())
        if e.mimeData().hasUrls():
            e.acceptProposedAction()
            self.urlsDroppedSignal.emit(e.mimeData().urls())

    def setModel(self,model):
        QtWidgets.QTableView.setModel(self,model)
        self.model().dataChanged.connect(self.drawButtons)

    @QtCore.Slot()
    def drawButtons(self):
        for i in range(self.model().rowCount()):
            if self.model()._editable:
                viewButton = QtWidgets.QPushButton("View/Edit")
            else:
                viewButton = QtWidgets.QPushButton("View")
            #viewButton.setAutoFillBackground(True)
            self.setIndexWidget(self.model().index(i,self.model().columnCount(QtCore.QModelIndex())-1), viewButton)
            viewButton.clicked.connect(functools.partial(self.showSignal.emit,i))

class SequenceTable(QtWidgets.QWidget):

    loadFromSeqFileSignal = QtCore.Signal()
    loadFromCoorFileSignal = QtCore.Signal()
    loadFromTextSignal = QtCore.Signal()
    loadDataSignal = QtCore.Signal(tuple)
    urlsDroppedSignal = QtCore.Signal(list)
    textDroppedSignal = QtCore.Signal(str)

    def __init__(self,parent=None):
        QtWidgets.QWidget.__init__(self,parent)
        layout = QtWidgets.QVBoxLayout()
        self.setLayout(layout)
        self.table = SequenceTableView()
        layout.addWidget(self.table)
        buttons = QtWidgets.QWidget()
        buttonLayout = QtWidgets.QHBoxLayout()
        buttons.setLayout(buttonLayout)
        plusButton = SquareButton("+")
        self.minusButton = SquareButton("-")
        buttonLayout.addWidget(plusButton)
        buttonLayout.addWidget(self.minusButton)
        buttonLayout.addStretch(3)
        layout.addWidget(buttons)
        buttonLayout.setMargin(0)
        buttonLayout.setSpacing(0)
        buttonLayout.setContentsMargins(0,0,0,0)
        layout.setContentsMargins(0,0,0,0)

        plusButton.setFixedWidth(24)
        self.minusButton.setFixedWidth(24)
        plusButton.setFixedHeight(24)
        self.minusButton.setFixedHeight(24)

        self.setStyleSheet("SquareButton{font-size: 18px;font-family: Arial;} SequenceTableView > QToolTip{font-family:courier;}")

        self.minusButton.setEnabled(False)

        self.minusButton.clicked.connect(self.checkDelete)

        plusMenu = QtWidgets.QMenu(self)
        loadFromSeq = plusMenu.addAction("Load sequence file")
        loadFromCoor = plusMenu.addAction("Load from coordinate file")
        loadFromText = plusMenu.addAction("Enter text")
        loadFromText.triggered.connect(self.insertFromText)
        plusButton.setMenu(plusMenu)
        plusButton.setPopupMode(QtWidgets.QToolButton.InstantPopup)
        loadFromSeq.triggered.connect(self.loadFromSeqFileSignal.emit)
        loadFromCoor.triggered.connect(self.loadFromCoorFileSignal.emit)
        self.table.urlsDroppedSignal.connect(self.urlsDroppedSignal.emit)
        self.table.textDroppedSignal.connect(self.textDroppedSignal.emit)

        @QtCore.Slot('QPos')
        def onCustomContextMenu(pos):
            index = self.table.indexAt(pos)
            if index.isValid():
                contextMenu = QtWidgets.QMenu(self.table)
                indices = (" ").join([str(z+1) for z in sorted([ x.row() for x in self.table.selectionModel().selectedRows() ])])
                if len(indices)>1:
                    r = "rows"
                else:
                    r = "row"
                if len(indices)==1:
                    editAct = contextMenu.addAction("Edit number of copies")
                    editAct.triggered.connect(self.doEditCopies)
                    editNameAct = contextMenu.addAction("Edit name")
                    editNameAct.triggered.connect(self.doEditName)
                    editDescAct = contextMenu.addAction("Edit description")
                    editDescAct.triggered.connect(self.doEditDesc)
                    editTypeAct = contextMenu.addAction("Edit sequence type")
                    editTypeAct.triggered.connect(self.doEditType)
                delAct = contextMenu.addAction("Delete "+r+" "+indices)
                delAct.triggered.connect(self.checkDelete)
                contextMenu.exec_(self.table.viewport().mapToGlobal(pos))

        self.table.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.table.customContextMenuRequested.connect(onCustomContextMenu)
        self.table.deleteSignal.connect(self.checkDelete)

    @QtCore.Slot()
    def insertFromText(self):
        win = QtWidgets.QDialog()
        layout = QtWidgets.QVBoxLayout()
        nameEdit = QtWidgets.QLineEdit()
        nameLabel = QtWidgets.QLabel("Name")
        descrEdit = QtWidgets.QLineEdit()
        descrLabel = QtWidgets.QLabel("Description")
        tb = QtWidgets.QTextEdit()
        copiesSpin = QtWidgets.QSpinBox()
        copiesSpin.setMinimum(0)
        copiesSpin.setMaximum((2<<30)-1)
        copiesSpin.setValue(1)
        copiesLabel = QtWidgets.QLabel("Number of copies")
        typeCombo = QtWidgets.QComboBox()
        typeCombo.addItem("PROTEIN")
        typeCombo.addItem("RNA")
        typeCombo.addItem("DNA")
        typeCombo.setCurrentIndex(0)
        typeLabel = QtWidgets.QLabel("Type")
        tb.setStyleSheet("QTextEdit{font-family:courier;}")
        tb.setTextInteractionFlags(QtCore.Qt.TextEditorInteraction)
        nameLine = QtWidgets.QWidget()
        copiesLine = QtWidgets.QWidget()
        descrLine = QtWidgets.QWidget()
        typeLine = QtWidgets.QWidget()
        nameLayout = QtWidgets.QHBoxLayout()
        copiesLayout = QtWidgets.QHBoxLayout()
        descrLayout = QtWidgets.QHBoxLayout()
        typeLayout = QtWidgets.QHBoxLayout()
        nameLine.setLayout(nameLayout)
        copiesLine.setLayout(copiesLayout)
        descrLine.setLayout(descrLayout)
        typeLine.setLayout(typeLayout)
        nameLayout.addWidget(nameLabel)
        nameLayout.addWidget(nameEdit)
        copiesLayout.addWidget(copiesLabel)
        copiesLayout.addWidget(copiesSpin)
        descrLayout.addWidget(descrLabel)
        descrLayout.addWidget(descrEdit)
        typeLayout.addWidget(typeLabel)
        typeLayout.addWidget(typeCombo)
        nameLayout.setContentsMargins(0,0,0,0)
        copiesLayout.setContentsMargins(0,0,0,0)
        descrLayout.setContentsMargins(0,0,0,0)
        typeLayout.setContentsMargins(0,0,0,0)
        layout.addWidget(nameLine)
        layout.addWidget(copiesLine)
        layout.addWidget(descrLine)
        layout.addWidget(typeLine)
        textLabel = QtWidgets.QLabel("Sequence: <i>(plain text, no header)</i>")
        layout.addWidget(textLabel)
        layout.addWidget(tb)
        win.setLayout(layout)
        dbb = QtWidgets.QDialogButtonBox()
        layout.addWidget(dbb)
        closeButton = dbb.addButton(QtWidgets.QDialogButtonBox.Save)
        closeButton.setEnabled(False)
        @QtCore.Slot()
        def sequenceChangedHandler():
            if len(tb.toPlainText())>0:
                closeButton.setEnabled(True)
            else:
                closeButton.setEnabled(False)
        tb.textChanged.connect(sequenceChangedHandler)
        discardButton = dbb.addButton(QtWidgets.QDialogButtonBox.Discard)
        closeButton.clicked.connect(win.accept)
        discardButton.clicked.connect(win.reject)
        rv = win.exec_()
        if rv == QtWidgets.QDialog.Accepted:
            name = nameEdit.text()
            copies = copiesSpin.value()
            desc = descrEdit.text()
            seqType = typeCombo.currentText()
            seq = "".join(str(tb.toPlainText()).split())
            self.table.model().addItem((name,copies,desc,seqType,seq))
            #self.loadDataSignal.emit((name,copies,desc,seqType,seq))

    @QtCore.Slot('QItemSelection','QItemSelection')
    def enableMinusDepOnSelection(self,sel,desel):
       if len(self.table.selectionModel().selectedRows()) == 0:
           self.minusButton.setEnabled(False)
       else:
           self.minusButton.setEnabled(True)

    def setModel(self,model):
        self.table.setModel(model)
        self.table.selectionModel().selectionChanged.connect(self.enableMinusDepOnSelection)

    @QtCore.Slot()
    def doEditType(self):
        if len(self.table.selectionModel().selectedRows())==1:
            idx = self.table.selectionModel().selectedRows()[0]
            editIdx = idx.sibling(idx.row(),3)
            self.table.edit(editIdx)

    @QtCore.Slot()
    def doEditName(self):
        if len(self.table.selectionModel().selectedRows())==1:
            idx = self.table.selectionModel().selectedRows()[0]
            editIdx = idx.sibling(idx.row(),0)
            self.table.edit(editIdx)

    @QtCore.Slot()
    def doEditDesc(self):
        if len(self.table.selectionModel().selectedRows())==1:
            idx = self.table.selectionModel().selectedRows()[0]
            editIdx = idx.sibling(idx.row(),2)
            self.table.edit(editIdx)

    @QtCore.Slot()
    def doEditCopies(self):
        if len(self.table.selectionModel().selectedRows())==1:
            idx = self.table.selectionModel().selectedRows()[0]
            editIdx = idx.sibling(idx.row(),1)
            self.table.edit(editIdx)

    @QtCore.Slot()
    def checkDelete(self):
        if len(self.table.selectionModel().selectedRows())>1:
             text = "Really delete sequences:<br/><br/>"
        else:
             text = "Really delete sequence:<br/><br/>"
        for x in self.table.selectionModel().selectedRows():
            if len(self.table.model().sequenceData(x.row()))>20:
                text += str(x.row()+1)+" <tt>"+self.table.model().sequenceData(x.row())[:20]+"</tt>...<br/>"
            else:
                text += str(x.row()+1)+" <tt>"+self.table.model().sequenceData(x.row())+"</tt><br/>"
        text += "<br/><br/>?"
        answer = QtWidgets.QMessageBox.question(self,"Delete",text, QtWidgets.QMessageBox.No|QtWidgets.QMessageBox.Yes)
        if answer == QtWidgets.QMessageBox.Yes:
#FIXME - Need to do something different ? Should I in fact be just emitting signal and leaving remove elsewhere?
            self.table.model().removeRows(self.table.selectionModel().selectedRows())
            self.minusButton.setEnabled(False)

class SequenceTableDialog(QtWidgets.QDialog):

    loadFromSeqFileSignal = QtCore.Signal()
    loadFromCoorFileSignal = QtCore.Signal()
    loadFromTextSignal = QtCore.Signal()
    loadDataSignal = QtCore.Signal(tuple)
    urlsDroppedSignal = QtCore.Signal(list)
    textDroppedSignal = QtCore.Signal(str)

    def keyPressEvent(self,e):
#FIXME - Never get here ... Lack of focus?
        print(e.key())
        if e.key()==QtCore.Qt.Key_Escape:
            self.reject()

    def __init__(self,parent=None):
        QtWidgets.QDialog.__init__(self,parent)
        st = SequenceTable(self)
        layout = QtWidgets.QVBoxLayout(self)
        self.setLayout(layout)
        layout.addWidget(st)
        st.loadFromSeqFileSignal.connect(self.loadFromSeqFileSignal.emit)
        st.loadFromCoorFileSignal.connect(self.loadFromCoorFileSignal.emit)
        st.loadFromTextSignal.connect(self.loadFromTextSignal.emit)
        st.loadDataSignal.connect(self.loadDataSignal.emit)
        st.urlsDroppedSignal.connect(self.urlsDroppedSignal.emit)
        st.textDroppedSignal.connect(self.textDroppedSignal.emit)
        dbb = QtWidgets.QDialogButtonBox()
        layout.addWidget(dbb)
        closeButton = dbb.addButton(QtWidgets.QDialogButtonBox.Save)
        discardButton = dbb.addButton(QtWidgets.QDialogButtonBox.Discard)
        closeButton.clicked.connect(self.accept)
        discardButton.clicked.connect(self.reject)
        layout.setContentsMargins(0,0,0,0)

if __name__ == "__main__":
    import sip
    sip.setdestroyonexit(False)

    app = QtWidgets.QApplication(sys.argv)

    win = SequenceTable()
    win.show()
    win.raise_()

    model = SequenceModel()
    win.setModel(model)

    def addRandomSequence():
        import random
        import string
        name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
        copies = random.randint(1,20)
        desc = "This sequence is "+''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(10))
        seqType = ("PROTEIN","RNA","DNA")[random.randint(0,2)]
        lenSeq = random.randint(1000,2000)
        if seqType == "PROTEIN":
            seq = ''.join(random.choice("ARNDBCEQZGHILKMFPSTWYZ") for _ in range(lenSeq))
        elif seqType == "RNA":
            seq = ''.join(random.choice("ACGTIU") for _ in range(lenSeq))
        else:
            seq = ''.join(random.choice("ACGT") for _ in range(lenSeq))
        model.addItem((name,copies,desc,seqType,seq))

    win.resize(700,200)

    for i in range(10):
        addRandomSequence()

    sys.exit(app.exec_())
