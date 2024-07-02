from __future__ import print_function

#TODO:
"""
   get/set data from "afar"
"""

import sys

from PySide6.QtCore import Slot, Signal
from PySide6.QtCore import QAbstractTableModel, QRegularExpression, Qt
from PySide6.QtWidgets import QApplication, QWidget, QHBoxLayout, QVBoxLayout, QTableView, QAbstractItemView, QHeaderView, QStyledItemDelegate, QSpinBox, QLineEdit, QToolButton, QSizePolicy, QStyleOptionToolButton, QStyle, QMessageBox, QMenu
from PySide6.QtGui import QRegularExpressionValidator, QPainter

class SquareButton(QToolButton):

    def __init__(self,text,parent=None):
        QToolButton.__init__(self,parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        self.setText(text)

    def paintEvent(self,_paintEvent ):
        styleOption = QStyleOptionToolButton()
        QToolButton.initStyleOption(self,styleOption)
        styleOption.features &= ~QStyleOptionToolButton.HasMenu

        painter = QPainter(self)
        self.style().drawComplexControl( QStyle.CC_ToolButton, styleOption, painter, self )

class AtomSelectionListModel(QAbstractTableModel):

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

    def addItem(self,item):
        self.beginResetModel()
        self._data.append(item)
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.endResetModel()
        self.dataChanged.emit(tl,br)

    def setEditable(self,val):
        self.beginResetModel()
        self._editable = val
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.endResetModel()
        self.dataChanged.emit(tl,br)

    def flags(self,index):
        if not self._editable:
            return QAbstractTableModel.flags(self,index) &~ Qt.ItemIsEditable
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.isValid():
            return QAbstractTableModel.flags(self,index) | Qt.ItemIsEditable | Qt.ItemIsEnabled | Qt.ItemIsSelectable;
        return Qt.NoItemFlags

    def headerData(self,section, orientation, role):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            if section == 0:
                return "Chain"
            elif section == 1:
                return "From residue"
            elif section == 2:
                return "To residue"
            elif section == 3:
                return "Specific atoms"
        else:
            return QAbstractTableModel.headerData(self,section, orientation, role)

    def data(self,index,role=Qt.DisplayRole):
        if role == Qt.DisplayRole or role == Qt.EditRole:
            return self._data[index.row()][index.column()]
        return None

    def setData(self,index,val,role):
        if role == Qt.EditRole:
            self._data[index.row()][index.column()] = val
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.dataChanged.emit(tl,br)
        return QAbstractTableModel.setData(self,index,val,role)

    def columnCount(self,parent):
        return 4

    def rowCount(self,parent=None):
        return len(self._data)

    def __init__(self,parent=None):
        QAbstractTableModel.__init__(self,parent)
        self._data = [["","","",""]]
        self._editable = True

    def selectionData(self):
        selData = []
        for dat in self._data:
            if dat[0] or dat[1] or dat[2] or dat[3]:
                selData.append({"chainId":dat[0],"firstRes":dat[1],"lastRes":dat[2],"atoms":dat[3]})
        return selData

    def setSelectionData(self,sels):
        self.beginResetModel()

        self._data = []
        for sel0 in sels:
            sel = sel0.get()
            item = [sel["chainId"],sel["firstRes"],sel["lastRes"],sel["atoms"]]
            self._data.append(item)

        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.endResetModel()
        self.dataChanged.emit(tl,br)

class RigidBodySelectionListModel(QAbstractTableModel):

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

    def addItem(self,item):
        self.beginResetModel()
        self._data.append(item)
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.endResetModel()
        self.dataChanged.emit(tl,br)

    def setEditable(self,val):
        self.beginResetModel()
        self._editable = val
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.endResetModel()
        self.dataChanged.emit(tl,br)

    def flags(self,index):
        if not self._editable:
            return QAbstractTableModel.flags(self,index) &~ Qt.ItemIsEditable
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.isValid():
            return QAbstractTableModel.flags(self,index) | Qt.ItemIsEditable | Qt.ItemIsEnabled | Qt.ItemIsSelectable;
        return Qt.NoItemFlags

    def headerData(self,section, orientation, role):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            if section == 0:
              return "Group ID"
            elif section == 1:
                return "Chain"
            elif section == 2:
                return "From residue"
            elif section == 3:
                return "To residue"
        else:
            return QAbstractTableModel.headerData(self,section, orientation, role)

    def data(self,index,role=Qt.DisplayRole):
        if role == Qt.DisplayRole or role == Qt.EditRole:
            return self._data[index.row()][index.column()]
        return None

    def setData(self,index,val,role):
        if role == Qt.EditRole:
            self._data[index.row()][index.column()] = val
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.dataChanged.emit(tl,br)
        return QAbstractTableModel.setData(self,index,val,role)

    def columnCount(self,parent):
        return 4

    def rowCount(self,parent=None):
        return len(self._data)

    def __init__(self,parent=None):
        QAbstractTableModel.__init__(self,parent)
        self._data = [["","","",""]]
        self._editable = True

    def selectionData(self):
        selData = []
        for dat in self._data:
            if dat[0] or dat[1] or dat[2] or dat[3]:
                selData.append({"groupId":dat[0],"chainId":dat[1],"firstRes":dat[2],"lastRes":dat[3]})
        return selData

    def setSelectionData(self,sels):
        self.beginResetModel()

        self._data = []
        for sel0 in sels:
            sel = sel0.get()
            item = [sel["groupId"],sel["chainId"],sel["firstRes"],sel["lastRes"]]
            self._data.append(item)

        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.endResetModel()
        self.dataChanged.emit(tl,br)

class OccupancySelectionListModel(QAbstractTableModel):

    def removeRows(self,idxs):
        self.beginResetModel()
        indices = [ x.row() for x in idxs ]
        dataNew = []
        for i,row in zip(range(len(self._data)),self._data):
            if not i in indices:
                dataNew.append(row)
        self._data = dataNew
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,5)
        self.endResetModel()
        self.dataChanged.emit(tl,br)

    def addItem(self,item):
        self.beginResetModel()
        self._data.append(item)
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,5)
        self.endResetModel()
        self.dataChanged.emit(tl,br)

    def setEditable(self,val):
        self.beginResetModel()
        self._editable = val
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,5)
        self.endResetModel()
        self.dataChanged.emit(tl,br)

    def flags(self,index):
        if not self._editable:
            return QAbstractTableModel.flags(self,index) &~ Qt.ItemIsEditable
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.isValid():
            return QAbstractTableModel.flags(self,index) | Qt.ItemIsEditable | Qt.ItemIsEnabled | Qt.ItemIsSelectable;
        return Qt.NoItemFlags

    def headerData(self,section, orientation, role):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            if section == 0:
              return "Group ID"
            elif section == 1:
                return "Chain(s)"
            elif section == 2:
                return "From residue"
            elif section == 3:
                return "To residue"
            elif section == 4:
                return "(Atom)"
            elif section == 5:
                return "(Alt)"
        else:
            return QAbstractTableModel.headerData(self,section, orientation, role)

    def data(self,index,role=Qt.DisplayRole):
        if role == Qt.DisplayRole or role == Qt.EditRole:
            #print(self._data)
            return self._data[index.row()][index.column()]
        return None

    def setData(self,index,val,role):
        print("setData")
        if role == Qt.EditRole:
            self._data[index.row()][index.column()] = val
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,5)
        self.dataChanged.emit(tl,br)
        return QAbstractTableModel.setData(self,index,val,role)

    def columnCount(self,parent):
        return 6

    def rowCount(self,parent=None):
        return len(self._data)

    def __init__(self,parent=None):
        QAbstractTableModel.__init__(self,parent)
        self._data = [["","","","","",""]]
        self._editable = True

    def selectionData(self):
        selData = []
        for dat in self._data:
            if dat[0] or dat[1] or dat[2] or dat[3] or dat[4] or dat[5]:
                selData.append({"groupId":dat[0],"chainIds":dat[1],"firstRes":dat[2],"lastRes":dat[3],"atoms":dat[4],"alt":dat[5]})
        return selData

    def setSelectionData(self,sels):
        self.beginResetModel()

        self._data = []
        for sel0 in sels:
            sel = sel0.get()
            item = [sel["groupId"],sel["chainIds"],sel["firstRes"],sel["lastRes"],sel["atoms"],sel["alt"]]
            self._data.append(item)

        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,5)
        self.endResetModel()
        self.dataChanged.emit(tl,br)

class OccupancyRelationListModel(QAbstractTableModel):

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

    def addItem(self,item):
        self.beginResetModel()
        self._data.append(item)
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.endResetModel()
        self.dataChanged.emit(tl,br)

    def setEditable(self,val):
        self.beginResetModel()
        self._editable = val
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.endResetModel()
        self.dataChanged.emit(tl,br)

    def flags(self,index):
        if not self._editable:
            return QAbstractTableModel.flags(self,index) &~ Qt.ItemIsEditable
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.isValid():
            return QAbstractTableModel.flags(self,index) | Qt.ItemIsEditable | Qt.ItemIsEnabled | Qt.ItemIsSelectable;
        return Qt.NoItemFlags

    def headerData(self,section, orientation, role):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            if section == 0:
              return "Group IDs (space seperated list)"
        else:
            return QAbstractTableModel.headerData(self,section, orientation, role)

    def data(self,index,role=Qt.DisplayRole):
        if role == Qt.DisplayRole or role == Qt.EditRole:
            return self._data[index.row()][index.column()]
        return None

    def setData(self,index,val,role):
        if role == Qt.EditRole:
            self._data[index.row()][index.column()] = val
        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.dataChanged.emit(tl,br)
        return QAbstractTableModel.setData(self,index,val,role)

    def columnCount(self,parent):
        return 1

    def rowCount(self,parent=None):
        return len(self._data)

    def __init__(self,parent=None):
        QAbstractTableModel.__init__(self,parent)
        self._data = [[""]]
        self._editable = True

    def selectionData(self):
        selData = []
        for dat in self._data:
            if dat[0]:
               selData.append({"groupIds":dat[0]})
        return selData

    def setSelectionData(self,sels):
        self.beginResetModel()

        self._data = []
        for sel0 in sels:
            sel = sel0.get()
            item = [sel["groupIds"]]
            self._data.append(item)

        tl = self.index(0,0)
        br = self.index(self.rowCount()-1,3)
        self.endResetModel()
        self.dataChanged.emit(tl,br)


class NoSpaceDelegate(QStyledItemDelegate):

    def __init__(self,parent=None):
        QStyledItemDelegate.__init__(self,parent)

    def createEditor(self,parent, option, index):
        editor = QLineEdit(parent)
        validator = QRegularExpressionValidator(QRegularExpression("\\S+"), self )
        editor.setValidator(validator)
        editor.setFrame(False)
        return editor

    def setEditorData(self,editor, index):
        value = index.data(Qt.EditRole)
        editor.setText(value)

    def setModelData(self,editor, model, index):
        model.setData(index, editor.text(), Qt.EditRole)

class AtomSelectionListTableView(QTableView):

    deleteSignal = Signal()

    def __init__(self,parent=None):
        QTableView.__init__(self,parent)
        self.setSelectionBehavior(QAbstractItemView.SelectRows)
        try:
            self.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        except:
            self.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.setEditTriggers(QAbstractItemView.AllEditTriggers)
#        chainDelegate = NoSpaceDelegate(self)
#        self.setItemDelegateForColumn(0,chainDelegate)
#        spinDelegateFrom = NoSpaceDelegate(self)
#        self.setItemDelegateForColumn(1,spinDelegateFrom)
#        spinDelegateTo = NoSpaceDelegate(self)
#        self.setItemDelegateForColumn(2,spinDelegateTo)

class MultiAtomSelection(QWidget):

    def __init__(self,parent=None):
        QWidget. __init__(self,parent)
        layout = QVBoxLayout()
        self.setLayout(layout)
        self.table = AtomSelectionListTableView()
        self.table.verticalHeader().hide()
        layout.addWidget(self.table)
        buttons = QWidget()
        buttonLayout = QHBoxLayout()
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

        self.minusButton.setEnabled(False)
        self.minusButton.clicked.connect(self.checkDelete)

        self.setStyleSheet("SquareButton{font-size: 18px;font-family: Arial;} AtomSelectionListTableView > QToolTip{font-family:courier;}")

        plusButton.clicked.connect(self.addRow)

        def onCustomContextMenu(pos):
            index = self.table.indexAt(pos)
            if index.isValid():
                contextMenu = QMenu(self.table)
                indices = (" ").join([str(z+1) for z in sorted([ x.row() for x in self.table.selectionModel().selectedRows() ])])
                if len(indices)>1:
                    r = "rows"
                else:
                    r = "row"

                delAct = contextMenu.addAction("Delete "+r+" "+indices)
                delAct.triggered.connect(self.checkDelete)
                contextMenu.exec_(self.table.viewport().mapToGlobal(pos))

        self.table.setContextMenuPolicy(Qt.CustomContextMenu)
        self.table.customContextMenuRequested.connect(onCustomContextMenu)
        self.table.deleteSignal.connect(self.checkDelete)

        self.setAttribute(Qt.WA_DeleteOnClose)

    def addRow(self):
        self.table.selectionModel().clear()
        self.table.model().addItem(["1"]+[""]*(self.table.model().columnCount(None)-1))
        self.minusButton.setEnabled(False)

    def checkDelete(self):
        if len(self.table.selectionModel().selectedRows())>1:
             text = "Really delete selection:<br/><br/>"
        else:
             text = "Really delete selection:<br/><br/>"
        for x in self.table.selectionModel().selectedRows():
            text += str(x.row()+1)+"<br/>"
        text += "<br/><br/>?"
        answer = QMessageBox.question(self,"Delete",text, QMessageBox.No|QMessageBox.Yes)
        if answer == QMessageBox.Yes:
#FIXME - Need to do something different ? Should I in fact be just emitting signal and leaving remove elsewhere?
            self.table.model().removeRows(self.table.selectionModel().selectedRows())
            self.minusButton.setEnabled(False)

    def enableMinusDepOnSelection(self,sel,desel):
        print(len(self.table.selectionModel().selectedRows()))
        if len(self.table.selectionModel().selectedRows()) == 0:
            self.minusButton.setEnabled(False)
        else:
            self.minusButton.setEnabled(True)

    def setModel(self,model):
        self.table.setModel(model)
        self.table.selectionModel().selectionChanged.connect(self.enableMinusDepOnSelection)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = MultiAtomSelection()
    win.show()
    win.raise_()
    model = AtomSelectionListModel()
    win.setModel(model)
    sys.exit(app.exec_())
