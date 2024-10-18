from __future__ import print_function

# Some classes to help handle dynamic GUI buttons

from PySide2 import QtCore,QtGui, QtWidgets

import math
import functools
# -------------------------------------------------------------
class ChoiceButtons(QtWidgets.QWidget):
#
# A class to handle a variable number of choice buttons, with annotations
# Buttons will be displayed vertically, with their name and a label
#   followed by an optional line of further annotation
#
# Usage (eg within a class):
#  1) Construction and display
#    # Create object ready for use: doesn't do much, no buttons
#    self.buttonobject = ChoiceButtons()
#    # put into place for display
#    currentframe.addWidget(self.buttonobject)
#    # connect to click process routine (see (2) below)
#    self.buttonobject.clickedSignal.connect(self.buttonClicked)
#
#  2) create function to process button clicks
#    def buttonClicked(self):
#      s = self.buttonobject.selected   # get name string from button click
#      storesomething = s   #  do something with string, eg store, or act on it
#
#  3) when required, set display dynamically
#    title = "A title for the button set"
#    # A list of names (strings) for each button:
#    #  these names will be returned to the buttonClicked function
#    choices = [button1, button2, button3, ...]
#    # Optional list of tags to display on the same line as the buttons
#    tags = [tag1, tag2, tag3, ...]  # same number as choices
#    # Optional list of further annotations for each button
#    notes = [note1, note2, note3, ...]  # same number as choices
#    # set up display
#    self.buttonobject.setChoices(title, buttons, labels, notes)
#

    clickedSignal = QtCore.Signal(str)

    def __init__(self,parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.selected = ""
        layout = QtWidgets.QVBoxLayout()
        self.setLayout(layout)
        layout.setSpacing(0)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def setChoices(self, title, choices, tags=None, notes=None, subtitle=None, singleChoice=True):
        # title       heading for the pane (bold)
        # choices     list of button names
        # tags        list of labels for button names, displayed on same line
        # notes       more information about the button, next line(s)

        tags_ = ['' for i in range(len(choices))]
        if tags_ is not None and (len(choices) == len(tags)):
            tags_ = tags
        notes_ = notes       # could be None
        if notes is not None:
            if (len(notes) == 0 or (len(notes) == 1 and notes[0] == '')):
                notes_ = None

        layout = self.layout()
        self.clearLayout(layout)

        if title is not None and len(title) > 0:
            # Header title
            layout.addWidget(QtWidgets.QLabel("<b>"+title+"</b>"))
        if subtitle is not None and len(subtitle) > 0:
            # Header title
            layout.addWidget(QtWidgets.QLabel("<i>"+subtitle+"</i>"))
        

        # draw each choice + labels etc
        for i in range(len(choices)):
            # choice + tag on same line
            linelayout = QtWidgets.QHBoxLayout()
            linelayout.setSpacing(0)
            label = QtWidgets.QLabel('>> ')
            linelayout.addWidget(label)
            c = str(choices[i])  # the choice
            if singleChoice:
                button = QtWidgets.QRadioButton(str(c))
            else:
                button = QtWidgets.QCheckBox(str(c))
            button.setMinimumWidth(80)
            linelayout.addWidget(button)
            button.clicked.connect(functools.partial(self.setSelected, c))
            button.clicked.connect(functools.partial(self.clickedSignal.emit, c))
            linelayout.setStretch(1,5)
            if tags_[i] != '':
                s = str(tags_[i])
                label = QtWidgets.QLabel(s)
                linelayout.addWidget(label)
                linelayout.setStretch(2,5)

            layout.addLayout(linelayout)

            # additional info if present
            if not notes_:
                continue
            if notes_ is not None or len(notes_) < i+1:
                if notes_[i] != '':
                    s = "    "+str(notes_[i])
                    layout.addWidget(QtWidgets.QLabel(s))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    @QtCore.Slot(str)
    def setSelected(self,s):
        self.selected = s

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def clearLayout(self, layout):
        # Clear everything in a layout, from stack overflow
        if layout is not None:
            while layout.count():
                child = layout.takeAt(0)
                if child.widget() is not None:
                    child.widget().deleteLater()
                elif child.layout() is not None:
                    self.clearLayout(child.layout())

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def addOtherText(self, header, textlist):
        # List of text strings to be dispaly below buttons
        layout = self.layout()
        layout.addWidget(QtWidgets.QLabel('<br/>'))
        if header is not None:
            label = QtWidgets.QLabel(str(header))
            label.setStyleSheet("color:SteelBlue; font-weight:bold;")
            layout.addWidget(label)

        for text in textlist:
            label = QtWidgets.QLabel(str(text))
            label.setStyleSheet("color:SteelBlue;")
            layout.addWidget(label)
            
#-------------------------------------------------------------------
class MyMessageBox:
    # A simple warning message box
    # Exits from displayText are: 0 for Open, 1 for Abort
    def __init__(self):
        self.msgBox = QtWidgets.QMessageBox()
        self.msgBox.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                  QtWidgets.QSizePolicy.Expanding)
        # default buttons
        buttonOpen = QtWidgets.QMessageBox.Open
        buttonAbort = QtWidgets.QMessageBox.Abort
        self.buttonlist = [buttonOpen, buttonAbort]
        self.default = 0
        self.icon = None

    def singleButton(self):
        buttonClose = QtWidgets.QMessageBox.Close
        self.buttonlist = [buttonClose]

    def setDefault(self, idx):
        # which button is default
        j = max(min(idx,len(self.buttonlist)), 0)
        self.default = j

    def setInformativeText(self, inftext):
        self.msgBox.setInformativeText(inftext)

    def warning(self):
        self.icon = QtWidgets.QMessageBox.Warning
        
    def displayText(self, text, mybuttons=None):
        if self.icon is not None:
            self.msgBox.setIcon(self.icon)
        self.msgBox.setText(text)
        # standard buttons
        if mybuttons is not None:
            self.buttonlist = mybuttons
        buttons = 0
        for b in self.buttonlist: buttons |= b
        self.msgBox.setStandardButtons(buttons)
        self.msgBox.setDefaultButton(self.buttonlist[self.default]);
        ret = self.msgBox.exec_();
        idx = self.buttonlist.index(ret)
        return idx
