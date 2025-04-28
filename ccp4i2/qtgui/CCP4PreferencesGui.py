"""
Copyright (C) 2011 University of York
Liz Potterton - create class to maintain GUIPreferences - Sept 2011
"""

from PySide2 import QtCore, QtWidgets

from ..core.CCP4Modules import PREFERENCES, TASKMANAGER
from ..utils.QApp import QTAPPLICATION


class CPreferencesWindow(QtWidgets.QDialog):

  def __init__(self,parent):
    QtWidgets.QDialog.__init__(self,parent)
    self.setLayout(QtWidgets.QVBoxLayout())

    widgetClass = TASKMANAGER().getTaskWidgetClass('guipreferences')
    widget = widgetClass(self,layoutMode='TAB')
    prefContainer = PREFERENCES()
    widget.setContainer(prefContainer)
    widget.draw()
    self.layout().addWidget(widget)
    self.buttons = QtWidgets.QDialogButtonBox(self)
    but = self.buttons.addButton(QtWidgets.QDialogButtonBox.Apply)
    but.clicked.connect(self.doApply)
    but = self.buttons.addButton(QtWidgets.QDialogButtonBox.Close)
    but.clicked.connect(self.close)
    line = QtWidgets.QHBoxLayout()
    line.addStretch(0.2)
    line.addWidget(self.buttons)
    line.addStretch(0.2)
    self.layout().addLayout(line)

    paddingAllowance = 24
    maxWidth = 1600
    # Remove the hardwired layout mostly for the Prefs GUI.
    ss = '''
                         QLabel#errorMessage  { color : red; max-width: ''' + str(maxWidth-paddingAllowance) + '''px; }
                         CComplexLineWidget  {  max-width: ''' + str(maxWidth-paddingAllowance) + '''px; }
                         CTabTaskWidget  {  max-width: ''' + str(maxWidth-paddingAllowance) + '''px; min-width: 600px }
                         CTextEdit {  max-width: ''' + str(maxWidth-paddingAllowance) + '''px; min-width: 60px }
                         CStackedWidget {  max-width: ''' + str(maxWidth-paddingAllowance) + '''px; min-width: 60px }
                         QFrame#taskLine {  max-width: ''' + str(maxWidth-paddingAllowance) + '''px; }
                         CProjectViewer  {  max-width: '''+str(maxWidth)+ '''px;  min-width: 60px }
                         CReportView {  max-width: '''+ str(maxWidth)+'''px;  min-width: 60px }
                         QScrollArea#messageScrollArea  {  max-width: ''' + str(maxWidth-paddingAllowance) + '''; min-width: 60px }
    '''
    self.setStyleSheet(ss)

  @QtCore.Slot()
  def doApply(self):
    PREFERENCES().save()
    try:
      PREFERENCES().EXEPATHLIST.setupExeLookup()
    except:
      pass

    try:
      QTAPPLICATION().setNamedStyle(PREFERENCES().WINDOWS_STYLE)
    except:
      QTAPPLICATION().setNamedStyle(PREFERENCES().WINDOWS_STLYE)

    from ..qtgui import CCP4ProjectViewer, CCP4StyleSheet, CCP4WebView
    CCP4StyleSheet.setStyleSheet()
    CCP4WebView.setGlobalSettings()

    for pV in CCP4ProjectViewer.CProjectViewer.Instances:
      try:
        pV.toolBar('project').setToolButtonStyle(int(PREFERENCES().TOOLBARBUTTONSSTYLE))
        pV.resizeColumnToContents(0)
      except:
        pass

    QTAPPLICATION().prefsChanged.emit()

  def close(self):
    self.doApply()
    QtWidgets.QDialog.close(self)
