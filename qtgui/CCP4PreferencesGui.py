"""
     qtgui/CCP4Peferences.py: CCP4 Gui Project
     Copyright (C) 2011 University of York

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
    Liz Potterton - create class to maintain GUIPreferences - Sept 2011
'''
from PySide2 import QtCore,QtGui, QtWidgets


class CPreferencesWindow(QtWidgets.QDialog):

  def __init__(self,parent):
    QtWidgets.QDialog.__init__(self,parent)
    self.setLayout(QtWidgets.QVBoxLayout())

    from core import CCP4Modules

    #print 'CPreferencesWindow.__init__',CCP4Modules.QTAPPLICATION().
   
    widgetClass = CCP4Modules.TASKMANAGER().getTaskWidgetClass('guipreferences')
    widget = widgetClass(self,layoutMode='TAB')
    prefContainer = CCP4Modules.PREFERENCES()
    widget.setContainer(prefContainer)
    widget.draw()
    self.layout().addWidget(widget)
    #except:
    #  pass
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
#Remove the hardwired layout mostly for the Prefs GUI.
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
    from core import CCP4Preferences,CCP4Modules
    from qtgui import CCP4ProjectManagerGui
    CCP4Preferences.CPreferences.insts.save()
    try:
      CCP4Preferences.CPreferences.insts.EXEPATHLIST.setupExeLookup()
    except:
      pass

    try:
      CCP4Modules.QTAPPLICATION().setNamedStyle(CCP4Preferences.CPreferences.insts.WINDOWS_STYLE)
    except:
      CCP4Modules.QTAPPLICATION().setNamedStyle(CCP4Preferences.CPreferences.insts.WINDOWS_STLYE)
    
    from qtgui import CCP4StyleSheet,CCP4WebView,CCP4ProjectViewer
    CCP4StyleSheet.setStyleSheet()
    CCP4WebView.setGlobalSettings()

    for pV in CCP4ProjectViewer.CProjectViewer.Instances:
      try:
        pV.toolBar('project').setToolButtonStyle(int(CCP4Preferences.CPreferences.insts.TOOLBARBUTTONSSTYLE))
        pV.resizeColumnToContents(0)
      except:
        pass

    CCP4Modules.QTAPPLICATION().prefsChanged.emit()

  def close(self):
    self.doApply()
    QtWidgets.QDialog.close(self)

