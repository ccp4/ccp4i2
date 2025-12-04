from __future__ import print_function

"""
     tasks/guipreferences.py: CCP4 GUI Project
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

"""
     Liz Potterton Sept 2011 - Create a task window for GUI preferences
"""
import sys
from baselayer import QtGui, QtWidgets,QtCore

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class CGuiPreferences(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'guipreferences'
  TASKVERSION = 0.0
  TASKMODULE='preferences'
  TASKTITLE='CCP4i2 Preferences'


  def drawContents(self):
    
    from core import CCP4Modules
    from core.CCP4Config import DEVELOPER

    self.container.WINDOWS_STYLE.setQualifier('enumerators',CCP4Modules.QTAPPLICATION().getStyleKeys())

    self.openFolder(title='Task interface')


    if DEVELOPER():
        self.createLine( [ 'label' , 'Use windows style', 'widget' , 'WINDOWS_STYLE'] )

    self.createLine( [ 'label' , 'Use font size',
                      'widget','GUI_FONT_SIZE',
                      'label','and report font size',
                      'widget','REPORT_FONT_SIZE'] )

    self.createLine( [ 'label' , 'Invalid/missing data highlighted by colouring ',
                       'widget', 'INVALID_FRAME_MODE',
                       'label' ,'width',
                       'widget', 'INVALID_FRAME_WIDTH',
                       'label' , 'and colour',
                       'widget', 'INVALID_FRAME_COLOUR' ], toggle = ['INVALID_FRAME_MODE', 'open', [ 1 ] ] )
    
    self.createLine( [ 'label' , 'Invalid/missing data highlighted by colouring ',
                       'widget', 'INVALID_FRAME_MODE',
                       'label' , 'and colour',
                       'widget', 'INVALID_FRAME_COLOUR' ], toggle = ['INVALID_FRAME_MODE', 'open', [ 0 ] ]  )
    
                      
    self.createLine( [ 'label', 'Show toolbar buttons as' , 'widget' , '-guiMode', 'radio', 'TOOLBARBUTTONSSTYLE'] )
    
    
    #self.setMenuText('TASK_WINDOW_LAYOUT',{ 'TAB':'tab frames','FOLDER':'folders'})
    # Show REFINE_MODE options as radio buttons with one button per line
    #self.createLine( ['label','Task windows laid out as',
    #                  'widget','-guiMode','radio','TASK_WINDOW_LAYOUT'])
    self.createLine( ['widget','HD_ICONS',
                      'label','Use HD Icons'])
    self.createLine( ['widget','COMPACT_TASK_MENU',
                      'label','Compact task menus'])
    self.createLine( ['widget','TABLES_ALTERNATING_COLOR',
                      'label','Tables with alternating coloured rows'])

    self.createLine( ['widget','AUTO_INFO_ON_FILE_IMPORT',
                       'label' , 'Query provenance of imported files' ] )
    self.createLine( ['label','Maximum number of lines in file to display in viewer',
                      'widget','TEXT_VIEW_LINE_LIMIT' ] )

    
    self.createLine( ['widget','EXTERNAL_FILES_IN_EXTERNAL_BROWSER',
                      'label' , 'Display non-local files in your usual web browser' ] )

    self.createLine( ['widget','EXTERNAL_FILES_IN_IFRAME',
                      'label' , 'Display custom html inside ccp4i2 reports' ] )

    deleteLineWidget = self.createLine( ['widget','DELETE_INTERACTIVE_JOBS','label','Delete interactive jobs (such as Coot) that have no output files' ] )
    deleteCB = deleteLineWidget.findChild(QtWidgets.QCheckBox)

    warnLineWidget = self.createLine( ['label', '&nbsp;','widget','SHOW_DELETE_INTERACTIVE_JOBS','label','Show warning when deleting interactive jobs' ] )
    warnCB = warnLineWidget.findChild(QtWidgets.QCheckBox)
    warnLabels = warnLineWidget.findChildren(QtWidgets.QLabel)
    warnLabel = None
    if len(warnLabels)>1:
        warnLabel = warnLabels[1]

    self.createLine( ['widget','JOB_LIST_DATE_TIME', 'label','Show full date and time in job list (requires restart)'])

    #I'd really like to either:
    #have access to the real QCheckBox/QLabel, i.e. let me write Qt
    #or to be able to pass an "enabled" dependency to createLine
    #or have a dependency in def.xml.

    if hasattr(deleteCB,"stateChanged") and hasattr(deleteCB.stateChanged,"connect") and hasattr(warnCB,"setEnabled"):
        if deleteCB.isChecked():
            warnCB.setEnabled(True)
        else:
            warnCB.setEnabled(False)
        deleteCB.stateChanged.connect(warnCB.setEnabled)

    if hasattr(deleteCB,"stateChanged") and hasattr(deleteCB.stateChanged,"connect") and hasattr(warnLabel,"setEnabled"):
        if deleteCB.isChecked():
            warnLabel.setEnabled(True)
        else:
            print("Disabling warnLabel")
            warnLabel.setEnabled(False)
        deleteCB.stateChanged.connect(warnLabel.setEnabled)

    #self.openFolder(title='Test alteratives')
    
    
    self.openFolder(title='Other software')
    
    self.createLine( ['label' , 'CCP4mg executable full path', 'widget','-jobCombo',False,'CCP4MG_EXECUTABLE' ] )
    self.createLine( ['label' , 'Coot executable full path', 'widget','-jobCombo',False,'COOT_EXECUTABLE' ] )
    self.createLine( ['label' , 'Directory containing Shelx programs', 'widget','-jobCombo',False,'SHELXDIR' ] )
    self.createLine( ['label' , 'Directory containing DIALS programs', 'widget','-jobCombo',False,'DIALSDIR' ] )
    self.createLine(['label', 'Directory containing BUSTER (GPL)', 'widget', '-jobCombo', False, 'BUSTERDIR'])
    self.createLine( ['advice', 'Enter full path to any program executable that is not the standard release version'])
    self.createLine( ['widget' , 'EXEPATHLIST' ] )

    
    self.openFolder(title='Advanced')
    self.createLine( ['label' ,'The following polling option may be necessary to enable report real time updates on networked systems' ])
    self.createLine( ['widget','FILESYSTEMWATCHERPOLLER','label','Use file system watcher polling mechanism (restart gui)'])
    self.createLine( ['widget','BZR_DOWNLOAD','label','Enable update of CCP4i2 from code repository (NOT normally recommended!)'])
    self.createLine( ['widget','DBLOCAL_QUIT_RUNNING','label','Allow quitting of ccp4i2 when jobs are running in "Local database mode" (NOT normally recommended!)'])
    self.createLine( ['advice','Deleting temporary files' ])
    self.createLine( ['widget','RETAIN_DIAGNOSTIC_FILES','label','Retain diagnostic files after job run'])
    self.createLine( ['widget','CLEANUP_ON_EXIT','label','Delete temporary files on CCP4i2 closedown'])
    self.createLine( ['widget','NATIVEFILEBROWSER','label','Use system native file browser'])
    self.createLine( ['widget','SHOW_TASK_MODE_BUTTONS','label','Show module buttons in task list (project must be closed and re-opened)'])
    self.createLine( ['widget','RESTORE_TO_TASKLIST','label','Show the task list on startup instead of the last open job (requires restart)'])
    self.createLine( ['widget','AUTO_UPDATE_REPORT80','label','Automatically update reports from old jobs to be compatible with CCP4 8.0'])
    self.createLine( ['widget','SHOW_WRAPPERS','label','Show testing/development tasks'])
    self.createLine( ['widget','DISABLE_WEBGL','label','Disable 3D views in reports (requires restart)'])
    self.createLine( [ 'label' , 'Invalid/missing data uses ',
                       'widget', 'INVALID_FRAME_MODE',
                       'label' , 'style' ], toggle = ['INVALID_FRAME_MODE', 'open', [ 2 ] ]  )

