from __future__ import print_function

"""
     CCP4StyleSheet.py: CCP4 GUI Project
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
     Liz Potterton Dec 2011 - Simple Qt Style Sheet
"""

import os
import sys
from core import CCP4Modules, CCP4Utils

from PySide2.QtGui import QPalette
from PySide2.QtWidgets import QLabel,QApplication

HIGHLIGHTCOLOUR = '#87C8D5'
LOWLIGHTCOLOUR ='#C1D8E0'
#LOWLIGHTCOLOUR = "#E5E5E5"
DARKCOLOUR = '#396169'
RUNCOLOUR = '#FFD777'

def setStyleSheet(app=None):
    if app is None:
        app = CCP4Modules.QTAPPLICATION()


    label = QLabel("Am I in the dark?")
    text_hsv_value = label.palette().color(QPalette.WindowText).value()
    bg_hsv_value = label.palette().color(QPalette.Background).value()
    isDarkMode = text_hsv_value > bg_hsv_value

    print("########################################")
    print("isDarkMode?",isDarkMode)
    print("########################################")
        
    preferences = CCP4Modules.PREFERENCES()
    #print 'setStyleSheet',preferences.INVALID_FRAME_WIDTH,preferences.INVALID_FRAME_COLOUR
    qticonsDir = os.path.join(CCP4Utils.getCCP4I2Dir(),'qticons').replace(os.path.sep,'/')
    # paddingAllowance is subtracted from complex widget width so they can be accomodated in the (notionally) 600px width
    # task input widget.  This is dependent on the padding taken up by some styling (eg maodern)
    # The CTabTaskWidget is actually >600 to allow for the padding around the various tab/scroll frames.
    paddingAllowance = 24
    screenWidth,screenHeight = app.screenSize()
    try:
        warningFrameWidth = min(int((preferences.INVALID_FRAME_WIDTH)/2),1)
    except:
        warningFrameWidth = 1
    print('StyleSheet screen size:',screenWidth,screenHeight)
    #print 'StyleSheet',qticonsDir
    if preferences.INVALID_FRAME_MODE == 0:
        styleSheet = '''QFrame [isValid="false"] { background-color: '''+str(preferences.INVALID_FRAME_COLOUR)+'''; }'''
    elif preferences.INVALID_FRAME_MODE == 1:
        styleSheet = '''QFrame [isValid="false"] { border: '''+str(preferences.INVALID_FRAME_WIDTH)+'''px solid '''+str(preferences.INVALID_FRAME_COLOUR)+'''; }'''
    else :
        paddingAllowance = 32
        #styleSheet = '''QFrame [isValid="false"] { background-image: url(qticons/incomplete.png); }'''
        styleSheet = '''QFrame { padding-left: 2px; }
                        QFrame [isValid="false"] { background: qlineargradient( x1:0.15 y1:0, x2:1 y2:0, stop:0 transparent, stop:1 #f55); border:0px; }
                        QFrame [isValid="false"] > QLineEdit { background: lemonchiffon; }
                        QFrame#tasksubframe {padding-top: 2px; padding-bottom: 2px; border:2px ridge; border-color: grey; border-radius:6px;} '''
        #51, 51, 51 40%
    styleSheet  += '''* { font-size: '''+str(preferences.GUI_FONT_SIZE)+'''px; }'''

    if not isDarkMode:
        if preferences.INVALID_FRAME_MODE == 2 :
            styleSheet += '''QFrame#jobHeaderFrame { background-image: url("'''+qticonsDir+'''/backgroundJobTitle.png"); padding-left: 10px;
                             padding-bottom: 5px; padding-top: 4px; padding-right: 10px; border: 1px solid; border-radius: 0px; border-color: #ec9; margin-bottom: 8px;}'''
        else :
            styleSheet += '''QFrame#followFrom { border: 2px solid '''+HIGHLIGHTCOLOUR+'''; }'''

    styleSheet +=     '''QFrame [hasWarning="true"] { border: '''+str(warningFrameWidth)+'''px solid '''+str(preferences.INVALID_FRAME_COLOUR)+'''; }
                         QFrame#highlight { border: 2px solid '''+HIGHLIGHTCOLOUR+'''; }
                         QStatusBar#statusWidget { border: 2px solid '''+HIGHLIGHTCOLOUR+'''; }
                         QLabel#emphasise { font-size: 14px; font-weight: bold }
                         QLabel#subtitle { color: #369; font-weight: bold; margin-top:12px; margin-bottom:4px; }
                         QLabel#italic { font-style: italic; padding-top: 4px; }
                         QLabel#bold { font-weight: bold }
                         QLabel#warning { font-weight: bold; color : red;}
                         QFrame#reportLinkFrame {  background-color: white }
                         QFrame#reportLinkFrame.QPushButton { foreground-color : blue; }
                         QTextEdit#readOnly { font-style: italic; }
                         QLabel#subtitle QToolTip { min-width: 300px; }
                         QFrame#taskTreeModule QLabel#title {  font-size: 14px; font-weight: bold }
                         QFrame#taskpipe  {  }
                         QFrame#tasktool  {  }
                         QFrame#taskpipe QLabel#icon { }
                         QFrame#tasktool QLabel#icon {  }
                         QFrame#taskpipe QLabel#title {  font-size: 14px; font-weight: bold }
                         QFrame#tasktool QLabel#title {  font-size: 14px; font-weight: bold }
                         QFrame#taskpipe QLabel#description {  font-size: 12px;  font-style: italic; }
                         QFrame#tasktool QLabel#description {  font-size: 12px;  font-style: italic; }
                         QFrame#taskpipe QLabel#title2 {  font-size: 14px; font-weight: bold }
                         QFrame#tasktool QLabel#title2 {  font-size: 14px; font-weight: bold }
                         QFrame#taskpipe QLabel#description2 {  font-size: 12px;  font-style: italic; }
                         QFrame#tasktool QLabel#description2 {  font-size: 12px;  font-style: italic; }

                         QLabel#bibreference_authorList { }
                         QLabel#bibreference_source {font-style: italic; }
                         QLabel#bibreference_selection { font-weight: bold; }
                         QCheckBox#bibreference_selection { font-weight: bold; }
                         QTextEdit#messageBox {  font: '''+str(preferences.GUI_FONT_SIZE)+'''px "Courier"; }
                         QLabel#errorMessage  { color : red; max-width: ''' + str(screenWidth-paddingAllowance) + '''px; }
                         CComplexLineWidget  {  max-width: ''' + str(screenWidth-paddingAllowance) + '''px; }
                         CTabTaskWidget  {  max-width: ''' + str(screenWidth-paddingAllowance) + '''px; min-width: 400px }
                         CTextEdit {  max-width: ''' + str(screenWidth-paddingAllowance) + '''px; min-width: 400px }
                         CStackedWidget {  max-width: ''' + str(screenWidth-paddingAllowance) + '''px; min-width: 400px }
                         QFrame#taskLine {  max-width: ''' + str(screenWidth-paddingAllowance) + '''px; }
                         CReportView {  max-width: '''+ str(screenWidth-50)+'''px;  min-width: 600px }
                         QLabel#jobListLabelRed  { color : red; }
                         QScrollArea#messageScrollArea  { max-width: ''' + str(screenWidth-paddingAllowance) + '''px; min-width: 600px }
                         QToolBar::separator { background: qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 #ccc, stop:1 #888); border: 1px solid #777; width: 4px; margin-top: 2px; margin-bottom: 2px; border-radius: 2px; }
                         QToolButton:pressed { background-color: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, stop: 0 #dadbde, stop: 1 #f6f7fa); }
                          '''

    if not isDarkMode :
       styleSheet +=    '''
                         QTreeView#taskTree {selection-background-color: '''+LOWLIGHTCOLOUR+''' ;  selection-color : black; }
                         QTreeView#projectWidget {selection-background-color: '''+LOWLIGHTCOLOUR+''' ;  selection-color : black; }
                         QTreeView#projectWidget::item:selected {background-color : '''+LOWLIGHTCOLOUR+''' ;  selection-color : black; }
                         QTreeView#projectWidget::branch:selected:!has-children  {background-color : '''+LOWLIGHTCOLOUR+''' ; }
                         QLabel#jobTitle { background-color:'''+LOWLIGHTCOLOUR+''';
                                          font-size: 14px; font-style: italic; font-weight: bold }
                         QLabel#jobStatus { background-color:'''+LOWLIGHTCOLOUR+''';
                                          font-size: 14px; font-style: italic; font-weight: bold }
                          '''

    if sys.platform == 'darwin':
        #Is it just me, or is Qt4's "grey-out" colour default a bit too subtle. (Qt5 seems better)
        #And maybe this could be for all widgets, not just QLabel?
        styleSheet += '''
                         QLabel:disabled { color: rgb(150, 150, 150); }
                         QTreeView:disabled{ color: rgb(150, 150, 150); }
                         '''
        if isDarkMode :
            styleSheet += """
QTabBar::tab {
    background: rgb(30, 30, 30);
    border: 1px solid #C4C4C3;
    border-bottom-color: #C2C7CB; /* same as the pane color */
    border-top-left-radius: 4px;
    border-top-right-radius: 4px;
    border-bottom-left-radius: 4px;
    border-bottom-right-radius: 4px;
    min-width: 8ex;
    padding: 2px;
    color: rgb(255, 255, 255);
}

QTabBar::tab:selected, QTabBar::tab:hover {
    background: rgb(10, 10, 10);
}

QTabBar::tab:selected {
    background: rgb(55, 55, 55);
    border-color: rgb(255, 255, 255);
    color: rgb(255, 255, 255);
}

QTabBar::tab:!selected {
    border-color: rgb(0, 0, 0);
    margin-top: 2px; /* make non-selected tabs look smaller */
}

QTabBar::close-button {
     image: url("""+qticonsDir+"""/tab-close.png);
     subcontrol-position: left;
 }
 QTabBar::close-button:hover {      image: url("""+qticonsDir+"""/tab-close-highlight.png)  }
            """
        else:
            styleSheet += """
QTabBar::tab {
    background: rgb(230, 230, 230);
    border: 1px solid #C4C4C3;
    border-bottom-color: #C2C7CB; /* same as the pane color */
    border-top-left-radius: 4px;
    border-top-right-radius: 4px;
    border-bottom-left-radius: 4px;
    border-bottom-right-radius: 4px;
    min-width: 8ex;
    padding: 2px;
}

QTabBar::tab:selected, QTabBar::tab:hover {
    background: rgb(210, 210, 210);
}

QTabBar::tab:selected {
    background: rgb(255, 255, 255);
    border-color: rgb(255, 255, 255);
    color: rgb(0, 0, 0);
}

QTabBar::tab:!selected {
    margin-top: 2px; /* make non-selected tabs look smaller */
}

QTabBar::close-button {
     image: url("""+qticonsDir+"""/tab-close.png);
     subcontrol-position: left;
 }
 QTabBar::close-button:hover {      image: url("""+qticonsDir+"""/tab-close-highlight.png)  }
            """

    app.setStyleSheet( styleSheet )
