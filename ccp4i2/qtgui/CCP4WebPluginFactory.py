from __future__ import print_function

"""
     CCP4WebPluginFactory.py: CCP4 GUI Project
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

"""
     Liz Potterton Jan 2010 - Simple demo CCP4WebPluginFactory and a test plugin
"""

import sys,os, functools
from PySide2 import QtWebEngine, QtWebEngineWidgets, QtGui, QtWidgets, QtCore
from core.CCP4ErrorHandling import *
from core.CCP4Config import DEVELOPER
from qtgui import CCP4ReportWidgets


##@package CCP4WebPluginFactory  (QtWebKit) Demo of webkit plgin widgets
# A trivial button to test PyQt/Javascript interface

class ccp4_test_plugin(QtWidgets.QWidget):
   def __init__(self,text='',parent=None,argumentNames=[],argumentValues=[]):
     QtWidgets.QWidget. __init__(self,parent)
     #for ii in range(len(argumentNames)): print 'ccp4_test_plugin',argumentNames[ii],argumentValues[ii]
     layout = QtWidgets.QGridLayout()
     self.label = QtWidgets.QLabel(text,self)
     self.label.setMinimumWidth(300)
     layout.addWidget(self.label,0,0)
     self.image = QtWidgets.QLabel(self)
     pixmap = QtGui.QPixmap(os.path.join(os.environ['CCP4I2'],'test','data','ramachandran.png')).scaled(200,200)
     self.image.setPixmap(pixmap)
     layout.addWidget(self.image,1,1)
     self.setLayout(layout)

   def handleImageClick(self):
     #print 'handleImageClick'
     pass

   @QtCore.Slot()
   def changeColour(self):
     #print 'changeColour'
     palette = QtGui.QPalette()
     palette.setColor(QtGui.QPalette.Normal,QtGui.QPalette.Window,QtGui.QColor(QtCore.Qt.cyan))
     self.setPalette(palette)
     
   @QtCore.Slot()
   def clearText(self):
     #print 'clearText'    
     self.label.setText('')

   @QtCore.Slot(str)
   def setLabelText(self, s):
     #print 'setLabelText',s
     self.label.setText(s)
     
   @QtCore.Slot(result=str)
   def value(self):
     #print 'value',self.label.text()
     return self.label.text()

   @QtCore.Slot()
   def update(self):
      # This is bad - assumes existance of another specific element on the html page
      # But does demo that can access html page
      frame = self.parent().page().mainFrame()
      qv = frame.evaluateJavaScript("test_text_input.value")
      #print 'ccp4_test_plugin.update qv',qv,qv.type(),qv
      self.label.setText(qv)
      #text_input = self.parent().page().mainFrame().findAllElements("test_text_input") - qt4.6
      #print "text_input",text_input.evaluateJavaScript("this.value")
      #self.setText(text_input.evaluateJavaScript("this.value"))
