from __future__ import print_function


"""
     CCP4ProjectViewerDemo.py: CCP4 GUI Project
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
   Liz Potterton Feb 2010 - Copied from earlier review.py used for Developers meeting demo
"""

##@package CCP4ProjectViewer (QtGui) Browser plugin for Project View demo                           
from PySide2 import QtGui, QtWidgets,QtCore
from qtgui import CCP4AbstractViewer,CCP4ProjectWidget



#------------------------------------------------------------------------------------------------------
def mimeType():
#------------------------------------------------------------------------------------------------------
#FIXME
      return None
      from qtcore import CCP4CustomMimeTypes            
      mimeType = CCP4CustomMimeTypes.CMimeType()
      mimeType.name = "application/ccp4-project"
      mimeType.description = "Review CCP4 project"
      mimeType.fileExtensions = ['CCP4project','ccp4_project']
      mimeType.viewers = [CProjectViewerDemo]
      return mimeType

class CCProjectViewerLayout( QtWidgets.QVBoxLayout):
  def __init__(self):
    QtWidgets.QVBoxLayout.__init__(self)
    self.setSizeConstraint(self.SetMaximumSize)

  def maximumSize(self):
    width,height = self.parentWidget().projectWidget.Size()
    print('CCProjectViewerLayout.maximumSize', width,height)
    return QtCore.QSize(width,height)

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
class CProjectViewerDemo(CCP4AbstractViewer.CAbstractViewer):
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

  MARGIN = 2
  
  def __init__(self,parent=None):
      CCP4AbstractViewer.CAbstractViewer.__init__(self,parent)
      self.projectWidget = None


#------------------------------------------------------------------------------------------------------
  def open(self,filename=''):
#------------------------------------------------------------------------------------------------------
    import os
    self.filename = filename
    self.setObjectName(os.path.splitext(os.path.basename(filename))[-1])
    layout = QtWidgets.QVBoxLayout()
    #layout = CCProjectViewerLayout()
    layout.setSpacing(CProjectViewerDemo.MARGIN)
    layout.setContentsMargins(CProjectViewerDemo.MARGIN,CProjectViewerDemo.MARGIN,CProjectViewerDemo.MARGIN,CProjectViewerDemo.MARGIN)
    layout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)
    self.projectWidget = CCP4ProjectWidget.CProjectWidget(self)
    self.projectWidget.setProjectData(filename)
    layout.addWidget(self.projectWidget)
    self.projectWidget.show()
    self.setLayout(layout)
    self.show()


#------------------------------------------------------------------------------------------------------
  def Size(self):
#------------------------------------------------------------------------------------------------------
    width,height = self.projectWidget.Size()
    '''
    Qt 4.6 will do this
    width = width + 2*self.layout().spacing() + self.layout().contentsMargins().left() + \
                                                self.layout().contentsMargins().right()
    height = height +  2*self.layout().spacing() + self.layout().contentsMargins().top() + \
                                                   self.layout().contentsMargins().bottom()
    '''
    width = width + 4* CProjectViewerDemo.MARGIN
    height = height +  4* CProjectViewerDemo.MARGIN
    #print 'CProjectViewer.Size',width,height
    return QtCore.QSize(width,height)

