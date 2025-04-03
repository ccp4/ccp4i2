"""
Copyright (C) 2009-2010 University of York
Liz Potterton Jan 2010 - Create CCP4AbstractViewer
"""

##@package CCP4ImageViewer  (QtGui) Web browser widget for image files

import os

from PySide2 import QtCore, QtGui, QtWidgets

from . import CCP4AbstractViewer
from ..core.CCP4ErrorHandling import CException, Severity


class CImageViewer(CCP4AbstractViewer.CAbstractViewer):

  ERROR_CODES = {}
  ERROR_CODES.update(CCP4AbstractViewer.CAbstractViewer.ERROR_CODES)
  ERROR_CODES.update( { 101 : { 'severity' : Severity.ERROR,
                                'description' : 'Image file apparently zero size' }
                         } )

  def __init__(self,parent,fileName=''):
      CCP4AbstractViewer.CAbstractViewer.__init__(self,parent)
      layout = QtWidgets.QVBoxLayout()
      self.label = QtWidgets.QLabel(self)
      layout.addWidget(self.label)
      self.setLayout(layout)
      self.fileName = None
      self.pixmap = QtGui.QPixmap()
      if fileName: self.open(fileName)
      self.isZoomed = 0

  def fitToWindow(self):
        #print 'mgImageViewerfitToWindow'
        #if self.pixmap.width()>self.width() or self.pixmap.height()>self.height():
        #self.isZoomed = 1 - self.isZoomed
        if not self.isZoomed:
          self.label.setPixmap(self.pixmap.scaled(self.size(),QtCore.Qt.KeepAspectRatio,QtCore.Qt.SmoothTransformation))
        else:
          self.label.setPixmap(self.pixmap)

  def mousePressEvent(self,event):
    if self.filename != '' and self.pixmap.width()>0 and self.pixmap.height()>0:
        try:
            drag = QtGui.QDrag(self)
            mimeData = QtCore.QMimeData()
            mimeData.setText(self.filename)
            #mimeData.setImageData(QtGui.QImage(self.filename))
            mimeData.setUrls([QtCore.QUrl.fromLocalFile(self.filename)])
            drag.setMimeData(mimeData)
            drag.setPixmap(self.pixmap.scaledToHeight(150))
            drag.start()
        except:
            print("Drag and drop failed")

  def open(self,fileName):
      try:
        self.pixmap = QtGui.QPixmap(fileName)
      except:
        raise CException(self.__class__,1,fileName)

      #print "mgImageViewer pixmap",self.pixmap,self.pixmap.width(),self.pixmap.height()

      w = int(self.pixmap.width())
      if w<=0:
        raise CException(self.__class__,101,fileName)
      
      self.label.setPixmap(self.pixmap)
      CCP4AbstractViewer.CAbstractViewer.open(self,fileName)

      self.label.show()
      #self.setWidget(self.label)

      return 0

  def mgSize(self):
    return self.pixmap.size()

  def Print(self,painter):
    painter.drawPixmap(0,0,self.pixmap)

  def Save(self,fileName):
    pass

  def isPrintable(self):
    return 1

  def isSaveable(self):
    return 1

  def isScaleable(self):
    return 1

  def scale(self,value):
    print('ImageViewer.scale',value)
    if not self.pixmap or not self.fileName: return
    self.label.setPixmap(self.pixmap.scaled(self.pixmap.width()*value,self.pixmap.height()*value,QtCore.Qt.KeepAspectRatio,QtCore.Qt.SmoothTransformation))


class CAnimationViewer(CCP4AbstractViewer.CAbstractViewer):

# Use QMovie to display gif or png animations
# Beware does NOT display mp4 etc. movies
  ERROR_CODES = {}
  ERROR_CODES.update(CCP4AbstractViewer.CAbstractViewer.ERROR_CODES)
  ERROR_CODES.update( { 101 : { 'severity' : Severity.ERROR,
                                'description' : 'Apparently invalid movie file' }
                         } )

  def __init__(self,parent,fileName=''):
      CCP4AbstractViewer.CAbstractViewer.__init__(self,parent)
      self.fileName=''
      layout = QtWidgets.QVBoxLayout()
      self.label = QtWidgets.QLabel(self)
      layout.addWidget(self.label)
      self.setLayout(layout)
      if fileName: self.open(fileName)

  def open(self,fileName):
    try:
      self.movie = QtGui.QMovie(fileName,'',self)
    except:
      raise CException(self.__class__,1,fileName)

    if not self.movie.isValid():
      raise CException(self.__class__,101,filenName)
  
    self.fileName = os.path.abspath(fileName)
    self.setObjectName(os.path.split(self.fileName)[-1])

    self.label.setMovie(self.movie)
    self.movie.start()

    return 0

  def mgSize(self):
    return self.movie.currentPixmap().size()
