from __future__ import print_function


"""
     CCP4TextViewer.py: CCP4 GUI Project
     Copyright (C) 2009-2010 University of York

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
     Liz Potterton Jan 2010 - Create CCP4AbstractViewer
"""
##@package CCP4TextViewer (QtGui) Web browser plugin to view text files and coordinate files
from PySide2 import QtGui, QtWidgets,QtCore
from qtgui import CCP4AbstractViewer
from core.CCP4ErrorHandling import *
from core.CCP4Modules import MIMETYPESHANDLER,QTAPPLICATION,PREFERENCES

class CTextBrowser(QtWidgets.QPlainTextEdit):

  mouseRelease = QtCore.Signal(object)
  mouseDoubleClick = QtCore.Signal(object)

  def wheelEvent(self,e):
      if ((e.modifiers() & QtCore.Qt.ControlModifier) or (e.modifiers() & QtCore.Qt.MetaModifier) and hasattr(self,"parent")) and hasattr(self.parent(),"setZoomFactor"):
          if e.delta() > 0:
              self.parent().setZoomFactor(self.parent().zoomFactor()*1.2)
          else:
              self.parent().setZoomFactor(self.parent().zoomFactor()/1.2)
      else:
          QtWidgets.QPlainTextEdit.wheelEvent(self,e)

  def __init__(self,parent=None):
    QtWidgets.QPlainTextEdit.__init__(self,parent)

  def mouseReleaseEvent(self,event):    
    self.mouseRelease.emit(event.globalPos())
    event.accept()
    
  def mouseDoubleClickEvent(self,event):    
    self.mouseDoubleClick.emit(event.globalPos())
    event.accept()


#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
class CTextViewer(CCP4AbstractViewer.CAbstractViewer):
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

  lineClicked = QtCore.Signal(str)
  wordClicked = QtCore.Signal(str)
  lineDoubleClicked = QtCore.Signal(str)
  wordDoubleClicked = QtCore.Signal(str)

  MENUTEXT = 'Display text'

  ERROR_CODES = {}
  ERROR_CODES.update(CCP4AbstractViewer.CAbstractViewer.ERROR_CODES)
  ERROR_CODES.update( { 101 : { 'severity' : SEVERITY_ERROR,
                                'description' : 'Error loading text into document' }
                        } )

# Display plain text using the QTextBrowser widget.
# Could alternatively us the QWebView which seems to handle text OK
    
  def eventFilter(self, watched, event):
    if event.type() == QtCore.QEvent.ShortcutOverride:
        keyEvent = event
        # Ignore only the Ctrl + V shortcut override, you can customize check for your needs
        if ((keyEvent.modifiers() & QtCore.Qt.ControlModifier) == QtCore.Qt.ControlModifier) and keyEvent.key() == QtCore.Qt.Key_Plus:
            pointSize = watched.font().pointSize()
            pixelSize = watched.font().pixelSize()
            if pointSize > 0:
                watched.setStyleSheet("QPlainTextEdit {  font-family: Courier ;font-size: "+str(pointSize+1)+"pt;}");
            elif pixelSize > 0:
                watched.setStyleSheet("QPlainTextEdit {  font-family: Courier ;font-size: "+str(pixelSize+1)+"px;}");
            else:
                watched.setStyleSheet("QPlainTextEdit {  font-family: Courier ;font-size: 13px;}");
            event.accept()
            return True
        if ((keyEvent.modifiers() & QtCore.Qt.ControlModifier) == QtCore.Qt.ControlModifier) and keyEvent.key() == QtCore.Qt.Key_Minus:
            pointSize = watched.font().pointSize()
            pixelSize = watched.font().pixelSize()
            if pointSize > 0:
                watched.setStyleSheet("QPlainTextEdit {  font-family: Courier ;font-size: "+str(pointSize-1)+"pt;}");
            elif pixelSize > 0:
                watched.setStyleSheet("QPlainTextEdit {  font-family: Courier ;font-size: "+str(pixelSize-1)+"px;}");
            else:
                watched.setStyleSheet("QPlainTextEdit {  font-family: Courier ;font-size: 13px;}");
            event.accept()
            return True
    return CCP4AbstractViewer.CAbstractViewer.eventFilter(self,watched, event);

#-------------------------------------------------------------------
  def __init__(self,parent,fileName=None):
#-------------------------------------------------------------------
      #print 'TextViewer.init',self,fileName
      CCP4AbstractViewer.CAbstractViewer.__init__(self,parent)
      layout = QtWidgets.QVBoxLayout()
      layout.setSpacing(0)
      layout.setContentsMargins(0,0,0,0)
      #self.viewer = QtWidgets.QTextBrowser(self)
      self.viewer = CTextBrowser(self)
      self.viewer.setReadOnly(True)
      self.setFocusProxy(self.viewer)
      self.viewer.installEventFilter(self)
      self.viewer.mouseRelease.connect(self.handleMouseClick)
      self.viewer.mouseDoubleClick.connect(self.handleMouseDoubleClick)
      layout.addWidget(self.viewer)
      
      self.fixedWidthFont = False
      self.document = None
      self.setLayout(layout)
      if fileName: self.open(fileName)
      self._zoomFactor = 1.0

  @QtCore.Slot('QPos')
  def handleMouseClick(self,pos):
    print('handleMouseClick',self.verticalScrollBar().value())
    text_cursor = self.viewer.cursorForPosition(pos)
    text_position = text_cursor.position()
    text_cursor.select(QtGui.QTextCursor.LineUnderCursor)
    line = str(text_cursor.selectedText())
    self.lineClicked.emit(line)
  
    word = self.getWord(text_position)
    self.wordClicked.emit(word)

  @QtCore.Slot('QPos')
  def handleMouseDoubleClick(self,pos):
    text_cursor = self.viewer.cursorForPosition(pos)
    text_position = text_cursor.position()
    text_cursor.select(QtGui.QTextCursor.LineUnderCursor)
    line = str(text_cursor.selectedText())
    self.lineDoubleClicked.emit(line)
    
    word = self.getWord(text_position)
    self.wordDoubleClicked.emit(word)
    
  def getWord(self,position):
    if not self.document: return ''
    text = str(self.document.toPlainText())
    fpos = position-1
    lpos = position
    while fpos>=0 and not [' ','\n'].count(text[fpos]):
      fpos = fpos -1
    while lpos<len(text) and not [' ','\n'].count(text[lpos]):
      lpos = lpos +1
      
    #print 'mgTextViewer.getWord',text[position-5:position+5],fpos,lpos,text[fpos:lpos]
    return text[fpos+1:lpos]


#-------------------------------------------------------------------
  def findText(self,subString='',direction=1,caseSensitive=0,wrapAroundDocument=0):
#-------------------------------------------------------------------
    #print 'TextViewer.findText',subString
    # wrapAroundDocument not supported by qt widget - here for compatibility
    flag = 0
    found = False
    if direction<0:
        if caseSensitive:
          flag = QtGui.QTextDocument.FindBackward | QtGui.QTextDocument.FindCaseSensitively
        else:
          flag = QtGui.QTextDocument.FindBackward
    elif caseSensitive :
          flag = QtGui.QTextDocument.FindCaseSensitively
      #print 'mgTextViewer.find flag',flag
    if flag:
        found = self.viewer.find(subString,flag)
    else:
        found = self.viewer.find(subString)
    #print 'findText.found',found
    if not found and wrapAroundDocument:
      # Try going via the QTextDocument to search from the top
      if flag:
          cursor = self.viewer.document().find(subString,flag)
      else:
          cursor = self.viewer.document().find(subString)
      #print 'findText.cursor',cursor,cursor.position()
      if cursor.position()>=0:
        self.viewer.setTextCursor(cursor)
        found = True
    if found: self.viewer.ensureCursorVisible()
    return found

#-------------------------------------------------------------------
  def highlightText(self):
#-------------------------------------------------------------------
    pass

  def goToTop(self):
    # This did not work
    #self.viewer.verticalScrollBar().triggerAction(QtWidgets.QScrollBar.SliderToMinimum)
    self.viewer.moveCursor (QtGui.QTextCursor.Start)
    self.viewer.ensureCursorVisible() ;
    #print 'goToTop',self.viewer.verticalScrollBar().value(),self.viewer.verticalScrollBar().minimum(),self.viewer.verticalScrollBar().maximum()
           
#-------------------------------------------------------------------
  def open(self,fileName=None,lineLimit=None,**kw):
#-------------------------------------------------------------------     
      import os
      try:
        f = open(fileName,'r')
      except:
        raise CException(self.__class__,1,fileName)
      if lineLimit is None or lineLimit is NotImplemented:
        lineLimit = PREFERENCES().TEXT_VIEW_LINE_LIMIT
      if lineLimit is not None and lineLimit is not NotImplemented:
        #print 'reading limited lines',lineLimit
        text = ''
        ii = 1
        latest = f.readline()
        while ii < lineLimit and latest != '':
          text = text + latest
          latest = f.readline()
          ii += 1
        if ii>=lineLimit:
          text = text + '\nDISPLAY LIMITED TO '+str(  lineLimit)+' LINES\nThe limit can be changed in Preferences\n' 
      else:
        try:
          text = f.read()
          f.close()
        except:
          raise CException(self.__class__,2,fileName)
      
      CCP4AbstractViewer.CAbstractViewer.open(self,fileName)
      self.fixedWidthFont = kw.get('fixedWidthFont',None)
      format = MIMETYPESHANDLER().formatFromFileExt(fileName=fileName)
      if self.fixedWidthFont is None:
        if format is not None:
          self.fixedWidthFont =  MIMETYPESHANDLER().getMimeTypeInfo(format,'fixedWidthFont')
        else:
          self.fixedWidthFont = False
      print('CTextViewer.open',self.fileName,format,self.fixedWidthFont)
      rv = self.loadText(text)
      self.watchFile()
      self.goToTop()
      return rv

  def reload(self):
    scrollPosition = self.viewer.verticalScrollBar().value()
    self.viewer.clear()
    try:
        f = open(self.fileName,'r')
    except:
        raise CException(self.__class__,1,'filename:'+str(self.fileName))
    try:
        text = f.read()
        f.close()
    except:
        raise CException(self.__class__,2,'filename:'+str(self.fileName))
    rv = self.loadText(text)
    self.viewer.verticalScrollBar().setValue(scrollPosition)
    return rv
      
#-------------------------------------------------------------------
  def loadText(self,text,html=0,fixedWidthFont=None):
#-------------------------------------------------------------------
    """
    self.document=QtGui.QTextDocument(self.viewer)
    try:
      if html:
        self.document.setHtml(text)
      else:
        self.document.setPlainText(text)
    except:
      raise CException(self.__class__,101,'filename:'+str(self.fileName))
    if fixedWidthFont is not None: self.fixedWidthFont = fixedWidthFont
    if self.fixedWidthFont: self.setFont(style='fixed_width')
    self.viewer.setDocument(self.document)
    """
    self.viewer.setPlainText(text)
    self.viewer.setStyleSheet("QPlainTextEdit { font-family: Courier ;}")

    return 0


  def Print(self,painter):
      text = ''
      try:
        f = open(self.fileName,'r')
        text = f.read()
        f.close()
      except:
        return 1
      painter.drawText(0,0,text)

  def setZoomFactor(self,val):
      try:
          fn = QtGui.QFont(self.parent().font())
          fnd = self.document.defaultFont()
          fn.setFamily(fnd.family())
      except:
          return
      if fn.pointSize() > -1:
          if val*fn.pointSize() > 1.0:
              fn.setPointSize(int(val*fn.pointSize()))
          else:
              return
      elif fn.pixelSize():
          if val*fn.pixelSize() > 1.0:
              fn.setPixelSize(int(val*fn.pixelSize()))
          else:
              return
      else:
          return
      self._zoomFactor = val
      self.document.setDefaultFont(fn)
    
  def zoomFactor(self):
      return self._zoomFactor
    
#-------------------------------------------------------------------
  def Save(self,fileName):
#-------------------------------------------------------------------
    #print 'mgTextViewer.Save',self,self.filename,filename
    from core import CCP4Utils
    CCP4Utils.saveFile(fileName,str(self.document.toPlainText()))

#-------------------------------------------------------------------
  def isPrintable(self):
#-------------------------------------------------------------------
    return 1

#-------------------------------------------------------------------
  def isSaveable(self):
#-------------------------------------------------------------------
    return 1
  
#-------------------------------------------------------------------
  def isSearchable(self):
#-------------------------------------------------------------------
    return 1

#-------------------------------------------------------------------
  def setFont(self,font=None,style=''):
#-------------------------------------------------------------------
    #print 'TextViewer.setFont',self.document,style=='fixed_width'
    if self.document is None: return 
    self._zoomFactor = 1.0
    if font is not None:
      self.document.setDefaultFont(font)
    elif style=='fixed_width':
     self.viewer.setStyleSheet("QPlainTextEdit { font-family: Courier ;}")
    else:
      self.viewer.setStyleSheet("")
    return 0


#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
class CCoordsViewer(CTextViewer):
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

  MENUTEXT='Display as text'

  def __init__(self,parent,fileName=''):
      CTextViewer.__init__(self,parent)
      if fileName: self.open(fileName)

#-------------------------------------------------------------------
  def loadText(self,text,html=0):
#-------------------------------------------------------------------
    CTextViewer.loadText(self,text=text,html=html)
    self.setFont(style='fixed_width')
    


'''
#-------------------------------------------------------------------
  def open(self,filename):
#-------------------------------------------------------------------
      # This method only needed if want ot mark up thr PDB file
      status,text = self.markupPDB(filename)
      if not status:
        self.filename = os.path.abspath(filename)
        self.setObjectName(os.path.split(self.filename)[-1])
        return self.loadText(text,html=1)
      else:
        return status
      font = QtGui.QFont()
      font.setFixedPitch(1)
      self.viewer.setCurrentFont(font)


#-------------------------------------------------------------------
  def markupPDB(self,filename):
#-------------------------------------------------------------------
      # Add markup to a PDB file - very slow!!
      output = ''
      n = 0
      try:
          f = open(filename)
          lines = f.readlines()
          f.close()
      except:
          return [1,'Error reading file']
      try:
          for line in lines:
            if line[0:4] == 'ATOM':
              n = n + 1
              #output = output + '<br/><a href = "ccp4:///atom/'+str(n)+'">'+line.strip('\n')+'</a>\n'
              output = output + '<br/>' +line.strip('\n') +   '<a href = "ccp4:///atom/'+str(n)+'">go to</a>\n'
            else:
              output = output +'<br/>' + line
          print 'markupPDB',output
      except:
          return [1,'Error marking up file']

      output = '<html>\n' + output + '\n<html>i\n'
      return [0,output]
'''


#-------------------------------------------------------------------
class CMtzHeaderViewer(CTextViewer):
#-------------------------------------------------------------------

  MENUTEXT='Display as text'

  def __init__(self,parent,fileName=''):
    CTextViewer.__init__(self,parent)
    if fileName: self.open(fileName)

  def open(self,fileName=None,**kw):
    self.loadText('Awaiting nicer MTZ dump for '+fileName)
    import os
    self.fileName = os.path.abspath(fileName)
    self.setObjectName(os.path.split(fileName)[-1])



class CScriptViewerLogger(QtCore.QObject):

  newLines = QtCore.Signal(list)

  def __init__(self,parent):
    QtCore.QObject.__init__(self,parent)
    self.log = []
  def write(self, data):
    #import sys
    #print>>sys.__stdout__, 'logging',data
    self.newLines.emit([data])
    self.log.append(data)

    
class CScriptViewer(CTextViewer):

  def __init__(self,parent,fileName=None):
    from qtgui import CCP4Widgets
    CTextViewer.__init__(self,parent)
    self.layout().insertWidget(0,CCP4Widgets.CItalicLabel('Enter script:'))
    self.viewer.setReadOnly(False)
    self.stdoutViewer = QtWidgets.QTextEdit(self)
    self.stdoutViewer.setReadOnly(True)
    self.layout().addWidget(CCP4Widgets.CItalicLabel('Standard output and errors:'))
    self.layout().addWidget(self.stdoutViewer)
    if fileName is not None: self.open(fileName)
    
  def open(self,fileName=None,**kw):
    if fileName is not None: CTextViewer.open(self,fileName=fileName)

  def isRunable(self):
    return 1

  def run(self):
    script = self.viewer.toPlainText().__str__()
    #print '*****CScriptViewer.run',script
    if len(script)>0:
      import sys,traceback
      logger = CScriptViewerLogger(self)
      logger.newLines.connect(self.appendStdoutViewer)
      sys.stdout = logger
      sys.stderr = logger
      try:
        exec(script)
      except Exception as e:
        self.appendStdoutViewer(e.__str__().split('\n'))
        exc_type,exc_value,exc_tb = sys.exc_info()[:3]
        self.appendStdoutViewer(traceback.format_tb(exc_tb))
      sys.stdout = sys.__stdout__
      sys.stderr = sys.__stderr__
      #self.appendStdoutViewer(logger)

    
  @QtCore.Slot(list)
  def appendStdoutViewer(self,lines):
    #print 'appendStdoutViewer',logger.log
    if len(lines)<=0: return
    self.stdoutViewer.setReadOnly(False)
    text = self.stdoutViewer.toPlainText().__str__()
    if len(text)>0 and text[-1]!='\n': text = text + '\n'
    for line in lines:
      if line != '\n': text = text + line + '\n'
    self.stdoutViewer.setPlainText(text)
    self.stdoutViewer.setReadOnly(True)
