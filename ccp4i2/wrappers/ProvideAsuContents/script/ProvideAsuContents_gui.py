from __future__ import print_function

"""
    wrappers/ProvideSequence/script/ProvideSequence_gui.py
    Martin Noble
    """

from PySide2 import QtGui, QtWidgets,QtCore,QtWebEngine, QtWebEngineWidgets
from qtgui.CCP4TaskWidget import CTaskWidget
from qtgui import CCP4Widgets
import os
import functools

class CTaskProvideAsuContents(CTaskWidget):
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'ProvideAsuContents'
    TASKVERSION = 0.0
    TASKMODULE='data_entry'
    TASKTITLE="Define AU contents"
    WHATNEXT = []
    DESCRIPTION = '''Define AU contents, estimate molecular weight, Matthews probability from list of sequences'''
    def isValid(self):
        inp = self.container.inputData
        requires = [inp.ASU_CONTENT,inp.HKLIN]
    
        invalidElements = super(CTaskProvideAsuContents, self).isValid()
        if self.matthewsInvalid:
            for require in requires:
                invalidElements.append(require)
            
        return invalidElements

    def drawContents(self):
        
      self.setProgramHelpFile('ProvideAsuContents')
      folder = self.openFolder(folderFunction='inputData',title='Enter sequences',followFrom=False)

      self.createLine(['advice','Optionally load existing AU content file to edit'])
      self.createLine(['widget','ASUCONTENTIN'])
      self.container.inputData.ASUCONTENTIN.dataChanged.connect(self.loadInputContent)
      self.container.inputData.ASUCONTENTIN.dataChanged.connect(functools.partial(self.loadInputContent,False))
      self.createLine(['subtitle','Specify the protein/nucleic acid sequences in the crystal','Provide a list of sequences, number of copies of each of them and experimental data. When your assumed composition leads to a sensible Matthews volume (shown in bottom panel), Run this task to commit the result'])
      self.createLine(['advice','Solvent content analysis will be done if you provide experimental data. Run the task to commit to database'])
      self.openSubFrame(frame=[True])
      self.createLine(['widget','ASU_CONTENT'])
      self.closeSubFrame()
      self.getWidget('ASUCONTENTIN').hideDefineWidget()
      self.getWidget('ASU_CONTENT').setMinimumHeight(200)

      self.createLine( ['subtitle' , 'Solvent content analysis', 'Select a reflection file containing the cell parameters' ] )
      self.openSubFrame(frame=[True])
      self.createLine( [ 'widget', 'HKLIN' ] )
      #self.resultWidget = QtWidgets.QTextEdit(self)
      self.resultWidget = QtWebEngineWidgets.QWebEngineView(self)
      self.resultWidget.setMinimumHeight(100)
      #self.resultWidget.setReadOnly(True)
      self.widget.currentFolderLayout.addWidget(self.resultWidget)
      self.closeSubFrame()

      self.container.inputData.ASUCONTENTIN.dataChanged.connect(self.runMatthews)
      self.container.inputData.ASU_CONTENT.dataChanged.connect(self.runMatthews)
      self.container.inputData.HKLIN.dataChanged.connect(self.runMatthews)
      self.runMatthews()

    @QtCore.Slot()
    def runMatthews(self):
      text = '<html>'
      text += '<style>'
      text += 'body {font-family: %s}' % ("Arial")
      text += '</style>'
      text += '<body>'

      #self.resultWidget.clear()
      self.resultWidget.setHtml("")
      self.matthewsInvalid = False
      if not self.container.inputData.HKLIN.isSet() or len(self.container.inputData.ASU_CONTENT)==0:
          text += "<p>Solvent content analysis will appear here when there is a valid sequence list and reflection file above.</p><p>When your assumed composition leads to a sensible Matthews volume, Run this task to commit the result.</p>"

      totWeight = 0.0

      polymerMode = ""
      text += "<table><tr><th style=\"text-align:left\">Name</th><th>Number of copies</th><th>Molecular weight</th></tr>"
      for seqObj in self.container.inputData.ASU_CONTENT:
          if seqObj.nCopies > 0:
              if seqObj.polymerType == "PROTEIN":
                  if polymerMode == "D":
                      polymerMode = "C"
                  elif polymerMode == "":
                      polymerMode = "P"
              if seqObj.polymerType in ["DNA","RNA"]:
                  if polymerMode == "P":
                      polymerMode = "C"
                  elif polymerMode == "":
                      polymerMode = "D"
          totWeight = totWeight + seqObj.molecularWeight(seqObj.polymerType)
          text = text + '<tr><td>{0}</td><td>{1} </td><td>{2:.2f}</td></tr>'.format(seqObj.name,seqObj.nCopies,float(seqObj.molecularWeight(seqObj.polymerType))) + '\n'
      text += "</table><br/>"

      if totWeight < 1e-6:
          self.resultWidget.setHtml('<b style="color: red;">Matthew\'s analysis suggests the current composition and cell volume are incompatible.</b>.')
          return

      text +=  '<br>Total sequence weight: {0:.2f}</br>'.format(float(totWeight))

      if not self.container.inputData.HKLIN.isSet() or len(self.container.inputData.ASU_CONTENT)==0:
          text += "</body></html>"
          self.resultWidget.setHtml(text)
          return

      rv = self.container.inputData.HKLIN.fileContent.matthewsCoeff(molWt=totWeight,polymerMode=polymerMode)
      vol = rv.get('cell_volume','Unkown')
      nmol=[]
      solv = []
      matt=[]
      prob=[]
      if vol == 'Unkown':
          text = text + '<p>Cell volume = Unknown</p>'
      else:
          text = text + '<p>Cell volume = {0:.1f}<p>\n'.format(float(vol)) 
      headText = ""
      if len(rv.get('results',[])) == 0:
          headText = '<b style="color: red;">Matthews analysis suggests the current composition and cell volume are incompatible.</b>'
          self.matthewsInvalid = True
      else:
          headText = headText +    '<table><tr><th>  Nmol  </th><th>  %solvent  </th><th>  Matthews  </th><th>  prob(Matthews)  </th></tr>\n'
      for result in rv.get('results',[]):
          headText = headText + '<tr><td>  {0}  </td><td>  {1:.2f}  </td><td>  {2:.2f}  </td><td>  {3:.2f}  </td></tr>'.format(result.get('nmol_in_asu'),result.get('percent_solvent'),result.get('matth_coef'),result.get('prob_matth')) + '\n'
          nmol.append(result.get('nmol_in_asu'))
          solv.append(result.get('percent_solvent'))
          matt.append(result.get('matth_coef'))
          prob.append(result.get('prob_matth'))
      headText = headText +    '</table><br/>'
      text = headText + text + '</body></html>'
      
      self.resultWidget.setHtml(text)

    @QtCore.Slot(bool)
    def loadInputContent(self,loadFile=True):
      print('CTaskProvideAsuContents.loadInputContent',self.container.inputData.ASUCONTENTIN)
      if self.container.inputData.ASUCONTENTIN.exists():
        if loadFile: self.container.inputData.ASUCONTENTIN.loadFile()
        self.container.inputData.ASU_CONTENT.blockSignals(True)
        self.container.inputData.ASU_CONTENT.set(self.container.inputData.ASUCONTENTIN.fileContent.seqList.get())
        self.getWidget('ASU_CONTENT').updateViewFromModel()
        self.container.inputData.ASU_CONTENT.blockSignals(False)
#This forces first line to be selected and hence updates other parts of GUI.
        try:
            self.getWidget('ASU_CONTENT').listWidget.currentRowChanged.emit(-1)
        except:
            pass
        try:
            self.getWidget('ASU_CONTENT').listWidget.currentRowChanged.emit(0)
        except:
            pass
        self.getWidget('ASU_CONTENT').validate()
        """
        if len(self.getWidget('ASU_CONTENT').model)<2:
           self.getWidget('ASU_CONTENT').setListVisible(0)
        else:
           self.getWidget('ASU_CONTENT').setListVisible(1)
        """

    def fix(self):
      #print 'CTaskProvideAsuContents.fix'
      
      report =  CTaskWidget.fix(self)
      for seqObj in self.container.inputData.ASU_CONTENT:
        #print 'CTaskProvideAsuContents.fix seq',str(seqObj.sequence)
        clean,err = seqObj.cleanupSequence(str(seqObj.sequence))
        #print 'CTaskProvideAsuContents.fix clean',clean
        seqObj.sequence.set(clean)
      return report
