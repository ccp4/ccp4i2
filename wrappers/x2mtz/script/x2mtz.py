"""
Copyright (C) 2015 STFC

Base class for importing reflection data provides processOutputFiles() method to split mtz
by either taking columns specified by HKLIN_OBS_COLUMNS and HKLIN_FREER_COLUMN or by finding best choice
automatically
"""

from lxml import etree

from core.CCP4PluginScript import CPluginScript
from core import CCP4Utils
from core.CCP4ErrorHandling import *


class x2mtz(CPluginScript):
    TASKNAME = 'x2mtz'
    MAINTAINER = 'liz.potterton@york.ac.uk'
    TASKMODULE = 'test'
    ERROR_CODES = {
      301: {'description': 'Input data file not found'},
      302: {'description': 'Failed automatic search for best reflection and freer data in converted file'}
    }

    def processOutputFiles(self):
      outputData = self.container.outputData
      inputData = self.container.inputData
      hklout = outputData.HKLOUT
      self.freeRcolumnLabel = None
      #print 'x2mtz.process',hklout
      if not hklout.exists():
        self.appendErrorReport(301)
        return CPluginScript.FAILED

      self.x2mtzXML = etree.Element('X2MTZ')

      iBestObs = -1
      iFree = -1
      if not inputData.HKLIN_OBS_COLUMNS.isSet() or not inputData.HKLIN_FREER_COLUMN.isSet():
        #Try to figure the 'best' obs data and freer data in the input file
        columnGroups = hklout.fileContent.getColumnGroups()
        for ii, columnGroup in enumerate(columnGroups):
          if columnGroup.columnGroupType == 'FreeR':
            iFree = ii
          elif columnGroup.columnGroupType == 'Obs':
            if iBestObs < 0 or columnGroup.contentFlag < columnGroups[iBestObs].contentFlag:
              iBestObs = ii
        # We must have Obs data, even if no FreeR data
        if iBestObs<0:
          self.appendErrorReport(302)
          return CPluginScript.FAILED

      mtzOut = []
      progCol = []
      if inputData.HKLIN_OBS_COLUMNS.isSet():
        mtzOut.append('OBSOUT')
        outputData.OBSOUT.contentFlag.set(inputData.HKLIN_OBS_CONTENT_FLAG)
        outputData.OBSOUT.annotation.set(outputData.OBSOUT.qualifiers('guiLabel')+' from '+inputData.HKLIN.stripedName())
        progCol.append(inputData.HKLIN_OBS_COLUMNS.__str__())
      elif iBestObs>=0:
        print ('Setting content flag',columnGroups[iBestObs].contentFlag,self.container.outputData.OBSOUT.columnNames(True))
        outputData.OBSOUT.contentFlag.set(columnGroups[iBestObs].contentFlag)
        outputData.OBSOUT.annotation.set(outputData.OBSOUT.qualifiers('guiLabel')+' from '+inputData.HKLIN.stripedName())
        mtzOut.append('OBSOUT')
        progCol.append(str(columnGroups[iBestObs].columnList[0].columnLabel))
        for col in columnGroups[iBestObs].columnList[1:]: progCol[-1] += ','+str(col.columnLabel)
                
      if inputData.HKLIN_FREER_COLUMN.isSet(allowDefault=False):
        mtzOut.append('FREEOUT')
        progCol.append(inputData.HKLIN_FREER_COLUMN.__str__())
        outputData.FREEOUT.annotation.set(outputData.FREEOUT.qualifiers('guiLabel')+' from '+inputData.HKLIN.stripedName())
      elif iFree>=0:
        outputData.FREEOUT.annotation.set(outputData.FREEOUT.qualifiers('guiLabel')+' from '+inputData.HKLIN.stripedName())
        mtzOut.append('FREEOUT')
        self.freeRcolumnLabel = str(columnGroups[iFree].columnList[0].columnLabel)
        progCol.append(self.freeRcolumnLabel)
        
      err = self.splitHklout(infile=str(hklout),programColumnNames=progCol,miniMtzsOut=mtzOut)

      # input names
      self.addElement(self.x2mtzXML, 'inputcolumnames',  progCol[0])

      # Get output column names for main file,
      #  cf CCP4PluginScript.splitHklout            
      mtzName = mtzOut[0]  # OBSOUT
      obj = self.container.outputData.get(mtzName)
      outputcolnames = obj.columnNames(True)  # as string
      self.addElement(self.x2mtzXML, 'outputcolumnnames', outputcolnames)
      
      with open (self.makeFileName('PROGRAMXML'),"w") as outputXML:
          CCP4Utils.writeXML(outputXML,etree.tostring(self.x2mtzXML,pretty_print=True))

      if err.maxSeverity()>SEVERITY_WARNING:
        print('ERROR in splitHklout')
        print(err.report())
        return CPluginScript.FAILED
      else:
        return CPluginScript.SUCCEEDED

    def addElement(self, containerXML, elementname, elementtext):
        #print 'addElement', elementname, type(elementtext), elementtext 
        e2 = etree.Element(elementname)
        e2.text = elementtext
        containerXML.append(e2)
