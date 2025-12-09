from __future__ import print_function

"""
     cif2mtz.py: CCP4 GUI Project
     Copyright (C) 2012 STFC
"""

import os
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.wrappers.x2mtz.script import x2mtz

class cif2mtz(x2mtz.x2mtz):

    TASKMODULE = 'test'      # Where this plugin will appear on the gui
    TASKTITLE = 'Import mmCIF reflection file' # A short title for gui menu
    TASKNAME = 'cif2mtz'   # Task name - should be same as class name
    TASKVERSION= 0.0               # Version of this plugin

    # used by the base class startProcess()
    TASKCOMMAND = 'cif2mtz'   # The command to run the executable
    # used by the base class makeCommandAndScript()
    COMLINETEMPLATE = None 
    COMTEMPLATE = None

    ERROR_CODES = { 301 : { 'description' : 'No output file found after cif2mtz conversion' },
                    302 : { 'description' : 'Output from cif2mtz does not contain recognised reflection or FreeR set data' }
                    }

    def makeCommandAndScript(self):

      inputData = self.container.inputData

      self.appendCommandLine(['HKLIN',inputData.HKLIN.fullPath])
      self.appendCommandLine(['HKLOUT',self.container.outputData.HKLOUT.fullPath])
      print('cif2mtz.makeCommandAndScript',inputData.SPACEGROUPCELL.spaceGroup,inputData.SPACEGROUPCELL.spaceGroup.isSet())
      if inputData.SPACEGROUPCELL.spaceGroup.isSet():
          self.appendCommandScript("SYMM '%s'"%str(inputData.SPACEGROUPCELL.spaceGroup))
      if inputData.SPACEGROUPCELL.cell.isSet():
          # angles default to 90 etc..
          inputData.SPACEGROUPCELL.cell.set(inputData.SPACEGROUPCELL.cell.fix(inputData.SPACEGROUPCELL.cell.get()))
          self.appendCommandScript("CELL %f %f %f %f %f %f"%
              (inputData.SPACEGROUPCELL.cell.a.__float__(),inputData.SPACEGROUPCELL.cell.b.__float__(),inputData.SPACEGROUPCELL.cell.c.__float__(),
              inputData.SPACEGROUPCELL.cell.alpha.__float__(),inputData.SPACEGROUPCELL.cell.beta.__float__(),inputData.SPACEGROUPCELL.cell.gamma))
      if inputData.CRYSTALNAME.isSet():
          self.appendCommandScript("XNAME '%s'"%str(inputData.CRYSTALNAME))
      if inputData.DATASETNAME.isSet():
          self.appendCommandScript("DNAME '%s'"%str(inputData.DATASETNAME))
      if inputData.MMCIF_SELECTED_BLOCK.isSet():
          self.appendCommandScript("BLOCK '%s'"%str(inputData.MMCIF_SELECTED_BLOCK))
      self.appendCommandScript('END')

      return 0

    '''
    def processOutputFiles(self):
      if not self.container.outputData.HKLOUT.exists():
        self.appendErrorReport(301)
        return CPluginScript.FAILED
      self.container.outputData.HKLOUT.annotation = self.container.outputData.HKLOUT.qualifiers('guiLabel')+' from ' + self.container.inputData.HKLIN.stripedName()

      if not self.container.controlParameters.SPLITMTZ: 
        return CPluginScript.SUCCEEDED
      
      columnGroups = self.container.outputData.HKLOUT.fileContent.getColumnGroups()
      iBestObs = -1
      iFree = -1
      for ii in range(len(columnGroups)):
        if columnGroups[ii].columnGroupType == 'FreeR':
          iFree = ii
        elif columnGroups[ii].columnGroupType == 'Obs':
          if iBestObs<0 or columnGroups[ii].contentFlag<columnGroups[iBestObs].contentFlag:
            iBestObs = ii
      #print 'processOutputFiles columnGroup indices',iFree,iBestObs
      if iBestObs<0 and iFree < 0:
        self.appendErrorReport(302)
        return CPluginScript.FAILED

      mtzOut = []
      progCol = []
      if iBestObs>=0:
          #print 'Setting content flag',columnGroups[iBestObs].contentFlag,self.container.outputData.OBSOUT.columnNames(True)
          self.container.outputData.OBSOUT.contentFlag = columnGroups[iBestObs].contentFlag
          self.container.outputData.OBSOUT.annotation = self.container.outputData.OBSOUT.qualifiers('guiLabel')+' from '+self.container.inputData.HKLIN.stripedName()
          mtzOut.append('OBSOUT')
          progCol.append(str(columnGroups[iBestObs].columnList[0].columnLabel))
          for col in columnGroups[iBestObs].columnList[1:]: progCol[-1] += ','+str(col.columnLabel)
      if iFree>=0:
          self.container.outputData.FREEOUT.annotation = self.container.outputData.FREEOUT.qualifiers('guiLabel')+' from '+self.container.inputData.HKLIN.stripedName()
          mtzOut.append('FREEOUT')
          progCol.append(str(columnGroups[iFree].columnList[0].columnLabel))
          
      #print 'processOutputFiles',mtzOut,progCol
      self.splitHklout(infile=str(self.container.outputData.HKLOUT),programColumnNames=progCol,miniMtzsOut=mtzOut)

      return CPluginScript.SUCCEEDED
    '''

