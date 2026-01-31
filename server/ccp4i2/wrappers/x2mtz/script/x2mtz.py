'''
Base class for importing reflection data provides processOutputFiles() method to split mtz
by either taking columns specified by HKLIN_OBS_COLUMNS and HKLIN_FREER_COLUMN or by finding best choice
automatically
'''

from lxml import etree

from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4ErrorHandling import SEVERITY_WARNING
from ccp4i2.core.CCP4PluginScript import CPluginScript


class x2mtz(CPluginScript):

    TASKNAME = 'x2mtz'
    MAINTAINER = 'liz.potterton@york.ac.uk'
    TASKMODULE = 'test'
    ERROR_CODES = { 301 : { 'description' : 'Input data file not found' },
                    302 : { 'description' : 'Failed automatic search for best reflection and freer data in converted file' }
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

      #print '\nx2mtz content', inputData.HKLIN.getFileContent()

      iBestObs = -1
      iFree = -1
      #print 'inputData.HKLIN_OBS_COLUMNS', inputData.HKLIN_OBS_COLUMNS
      if not inputData.HKLIN_OBS_COLUMNS.isSet() or not inputData.HKLIN_FREER_COLUMN.isSet():
        #Try to figure the 'best' obs data and freer data in the input file
        columnGroups = hklout.fileContent.getColumnGroups()
        #        for cg in columnGroups:
        #           print 'processOutputFiles',cg.get()
        #           print "* type",cg.columnGroupType
        for ii in range(len(columnGroups)):
          if columnGroups[ii].columnGroupType == 'FreeR':
            iFree = ii
          elif columnGroups[ii].columnGroupType == 'Obs':
            #print "@cg ", columnGroups[ii].contentFlag, columnGroups[iBestObs].contentFlag
            if iBestObs<0 or columnGroups[ii].contentFlag<columnGroups[iBestObs].contentFlag:
              #print "select",ii
              iBestObs = ii
        #print 'x2mtz.processOutputFiles best columnGroup indices',iFree,iBestObs
        #if iBestObs<0 and iFree < 0:
        #  We must have Obs data, even if no FreeR data
        if iBestObs<0:
          self.appendErrorReport(302)
          return CPluginScript.FAILED

      #print ('x2mtz.processOutputFiles HKLIN_OBS_COLUMNS',inputData.HKLIN_OBS_COLUMNS,type(inputData.HKLIN_OBS_COLUMNS),inputData.HKLIN_OBS_COLUMNS.isSet(),inputData.HKLIN_FREER_COLUMN.isSet())

      mtzOut = []
      progCol = []
      if inputData.HKLIN_OBS_COLUMNS.isSet():
        mtzOut.append('OBSOUT')
        # Get content flag - use provided value or derive from column groups
        obs_content_flag = inputData.HKLIN_OBS_CONTENT_FLAG
        # Extract int value if it's a CInt
        while hasattr(obs_content_flag, 'value'):
            obs_content_flag = obs_content_flag.value
        obs_content_flag = int(obs_content_flag) if obs_content_flag else 0

        # If content flag is 0 or unset, try to derive from column groups
        if obs_content_flag == 0:
            col_names = inputData.HKLIN_OBS_COLUMNS.__str__().split(',')
            for cg in columnGroups:
                if cg.columnGroupType == 'Obs':
                    cg_labels = [str(c.columnLabel) for c in cg.columnList]
                    if cg_labels == col_names:
                        cf = cg.contentFlag
                        while hasattr(cf, 'value'):
                            cf = cf.value
                        obs_content_flag = int(cf) if cf else 0
                        break

        outputData.OBSOUT.contentFlag.set(obs_content_flag)
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
        
      #print ('x2mtz to splitHklout',mtzOut,progCol)
      err = self.splitHklout(infile=str(hklout),programColumnNames=progCol,miniMtzsOut=mtzOut)

      # input names
      self.addElement(self.x2mtzXML, 'inputcolumnames',  progCol[0])

      # Get output column names for main file,
      #  cf CCP4PluginScript.splitHklout            
      mtzName = mtzOut[0]  # OBSOUT
      #print("OP things", mtzName, self.container.outputData)
      obj = self.container.outputData.get(mtzName)
      outputcolnames = obj.columnNames(True)  # as string
      #print('outputcolnames', outputcolnames)
      self.addElement(self.x2mtzXML, 'outputcolumnnames', outputcolnames)
      
      with open (self.makeFileName('PROGRAMXML'),"w") as outputXML:
          #print("*x2mtz write XML to ",outputXML)
          CCP4Utils.writeXML(outputXML,etree.tostring(self.x2mtzXML,pretty_print=True))

      #print( 'x2mtz splitHklout err',err)
      if err.maxSeverity()>SEVERITY_WARNING:
        print('ERROR in splitHklout')
        print(err.report())
        return CPluginScript.FAILED
      else:
        return CPluginScript.SUCCEEDED
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    def addElement(self, containerXML, elementname, elementtext):
        #print 'addElement', elementname, type(elementtext), elementtext 
        e2 = etree.Element(elementname)
        e2.text = elementtext
        containerXML.append(e2)
