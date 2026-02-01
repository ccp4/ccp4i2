import json
import os
import sys
from collections import OrderedDict

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling, CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class mrbump_model_prep(CPluginScript):

    TASKNAME = 'mrbump_model_prep'
    TASKCOMMAND = 'mrbump'
    TASKTITLE = 'MrBUMP model preparation'
    MAINTAINER = 'stuart.mcnicholas@york.ac.uk'
    WHATNEXT = []
    ERROR_CODES = { 301 : { 'description' : 'Failed somehow in MrBUMP model prep' }
                    }

    def writeProgramXML(self,logDir):

        MRBUMPFIELDS = ['chainSource', 'coverage', 'eLLG', 'evalue', 'experiment', 'mgName', 'modelName', 'modelPDBfile', 'rank', 'resolution', 'score', 'seqID', 'source', 'sourceChainID', 'tarEnd', 'tarGroupEnd', 'tarGroupStart', 'tarStart', 'type']

        from mrbump.output import modelout
        mjson=modelout.Json()
        modelsJsonFile=os.path.join(logDir, "models.json")
        if os.path.isfile(modelsJsonFile):
             model_dict=mjson.readJson(modelsJsonFile)
             xmlroot = etree.Element('mrbump_model_prep')

             for k,v in model_dict.items():
                 model = etree.SubElement(xmlroot,"model")
                 for field in MRBUMPFIELDS:
                     ele = etree.SubElement(model,field)
                     ele.text = str(getattr(v,field))

             with open(modelsJsonFile) as f:
                 j = json.load(f, object_pairs_hook=OrderedDict)

             for k,v in j.items():
                 theKey = k
                 break
                 
             bestFile = model_dict[theKey].modelPDBfile
             bestModel = etree.SubElement(xmlroot,"bestModel")
             bestModel.text = str(bestFile)

             with open(str(self.makeFileName('PROGRAMXML')), 'w') as ostream:
                 CCP4Utils.writeXML(ostream,etree.tostring(xmlroot,pretty_print=True))

    def findOutputFileFromLog(self,logDir):

        from mrbump.output import modelout
        mjson=modelout.Json()
        modelsJsonFile=os.path.join(logDir, "models.json")
        if os.path.isfile(modelsJsonFile):
             model_dict=mjson.readJson(modelsJsonFile)
             with open(modelsJsonFile) as f:
                 j = json.load(f, object_pairs_hook=OrderedDict)

             for k,v in j.items():
                 theKey = k
                 break
                 
             newFile = model_dict[theKey].modelPDBfile
             return newFile
        else:
                 sys.write("Error: Can't find MrBUMP models json file:\n %s\n" % modelsJsonFile)

    def processInputFiles(self):
        from ccp4i2.core import CCP4XtalData
        error = None
        self.hklin = None
        dataObjects = []
        if not self.container.inputData.F_SIGF.isSet():
            return CPluginScript.SUCCEEDED
        #Append Observation with representation dependent on whether we are detwining on Is or not
        dataObjects += [['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]]
        #Include FreeRflag if called for
        if self.container.inputData.FREERFLAG.isSet():
            dataObjects += ['FREERFLAG']
        self.hklin,error = self.makeHklin(dataObjects)
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        else:
            return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):

      inp = self.container.inputData

      keyin = "GESMAX 1\n" 
      keyin += "PICKLE False\n" 
      keyin += "MDLC True\n" 
      keyin += "MDLM False\n" 
      keyin += "MDLS False\n" 

      if inp.SEARCH_AFDB and inp.AFDBLEVEL != 'None':
          keyin += "AFLEVEL "+str(inp.AFDBLEVEL)+"\n" 
      else:
          keyin += "RLEVEL all\n" 

      if inp.SEARCH_PDB and inp.REDUNDANCYLEVEL:
          if str(inp.REDUNDANCYLEVEL) == '110':
             keyin += "RLEVEL ALL\n" 
          elif str(inp.REDUNDANCYLEVEL) == '100':
             keyin += "RLEVEL 100\n" 
          elif str(inp.REDUNDANCYLEVEL) == '95':
             keyin += "RLEVEL 95\n" 
          elif str(inp.REDUNDANCYLEVEL) == '90':
             keyin += "RLEVEL 90\n" 
          elif str(inp.REDUNDANCYLEVEL) == '70':
             keyin += "RLEVEL 70\n" 
          elif str(inp.REDUNDANCYLEVEL) == '50':
             keyin += "RLEVEL 50\n" 
      elif inp.SEARCH_PDB:
         keyin += "RLEVEL 95\n" 

      if inp.MRMAX:
          keyin += "MRNUM %s\n" % str(inp.MRMAX)
      else:
          keyin += "MRNUM 10\n" 

      keyin += "EXIT True\n" 

      keyfile = self.workDirectory / "keywords.txt"
      with keyfile.open("w") as kf:
          kf.write(keyin)

      if self.container.inputData.F_SIGF.isSet():
          self.appendCommandLine( [ 'HKLIN', self.hklin ] )
      
      seqFile = self.workDirectory / 'SEQIN.fasta'
      inp.ASUIN.writeFasta(fileName=seqFile)
      self.appendCommandLine( [ 'SEQIN', seqFile ] )
      self.appendCommandLine( [ 'KEYIN', keyfile ] )

      return CPluginScript.SUCCEEDED

    def processOutputFiles(self):

        log_dir = str(self.workDirectory / "search_mrbump_1" / "logs")
        self.writeProgramXML(log_dir)
        xyzout = self.findOutputFileFromLog(log_dir)
        if os.path.exists(xyzout):
            self.container.outputData.XYZOUT.set(xyzout)
            self.container.outputData.XYZOUT.annotation.set('Model from MrBump model search')

        return CPluginScript.SUCCEEDED
